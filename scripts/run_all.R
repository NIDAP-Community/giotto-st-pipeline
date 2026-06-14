#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cli)
  library(rlang)
  library(jsonlite)
  library(Matrix)
})

# ---- helpers ---------------------------------------------------------------

timestamp_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_run_parameters <- function(output_dir, params) {
  meta_dir <- file.path(output_dir, "metadata")
  ensure_dir(meta_dir)
  path <- file.path(meta_dir, "run_parameters.json")
  jsonlite::write_json(params, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}

write_provenance <- function(output_dir, provenance) {
  meta_dir <- file.path(output_dir, "metadata")
  ensure_dir(meta_dir)
  path <- file.path(meta_dir, "provenance.json")
  jsonlite::write_json(provenance, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}

nonempty_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(NULL)
  }
  value
}

command_output <- function(command, args = character(), wd = NULL) {
  if (is.null(Sys.which(command)) || !nzchar(Sys.which(command))) {
    return(NULL)
  }
  old_wd <- NULL
  if (!is.null(wd)) {
    if (!dir.exists(wd)) {
      return(NULL)
    }
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(wd)
  }
  out <- tryCatch(
    suppressWarnings(system2(command, args, stdout = TRUE, stderr = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out) || length(out) == 0 || !is.null(attr(out, "status"))) {
    return(NULL)
  }
  paste(out, collapse = "\n")
}

sha256_file <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }

  sha <- command_output("sha256sum", path)
  if (!is.null(sha)) {
    return(strsplit(sha, "\\s+")[[1]][1])
  }

  sha <- command_output("shasum", c("-a", "256", path))
  if (!is.null(sha)) {
    return(strsplit(sha, "\\s+")[[1]][1])
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }

  NULL
}

read_pipeline_version <- function() {
  repo_root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  version_file <- file.path(repo_root, "VERSION")
  if (file.exists(version_file)) {
    return(trimws(readLines(version_file, n = 1, warn = FALSE)))
  }
  nonempty_env("GIOTTO_PIPELINE_VERSION") %||% "unknown"
}

collect_provenance <- function() {
  repo_root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  git_commit <- command_output("git", c("-C", repo_root, "rev-parse", "HEAD"))
  git_status <- command_output("git", c("-C", repo_root, "status", "--porcelain"))
  renv_lock_path <- file.path(repo_root, "renv.lock")

  source_dirty <- NULL
  if (!is.null(git_status)) {
    source_dirty <- nzchar(git_status)
  }

  pipeline_version <- read_pipeline_version()
  source_commit <- nonempty_env("GIOTTO_PIPELINE_SOURCE_COMMIT") %||% git_commit %||% "unknown"
  lock_sha <- nonempty_env("GIOTTO_PIPELINE_RENV_LOCK_SHA256") %||% sha256_file(renv_lock_path)

  # Warn on unresolved provenance fields
  if (identical(pipeline_version, "unknown")) {
    cli::cli_warn("Pipeline version could not be determined (no VERSION file or GIOTTO_PIPELINE_VERSION env var).")
  }
  if (identical(source_commit, "unknown")) {
    cli::cli_warn("Source commit unknown: not a git repo and GIOTTO_PIPELINE_SOURCE_COMMIT not set.")
  }
  if (is.null(lock_sha)) {
    cli::cli_warn("renv.lock SHA-256 could not be computed (file missing or no hash tool available).")
  }

  list(
    pipeline = list(
      name = "giotto-st-pipeline",
      version = pipeline_version,
      source_commit = source_commit,
      source_dirty = source_dirty,
      renv_lock_sha256 = lock_sha,
      renv_lock_path = if (file.exists(renv_lock_path)) renv_lock_path else NULL
    ),
    container = list(
      runtime = nonempty_env("GIOTTO_PIPELINE_CONTAINER_RUNTIME"),
      image = nonempty_env("GIOTTO_PIPELINE_CONTAINER_IMAGE"),
      image_digest = nonempty_env("GIOTTO_PIPELINE_CONTAINER_DIGEST"),
      sif_path = nonempty_env("GIOTTO_PIPELINE_SIF_PATH"),
      sif_sha256 = nonempty_env("GIOTTO_PIPELINE_SIF_SHA256")
    ),
    runtime = list(
      r_version = R.version.string,
      platform = R.version$platform,
      timestamp_utc = timestamp_utc()
    )
  )
}

abort_missing <- function(flag, path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(invisible(TRUE))
  if (!file.exists(path) && !dir.exists(path)) {
    cli::cli_abort("{.var {flag}} not found: {.path {path}}")
  }
  invisible(TRUE)
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", cmd_args[grep(file_flag, cmd_args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path)))
  }
  getwd()
}

source_helper <- function(rel_path) {
  base_dir <- get_script_dir()
  candidate <- normalizePath(file.path(base_dir, "..", rel_path), mustWork = FALSE)
  if (!file.exists(candidate)) {
    cli::cli_abort("Required helper not found: {.path {rel_path}} (looked in {.path {candidate}})")
  }
  source(candidate, local = FALSE)
  invisible(candidate)
}

split_prefixes <- function(value) {
  if (is.null(value)) {
    return(character())
  }
  chars <- as.character(value)
  chars <- chars[!is.na(chars)]
  if (length(chars) == 0) {
    return(character())
  }
  pieces <- unlist(strsplit(chars, ",", fixed = FALSE), use.names = FALSE)
  pieces <- trimws(pieces)
  pieces <- pieces[nzchar(pieces)]
  unique(pieces)
}

load_support_files <- function() {
  source_helper("R/ingest_xenium.R")
  source_helper("R/ingest_visium_hd.R")
  source_helper("R/utils_sparse.R")
  source_helper("R/ingest_h5ad.R")
  source_helper("R/pipeline_basic.R")
  invisible(TRUE)
}

normalize_params <- function(params) {
  params$input_format <- tolower(params$input_format %||% "auto")
  params$stage <- tolower(params$stage %||% "all")
  params$output_dir <- params$output_dir %||% "/output"
  params$project_id <- params$project_id %||% "giotto-st"
  params$dry_run <- isTRUE(params$dry_run)
  params$verbose <- isTRUE(params$verbose)
  params$cores <- as.integer(params$cores %||% 4)
  if (is.na(params$cores) || params$cores < 1) {
    cli::cli_abort("Invalid {.var --cores} value: {params$cores}")
  }

  params$spatial_point_size <- suppressWarnings(as.numeric(params$spatial_point_size %||% 2.25))
  if (is.na(params$spatial_point_size) || params$spatial_point_size <= 0) {
    cli::cli_abort("Invalid {.var --spatial_point_size} value: {params$spatial_point_size}. Use a positive number.")
  }
  params$umap_point_size <- suppressWarnings(as.numeric(params$umap_point_size %||% 1.5))
  if (is.na(params$umap_point_size) || params$umap_point_size <= 0) {
    cli::cli_abort("Invalid {.var --umap_point_size} value: {params$umap_point_size}. Use a positive number.")
  }
  params$spatial_legend_text <- suppressWarnings(as.numeric(params$spatial_legend_text %||% 12))
  if (is.na(params$spatial_legend_text) || params$spatial_legend_text <= 0) {
    cli::cli_abort("Invalid {.var --spatial_legend_text} value: {params$spatial_legend_text}. Use a positive number.")
  }
  params$spatial_legend_symbol_size <- suppressWarnings(as.numeric(params$spatial_legend_symbol_size %||% 1.4))
  if (is.na(params$spatial_legend_symbol_size) || params$spatial_legend_symbol_size <= 0) {
    cli::cli_abort("Invalid {.var --spatial_legend_symbol_size} value: {params$spatial_legend_symbol_size}. Use a positive number.")
  }
  params$spatial_axis_text <- suppressWarnings(as.numeric(params$spatial_axis_text %||% 12))
  if (is.na(params$spatial_axis_text) || params$spatial_axis_text <= 0) {
    cli::cli_abort("Invalid {.var --spatial_axis_text} value: {params$spatial_axis_text}. Use a positive number.")
  }
  params$spatial_axis_title <- suppressWarnings(as.numeric(params$spatial_axis_title %||% 12))
  if (is.na(params$spatial_axis_title) || params$spatial_axis_title <= 0) {
    cli::cli_abort("Invalid {.var --spatial_axis_title} value: {params$spatial_axis_title}. Use a positive number.")
  }
  params$pca_dims <- suppressWarnings(as.integer(params$pca_dims %||% 10))
  if (is.na(params$pca_dims) || params$pca_dims < 2) {
    cli::cli_abort("Invalid {.var --pca_dims} value: {params$pca_dims}. Use an integer of 2 or greater.")
  }
  params$neighbor_k <- suppressWarnings(as.integer(params$neighbor_k %||% 20))
  if (is.na(params$neighbor_k) || params$neighbor_k < 1) {
    cli::cli_abort("Invalid {.var --neighbor_k} value: {params$neighbor_k}. Use a positive integer.")
  }
  params$cluster_resolution <- suppressWarnings(as.numeric(params$cluster_resolution %||% 0.4))
  if (is.na(params$cluster_resolution) || params$cluster_resolution <= 0) {
    cli::cli_abort("Invalid {.var --cluster_resolution} value: {params$cluster_resolution}. Use a positive number.")
  }

  valid_stages <- c("all", "validate", "ingest", "qc", "analyze", "export")
  if (!params$stage %in% valid_stages) {
    cli::cli_abort("Invalid {.var --stage} value: {params$stage}. Use one of {paste(valid_stages, collapse = ', ')}")
  }

  params$seed <- as.integer(params$seed %||% 1)
  if (is.na(params$seed)) {
    cli::cli_abort("Invalid {.var --seed} value")
  }

  if (!is.null(params$max_cells)) {
    if (is.character(params$max_cells) && tolower(params$max_cells) == "none") {
      params$max_cells <- NA_integer_
    } else {
      params$max_cells <- suppressWarnings(as.integer(params$max_cells))
      if (is.na(params$max_cells) || params$max_cells < 1) {
        params$max_cells <- NA_integer_
      }
    }
  } else {
    params$max_cells <- NA_integer_
  }

  if (!is.null(params$min_genes_per_cell)) {
    params$min_genes_per_cell <- suppressWarnings(as.integer(params$min_genes_per_cell))
    if (is.na(params$min_genes_per_cell) || params$min_genes_per_cell < 1) {
      params$min_genes_per_cell <- NA_integer_
    }
  } else {
    params$min_genes_per_cell <- NA_integer_
  }

  if (!is.null(params$min_total_expr_per_cell)) {
    params$min_total_expr_per_cell <- suppressWarnings(as.integer(params$min_total_expr_per_cell))
    if (is.na(params$min_total_expr_per_cell) || params$min_total_expr_per_cell < 1) {
      params$min_total_expr_per_cell <- NA_integer_
    }
  } else {
    params$min_total_expr_per_cell <- NA_integer_
  }

  if (!is.null(params$max_mito_pct)) {
    params$max_mito_pct <- suppressWarnings(as.numeric(params$max_mito_pct))
    if (is.na(params$max_mito_pct) || params$max_mito_pct < 0) {
      params$max_mito_pct <- NA_real_
    }
  } else {
    params$max_mito_pct <- NA_real_
  }

  raw_prefix <- params$mito_gene_prefixes
  if (is.null(raw_prefix) || (length(raw_prefix) == 1 && is.na(raw_prefix))) {
    prefixes <- c("MT-")
  } else {
    prefixes <- split_prefixes(raw_prefix)
    if (length(prefixes) == 1 && toupper(prefixes) == "NONE") {
      prefixes <- character()
    }
  }
  prefixes <- toupper(prefixes)
  params$mito_gene_prefixes <- prefixes

  if (!is.null(params$python_path) && nzchar(params$python_path) && tolower(params$python_path) != "none") {
    if (!file.exists(params$python_path)) {
      cli::cli_abort("Provided {.var --python_path} does not exist: {.path {params$python_path}}")
    }
    params$python_path <- normalizePath(params$python_path, mustWork = TRUE)
  } else {
    params$python_path <- NULL
  }

  if (!is.null(params$input_dir) && nzchar(params$input_dir)) {
    if (dir.exists(params$input_dir)) {
      params$input_dir <- normalizePath(params$input_dir, mustWork = TRUE)
    } else if (file.exists(params$input_dir)) {
      params$input_dir <- normalizePath(params$input_dir, mustWork = TRUE)
    }
  }

  if (!is.null(params$input_path) && nzchar(params$input_path)) {
    if (!file.exists(params$input_path)) {
      cli::cli_abort("h5ad input file not found: {.path {params$input_path}}")
    }
    params$input_path <- normalizePath(params$input_path, mustWork = TRUE)
  }

  if (!is.null(params$input_object) && nzchar(params$input_object) && tolower(params$input_object) != "none") {
    if (!file.exists(params$input_object)) {
      cli::cli_abort("Input object not found: {.path {params$input_object}}")
    }
    params$input_object <- normalizePath(params$input_object, mustWork = TRUE)
  } else {
    params$input_object <- NULL
  }

  params
}

detect_input_format <- function(params) {
  fmt <- params$input_format
  if (!identical(fmt, "auto")) {
    return(fmt)
  }

  if (!is.null(params$input_path) && nzchar(params$input_path)) {
    ext <- tolower(tools::file_ext(params$input_path))
    if (identical(ext, "h5ad")) {
      return("h5ad")
    }
  }

  in_dir <- params$input_dir
  if (!is.null(in_dir) && dir.exists(in_dir)) {
    in_dir_norm <- normalizePath(in_dir, mustWork = FALSE)
    if (file.exists(file.path(in_dir_norm, "cell_feature_matrix.h5"))) {
      return("xenium")
    }
    subdirs <- list.dirs(in_dir_norm, recursive = FALSE, full.names = TRUE)
    if (any(file.exists(file.path(subdirs, "cell_feature_matrix.h5")))) {
      return("xenium")
    }

    visium_markers <- c(
      file.path(in_dir_norm, "filtered_feature_bc_matrix.h5"),
      file.path(in_dir_norm, "filtered_feature_bc_matrix"),
      file.path(in_dir_norm, "raw_feature_bc_matrix.h5"),
      file.path(in_dir_norm, "raw_feature_bc_matrix")
    )
    spatial_dir <- file.path(in_dir_norm, "spatial")
    if (any(file.exists(visium_markers)) || any(dir.exists(visium_markers))) {
      if (dir.exists(spatial_dir)) {
        return("visium")
      }
    }
    visium_subdirs <- subdirs[vapply(subdirs, function(sd) {
      any(file.exists(file.path(sd, c("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")))) ||
        any(dir.exists(file.path(sd, c("filtered_feature_bc_matrix", "raw_feature_bc_matrix"))))
    }, logical(1))]
    if (length(visium_subdirs) > 0) {
      has_spatial <- vapply(visium_subdirs, function(sd) dir.exists(file.path(sd, "spatial")), logical(1))
      if (any(has_spatial)) {
        return("visium")
      }
    }
  }

  if (!is.null(in_dir) && file.exists(in_dir)) {
    ext <- tolower(tools::file_ext(in_dir))
    if (identical(ext, "h5ad")) {
      return("h5ad")
    }
  }

  fmt
}

# ---- CLI parsing ------------------------------------------------------------

parse_args <- function() {
  # Use optparse/argparse if available; otherwise fall back to a minimal parser.
  # Copilot should prefer optparse for consistency with many R CLIs.
  if (requireNamespace("optparse", quietly = TRUE)) {
    optparse <- asNamespace("optparse")
    cli_tokens <- commandArgs(trailingOnly = TRUE)
    specified_options <- unique(sub("=.*$", "", sub("^--", "", cli_tokens[grepl("^--", cli_tokens)])))

    option_list <- list(
      optparse$make_option(c("--config"), type = "character", default = NULL,
                           help = "Path to YAML/JSON config file. CLI flags override config."),
      optparse$make_option(c("--version"), action = "store_true", default = FALSE,
                           help = "Print pipeline version and exit."),
      optparse$make_option(c("--stage"), type = "character", default = "all",
                           help = "Workflow stage: all|validate|ingest|qc|analyze|export (default: all)."),
      optparse$make_option(c("--input_format"), type = "character", default = "auto",
               help = "Input format: auto|xenium|visium|h5ad|matrix. 'matrix' = explicit expr+spatial(+meta)."),
      optparse$make_option(c("--input_dir"), type = "character", default = NULL,
               help = "Standardized ST input directory (recommended for general support)."),
      optparse$make_option(c("--input_path"), type = "character", default = NULL,
               help = "Path to a single input file (e.g., .h5ad)."),
      optparse$make_option(c("--input_object"), type = "character", default = NULL,
                           help = "Path to an existing Giotto RDS object for qc|analyze|export stages."),
      optparse$make_option(c("--expr"), type = "character", default = NULL,
                           help = "Expression matrix file (genes x spots/cells). Used for input_format=matrix."),
      optparse$make_option(c("--spatial"), type = "character", default = NULL,
                           help = "Spatial coordinates file. Used for input_format=matrix."),
      optparse$make_option(c("--meta"), type = "character", default = NULL,
                           help = "Optional metadata file (spot/cell annotations)."),
      optparse$make_option(c("--output_dir"), type = "character", default = "/output",
                           help = "Base output directory (default: /output)."),
      optparse$make_option(c("--project_id"), type = "character", default = "giotto-st",
                           help = "Short identifier used in output filenames."),
      optparse$make_option(c("--python_path"), type = "character", default = NULL,
               help = "Path to Python binary for Giotto (optional)."),
      optparse$make_option(c("--cores"), type = "integer", default = 4,
               help = "Number of cores to use (default: 4)."),
      optparse$make_option(c("--spatial_point_size"), type = "double", default = 2.25,
                           help = "Dot size for exported spatial map (default: 2.25)."),
      optparse$make_option(c("--umap_point_size"), type = "double", default = 1.5,
                           help = "Dot size for exported UMAP plot (default: 1.5)."),
      optparse$make_option(c("--spatial_legend_text"), type = "double", default = 12,
                           help = "Legend text size for exported spatial map (default: 12)."),
      optparse$make_option(c("--spatial_legend_symbol_size"), type = "double", default = 1.4,
                           help = "Legend symbol size for exported spatial map (default: 1.4)."),
      optparse$make_option(c("--spatial_axis_text"), type = "double", default = 12,
                           help = "Axis tick text size for exported spatial map (default: 12)."),
      optparse$make_option(c("--spatial_axis_title"), type = "double", default = 12,
                           help = "Axis title size for exported spatial map (default: 12)."),
      optparse$make_option(c("--pca_dims"), type = "integer", default = 10,
                           help = "Maximum number of PCA dimensions to use for UMAP and clustering (default: 10)."),
      optparse$make_option(c("--neighbor_k"), type = "integer", default = 20,
                           help = "Nearest-neighbor graph k before clustering (default: 20; capped at cells - 1)."),
      optparse$make_option(c("--cluster_resolution"), type = "double", default = 0.4,
                           help = "Leiden clustering resolution (default: 0.4)."),
      optparse$make_option(c("--seed"), type = "integer", default = 1,
                           help = "Random seed for reproducibility."),
      optparse$make_option(c("--max_cells"), type = "integer", default = NA,
               help = "Optional cap on number of cells/spots to process (random downsample)."),
       optparse$make_option(c("--min_genes_per_cell"), type = "integer", default = NA,
          help = "Optional QC filter: drop cells with detected genes below this threshold."),
       optparse$make_option(c("--min_total_expr_per_cell"), type = "integer", default = NA,
          help = "Optional QC filter: drop cells with total expression below this threshold."),
       optparse$make_option(c("--max_mito_pct"), type = "double", default = NA,
          help = "Optional QC filter: drop cells whose mitochondrial fraction exceeds this percentage."),
        optparse$make_option(c("--mito_gene_prefixes"), type = "character", default = NA,
             help = "Comma-separated gene symbol prefixes to treat as mitochondrial (case-insensitive). Use 'none' to disable."),
      optparse$make_option(c("--dry_run"), action = "store_true", default = FALSE,
                           help = "Validate inputs and config, then exit without running analysis."),
      optparse$make_option(c("--verbose"), action = "store_true", default = FALSE,
                           help = "Enable verbose logging.")
    )

    parser <- optparse$OptionParser(option_list = option_list)
    args <- optparse$parse_args(parser)
    attr(args, "specified_options") <- specified_options

    return(args)
  }

  cli::cli_abort("Package {.pkg optparse} is required for scripts/run_all.R. Install it (local dev) or bake into container later.")
}

# ---- config loading ---------------------------------------------------------

read_config <- function(path) {
  if (is.null(path)) return(list())
  abort_missing("--config", path)

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yaml", "yml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      cli::cli_abort("Config is YAML but {.pkg yaml} is not installed.")
    }
    return(yaml::read_yaml(path))
  }

  if (ext %in% c("json")) {
    return(jsonlite::read_json(path, simplifyVector = TRUE))
  }

  cli::cli_abort("Unsupported config extension: {.val {ext}} (use .yml/.yaml or .json)")
}

merge_config <- function(cfg, args) {
  # CLI overrides config when CLI is non-NULL / non-empty.
  out <- cfg
  specified_options <- attr(args, "specified_options") %||% names(args)
  for (nm in names(args)) {
    if (!nm %in% specified_options) {
      next
    }
    val <- args[[nm]]
    if (is.logical(val)) {
      # for flags, always respect CLI explicit value
      out[[nm]] <- val
    } else if (!is.null(val) && !(is.character(val) && !nzchar(val))) {
      out[[nm]] <- val
    }
  }
  out
}

# ---- ingest (general ST) ----------------------------------------------------

create_giotto_object <- function(params, output_dir) {
  if (!requireNamespace("Giotto", quietly = TRUE)) {
    cli::cli_abort("{.pkg Giotto} is not installed. Install locally (dev) or bake into container later.")
  }

  fmt <- tolower(params$input_format %||% "auto")

  if (fmt == "xenium") {
    abort_missing("--input_dir", params$input_dir)
    return(ingest_xenium(
      input_dir = params$input_dir,
      output_dir = output_dir,
      project_id = params$project_id,
      python_path = params$python_path,
      cores = params$cores
    ))
  }

  if (fmt == "visium") {
    abort_missing("--input_dir", params$input_dir)
    return(ingest_visium_hd(
      input_dir = params$input_dir,
      output_dir = output_dir,
      project_id = params$project_id,
      python_path = params$python_path,
      cores = params$cores
    ))
  }

  if (fmt == "h5ad") {
    input_path <- params$input_path %||% params$input_dir
    abort_missing("--input_path", input_path)
    return(ingest_h5ad(
      input_path = input_path,
      output_dir = output_dir,
      project_id = params$project_id,
      python_path = params$python_path,
      cores = params$cores
    ))
  }

  if (fmt == "matrix") {
    abort_missing("--expr", params$expr)
    abort_missing("--spatial", params$spatial)
    if (!is.null(params$meta)) abort_missing("--meta", params$meta)
    cli::cli_abort("input_format=matrix is not implemented yet.")
  }

  if (fmt == "auto") {
    cli::cli_abort("input_format=auto could not resolve a supported ingest path. Specify {.var --input_format} explicitly.")
  }

  cli::cli_abort("Unsupported input_format: {.val {params$input_format}}")
}

# ---- analysis stages (skeleton) --------------------------------------------

run_pipeline <- function(gobj, params, stats) {
  run_basic_pipeline(
    gobj = gobj,
    stats = stats,
    output_dir = params$output_dir,
    project_id = params$project_id,
    cores = params$cores,
    spatial_point_size = params$spatial_point_size,
    umap_point_size = params$umap_point_size,
    spatial_legend_text = params$spatial_legend_text,
    spatial_legend_symbol_size = params$spatial_legend_symbol_size,
    spatial_axis_text = params$spatial_axis_text,
    spatial_axis_title = params$spatial_axis_title,
    pca_dims = params$pca_dims,
    neighbor_k = params$neighbor_k,
    cluster_resolution = params$cluster_resolution
  )
}

infer_project_id_from_object <- function(path) {
  name <- basename(path)
  name <- sub("\\.rds$", "", name, ignore.case = TRUE)
  name <- sub("_ingested_giotto$", "", name)
  name <- sub("_qc_giotto$", "", name)
  name <- sub("_analyzed_giotto$", "", name)
  name <- sub("_giotto_object$", "", name)
  name
}

resolve_project_id <- function(project_id, input_object = NULL) {
  if (!is.null(project_id) && nzchar(project_id) && !identical(project_id, "giotto-st")) {
    return(project_id)
  }
  if (!is.null(input_object) && nzchar(input_object)) {
    inferred <- infer_project_id_from_object(input_object)
    if (nzchar(inferred)) {
      return(inferred)
    }
  }
  project_id
}

load_stage_object <- function(path) {
  gobj <- readRDS(path)
  if (!inherits(gobj, "giotto")) {
    cli::cli_abort("Input object is not a Giotto object: {.path {path}}")
  }

  stats <- derive_giotto_stats(gobj)
  list(giotto = gobj, stats = stats)
}

write_stage_object <- function(gobj, output_dir, project_id, suffix) {
  object_path <- file.path(output_dir, "objects", paste0(project_id, "_", suffix, "_giotto.rds"))
  saveRDS(gobj, object_path)
  cli::cli_alert_success("Saved stage object: {.path {object_path}}")
  object_path
}

write_filter_summary <- function(output_dir, project_id, qc_filter) {
  if (length(qc_filter$details) == 0) {
    return(NULL)
  }

  qc_rows <- lapply(qc_filter$details, function(detail) {
    data.frame(
      type = detail$type %||% NA_character_,
      threshold = detail$threshold %||% NA_real_,
      total_before = detail$total_before %||% NA_real_,
      total_after = detail$total_after %||% NA_real_,
      removed = detail$removed %||% NA_real_,
      mito_prefixes = if (is.null(detail$prefixes) || length(detail$prefixes) == 0) NA_character_ else paste(detail$prefixes, collapse = ";"),
      mito_genes_detected = detail$mito_genes_detected %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })

  qc_table <- do.call(rbind, qc_rows)
  qc_summary_path <- file.path(output_dir, "metadata", paste0(project_id, "_filter_summary.csv"))
  utils::write.csv(qc_table, qc_summary_path, row.names = FALSE, quote = TRUE)
  qc_summary_path
}

perform_validation <- function(params, run_record) {
  if (identical(params$input_format, "xenium")) {
    layout <- validate_xenium_inputs(params$input_dir, params$project_id)
    params$project_id <- layout$project_id
    run_record$params <- params
    run_record$ingest <- list(
      run_dir = layout$run_dir,
      h5_path = layout$h5_path,
      cells_path = layout$cells_path
    )
    return(list(params = params, run_record = run_record, message = paste("Validated Xenium inputs in", layout$run_dir)))
  }

  if (identical(params$input_format, "visium")) {
    layout <- validate_visium_inputs(params$input_dir, params$project_id)
    params$project_id <- layout$project_id
    run_record$params <- params
    run_record$ingest <- list(
      run_dir = layout$run_dir,
      matrix = layout$matrix,
      spatial = layout$spatial,
      image = layout$image
    )
    return(list(params = params, run_record = run_record, message = paste("Validated Visium inputs in", layout$run_dir)))
  }

  if (identical(params$input_format, "h5ad")) {
    input_path <- params$input_path %||% params$input_dir
    layout <- validate_h5ad_inputs(input_path, params$project_id, params$python_path)
    params$project_id <- layout$project_id
    run_record$params <- params
    run_record$ingest <- list(
      h5ad_path = layout$path,
      n_genes = layout$n_genes,
      n_cells = layout$n_cells
    )
    return(list(params = params, run_record = run_record, message = paste("Validated h5ad inputs at", layout$path)))
  }

  cli::cli_abort("Validation is not implemented for input_format={.val {params$input_format}}")
}

maybe_downsample_giotto <- function(gobj, stats, max_cells) {
  if (is.na(max_cells) || max_cells < 1) {
    return(list(giotto = gobj, stats = stats, downsampled = FALSE, removed = 0L))
  }

  expr_raw <- methods::slot(gobj, "raw_exprs")
  current_cells <- stats$n_cells %||% ncol(expr_raw)
  if (current_cells <= max_cells) {
    return(list(giotto = gobj, stats = stats, downsampled = FALSE, removed = 0L))
  }

  cell_dt <- Giotto::pDataDT(gobj)
  if (!"cell_ID" %in% names(cell_dt)) {
    cli::cli_abort("Giotto object lacks cell_ID metadata; cannot downsample.")
  }

  keep_ids <- sample(cell_dt$cell_ID, max_cells)
  cli::cli_alert_info("Downsampling from {current_cells} to {length(keep_ids)} cells/spots")
  gobj_ds <- Giotto::subsetGiotto(gobj, cell_ids = keep_ids)

  stats$n_cells <- length(keep_ids)

  list(giotto = gobj_ds, stats = stats, downsampled = TRUE, removed = current_cells - length(keep_ids))
}

compute_cell_metrics <- function(expr_raw, mito_prefixes) {
  total_cells <- as.integer(ncol(expr_raw))
  gene_counts <- Matrix::colSums(expr_raw > 0)
  total_expr <- Matrix::colSums(expr_raw)

  mito_prefixes_upper <- unique(toupper(mito_prefixes))
  mito_prefixes_upper <- mito_prefixes_upper[nzchar(mito_prefixes_upper)]
  mito_mask <- rep(FALSE, nrow(expr_raw))
  mito_genes_detected <- 0L
  mito_pct <- rep(NA_real_, total_cells)

  if (length(mito_prefixes_upper) > 0) {
    gene_names_upper <- toupper(rownames(expr_raw) %||% rep("", nrow(expr_raw)))
    for (pref in mito_prefixes_upper) {
      mito_mask <- mito_mask | startsWith(gene_names_upper, pref)
    }
    mito_genes_detected <- sum(mito_mask)
    if (mito_genes_detected > 0) {
      mito_counts <- Matrix::colSums(expr_raw[mito_mask, , drop = FALSE])
      mito_pct <- ifelse(total_expr > 0, (mito_counts / total_expr) * 100, 0)
    } else {
      mito_pct <- rep(0, total_cells)
    }
  }

  list(
    gene_counts = as.numeric(gene_counts),
    total_expr = as.numeric(total_expr),
    mito_pct = as.numeric(mito_pct),
    mito_genes_detected = mito_genes_detected,
    mito_prefixes = mito_prefixes_upper,
    total_cells = total_cells
  )
}

apply_qc_filters <- function(gobj, stats, thresholds) {
  expr_raw <- methods::slot(gobj, "raw_exprs")
  metrics <- compute_cell_metrics(expr_raw, thresholds$mito_prefixes %||% character())
  if (metrics$total_cells == 0) {
    return(list(giotto = gobj, stats = stats, filtered = FALSE, removed = 0L, details = list(), summary = list()))
  }

    keep <- rep(TRUE, metrics$total_cells)
    details <- list()

    record_detail <- function(type, threshold, before, removed, extra = NULL) {
      detail <- list(
        type = type,
        threshold = threshold,
        total_before = as.integer(before),
        removed = as.integer(removed),
        total_after = as.integer(before - removed)
      )
    if (!is.null(extra)) {
      detail <- c(detail, extra)
    }
    details <<- append(details, list(detail))
  }

  if (!is.na(thresholds$min_genes_per_cell) && thresholds$min_genes_per_cell > 0) {
      before <- sum(keep)
    fail <- metrics$gene_counts < thresholds$min_genes_per_cell
      removed <- as.integer(sum(fail & keep))
      record_detail("min_genes_per_cell", thresholds$min_genes_per_cell, before, removed)
      if (removed > 0) {
        keep[fail & keep] <- FALSE
    }
  }

  if (!is.na(thresholds$min_total_expr_per_cell) && thresholds$min_total_expr_per_cell > 0) {
      before <- sum(keep)
    fail <- metrics$total_expr < thresholds$min_total_expr_per_cell
      removed <- as.integer(sum(fail & keep))
      record_detail("min_total_expr_per_cell", thresholds$min_total_expr_per_cell, before, removed)
      if (removed > 0) {
        keep[fail & keep] <- FALSE
    }
  }

  if (!is.na(thresholds$max_mito_pct) && thresholds$max_mito_pct >= 0 && any(!is.na(metrics$mito_pct))) {
      before <- sum(keep)
      fail <- metrics$mito_pct > thresholds$max_mito_pct
      fail_clean <- ifelse(is.na(fail), FALSE, fail)
      removed <- as.integer(sum(fail_clean & keep))
    record_detail(
      "max_mito_pct",
      thresholds$max_mito_pct,
        before,
        removed,
      extra = list(prefixes = metrics$mito_prefixes, mito_genes_detected = metrics$mito_genes_detected)
    )
    if (removed > 0) {
        idx <- which(fail_clean & keep)
        if (length(idx) > 0) {
          keep[idx] <- FALSE
        }
    }
  }

  total_removed <- as.integer(sum(!keep))
  if (total_removed == metrics$total_cells && total_removed > 0) {
    cli::cli_abort("QC filters removed all cells. Loosen thresholds and retry.")
  }

  if (total_removed > 0) {
    keep_ids <- colnames(expr_raw)[keep]
    cli::cli_alert_info(
      "QC filtering removed {total_removed} of {metrics$total_cells} cells. Remaining: {length(keep_ids)}"
    )
    gobj <- Giotto::subsetGiotto(gobj, cell_ids = keep_ids)
  }

  stats$n_cells <- as.integer(metrics$total_cells - total_removed)
  summary <- list(
    total_cells_before = metrics$total_cells,
    total_cells_after = stats$n_cells,
    mito_genes_detected = metrics$mito_genes_detected,
    mito_prefixes = metrics$mito_prefixes
  )
  stats$qc_overview <- summary

  list(
    giotto = gobj,
    stats = stats,
    filtered = total_removed > 0,
    removed = total_removed,
    details = details,
    summary = summary
  )
}

# ---- main ------------------------------------------------------------------

main <- function() {
  load_support_files()

  # Handle --version before full arg parsing
  cli_tokens <- commandArgs(trailingOnly = TRUE)
  if ("--version" %in% cli_tokens) {
    ver <- read_pipeline_version()
    cat(sprintf("giotto-st-pipeline %s\n", ver))
    return(invisible(0))
  }

  args <- parse_args()
  cfg <- read_config(args$config)
  params <- merge_config(cfg, args)
  params <- normalize_params(params)
  params$project_id <- resolve_project_id(params$project_id, params$input_object)
  if (params$stage %in% c("all", "validate", "ingest") || is.null(params$input_object)) {
    params$input_format <- detect_input_format(params)
  }

  if (identical(params$input_format, "xenium") && (is.null(params$input_dir) || !nzchar(params$input_dir))) {
    cli::cli_abort("{.var --input_dir} is required for input_format=xenium")
  }

  if (params$stage %in% c("analyze", "export") && is.null(params$input_object)) {
    cli::cli_abort("{.var --input_object} is required for stage={params$stage}")
  }

  output_dir <- params$output_dir
  ensure_dir(output_dir)

  # Standard output subdirs
  ensure_dir(file.path(output_dir, "metadata"))
  ensure_dir(file.path(output_dir, "objects"))
  ensure_dir(file.path(output_dir, "qc"))
  ensure_dir(file.path(output_dir, "plots"))
  ensure_dir(file.path(output_dir, "tables"))
  ensure_dir(file.path(output_dir, "logs"))

  set.seed(params$seed)

  provenance <- collect_provenance()
  provenance_path <- write_provenance(output_dir, provenance)

  # Record run parameters early so controlled failures still produce metadata.
  run_record <- list(
    tool = "giotto-st-pipeline",
    script = "scripts/run_all.R",
    started_utc = timestamp_utc(),
    params = params,
    stage = params$stage,
    provenance = provenance,
    session = list(
      r_version = R.version.string,
      platform = R.version$platform
    )
  )
  run_params_path <- write_run_parameters(output_dir, run_record)
  cli::cli_alert_info("Wrote provenance metadata: {.path {provenance_path}}")
  cli::cli_alert_info("Wrote run parameters: {.path {run_params_path}}")

  if (isTRUE(params$dry_run) || identical(params$stage, "validate")) {
    validation <- perform_validation(params, run_record)
    run_record <- validation$run_record
    params <- validation$params
    run_record$status <- if (isTRUE(params$dry_run)) "dry_run" else "validated"
    run_record$completed_utc <- timestamp_utc()
    write_run_parameters(output_dir, run_record)
    cli::cli_alert_success(validation$message)
    return(invisible(0))
  }

  ingest <- NULL
  if (is.null(params$input_object) || params$stage %in% c("all", "ingest")) {
    cli::cli_h1("Ingest")
    ingest <- create_giotto_object(params, output_dir)
    params$project_id <- ingest$project_id
    cli::cli_alert_info("Resolved project_id: {params$project_id}")
    run_record$params <- params
    run_record$ingest <- list(
      run_dir = ingest$run_dir,
      n_genes = ingest$stats$n_genes,
      n_cells = ingest$stats$n_cells
    )
    ingest_files <- ingest$files %||% list()
    ingest_files$run_dir <- NULL
    run_record$ingest <- c(run_record$ingest, ingest_files)
    write_run_parameters(output_dir, run_record)
  } else {
    cli::cli_h1("Load input object")
    loaded <- load_stage_object(params$input_object)
    params$project_id <- resolve_project_id(params$project_id, params$input_object)
    run_record$params <- params
    run_record$input_object <- params$input_object
    ingest <- list(
      giotto = loaded$giotto,
      stats = loaded$stats,
      project_id = params$project_id,
      run_dir = dirname(params$input_object),
      files = list(input_object = params$input_object)
    )
    write_run_parameters(output_dir, run_record)
  }

  if (identical(params$stage, "ingest")) {
    ingested_path <- write_stage_object(ingest$giotto, output_dir, params$project_id, "ingested")
    run_record$stage_outputs <- list(ingested_object = ingested_path)
    run_record$status <- "success"
    run_record$completed_utc <- timestamp_utc()
    write_run_parameters(output_dir, run_record)
    return(invisible(0))
  }

  if (params$stage %in% c("all", "qc")) {
    downsample <- maybe_downsample_giotto(ingest$giotto, ingest$stats, params$max_cells)
    ingest$giotto <- downsample$giotto
    ingest$stats <- downsample$stats
    if (isTRUE(downsample$downsampled)) {
      run_record$ingest$n_cells <- ingest$stats$n_cells
      run_record$ingest$downsample <- list(
        max_cells = params$max_cells,
        removed = downsample$removed
      )
      write_run_parameters(output_dir, run_record)
    }

    qc_thresholds <- list(
      min_genes_per_cell = params$min_genes_per_cell,
      min_total_expr_per_cell = params$min_total_expr_per_cell,
      max_mito_pct = params$max_mito_pct,
      mito_prefixes = params$mito_gene_prefixes
    )
    qc_filter <- apply_qc_filters(ingest$giotto, ingest$stats, qc_thresholds)
    ingest$giotto <- qc_filter$giotto
    ingest$stats <- qc_filter$stats
    run_record$ingest$n_cells <- ingest$stats$n_cells
    run_record$ingest$qc_filters <- qc_filter$details
    run_record$ingest$qc_overview <- qc_filter$summary
    qc_summary_path <- write_filter_summary(output_dir, params$project_id, qc_filter)
    if (!is.null(qc_summary_path)) {
      run_record$ingest$qc_filter_summary_path <- qc_summary_path
    }
    write_run_parameters(output_dir, run_record)

    if (identical(params$stage, "qc")) {
      qc_path <- write_stage_object(ingest$giotto, output_dir, params$project_id, "qc")
      run_record$stage_outputs <- list(qc_object = qc_path)
      run_record$status <- "success"
      run_record$completed_utc <- timestamp_utc()
      write_run_parameters(output_dir, run_record)
      return(invisible(0))
    }
  }

  if (params$stage %in% c("all", "analyze")) {
    cli::cli_h1("Run analysis")
    analysis <- run_analysis_pipeline(
      ingest$giotto,
      ingest$stats,
      cores = params$cores,
      pca_dims = params$pca_dims,
      neighbor_k = params$neighbor_k,
      cluster_resolution = params$cluster_resolution
    )
    ingest$giotto <- analysis$giotto
    ingest$stats <- analysis$stats
    run_record$pipeline <- list(cluster_column = analysis$cluster_column)
    write_run_parameters(output_dir, run_record)

    if (identical(params$stage, "analyze")) {
      analyzed_path <- write_stage_object(ingest$giotto, output_dir, params$project_id, "analyzed")
      run_record$stage_outputs <- list(analyzed_object = analyzed_path)
      run_record$status <- "success"
      run_record$completed_utc <- timestamp_utc()
      write_run_parameters(output_dir, run_record)
      return(invisible(0))
    }
  }

  cli::cli_h1("Export results")
  pipeline <- export_pipeline_outputs(
    gobj = ingest$giotto,
    stats = ingest$stats,
    output_dir = output_dir,
    project_id = params$project_id,
    cluster_column = run_record$pipeline$cluster_column %||% NULL,
    spatial_point_size = params$spatial_point_size,
    umap_point_size = params$umap_point_size,
    spatial_legend_text = params$spatial_legend_text,
    spatial_legend_symbol_size = params$spatial_legend_symbol_size,
    spatial_axis_text = params$spatial_axis_text,
    spatial_axis_title = params$spatial_axis_title
  )
  pipeline_outputs <- pipeline$outputs %||% list()

  run_record$pipeline <- c(run_record$pipeline %||% list(), list(outputs = pipeline_outputs))

  run_record$status <- "success"
  run_record$completed_utc <- timestamp_utc()
  write_run_parameters(output_dir, run_record)

  cli::cli_alert_success("Run complete.")
  invisible(0)
}

# Execute with explicit error handling
tryCatch(
  main(),
  error = function(e) {
    # Ensure a clean, informative CLI error
    msg <- conditionMessage(e)
    cli::cli_alert_danger("Run failed: {msg}")
    quit(status = 1)
  }
)
