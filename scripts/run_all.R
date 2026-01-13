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
  params$output_dir <- params$output_dir %||% "/output"
  params$cores <- as.integer(params$cores %||% 4)
  if (is.na(params$cores) || params$cores < 1) {
    cli::cli_abort("Invalid {.var --cores} value: {params$cores}")
  }

  params$seed <- as.integer(params$seed %||% 1)
  if (is.na(params$seed)) {
    cli::cli_abort("Invalid {.var --seed} value")
  }

  if (!is.null(params$max_cells)) {
    params$max_cells <- suppressWarnings(as.integer(params$max_cells))
    if (is.na(params$max_cells) || params$max_cells < 1) {
      params$max_cells <- NA_integer_
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

  if (!is.null(params$python_path) && nzchar(params$python_path)) {
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

    option_list <- list(
      optparse$make_option(c("--config"), type = "character", default = NULL,
                           help = "Path to YAML/JSON config file. CLI flags override config."),
      optparse$make_option(c("--input_format"), type = "character", default = "auto",
               help = "Input format: auto|xenium|visium|matrix. 'matrix' = explicit expr+spatial(+meta)."),
      optparse$make_option(c("--input_dir"), type = "character", default = NULL,
               help = "Standardized ST input directory (recommended for general support)."),
      optparse$make_option(c("--input_path"), type = "character", default = NULL,
               help = "Path to a single input file (e.g., .h5ad)."),
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
  for (nm in names(args)) {
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
    cores = params$cores
  )
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

  args <- parse_args()
  cfg <- read_config(args$config)
  params <- merge_config(cfg, args)
  params <- normalize_params(params)
  params$input_format <- detect_input_format(params)

  if (identical(params$input_format, "xenium") && (is.null(params$input_dir) || !nzchar(params$input_dir))) {
    cli::cli_abort("{.var --input_dir} is required for input_format=xenium")
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

  # Record run parameters early (mirrors multi-gene-correlations style)
  run_record <- list(
    tool = "giotto-st-pipeline",
    script = "scripts/run_all.R",
    started_utc = timestamp_utc(),
    params = params,
    session = list(
      r_version = R.version.string,
      platform = R.version$platform
    )
  )
  run_params_path <- write_run_parameters(output_dir, run_record)
  cli::cli_alert_info("Wrote run parameters: {.path {run_params_path}}")

  if (isTRUE(params$dry_run)) {
    if (identical(params$input_format, "xenium")) {
      layout <- validate_xenium_inputs(params$input_dir, params$project_id)
      params$project_id <- layout$project_id
      run_record$params <- params
      run_record$ingest <- list(
        run_dir = layout$run_dir,
        h5_path = layout$h5_path,
        cells_path = layout$cells_path
      )
      run_record$status <- "dry_run"
      run_record$completed_utc <- timestamp_utc()
      write_run_parameters(output_dir, run_record)
      cli::cli_alert_success("Dry run: Xenium inputs validated in {.path {layout$run_dir}}")
      return(invisible(0))
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
      run_record$status <- "dry_run"
      run_record$completed_utc <- timestamp_utc()
      write_run_parameters(output_dir, run_record)
      cli::cli_alert_success("Dry run: Visium inputs validated in {.path {layout$run_dir}}")
      return(invisible(0))
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
      run_record$status <- "dry_run"
      run_record$completed_utc <- timestamp_utc()
      write_run_parameters(output_dir, run_record)
      cli::cli_alert_success("Dry run: h5ad inputs validated at {.path {layout$path}}")
      return(invisible(0))
    }

    cli::cli_abort("Dry run is not implemented for input_format={.val {params$input_format}}")
  }

  # Ingest
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
    if (qc_filter$filtered || length(qc_filter$details) > 0) {
      run_record$ingest$n_cells <- ingest$stats$n_cells

      qc_summary_path <- NULL
      if (length(qc_filter$details) > 0) {
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
        qc_summary_path <- file.path(output_dir, "metadata", paste0(params$project_id, "_filter_summary.csv"))
        utils::write.csv(qc_table, qc_summary_path, row.names = FALSE, quote = TRUE)
      }

      run_record$ingest$qc_filters <- qc_filter$details
      run_record$ingest$qc_overview <- qc_filter$summary
      if (!is.null(qc_summary_path)) {
        run_record$ingest$qc_filter_summary_path <- qc_summary_path
      }
      write_run_parameters(output_dir, run_record)
  }

  # Run pipeline
  cli::cli_h1("Run pipeline")
  pipeline <- run_pipeline(ingest$giotto, params, ingest$stats)
  pipeline_outputs <- pipeline$outputs %||% list()

  run_record$pipeline <- list(
    cluster_column = pipeline$cluster_column,
    outputs = pipeline_outputs
  )

  # Save final object
  out_rds <- file.path(output_dir, "objects", paste0(params$project_id, "_giotto_object.rds"))
  saveRDS(pipeline$giotto, out_rds)
  cli::cli_alert_success("Saved Giotto object: {.path {out_rds}}")

  run_record$pipeline$outputs$giotto_object <- out_rds

  session_path <- file.path(output_dir, "metadata", "session_info.txt")
  writeLines(capture.output(sessionInfo()), session_path)
  cli::cli_alert_info("Wrote session info: {.path {session_path}}")

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
