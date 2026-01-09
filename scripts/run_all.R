#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cli)
  library(rlang)
  library(jsonlite)
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

load_support_files <- function() {
  source_helper("R/ingest_xenium.R")
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

  if (!is.null(params$python_path) && nzchar(params$python_path)) {
    if (!file.exists(params$python_path)) {
      cli::cli_abort("Provided {.var --python_path} does not exist: {.path {params$python_path}}")
    }
    params$python_path <- normalizePath(params$python_path, mustWork = TRUE)
  } else {
    params$python_path <- NULL
  }

  if (!is.null(params$input_dir) && nzchar(params$input_dir) && dir.exists(params$input_dir)) {
    params$input_dir <- normalizePath(params$input_dir, mustWork = TRUE)
  }

  params
}

detect_input_format <- function(params) {
  fmt <- params$input_format
  if (!identical(fmt, "auto")) {
    return(fmt)
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

  if (fmt == "matrix") {
    abort_missing("--expr", params$expr)
    abort_missing("--spatial", params$spatial)
    if (!is.null(params$meta)) abort_missing("--meta", params$meta)
    cli::cli_abort("input_format=matrix is not implemented yet.")
  }

  if (fmt == "visium") {
    abort_missing("--input_dir", params$input_dir)
    cli::cli_abort("input_format=visium is not implemented yet.")
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
    h5_path = ingest$files$h5_path,
    cells_path = ingest$files$cells_path,
    n_genes = ingest$stats$n_genes,
    n_cells = ingest$stats$n_cells
  )
  write_run_parameters(output_dir, run_record)

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
