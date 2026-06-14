#!/usr/bin/env Rscript
# validate_config.R — Pre-flight config checker with actionable suggestions
#
# Usage:
#   Rscript scripts/validate_config.R configs/my_run.yaml
#   Rscript scripts/validate_config.R --config configs/my_run.yaml
#
# Exits 0 if the config is valid, 1 if problems were found.
# Designed to be fast (no Giotto load, no analysis) and beginner-friendly.

suppressPackageStartupMessages({

  library(cli)
  library(rlang)
})

# ---- helpers ----------------------------------------------------------------

is_not_set <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(TRUE)
  }
  if (length(x) == 1 && is.na(x)) {
    return(TRUE)
  }
  if (is.character(x)) {
    return(!nzchar(trimws(x)) || tolower(trimws(x)) %in% c("none", "null", "na"))
  }
  FALSE
}

read_pipeline_version <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  version_file <- file.path(repo_root, "VERSION")
  if (file.exists(version_file)) {
    return(trimws(readLines(version_file, n = 1, warn = FALSE)))
  }

  "unknown"
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

# ---- config loading ---------------------------------------------------------

load_yaml_config <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("Config file not found: {.path {path}}")
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yaml", "yml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      cli::cli_abort("Config is YAML but {.pkg yaml} is not installed.")
    }
    return(yaml::read_yaml(path))
  }
  if (ext == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      cli::cli_abort("Config is JSON but {.pkg jsonlite} is not installed.")
    }
    return(jsonlite::read_json(path, simplifyVector = TRUE))
  }
  cli::cli_abort("Unsupported config format: {.val {ext}}. Use .yaml, .yml, or .json.")
}

# ---- validators -------------------------------------------------------------

check_required_fields <- function(cfg, issues) {
  if (is.null(cfg$input_format) || !nzchar(cfg$input_format %||% "")) {
    issues <- c(issues, list(list(
      severity = "error",
      field = "input_format",
      message = "Missing required field 'input_format'.",
      suggestion = "Set to one of: visium, xenium, h5ad, matrix, or auto."
    )))
  }

  has_input_dir <- !is_not_set(cfg$input_dir)
  has_input_path <- !is_not_set(cfg$input_path)
  has_input_object <- !is_not_set(cfg$input_object)


  if (!has_input_dir && !has_input_path && !has_input_object) {
    issues <- c(issues, list(list(
      severity = "error",
      field = "input_dir / input_path / input_object",
      message = "No input source specified.",
      suggestion = "Set 'input_dir' (for visium/xenium), 'input_path' (for h5ad), or 'input_object' (for resuming from checkpoint)."
    )))
  }

  if (is.null(cfg$output_dir) || !nzchar(cfg$output_dir %||% "")) {
    issues <- c(issues, list(list(
      severity = "warning",
      field = "output_dir",
      message = "No output_dir specified; will default to '/output'.",
      suggestion = "Set 'output_dir' to a writable path (e.g., 'results/my_sample' for local runs)."
    )))
  }

  if (is.null(cfg$project_id) || !nzchar(cfg$project_id %||% "")) {
    issues <- c(issues, list(list(
      severity = "warning",
      field = "project_id",
      message = "No project_id specified; will default to 'giotto-st'.",
      suggestion = "Set 'project_id' to a short, descriptive name (e.g., 'visium_sample_A')."
    )))
  }

  issues
}

check_input_format <- function(cfg, issues) {
  fmt <- tolower(cfg$input_format %||% "")
  valid_formats <- c("visium", "xenium", "h5ad", "matrix", "auto")

  if (nzchar(fmt) && !fmt %in% valid_formats) {
    issues <- c(issues, list(list(
      severity = "error",
      field = "input_format",
      message = sprintf("Invalid input_format: '%s'.", cfg$input_format),
      suggestion = sprintf("Use one of: %s.", paste(valid_formats, collapse = ", "))
    )))
  }

  issues
}

check_paths_exist <- function(cfg, issues) {
  # input_dir

  if (!is_not_set(cfg$input_dir)) {
    path <- cfg$input_dir
    # Skip container paths (starting with /)
    if (!startsWith(path, "/data") && !startsWith(path, "/output")) {
      if (!dir.exists(path) && !file.exists(path)) {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_dir",
          message = sprintf("input_dir does not exist: '%s'.", path),
          suggestion = "Check that the path is correct and accessible. For container runs, use /data."
        )))
      }
    }
  }

  # input_path
  if (!is_not_set(cfg$input_path)) {
    path <- cfg$input_path
    if (!startsWith(path, "/data") && !startsWith(path, "/output")) {
      if (!file.exists(path)) {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_path",
          message = sprintf("input_path file does not exist: '%s'.", path),
          suggestion = "Provide the full path to your .h5ad file."
        )))
      }
    }
  }

  # input_object
  if (!is_not_set(cfg$input_object)) {
    path <- cfg$input_object
    if (!startsWith(path, "/data") && !startsWith(path, "/output")) {
      if (!file.exists(path)) {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_object",
          message = sprintf("input_object file does not exist: '%s'.", path),
          suggestion = "Provide the path to a previously saved *_giotto.rds checkpoint file."
        )))
      }
    }
  }

  issues
}

check_format_vs_data <- function(cfg, issues) {
  fmt <- tolower(cfg$input_format %||% "auto")
  input_dir <- cfg$input_dir %||% ""

  # Only check if input_dir exists locally (skip container paths)
  if (!nzchar(input_dir) || startsWith(input_dir, "/data") || !dir.exists(input_dir)) {
    return(issues)
  }

  input_dir_norm <- normalizePath(input_dir, mustWork = FALSE)

  if (fmt == "visium") {
    visium_markers <- c(
      file.path(input_dir_norm, "filtered_feature_bc_matrix.h5"),
      file.path(input_dir_norm, "filtered_feature_bc_matrix"),
      file.path(input_dir_norm, "raw_feature_bc_matrix.h5"),
      file.path(input_dir_norm, "raw_feature_bc_matrix")
    )
    spatial_dir <- file.path(input_dir_norm, "spatial")
    has_matrix <- any(file.exists(visium_markers)) || any(dir.exists(visium_markers))
    has_spatial <- dir.exists(spatial_dir)

    if (!has_matrix) {
      # Check if it looks like Xenium instead
      if (file.exists(file.path(input_dir_norm, "cell_feature_matrix.h5"))) {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_format",
          message = "input_format is 'visium' but input_dir looks like a Xenium output.",
          suggestion = "Change input_format to 'xenium'."
        )))
      } else {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_dir",
          message = "input_format is 'visium' but no Spaceranger matrix files found.",
          suggestion = "Point input_dir to the Spaceranger 'outs/' directory (should contain filtered_feature_bc_matrix.h5 and spatial/)."
        )))
      }
    } else if (!has_spatial) {
      issues <- c(issues, list(list(
        severity = "warning",
        field = "input_dir",
        message = "Found expression matrix but no 'spatial/' subdirectory.",
        suggestion = "Ensure the Spaceranger outs/ directory contains a 'spatial/' folder with tissue positions."
      )))
    }
  }

  if (fmt == "xenium") {
    has_cfm <- file.exists(file.path(input_dir_norm, "cell_feature_matrix.h5"))
    # Check subdirs
    if (!has_cfm) {
      subdirs <- list.dirs(input_dir_norm, recursive = FALSE, full.names = TRUE)
      has_cfm <- any(file.exists(file.path(subdirs, "cell_feature_matrix.h5")))
    }
    if (!has_cfm) {
      # Check if it looks like Visium instead
      visium_markers <- c(
        file.path(input_dir_norm, "filtered_feature_bc_matrix.h5"),
        file.path(input_dir_norm, "raw_feature_bc_matrix.h5")
      )
      if (any(file.exists(visium_markers))) {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_format",
          message = "input_format is 'xenium' but input_dir looks like Visium/Spaceranger output.",
          suggestion = "Change input_format to 'visium'."
        )))
      } else {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_dir",
          message = "input_format is 'xenium' but no cell_feature_matrix.h5 found.",
          suggestion = "Point input_dir to the Xenium output directory (or its parent)."
        )))
      }
    }
  }

  if (fmt == "h5ad") {
    # If using input_dir, check for .h5ad files
    input_path <- cfg$input_path %||% ""
    if (!nzchar(input_path) && dir.exists(input_dir_norm)) {
      h5ad_files <- list.files(input_dir_norm, pattern = "\\.h5ad$", full.names = TRUE)
      if (length(h5ad_files) == 0) {
        issues <- c(issues, list(list(
          severity = "error",
          field = "input_dir",
          message = "input_format is 'h5ad' but no .h5ad files found in input_dir.",
          suggestion = "Use 'input_path' to specify the .h5ad file directly, or point input_dir to a directory containing one."
        )))
      }
    }
  }

  issues
}

check_numeric_params <- function(cfg, issues) {
  numeric_checks <- list(
    list(field = "cores", min = 1, type = "integer"),
    list(field = "pca_dims", min = 2, type = "integer"),
    list(field = "neighbor_k", min = 1, type = "integer"),
    list(field = "cluster_resolution", min = 0.001, type = "numeric"),
    list(field = "spatial_point_size", min = 0.01, type = "numeric"),
    list(field = "umap_point_size", min = 0.01, type = "numeric"),
    list(field = "spatial_legend_text", min = 0.01, type = "numeric"),
    list(field = "spatial_legend_symbol_size", min = 0.01, type = "numeric"),
    list(field = "spatial_axis_text", min = 0.01, type = "numeric"),
    list(field = "spatial_axis_title", min = 0.01, type = "numeric"),
    list(field = "max_cells", min = 1, type = "integer"),
    list(field = "min_genes_per_cell", min = 1, type = "integer"),
    list(field = "min_total_expr_per_cell", min = 1, type = "integer"),
    list(field = "max_mito_pct", min = 0, type = "numeric"),
    list(field = "seed", min = NA, type = "integer")
  )

  for (check in numeric_checks) {
    val <- cfg[[check$field]]
    if (is_not_set(val)) next
    num_val <- suppressWarnings(as.numeric(val))
    if (is.na(num_val)) {
      issues <- c(issues, list(list(
        severity = "error",
        field = check$field,
        message = sprintf("'%s' is not a valid number: '%s'.", check$field, as.character(val)),
        suggestion = sprintf("Provide a %s value.", check$type)
      )))
    } else if (!is.na(check$min) && num_val < check$min) {
      issues <- c(issues, list(list(
        severity = "error",
        field = check$field,
        message = sprintf("'%s' value %s is below minimum %s.", check$field, num_val, check$min),
        suggestion = sprintf("Use a value >= %s.", check$min)
      )))
    }
  }

  issues
}

check_stage <- function(cfg, issues) {
  stage <- tolower(cfg$stage %||% "all")
  valid_stages <- c("all", "validate", "ingest", "qc", "analyze", "export")
  if (!stage %in% valid_stages) {
    issues <- c(issues, list(list(
      severity = "error",
      field = "stage",
      message = sprintf("Invalid stage: '%s'.", cfg$stage),
      suggestion = sprintf("Use one of: %s.", paste(valid_stages, collapse = ", "))
    )))
  }

  # Check stage/input_object dependency
  if (stage %in% c("analyze", "export")) {
    if (is_not_set(cfg$input_object)) {
      issues <- c(issues, list(list(
        severity = "error",
        field = "input_object",
        message = sprintf("Stage '%s' requires 'input_object' (a saved Giotto .rds checkpoint).", stage),
        suggestion = "Provide 'input_object' path, or use stage 'all' to run the full pipeline."
      )))
    }
  }

  issues
}

check_mito_prefixes <- function(cfg, issues) {
  val <- cfg$mito_gene_prefixes
  if (is_not_set(val)) return(issues)
  if (is.character(val) && nzchar(val)) {
    pieces <- trimws(unlist(strsplit(val, ",")))
    pieces <- pieces[nzchar(pieces)]
    if (length(pieces) == 0) {
      issues <- c(issues, list(list(
        severity = "warning",
        field = "mito_gene_prefixes",
        message = "mito_gene_prefixes is set but empty after parsing.",
        suggestion = "Use 'MT-' for human, 'mt-' for mouse, or 'none' to disable mitochondrial filtering."
      )))
    }
  }
  issues
}

# ---- reporter ---------------------------------------------------------------

report_issues <- function(issues, config_path) {
  errors <- Filter(function(x) identical(x$severity, "error"), issues)
  warnings <- Filter(function(x) identical(x$severity, "warning"), issues)

  if (length(issues) == 0) {
    cli::cli_alert_success("Config is valid: {.path {config_path}}")
    return(invisible(0L))
  }

  if (length(errors) > 0) {
    cli::cli_h2("Errors ({length(errors)})")
    for (issue in errors) {
      cli::cli_alert_danger("{.field {issue$field}}: {issue$message}")
      cli::cli_alert_info("  {.emph Suggestion}: {issue$suggestion}")
    }
  }

  if (length(warnings) > 0) {
    cli::cli_h2("Warnings ({length(warnings)})")
    for (issue in warnings) {
      cli::cli_alert_warning("{.field {issue$field}}: {issue$message}")
      cli::cli_alert_info("  {.emph Suggestion}: {issue$suggestion}")
    }
  }

  cli::cli_rule()
  if (length(errors) > 0) {
    cli::cli_alert_danger("{length(errors)} error(s) found. Fix these before running the pipeline.")
    return(invisible(1L))
  }
  cli::cli_alert_warning("{length(warnings)} warning(s). Pipeline will run but review the suggestions above.")
  invisible(0L)
}

# ---- main -------------------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Handle --version

  if ("--version" %in% args) {
    ver <- read_pipeline_version()
    cat(sprintf("giotto-st-pipeline validate_config %s\n", ver))
    return(invisible(0L))
  }

  # Parse config path from args
  config_path <- NULL
  if ("--config" %in% args) {
    idx <- which(args == "--config")
    if (idx < length(args)) {
      config_path <- args[idx + 1]
    }
  } else if (length(args) >= 1 && !startsWith(args[1], "--")) {
    config_path <- args[1]
  }

  if (is.null(config_path)) {
    cli::cli_h1("Config Validator")
    cat("Usage:\n")
    cat("  Rscript scripts/validate_config.R configs/my_run.yaml\n")
    cat("  Rscript scripts/validate_config.R --config configs/my_run.yaml\n\n")
    cat("Validates a pipeline config file and reports issues with suggestions.\n")
    return(invisible(0L))
  }

  ver <- read_pipeline_version()
  cli::cli_h1("giotto-st-pipeline config validator (v{ver})")
  cli::cli_alert_info("Checking: {.path {config_path}}")

  cfg <- tryCatch(
    load_yaml_config(config_path),
    error = function(e) {
      cli::cli_alert_danger("Failed to read config: {conditionMessage(e)}")
      quit(status = 1)
    }
  )

  issues <- list()
  issues <- check_required_fields(cfg, issues)
  issues <- check_input_format(cfg, issues)
  issues <- check_paths_exist(cfg, issues)
  issues <- check_format_vs_data(cfg, issues)
  issues <- check_numeric_params(cfg, issues)
  issues <- check_stage(cfg, issues)
  issues <- check_mito_prefixes(cfg, issues)

  exit_code <- report_issues(issues, config_path)
  quit(status = exit_code)
}

main()
