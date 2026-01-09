suppressPackageStartupMessages({
  library(cli)
  library(data.table)
  library(Giotto)
})

# Resolve a Xenium run directory and project identifier.
resolve_xenium_run_dir <- function(input_dir, project_id = NULL) {
  if (is.null(input_dir) || !nzchar(input_dir)) {
    cli::cli_abort("{.var input_dir} is required for Xenium ingest.")
  }

  if (!dir.exists(input_dir)) {
    cli::cli_abort("Xenium input directory not found: {.path {input_dir}}")
  }

  normalized_input <- normalizePath(input_dir, mustWork = TRUE)
  candidate <- file.path(normalized_input, "cell_feature_matrix.h5")

  if (file.exists(candidate)) {
    chosen_dir <- normalized_input
  } else {
    subdirs <- list.dirs(normalized_input, recursive = FALSE, full.names = TRUE)
    subdirs <- subdirs[file.exists(file.path(subdirs, "cell_feature_matrix.h5"))]

    if (length(subdirs) == 0) {
      cli::cli_abort("No Xenium output directories with cell_feature_matrix.h5 found under {.path {normalized_input}}")
    }

    if (length(subdirs) == 1) {
      chosen_dir <- subdirs
    } else {
      if (!is.null(project_id) && nzchar(project_id) && project_id != "giotto-st") {
        matches <- subdirs[basename(subdirs) %in% c(project_id, paste0("output-", project_id))]
        if (length(matches) == 1) {
          chosen_dir <- matches
        } else {
          cli::cli_abort(
            "Multiple Xenium runs detected under {.path {normalized_input}}. Set {.var --input_dir} to a specific output-* directory or pass {.var --project_id} matching one of: {paste(basename(subdirs), collapse = ', ')}"
          )
        }
      } else {
        cli::cli_abort(
          "Multiple Xenium runs detected under {.path {normalized_input}}. Set {.var --input_dir} to a specific output-* directory or supply {.var --project_id} to disambiguate. Found: {paste(basename(subdirs), collapse = ', ')}"
        )
      }
    }
  }

  resolved_project_id <- project_id
  if (is.null(resolved_project_id) || !nzchar(resolved_project_id) || resolved_project_id == "giotto-st") {
    resolved_project_id <- sub("^output-", "", basename(chosen_dir))
  }

  list(run_dir = chosen_dir, project_id = resolved_project_id)
}

validate_xenium_inputs <- function(input_dir, project_id = NULL) {
  resolved <- resolve_xenium_run_dir(input_dir, project_id)

  h5_path <- file.path(resolved$run_dir, "cell_feature_matrix.h5")
  if (!file.exists(h5_path)) {
    cli::cli_abort("Missing cell_feature_matrix.h5 in {.path {resolved$run_dir}}")
  }

  cells_candidates <- file.path(resolved$run_dir, c("cells.csv.gz", "cells.csv"))
  cells_candidates <- cells_candidates[file.exists(cells_candidates)]
  if (length(cells_candidates) == 0) {
    cli::cli_abort("Missing cells.csv(.gz) in {.path {resolved$run_dir}}")
  }

  list(
    run_dir = resolved$run_dir,
    project_id = resolved$project_id,
    h5_path = h5_path,
    cells_path = cells_candidates[1]
  )
}

ingest_xenium <- function(input_dir, output_dir, project_id = NULL, python_path = NULL, cores = 4) {
  layout <- validate_xenium_inputs(input_dir, project_id)

  spatial_meta <- data.table::fread(layout$cells_path)
  required_cols <- c("cell_id", "x_centroid", "y_centroid")
  if (!all(required_cols %in% names(spatial_meta))) {
    missing <- setdiff(required_cols, names(spatial_meta))
    cli::cli_abort("cells table in {.path {layout$run_dir}} lacks required columns: {paste(missing, collapse = ', ')}")
  }

  expr_list <- Giotto::get10Xmatrix_h5(layout$h5_path)
  expr_mtx <- expr_list[["Gene Expression"]]
  if (is.null(expr_mtx)) {
    expr_mtx <- expr_list[[1L]]
    cli::cli_warn("Using first matrix in {.path {layout$h5_path}} because 'Gene Expression' slot was not found")
  }

  setnames(spatial_meta, old = required_cols, new = c("cell_ID", "sdimx", "sdimy"))

  spatial_meta[, sdimx := as.numeric(sdimx)]
  spatial_meta[, sdimy := as.numeric(sdimy)]
  if (anyNA(spatial_meta$sdimx) || anyNA(spatial_meta$sdimy)) {
    cli::cli_abort("Non-numeric spatial coordinates detected in {.path {layout$run_dir}}")
  }

  common_cells <- colnames(expr_mtx)[colnames(expr_mtx) %in% spatial_meta$cell_ID]
  if (length(common_cells) == 0) {
    cli::cli_abort("No overlapping cell IDs between expression matrix and metadata in {.path {layout$run_dir}}")
  }

  expr_mtx <- expr_mtx[, common_cells, drop = FALSE]
  spatial_meta <- spatial_meta[match(common_cells, spatial_meta$cell_ID)]

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }

  instr_args <- list(
    save_dir = plots_dir,
    save_plot = TRUE,
    show_plot = FALSE
  )
  if (!is.null(python_path) && nzchar(python_path)) {
    instr_args$python_path <- python_path
  }

  instructions <- do.call(Giotto::createGiottoInstructions, instr_args)

  gobj <- Giotto::createGiottoObject(
    raw_exprs = expr_mtx,
    spatial_locs = spatial_meta[, .(cell_ID, sdimx, sdimy)],
    instructions = instructions,
    cores = cores
  )

  list(
    giotto = gobj,
    project_id = layout$project_id,
    run_dir = layout$run_dir,
    stats = list(n_genes = nrow(expr_mtx), n_cells = length(common_cells)),
    files = layout
  )
}
