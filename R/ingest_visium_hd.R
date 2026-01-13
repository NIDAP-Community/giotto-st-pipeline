suppressPackageStartupMessages({
  library(cli)
  library(data.table)
  library(Giotto)
})

# Identify the Visium/Visium HD run directory to ingest.
resolve_visium_run_dir <- function(input_dir, project_id = NULL) {
  if (is.null(input_dir) || !nzchar(input_dir)) {
    cli::cli_abort("{.var input_dir} is required for Visium ingest.")
  }
  if (!dir.exists(input_dir)) {
    cli::cli_abort("Visium input directory not found: {.path {input_dir}}")
  }

  normalized_input <- normalizePath(input_dir, mustWork = TRUE)

  has_visium_structure <- function(path) {
    matrix_candidates <- c(
      file.path(path, "filtered_feature_bc_matrix.h5"),
      file.path(path, "raw_feature_bc_matrix.h5")
    )
    matrix_dirs <- c(
      file.path(path, "filtered_feature_bc_matrix"),
      file.path(path, "raw_feature_bc_matrix")
    )
    spatial_dir <- file.path(path, "spatial")
    spatial_zip <- Sys.glob(file.path(path, "*spatial*.zip"))

    matrix_present <- any(file.exists(matrix_candidates)) || any(dir.exists(matrix_dirs))
    spatial_present <- dir.exists(spatial_dir) || length(spatial_zip) > 0
    matrix_present && spatial_present
  }

  if (has_visium_structure(normalized_input)) {
    chosen_dir <- normalized_input
  } else {
    subdirs <- list.dirs(normalized_input, recursive = FALSE, full.names = TRUE)
    subdirs <- subdirs[vapply(subdirs, has_visium_structure, logical(1))]

    if (length(subdirs) == 0) {
      cli::cli_abort(
        "No Visium-like outputs found under {.path {normalized_input}}. Expect Spaceranger/Visium HD outs directory with expression + spatial files."
      )
    }

    if (length(subdirs) == 1) {
      chosen_dir <- subdirs
    } else {
      if (!is.null(project_id) && nzchar(project_id) && project_id != "giotto-st") {
        matches <- subdirs[basename(subdirs) %in% c(project_id, paste0(project_id, "_outs"))]
        if (length(matches) == 1) {
          chosen_dir <- matches
        } else {
          cli::cli_abort(
            "Multiple Visium-like directories detected. Provide a more specific {.var --input_dir} or set {.var --project_id} matching one of: {paste(basename(subdirs), collapse = ', ')}"
          )
        }
      } else {
        cli::cli_abort(
          "Multiple Visium-like directories detected. Provide a more specific {.var --input_dir} or supply {.var --project_id} to disambiguate. Found: {paste(basename(subdirs), collapse = ', ')}"
        )
      }
    }
  }

  resolved_project_id <- project_id
  if (is.null(resolved_project_id) || !nzchar(resolved_project_id) || resolved_project_id == "giotto-st") {
    resolved_project_id <- basename(chosen_dir)
    resolved_project_id <- gsub("[[:space:]]+", "_", resolved_project_id)
  }

  list(run_dir = normalizePath(chosen_dir, mustWork = TRUE), project_id = resolved_project_id)
}

# Locate the spatial metadata file, supporting CSV and Parquet outputs.
locate_visium_spatial_file <- function(spatial_dir) {
  if (!dir.exists(spatial_dir)) {
    cli::cli_abort("Visium spatial directory not found: {.path {spatial_dir}}")
  }

  priority <- c(
    "tissue_positions.parquet",
    "tissue_positions.parquet.gz",
    "tissue_positions_list.parquet",
    "tissue_positions.csv",
    "tissue_positions_list.csv",
    "tissue_positions_list.csv.gz"
  )

  candidates <- unique(c(priority, basename(Sys.glob(file.path(spatial_dir, "tissue_positions*")))))

  for (candidate in candidates) {
    path <- file.path(spatial_dir, candidate)
    if (file.exists(path)) {
      return(path)
    }
  }

  cli::cli_abort("No tissue position file detected under {.path {spatial_dir}}. Expected tissue_positions.parquet or tissue_positions_list.csv.")
}

# Read spatial metadata and harmonize columns for Giotto.
read_visium_spatial <- function(spatial_path) {
  ext <- tolower(tools::file_ext(spatial_path))
  is_parquet <- grepl("\\.parquet(\\.gz)?$", spatial_path, ignore.case = TRUE)

  if (is_parquet) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      cli::cli_abort(
        "Detected {.file {basename(spatial_path)}} but {.pkg arrow} is not installed. Install arrow via renv or provide a CSV tissue positions file."
      )
    }
    dt <- data.table::as.data.table(arrow::read_parquet(spatial_path))
  } else {
    dt <- data.table::fread(spatial_path)
  }

  if (!"barcode" %in% names(dt)) {
    possible <- grep("barcode", names(dt), value = TRUE)
    if (length(possible) == 1) {
      data.table::setnames(dt, possible, "barcode")
    } else {
      cli::cli_abort("Spatial metadata at {.path {spatial_path}} lacks a barcode column.")
    }
  }

  lower_names <- tolower(names(dt))
  data.table::setnames(dt, names(dt), lower_names)

  coord_pairs <- list(
    c("x_centroid", "y_centroid"),
    c("x", "y"),
    c("col_pxl", "row_pxl"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("sdimx", "sdimy")
  )

  coord_cols <- NULL
  for (pair in coord_pairs) {
    if (all(pair %in% lower_names)) {
      coord_cols <- pair
      break
    }
  }

  if (is.null(coord_cols)) {
    cli::cli_abort(
      "Unable to identify spatial coordinate columns in {.path {spatial_path}}. Checked for centroid or pixel columns."
    )
  }

  dt[, cell_id := barcode]
  dt[, sdimx := as.numeric(get(coord_cols[1]))]
  dt[, sdimy := as.numeric(get(coord_cols[2]))]

  if (anyNA(dt$sdimx) || anyNA(dt$sdimy)) {
    cli::cli_abort("Spatial coordinates contain non-numeric values in {.path {spatial_path}}")
  }

  data.table::setkey(dt, cell_id)

  dt
}

validate_visium_inputs <- function(input_dir, project_id = NULL) {
  resolved <- resolve_visium_run_dir(input_dir, project_id)
  run_dir <- resolved$run_dir

  spatial_dir <- file.path(run_dir, "spatial")
  spatial_zip <- NULL
  extracted_dir <- NULL
  if (!dir.exists(spatial_dir)) {
    zip_candidates <- Sys.glob(file.path(run_dir, "*spatial*.zip"))
    if (length(zip_candidates) == 0) {
      cli::cli_abort(
        "Visium spatial directory not found and no *spatial*.zip archive detected under {.path {run_dir}}."
      )
    }
    spatial_zip <- normalizePath(zip_candidates[1], mustWork = TRUE)
    extracted_dir <- tempfile("visium_spatial_")
    dir.create(extracted_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(spatial_zip, exdir = extracted_dir)
    candidate_spatial <- file.path(extracted_dir, "spatial")
    if (dir.exists(candidate_spatial)) {
      spatial_dir <- candidate_spatial
    } else {
      spatial_dir <- extracted_dir
    }
  }

  spatial_file <- locate_visium_spatial_file(spatial_dir)

  matrix_info <- NULL
  matrix_dirs <- c(
    file.path(run_dir, "filtered_feature_bc_matrix"),
    file.path(run_dir, "raw_feature_bc_matrix")
  )
  existing_dirs <- matrix_dirs[dir.exists(matrix_dirs)]
  if (length(existing_dirs) > 0) {
    matrix_info <- list(type = "mtx", path = normalizePath(existing_dirs[1], mustWork = TRUE))
  } else {
    matrix_paths <- c(
      file.path(run_dir, "filtered_feature_bc_matrix.h5"),
      file.path(run_dir, "raw_feature_bc_matrix.h5")
    )
    existing_h5 <- matrix_paths[file.exists(matrix_paths)]
    if (length(existing_h5) > 0) {
      matrix_info <- list(type = "h5", path = normalizePath(existing_h5[1], mustWork = TRUE))
    }
  }

  if (is.null(matrix_info)) {
    cli::cli_abort(
      "No filtered/raw feature matrix found under {.path {run_dir}}. Expected filtered_feature_bc_matrix(.h5)."
    )
  }

  images <- Sys.glob(file.path(spatial_dir, "*aligned*image.*"))
  if (length(images) == 0) {
    images <- Sys.glob(file.path(spatial_dir, "*tissue_image.*"))
  }
  if (length(images) > 0) {
    images <- normalizePath(images[1], mustWork = FALSE)
  } else {
    images <- NULL
  }

  list(
    run_dir = run_dir,
    project_id = resolved$project_id,
    matrix = matrix_info,
    spatial = list(
      path = normalizePath(spatial_file, mustWork = TRUE),
      source = if (!is.null(spatial_zip)) "zip" else "directory",
      zip_path = spatial_zip,
      extracted_dir = if (!is.null(extracted_dir)) normalizePath(spatial_dir, mustWork = TRUE) else spatial_dir
    ),
    image = images
  )
}

#' Ingest Visium or Visium HD outputs into a Giotto object.
ingest_visium_hd <- function(input_dir, output_dir, project_id = NULL, python_path = NULL, cores = 4) {
  layout <- validate_visium_inputs(input_dir, project_id)

  # Expression matrix
  if (identical(layout$matrix$type, "h5")) {
    expr_list <- Giotto::get10Xmatrix_h5(layout$matrix$path)
    expr_mtx <- expr_list[["Gene Expression"]]
    if (is.null(expr_mtx)) {
      expr_mtx <- expr_list[[1L]]
      cli::cli_warn(
        "Using first matrix in {.path {layout$matrix$path}} because 'Gene Expression' slot was not found."
      )
    }
  } else {
    expr_mtx <- Giotto::get10Xmatrix(layout$matrix$path)
  }

  spatial_dt <- read_visium_spatial(layout$spatial$path)

  common_cells <- intersect(colnames(expr_mtx), spatial_dt$cell_id)
  if (length(common_cells) == 0) {
    cli::cli_abort(
      "No overlapping barcodes between expression matrix and spatial metadata in {.path {layout$run_dir}}"
    )
  }

  expr_mtx <- expr_mtx[, common_cells, drop = FALSE]
  spatial_dt <- spatial_dt[common_cells]

  # Prepare instructions (save plots off by default; pipeline will handle plotting)
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

  spatial_locs <- spatial_dt[, .(cell_ID = cell_id, sdimx, sdimy)]
  cell_metadata <- data.table::copy(spatial_dt)
  cell_metadata[, `:=`(sdimx = NULL, sdimy = NULL)]
  data.table::setnames(cell_metadata, "cell_id", "cell_ID")

  giotto_images <- NULL
  if (!is.null(layout$image) && file.exists(layout$image)) {
    if (!requireNamespace("magick", quietly = TRUE)) {
      cli::cli_warn(
        "Visium image detected at {.path {layout$image}} but {.pkg magick} is missing; proceeding without spatial image."
      )
    } else {
      try({
        img <- magick::image_read(layout$image)
        giotto_img <- Giotto::createGiottoImage(
          gobject = NULL,
          spatial_locs = spatial_locs,
          mg_object = img,
          name = "image"
        )
        giotto_images <- list(image = giotto_img)
      }, silent = TRUE)
    }
  }

  gobj <- Giotto::createGiottoObject(
    raw_exprs = expr_mtx,
    spatial_locs = spatial_locs,
    instructions = instructions,
    cell_metadata = cell_metadata,
    images = giotto_images,
    cores = cores
  )

  stats <- list(n_genes = nrow(expr_mtx), n_cells = ncol(expr_mtx))

  list(
    giotto = gobj,
    project_id = layout$project_id,
    run_dir = layout$run_dir,
    stats = stats,
    files = list(
      matrix = layout$matrix,
      spatial = layout$spatial,
      image = layout$image
    )
  )
}
