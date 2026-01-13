suppressPackageStartupMessages({
  library(cli)
  library(data.table)
  library(Giotto)
  library(reticulate)
})

if (!exists("build_dgCMatrix_from_py", mode = "function")) {
  cli::cli_abort("Missing sparse utility: source R/utils_sparse.R before calling ingest_h5ad().")
}

normalize_py_path <- function(python_path) {
  if (is.null(python_path) || !nzchar(python_path)) {
    return(NULL)
  }
  normalizePath(python_path, mustWork = TRUE)
}

use_python_if_supplied <- function(python_path) {
  if (!is.null(python_path) && nzchar(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }
}

load_anndata_modules <- function(python_path = NULL) {
  use_python_if_supplied(python_path)

  if (!reticulate::py_module_available("anndata")) {
    cli::cli_abort(
      "Python module {.pkg anndata} is required to ingest h5ad files. Install it in your Giotto python environment or supply {.var --python_path}."
    )
  }
  if (!reticulate::py_module_available("scipy.sparse")) {
    cli::cli_abort(
      "Python module {.pkg scipy.sparse} is required to ingest h5ad files. Install it in your Giotto python environment or supply {.var --python_path}."
    )
  }

  list(
    anndata = reticulate::import("anndata", delay_load = FALSE),
    sparse = reticulate::import("scipy.sparse", delay_load = FALSE)
  )
}

read_h5ad_safe <- function(path, modules, backed = "r") {
  tryCatch(
    modules$anndata$read_h5ad(path, backed = backed),
    error = function(err) {
      if (!identical(backed, "r")) {
        cli::cli_abort("Failed to read {.path {path}} as h5ad: {err$message}")
      }
      cli::cli_alert_warning(
        "Backed read failed for {.path {path}} ({err$message}). Falling back to in-memory load."
      )
      modules$anndata$read_h5ad(path)
    }
  )
}

candidate_layer_names <- c("counts", "count", "data", "matrix", "raw")

extract_expression_matrix <- function(adata, sparse_module) {
  gene_ids_main <- reticulate::py_to_r(adata$var_names$tolist())
  raw_gene_ids <- NULL
  if (reticulate::py_has_attr(adata, "raw") && !is_py_none(adata$raw)) {
    raw_gene_ids <- reticulate::py_to_r(adata$raw$var_names$tolist())
  }

  try_source <- function(candidate, gene_ids, label) {
    if (is_py_none(candidate)) {
      return(NULL)
    }
    mat <- NULL
    if (inherits(candidate, "Matrix") || is.matrix(candidate) || inherits(candidate, "array")) {
      mat <- Matrix::Matrix(candidate, sparse = TRUE)
    } else {
      mat <- build_dgCMatrix_from_py(candidate, sparse_module)
    }
    if (is.null(mat)) {
      return(NULL)
    }
    dims <- dim(mat)
    if (is.null(dims) || any(dims == 0)) {
      return(NULL)
    }
    if (Matrix::nnzero(mat) == 0) {
      return(NULL)
    }
    gene_vec <- gene_ids
    if (is.null(gene_vec)) {
      gene_vec <- paste0("gene_", seq_len(dims[2]))
    } else if (length(gene_vec) != dims[2]) {
      gene_vec <- make.unique(as.character(gene_vec[seq_len(min(length(gene_vec), dims[2]))]))
      if (length(gene_vec) < dims[2]) {
        gene_vec <- c(gene_vec, paste0("gene_", seq((length(gene_vec) + 1), dims[2])))
      }
    }
    list(
      matrix = Matrix::t(mat),
      gene_ids = make.unique(as.character(gene_vec))[seq_len(dims[2])],
      source = label
    )
  }

  expr <- try_source(adata$X, gene_ids_main, "X")
  if (!is.null(expr)) {
    return(expr)
  }

  layer_items <- list()
  try({
    reticulate::iterate(adata$layers$items(), function(kv) {
      key_chr <- as.character(kv[[1]])
      layer_items[[key_chr]] <<- kv[[2]]
    })
  }, silent = TRUE)

  if (length(layer_items) > 0) {
    layer_keys_chr <- names(layer_items)
    layer_keys_lower <- tolower(layer_keys_chr)
    for (candidate in candidate_layer_names) {
      idx <- which(layer_keys_lower == candidate)
      if (length(idx) == 0) next
      key <- layer_keys_chr[[idx[[1]]]]
      py_layer <- layer_items[[key]]
      expr <- try_source(py_layer, gene_ids_main, sprintf("layers[%s]", key))
      if (!is.null(expr)) {
        return(expr)
      }
    }
    for (key in layer_keys_chr) {
      py_layer <- layer_items[[key]]
      expr <- try_source(py_layer, gene_ids_main, sprintf("layers[%s]", key))
      if (!is.null(expr)) {
        return(expr)
      }
    }
  }

  if (!is.null(raw_gene_ids)) {
    expr <- try_source(adata$raw$X, raw_gene_ids, "raw")
    if (!is.null(expr)) {
      return(expr)
    }
  }

  cli::cli_abort("Unable to extract an expression matrix from supplied h5ad. Checked X, layers, and raw slots.")
}

extract_spatial_from_obs <- function(obs_dt) {
  lower_names <- tolower(names(obs_dt))
  candidate_pairs <- list(
    c("imagecol", "imagerow"),
    c("array_col", "array_row"),
    c("x_centroid", "y_centroid"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("col", "row"),
    c("x", "y")
  )
  for (pair in candidate_pairs) {
    if (all(pair %in% lower_names)) {
      x_col <- names(obs_dt)[match(pair[1], lower_names)]
      y_col <- names(obs_dt)[match(pair[2], lower_names)]
      coords <- obs_dt[, .(sdimx = as.numeric(get(x_col)), sdimy = as.numeric(get(y_col)))]
      if (anyNA(coords$sdimx) || anyNA(coords$sdimy)) {
        next
      }
      return(coords)
    }
  }
  NULL
}

validate_h5ad_inputs <- function(input_path, project_id = NULL, python_path = NULL) {
  if (is.null(input_path) || !nzchar(input_path)) {
    cli::cli_abort("{.var input_path} is required for h5ad ingest.")
  }
  if (!file.exists(input_path)) {
    cli::cli_abort("h5ad file not found: {.path {input_path}}")
  }

  normalized_path <- normalizePath(input_path, mustWork = TRUE)
  modules <- load_anndata_modules(python_path)

  adata <- read_h5ad_safe(normalized_path, modules, backed = "r")
  on.exit({
    try({
      if (reticulate::py_has_attr(adata, "file")) {
        adata$file$close()
      }
    }, silent = TRUE)
  }, add = TRUE)

  shape <- reticulate::py_to_r(adata$shape)
  n_cells <- as.integer(shape[[1]])
  n_genes <- as.integer(shape[[2]])

  obs_df <- reticulate::py_to_r(adata$obs)
  obs_cols <- names(obs_df)
  obsm_keys <- reticulate::py_to_r(adata$obsm$keys())
  obsm_keys_lower <- tolower(obsm_keys)
  lower_obs <- tolower(obs_cols)

  spatial_available <- any(c("imagecol", "array_col", "x_centroid", "pxl_col_in_fullres", "col", "x") %in% lower_obs) &&
    any(c("imagerow", "array_row", "y_centroid", "pxl_row_in_fullres", "row", "y") %in% lower_obs)
  if (!spatial_available) {
    spatial_available <- "spatial" %in% obsm_keys_lower
  }

  if (!spatial_available) {
    cli::cli_abort(
      "No spatial coordinates detected in {.path {input_path}}. Provide columns such as imagecol/imagerow or obsm['spatial']."
    )
  }

  resolved_project_id <- project_id
  if (is.null(resolved_project_id) || !nzchar(resolved_project_id) || resolved_project_id == "giotto-st") {
    resolved_project_id <- tools::file_path_sans_ext(basename(normalized_path))
  }

  list(
    path = normalized_path,
    project_id = resolved_project_id,
    n_cells = n_cells,
    n_genes = n_genes
  )
}

ingest_h5ad <- function(input_path, output_dir, project_id = NULL, python_path = NULL, cores = 4) {
  layout <- validate_h5ad_inputs(input_path, project_id, python_path)

  modules <- load_anndata_modules(python_path)
  adata <- read_h5ad_safe(layout$path, modules, backed = NULL)

  on.exit({
    try({
      if (reticulate::py_has_attr(adata, "file")) {
        adata$file$close()
      }
    }, silent = TRUE)
  }, add = TRUE)

  sparse_module <- modules$sparse
  expr_info <- extract_expression_matrix(adata, sparse_module)
  expr_mtx <- expr_info$matrix

  cell_ids <- reticulate::py_to_r(adata$obs_names$tolist())
  cell_ids <- make.unique(as.character(cell_ids))
  gene_ids <- make.unique(as.character(expr_info$gene_ids))

  dimnames(expr_mtx) <- list(gene_ids, cell_ids)

  obs_df <- reticulate::py_to_r(adata$obs)
  obs_dt <- data.table::as.data.table(obs_df)
  obs_dt[, cell_ID := cell_ids]

  spatial_dt <- extract_spatial_from_obs(obs_dt)
  if (is.null(spatial_dt)) {
    obsm_keys <- reticulate::py_to_r(adata$obsm$keys())
    match_idx <- match("spatial", tolower(obsm_keys))
    if (!is.na(match_idx)) {
      key <- obsm_keys[[match_idx]]
      coords <- reticulate::py_to_r(adata$obsm[[key]])
      if (is.matrix(coords) && ncol(coords) >= 2) {
        spatial_dt <- data.table::data.table(sdimx = as.numeric(coords[, 1]), sdimy = as.numeric(coords[, 2]))
      }
    }
  }

  if (is.null(spatial_dt) || anyNA(spatial_dt$sdimx) || anyNA(spatial_dt$sdimy)) {
    cli::cli_abort("Unable to derive spatial coordinates from {.path {layout$path}}. Ensure imagecol/imagerow columns or obsm['spatial'] exist.")
  }

  spatial_dt[, cell_ID := cell_ids]
  spatial_dt <- spatial_dt[, .(cell_ID, sdimx, sdimy)]

  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }

  instr_args <- list(
    save_dir = plots_dir,
    save_plot = TRUE,
    show_plot = FALSE
  )
  resolved_python_path <- normalize_py_path(python_path)
  if (!is.null(resolved_python_path)) {
    instr_args$python_path <- resolved_python_path
  }
  instructions <- do.call(Giotto::createGiottoInstructions, instr_args)

  gobj <- Giotto::createGiottoObject(
    raw_exprs = expr_mtx,
    spatial_locs = spatial_dt,
    instructions = instructions,
    cell_metadata = obs_dt,
    cores = cores
  )

  list(
    giotto = gobj,
    project_id = layout$project_id,
    run_dir = dirname(layout$path),
    stats = list(
      n_genes = nrow(expr_mtx),
      n_cells = ncol(expr_mtx),
      expression_source = expr_info$source
    ),
    files = list(
      h5ad_path = layout$path,
      expression_source = expr_info$source
    )
  )
}
