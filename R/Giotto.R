target <- normalizePath("~/R/4.4.3_lib", mustWork = TRUE)
target
target_lib <- target
.libPaths(c(target, .libPaths()))
.libPaths()
python_path <- normalizePath("~/.local/share/r-miniconda/envs/giotto_env/bin/python", mustWork = TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

required_pkgs <- c("Giotto", "data.table", "Matrix", "magrittr", "spatstat.geom", "spatstat.core")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages(lib.loc = .libPaths())))
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(magrittr)
  library(spatstat.geom)
  library(spatstat.core)
  library(Giotto)
})

data_root <- "/data/STAG/data"
results_root <- file.path(getwd(), "giotto_runs")
if (!dir.exists(results_root)) dir.create(results_root, recursive = TRUE)

output_dirs <- list.dirs(data_root, recursive = FALSE, full.names = TRUE)
output_dirs <- output_dirs[grepl("^output-", basename(output_dirs))]
h5_paths <- file.path(output_dirs, "cell_feature_matrix.h5")
h5_paths <- h5_paths[file.exists(h5_paths)]

if (length(h5_paths) == 0) {
  stop("No cell_feature_matrix.h5 files found under ", data_root)
}

run_xenium_analysis <- function(h5_path) {
  base_dir <- dirname(h5_path)
  run_id <- basename(base_dir)
  message("Processing ", run_id)

  cells_file <- file.path(base_dir, "cells.csv.gz")
  if (!file.exists(cells_file)) {
    fallback <- file.path(base_dir, "cells.csv")
    if (file.exists(fallback)) {
      cells_file <- fallback
    } else {
      stop("Missing cells.csv(.gz) for ", run_id)
    }
  }

  spatial_meta <- fread(cells_file)
  required_cols <- c("cell_id", "x_centroid", "y_centroid")
  if (!all(required_cols %in% names(spatial_meta))) {
    stop("cells table for ", run_id, " lacks required columns: ",
         paste(setdiff(required_cols, names(spatial_meta)), collapse = ", "))
  }

  expr_list <- get10Xmatrix_h5(h5_path)
  expr_mtx <- expr_list[["Gene Expression"]]
  if (is.null(expr_mtx)) {
    expr_mtx <- expr_list[[1L]]
    warning("Using first available matrix for ", run_id)
  }

  setnames(spatial_meta, old = required_cols, new = c("cell_ID", "sdimx", "sdimy"))
  common_cells <- intersect(colnames(expr_mtx), spatial_meta$cell_ID)

  if (length(common_cells) == 0) {
    stop("No overlapping cells between expression matrix and metadata for ", run_id)
  }

  expr_mtx <- expr_mtx[, common_cells, drop = FALSE]
  spatial_meta <- spatial_meta[cell_ID %in% common_cells]

  run_dir <- file.path(results_root, run_id)
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)

  instr <- createGiottoInstructions(save_dir = run_dir,
                                    save_plot = TRUE,
                                    show_plot = FALSE,
                                    python_path = python_path)

  gobj <- createGiottoObject(raw_exprs = expr_mtx,
                             spatial_locs = spatial_meta[, .(cell_ID, sdimx, sdimy)],
                             instructions = instr,
                             cores = 4)

  gobj <- normalizeGiotto(gobj, scalefactor = 6000)
  gobj <- addStatistics(gobj)

  cell_meta <- pDataDT(gobj)
  if ("nr_genes" %in% names(cell_meta)) {
    cell_meta[is.na(nr_genes), nr_genes := 0]
  }
  if ("total_expr" %in% names(cell_meta)) {
    cell_meta[is.na(total_expr), total_expr := 0]
  }

  max_dim <- max(2, min(10, nrow(expr_mtx) - 1, length(common_cells) - 1))
  dims_to_use <- seq_len(max_dim)

  gobj <- runPCA(gobj,
                  expression_values = "normalized",
                  genes_to_use = NULL,
                  ncp = max_dim,
                  method = "factominer")
  gobj <- runUMAP(gobj, dimensions_to_use = dims_to_use)
  gobj <- createNearestNetwork(gobj, dimensions_to_use = dims_to_use, k = 20)

  gobj <- tryCatch({
    doLeidenCluster(gobj, resolution = 0.4)
  }, error = function(e) {
    message("Leiden clustering failed for ", run_id, ": ", conditionMessage(e))
    return(gobj)
  })

  cell_meta <- pDataDT(gobj)
  if ("leiden_clus" %in% names(cell_meta)) {
    cluster_column <- "leiden_clus"
  } else if ("cluster" %in% names(cell_meta)) {
    cluster_column <- "cluster"
  } else {
    cluster_column <- "nr_genes"
  }
  if (cluster_column %in% names(cell_meta) && anyNA(cell_meta[[cluster_column]])) {
    if (is.numeric(cell_meta[[cluster_column]])) {
      cell_meta[is.na(get(cluster_column)), (cluster_column) := 0]
    } else {
      cell_meta[is.na(get(cluster_column)), (cluster_column) := "unknown"]
    }
  }

  spatPlot(gobj, show_image = FALSE, cell_color = cluster_column)
  plotUMAP(gobj, cell_color = cluster_column, point_size = 1.5)

  saveRDS(gobj, file = file.path(run_dir, "giotto_object.rds"))
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))

  invisible(TRUE)
}

results <- lapply(h5_paths, function(path) {
  tryCatch(run_xenium_analysis(path),
           error = function(e) {
             message("Failed for ", basename(dirname(path)), ": ", conditionMessage(e))
             return(FALSE)
           })
})

names(results) <- basename(dirname(h5_paths))

completed <- names(which(vapply(results, function(x) isTRUE(x), logical(1))))
failed <- names(which(vapply(results, function(x) identical(x, FALSE), logical(1))))

if (length(completed) > 0) {
  message("Completed runs: ", paste(completed, collapse = ", "))
}
if (length(failed) > 0) {
  message("Failed runs: ", paste(failed, collapse = ", "))
}




