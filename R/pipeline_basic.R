suppressPackageStartupMessages({
  library(cli)
  library(data.table)
  library(Giotto)
})

# Run the core Giotto workflow shared between ST modalities.
run_basic_pipeline <- function(gobj, stats, output_dir, project_id, cores = 4) {
  if (is.null(stats) || is.null(stats$n_genes) || is.null(stats$n_cells)) {
    expr_raw <- methods::slot(gobj, "raw_exprs")
    stats <- list(n_genes = nrow(expr_raw), n_cells = ncol(expr_raw))
  }

  if (stats$n_genes < 2 || stats$n_cells < 2) {
    cli::cli_abort("At least two genes and two cells are required after ingest (found {stats$n_genes} genes, {stats$n_cells} cells)")
  }

  cli::cli_alert_info("Normalizing expression values")
  gobj <- Giotto::normalizeGiotto(gobj, scalefactor = 6000)

  cli::cli_alert_info("Adding expression statistics")
  gobj <- Giotto::addStatistics(gobj)

  cell_meta <- Giotto::pDataDT(gobj)
  if ("nr_genes" %in% names(cell_meta)) {
    cell_meta[is.na(nr_genes), nr_genes := 0]
  }
  if ("total_expr" %in% names(cell_meta)) {
    cell_meta[is.na(total_expr), total_expr := 0]
  }

  ncp <- min(10, stats$n_genes, stats$n_cells)
  ncp <- max(2, ncp)
  dims_to_use <- seq_len(ncp)

  cli::cli_alert_info("Running PCA ({length(dims_to_use)} components)")
  gobj <- Giotto::runPCA(
    gobj,
    expression_values = "normalized",
    ncp = length(dims_to_use),
    method = "factominer"
  )

  cli::cli_alert_info("Running UMAP")
  gobj <- Giotto::runUMAP(gobj, dimensions_to_use = dims_to_use)

  cli::cli_alert_info("Creating nearest-neighbour graph")
  k_neigh <- max(1, min(20, stats$n_cells - 1))
  gobj <- Giotto::createNearestNetwork(gobj, dimensions_to_use = dims_to_use, k = k_neigh)

  cli::cli_alert_info("Running Leiden clustering")
  gobj <- tryCatch(
    Giotto::doLeidenCluster(gobj, resolution = 0.4),
    error = function(e) {
      cli::cli_warn("Leiden clustering failed: {conditionMessage(e)}; continuing without clusters")
      gobj
    }
  )

  cell_meta <- Giotto::pDataDT(gobj)
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

  tables_dir <- file.path(output_dir, "tables")
  plots_dir <- file.path(output_dir, "plots")
  qc_dir <- file.path(output_dir, "qc")
  if (!dir.exists(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Basic QC summaries before downstream embedding
  qc_prefix <- file.path(qc_dir, paste0(project_id, "_"))
  nr_genes <- if ("nr_genes" %in% names(cell_meta)) cell_meta$nr_genes else rep(NA_real_, nrow(cell_meta))
  total_expr <- if ("total_expr" %in% names(cell_meta)) cell_meta$total_expr else rep(NA_real_, nrow(cell_meta))
  qc_metrics <- data.table::data.table(
    metric = c(
      "cells_total",
      "genes_total",
      "median_genes_per_cell",
      "iqr_genes_per_cell_lower",
      "iqr_genes_per_cell_upper",
      "median_total_expr_per_cell",
      "iqr_total_expr_per_cell_lower",
      "iqr_total_expr_per_cell_upper",
      "pct_cells_lt_200_genes",
      "pct_cells_lt_500_total_expr"
    ),
    value = c(
      stats$n_cells,
      stats$n_genes,
      stats::median(nr_genes, na.rm = TRUE),
      stats::quantile(nr_genes, 0.25, na.rm = TRUE),
      stats::quantile(nr_genes, 0.75, na.rm = TRUE),
      stats::median(total_expr, na.rm = TRUE),
      stats::quantile(total_expr, 0.25, na.rm = TRUE),
      stats::quantile(total_expr, 0.75, na.rm = TRUE),
      mean(nr_genes < 200, na.rm = TRUE) * 100,
      mean(total_expr < 500, na.rm = TRUE) * 100
    )
  )

  qc_metrics_path <- paste0(qc_prefix, "qc_metrics.csv")
  data.table::fwrite(qc_metrics, file = qc_metrics_path)

  qc_nr_genes_hist_path <- paste0(qc_prefix, "nr_genes_hist.png")
  grDevices::png(qc_nr_genes_hist_path, width = 1600, height = 1200, res = 200)
  hist(nr_genes, breaks = 50, col = "steelblue", border = "white",
       main = "Detected genes per cell", xlab = "Number of genes")
  grDevices::dev.off()

  qc_total_expr_hist_path <- paste0(qc_prefix, "total_expr_hist.png")
  grDevices::png(qc_total_expr_hist_path, width = 1600, height = 1200, res = 200)
  hist(total_expr, breaks = 50, col = "darkseagreen", border = "white",
       main = "Total expression counts per cell", xlab = "Total expression")
  grDevices::dev.off()

  qc_scatter_path <- paste0(qc_prefix, "genes_vs_expr.png")
  grDevices::png(qc_scatter_path, width = 1600, height = 1200, res = 200)
  graphics::plot(
    total_expr,
    nr_genes,
    pch = 19,
    col = grDevices::adjustcolor("#1f77b4", alpha.f = 0.35),
    cex = 0.6,
    main = "Genes vs total expression",
    xlab = "Total expression",
    ylab = "Genes detected"
  )
  graphics::grid()
  grDevices::dev.off()

  qc_summary_path <- paste0(qc_prefix, "qc_summary.txt")
  genes_q <- stats::quantile(nr_genes, c(0.25, 0.5, 0.75), na.rm = TRUE)
  expr_q <- stats::quantile(total_expr, c(0.25, 0.5, 0.75), na.rm = TRUE)
  low_gene_pct <- qc_metrics[metric == "pct_cells_lt_200_genes", value]
  low_expr_pct <- qc_metrics[metric == "pct_cells_lt_500_total_expr", value]
  qc_notes <- c(
    sprintf("Cells analysed: %d; genes quantified: %d.", stats$n_cells, stats$n_genes),
    sprintf(
      "Genes/cell median %.0f (IQR %.0f-%.0f).",
      genes_q[["50%"]], genes_q[["25%"]], genes_q[["75%"]]
    ),
    sprintf(
      "Total expression/cell median %.0f (IQR %.0f-%.0f).",
      expr_q[["50%"]], expr_q[["25%"]], expr_q[["75%"]]
    ),
    sprintf("%.1f%% of cells fall below 200 genes per cell.", low_gene_pct),
    sprintf("%.1f%% of cells fall below 500 total expression counts.", low_expr_pct)
  )
  if (!is.na(low_gene_pct) && low_gene_pct > 10) {
    qc_notes <- c(qc_notes, "Consider tightening filtering thresholds; low-complexity cells exceed 10%.")
  }
  if (!is.na(low_expr_pct) && low_expr_pct > 10) {
    qc_notes <- c(qc_notes, "Expression depth shows a sizeable low-count tail; review sequencing depth or cell filtering.")
  }
  writeLines(qc_notes, qc_summary_path)

  clusters_path <- file.path(tables_dir, "clusters.csv")
  cli::cli_alert_info("Writing cluster metadata table")
  data.table::fwrite(cell_meta, file = clusters_path)

  save_prefix <- project_id
  spatial_path <- file.path(plots_dir, paste0(save_prefix, "_spatial.png"))
  umap_path <- file.path(plots_dir, paste0(save_prefix, "_umap.png"))

  cli::cli_alert_info("Rendering spatial plot")
  Giotto::spatPlot2D(
    gobj,
    show_image = FALSE,
    show_plot = FALSE,
    cell_color = cluster_column,
    save_plot = TRUE,
    save_param = list(
      save_dir = plots_dir,
      save_name = paste0(save_prefix, "_spatial"),
      save_format = "png"
    )
  )

  cli::cli_alert_info("Rendering UMAP plot")
  Giotto::plotUMAP(
    gobj,
    cell_color = cluster_column,
    point_size = 1.5,
    show_plot = FALSE,
    save_plot = TRUE,
    save_param = list(
      save_dir = plots_dir,
      save_name = paste0(save_prefix, "_umap"),
      save_format = "png"
    )
  )

  list(
    giotto = gobj,
    cluster_column = cluster_column,
    outputs = list(
      clusters_table = clusters_path,
      spatial_plot = spatial_path,
      umap_plot = umap_path,
      qc_metrics = qc_metrics_path,
      qc_summary = qc_summary_path,
      qc_nr_genes_hist = qc_nr_genes_hist_path,
      qc_total_expr_hist = qc_total_expr_hist_path,
      qc_genes_vs_expr = qc_scatter_path
    )
  )
}
