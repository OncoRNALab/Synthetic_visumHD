#!/usr/bin/env Rscript
##
## Run SRTsim to generate biologically realistic spatial counts and coordinates.
## Exports: expression_matrix.tsv (spots x genes), spot_coordinates.tsv (spot_id, barcode, x, y).
## If SRTsim is not available, falls back to a simple spatial Negative Binomial simulator.
## All randomness controlled by --seed.
##

library(methods)
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(sub("^--seed=(.*)$", "\\1", args[grepl("^--seed=", args)]))
if (length(seed) == 0) seed <- 42
out_matrix <- sub("^--output-matrix=(.*)$", "\\1", args[grepl("^--output-matrix=", args)])
if (length(out_matrix) == 0) out_matrix <- "expression_matrix.tsv"
out_coords <- sub("^--output-coords=(.*)$", "\\1", args[grepl("^--output-coords=", args)])
if (length(out_coords) == 0) out_coords <- "spot_coordinates.tsv"
n_genes <- as.integer(sub("^--n-genes=(.*)$", "\\1", args[grepl("^--n-genes=", args)]))
if (length(n_genes) == 0) n_genes <- 500L
n_spots <- as.integer(sub("^--n-spots=(.*)$", "\\1", args[grepl("^--n-spots=", args)]))
if (length(n_spots) == 0) n_spots <- 1000L

set.seed(seed)

srtsim_available <- FALSE
if (requireNamespace("SRTsim", quietly = TRUE)) {
  tryCatch({
    suppressPackageStartupMessages(library(SRTsim, character.only = TRUE))
    srtsim_available <- TRUE
  }, error = function(e) NULL)
}

if (srtsim_available) {
  ## Reference-based: use toyData, fit, generate counts
  data(toyData, package = "SRTsim")
  toyCount <- toyData$toyCount
  toyInfo  <- toyData$toyInfo
  if (!all(c("x", "y") %in% colnames(toyInfo)))
    toyInfo <- data.frame(x = seq_len(ncol(toyCount)), y = 0, label = "spot")
  srt <- SRTsim::createSRT(count_in = toyCount, loc_in = toyInfo)
  srt <- SRTsim::srtsim_fit(srt, sim_schem = "tissue")
  set.seed(seed)
  srt <- SRTsim::srtsim_count(srt)
  count_mat <- as.matrix(SRTsim::simCounts(srt))
  coldat    <- as.data.frame(SRTsim::simcolData(srt))
  if (!"x" %in% colnames(coldat)) coldat$x <- seq_len(ncol(count_mat))
  if (!"y" %in% colnames(coldat)) coldat$y <- 0
  x <- coldat$x
  y <- coldat$y
  gene_ids <- rownames(count_mat)
  spot_ids <- paste0("spot_", seq_len(ncol(count_mat)))
} else {
  ## Fallback: simple spatial Negative Binomial simulator
  ## Grid of sqrt(n_spots) x sqrt(n_spots); NB counts with simple spatial correlation
  n_side <- max(1L, as.integer(sqrt(n_spots)))
  n_spots <- n_side * n_side
  x <- rep(seq_len(n_side), each = n_side)
  y <- rep(seq_len(n_side), times = n_side)
  spot_ids <- paste0("spot_", seq_len(n_spots))
  gene_ids <- paste0("gene_", seq_len(n_genes))
  mu_base  <- exp(rnorm(n_genes, mean = 0, sd = 1.5))
  size     <- 1 / (0.3 + runif(n_genes))
  count_mat <- matrix(0L, nrow = n_genes, ncol = n_spots)
  for (g in seq_len(n_genes)) {
    mu_spot <- mu_base[g] * (1 + 0.3 * sin(x / 2) * cos(y / 2))
    count_mat[g, ] <- rnbinom(n_spots, size = size[g], mu = pmax(0.1, mu_spot))
  }
  rownames(count_mat) <- gene_ids
}

## Spot coordinates: spot_id, barcode, x, y (36nt barcode per spot)
bases <- c("A", "C", "G", "T")
barcode_len <- 36L
barcodes <- character(n_spots)
for (i in seq_len(n_spots)) {
  repeat {
    bc <- paste(sample(bases, barcode_len, replace = TRUE), collapse = "")
    if (!bc %in% barcodes) break
  }
  barcodes[i] <- bc
}
coord_df <- data.frame(
  spot_id  = spot_ids,
  barcode  = barcodes,
  x        = x,
  y        = y,
  stringsAsFactors = FALSE
)
write.table(coord_df, out_coords, sep = "\t", row.names = FALSE, quote = FALSE)

## Expression matrix: spot_id as column names, gene_id as row names; first column = spot_id
mat_for_export <- as.data.frame(t(count_mat))
mat_for_export <- cbind(spot_id = spot_ids, mat_for_export)
write.table(mat_for_export, out_matrix, sep = "\t", row.names = FALSE, quote = FALSE)

quit(save = "no", status = 0)
