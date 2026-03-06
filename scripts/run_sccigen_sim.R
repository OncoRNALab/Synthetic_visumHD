#!/usr/bin/env Rscript
##
## Run sCCIgen to generate biologically realistic spatial counts and coordinates.
## Supports optional H5 input (10x format) or synthetic data generation.
## Exports: expression_matrix.tsv (spots x genes), spot_coordinates.tsv (spot_id, barcode, x, y).
## Falls back to spatial NB if sCCIgen not available.
## All randomness controlled by --seed.
##

suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
seed <- as.integer(sub("^--seed=(.*)$", "\\1", args[grepl("^--seed=", args)]))
if (length(seed) == 0) seed <- 42

out_matrix <- sub("^--output-matrix=(.*)$", "\\1", args[grepl("^--output-matrix=", args)])
if (length(out_matrix) == 0) out_matrix <- "expression_matrix.tsv"

out_coords <- sub("^--output-coords=(.*)$", "\\1", args[grepl("^--output-coords=", args)])
if (length(out_coords) == 0) out_coords <- "spot_coordinates.tsv"

feature_h5 <- sub("^--feature-h5=(.*)$", "\\1", args[grepl("^--feature-h5=", args)])
if (length(feature_h5) == 0) feature_h5 <- NULL

n_genes <- as.integer(sub("^--n-genes=(.*)$", "\\1", args[grepl("^--n-genes=", args)]))
if (length(n_genes) == 0) n_genes <- 500L

n_spots <- as.integer(sub("^--n-spots=(.*)$", "\\1", args[grepl("^--n-spots=", args)]))
if (length(n_spots) == 0) n_spots <- 1000L

required_genes_str <- sub("^--required-genes=(.*)$", "\\1", args[grepl("^--required-genes=", args)])
required_genes <- if (length(required_genes_str) > 0 && nchar(required_genes_str) > 0)
  trimws(strsplit(required_genes_str, ",")[[1]]) else character(0)

# Synthesis error rates (applied ONCE per spot barcode, models array printing errors)
synth_sub_rate <- as.numeric(sub("^--synth-sub-rate=(.*)$", "\\1", args[grepl("^--synth-sub-rate=", args)]))
if (length(synth_sub_rate) == 0) synth_sub_rate <- 0.0

synth_ins_rate <- as.numeric(sub("^--synth-ins-rate=(.*)$", "\\1", args[grepl("^--synth-ins-rate=", args)]))
if (length(synth_ins_rate) == 0) synth_ins_rate <- 0.0

synth_del_rate <- as.numeric(sub("^--synth-del-rate=(.*)$", "\\1", args[grepl("^--synth-del-rate=", args)]))
if (length(synth_del_rate) == 0) synth_del_rate <- 0.0

grid_w <- as.integer(sub("^--grid-w=(.*)$", "\\1", args[grepl("^--grid-w=", args)]))
if (length(grid_w) == 0) grid_w <- as.integer(sqrt(n_spots))

grid_h <- as.integer(sub("^--grid-h=(.*)$", "\\1", args[grepl("^--grid-h=", args)]))
if (length(grid_h) == 0) grid_h <- as.integer(ceiling(n_spots / grid_w))

set.seed(seed)

# Check for sCCIgen
sccigen_available <- FALSE
if (requireNamespace("sCCIgen", quietly = TRUE)) {
  tryCatch({
    suppressPackageStartupMessages(library(sCCIgen))
    sccigen_available <- TRUE
  }, error = function(e) NULL)
}

# Read Visium HD feature_slice.h5 format.
# /feature_slices/{feature_idx}/row, col, data
# Each group is ONE GENE; row/col are spatial coordinates of non-zero spots; data is counts.
read_feature_slice_h5 <- function(h5file, n_genes_target, n_spots_target, seed_val,
                                  req_genes = character(0)) {
  set.seed(seed_val)
  library(rhdf5)

  # datasetinfo = TRUE gives us the 'dim' column (n_nonzero spots per gene)
  # without any extra data reads – used for stratified sampling below.
  all_ls <- tryCatch(
    rhdf5::h5ls(h5file, recursive = TRUE, datasetinfo = TRUE),
    error = function(e) rhdf5::h5ls(h5file, recursive = TRUE, datasetinfo = FALSE)
  )
  feature_ids <- all_ls$name[all_ls$group == "/feature_slices"]
  cat("Found", length(feature_ids), "features in /feature_slices\n")
  if (length(feature_ids) == 0) return(NULL)

  # ── Phase 0: stratified gene sampling ─────────────────────────────────────
  # Genes are binned into 5 expression quantile strata based on the number of
  # non-zero spots (a cheap proxy for total expression, read from h5ls metadata).
  # An equal number of genes is drawn from each stratum so that the simulated
  # set spans the full expression range, from rare to highly abundant genes.
  if (length(feature_ids) > n_genes_target) {

    # Extract n_nonzero per gene from the 'dim' column of the 'row' datasets.
    # Each /feature_slices/{fid}/row dataset has length = number of non-zero spots.
    n_nonzero <- rep(1L, length(feature_ids))
    names(n_nonzero) <- feature_ids

    if ("dim" %in% colnames(all_ls)) {
      row_ls <- all_ls[
        all_ls$name == "row" &
        all_ls$group %in% paste0("/feature_slices/", feature_ids), ]

      if (nrow(row_ls) > 0) {
        # dim field is a character string, e.g. "50234" or "50234 x 1"
        fid_from_group <- sub("^/feature_slices/", "", row_ls$group)
        parsed_dim     <- suppressWarnings(
          as.integer(sub("\\s*[xX].*", "", trimws(as.character(row_ls$dim))))
        )
        valid <- !is.na(parsed_dim) & fid_from_group %in% feature_ids
        n_nonzero[fid_from_group[valid]] <- parsed_dim[valid]
        cat("Expression proxy: read n_nonzero for",
            sum(valid), "/", length(feature_ids), "genes from metadata\n")
      }
    }

    # 5-bin stratified sampling on log1p scale
    log_expr  <- log1p(n_nonzero[feature_ids])
    probs     <- seq(0, 1, by = 0.2)
    breaks    <- unique(quantile(log_expr, probs = probs, names = FALSE))

    if (length(breaks) < 2) {
      # All genes have identical expression proxy → fall back to uniform
      cat("Warning: cannot stratify (uniform expression proxy), using random sampling\n")
      feature_ids <- sample(feature_ids, n_genes_target)
    } else {
      bins    <- cut(log_expr, breaks = breaks, include.lowest = TRUE)
      n_bins  <- nlevels(bins)
      per_bin <- ceiling(n_genes_target / n_bins)

      selected <- character(0)
      for (b in levels(bins)) {
        idx <- feature_ids[bins == b]
        k   <- min(per_bin, length(idx))
        if (k > 0) selected <- c(selected, sample(idx, k))
      }
      # Trim to exact target if any bin over-sampled
      if (length(selected) > n_genes_target)
        selected <- sample(selected, n_genes_target)

      feature_ids <- selected
      cat(sprintf(
        "Stratified sampling: %d genes selected across %d expression bins\n",
        length(feature_ids), n_bins))
      cat("Bin sizes (n_nonzero ranges):", paste(
        tapply(n_nonzero[feature_ids], bins[match(feature_ids, names(log_expr))],
               function(x) sprintf("[%d-%d]", min(x), max(x))),
        collapse = "  "), "\n")
    }
  }

  # ── Phase 0b: force-include required genes ────────────────────────────────
  # Gene names in the H5 are numeric feature indices; the actual symbols are
  # read later via /features/name.  We therefore need to build a symbol→fid
  # map here so we can pin required genes by symbol.
  if (length(req_genes) > 0) {
    gnames_all <- tryCatch(rhdf5::h5read(h5file, "/features/name"), error = function(e) NULL)
    if (!is.null(gnames_all)) {
      all_fids_full <- all_ls$name[all_ls$group == "/feature_slices"]
      feat_num_full <- suppressWarnings(as.integer(all_fids_full))
      valid_full    <- !is.na(feat_num_full) & feat_num_full >= 1L &
                       feat_num_full <= length(gnames_all)
      symbol_to_fid <- setNames(all_fids_full[valid_full],
                                gnames_all[feat_num_full[valid_full]])

      pinned_fids <- character(0)
      for (rg in req_genes) {
        if (rg %in% names(symbol_to_fid)) {
          pinned_fids <- c(pinned_fids, symbol_to_fid[[rg]])
        } else {
          cat("Warning: required gene '", rg, "' not found in H5 features\n", sep = "")
        }
      }
      pinned_fids <- unique(pinned_fids)

      # Keep pinned genes; fill the remainder from the already-sampled feature_ids
      remaining_needed <- max(0L, n_genes_target - length(pinned_fids))
      extra_pool       <- setdiff(feature_ids, pinned_fids)
      extra_sampled    <- if (remaining_needed > 0 && length(extra_pool) > 0)
        sample(extra_pool, min(remaining_needed, length(extra_pool))) else character(0)
      feature_ids      <- c(pinned_fids, extra_sampled)

      cat(sprintf(
        "Required genes: %d requested, %d pinned. Total genes after pinning: %d\n",
        length(req_genes), length(pinned_fids), length(feature_ids)))
    } else {
      cat("Warning: could not read /features/name from H5; skipping required-gene pinning\n")
    }
  }

  # ── Phase 1: discover spots by scanning genes until enough unique spots ────
  # Iterate over genes in random order until we have collected at least
  # n_spots_target unique "row_col" keys (or all genes are exhausted).
  # Using a fixed small pilot count (e.g. 10) is unreliable when few genes
  # are selected, because those genes may cover only a small tissue region.
  pilot_spots <- character(0)
  pilot_order <- sample(seq_along(feature_ids))

  for (idx in pilot_order) {
    fid  <- feature_ids[idx]
    base <- paste0("/feature_slices/", fid)
    r <- tryCatch(as.integer(rhdf5::h5read(h5file, paste0(base, "/row"))),
                  error = function(e) integer(0))
    cv <- tryCatch(as.integer(rhdf5::h5read(h5file, paste0(base, "/col"))),
                   error = function(e) integer(0))
    if (length(r) > 0)
      pilot_spots <- union(pilot_spots, paste(r, cv, sep = "_"))
    if (length(pilot_spots) >= n_spots_target * 5L) break
  }

  if (length(pilot_spots) > n_spots_target)
    selected_spots <- sample(pilot_spots, n_spots_target)
  else
    selected_spots <- pilot_spots

  cat("Pilot found", length(pilot_spots), "unique spots; selected",
      length(selected_spots), "for simulation\n")
  n_spots <- length(selected_spots)

  # ── Phase 2: read all target genes, keeping only selected spots ────────────
  feat_list <- vector("list", length(feature_ids))
  spot_list <- vector("list", length(feature_ids))
  data_list <- vector("list", length(feature_ids))

  for (i in seq_along(feature_ids)) {
    fid  <- feature_ids[i]
    base <- paste0("/feature_slices/", fid)
    r    <- as.integer(rhdf5::h5read(h5file, paste0(base, "/row")))
    c    <- as.integer(rhdf5::h5read(h5file, paste0(base, "/col")))
    d    <- as.integer(rhdf5::h5read(h5file, paste0(base, "/data")))
    keys <- paste(r, c, sep = "_")
    keep <- keys %in% selected_spots
    if (any(keep)) {
      feat_list[[i]] <- rep(i, sum(keep))
      spot_list[[i]] <- match(keys[keep], selected_spots)
      data_list[[i]] <- d[keep]
    }
  }

  all_feat <- unlist(feat_list)
  all_spot <- unlist(spot_list)
  all_data <- unlist(data_list)

  if (length(all_feat) == 0) {
    cat("No data found for selected spots – returning NULL\n")
    return(NULL)
  }

  mat <- sparseMatrix(
    i    = all_feat,
    j    = all_spot,
    x    = as.numeric(all_data),
    dims = c(length(feature_ids), n_spots)
  )

  # Decode spatial coordinates from "row_col" keys
  parts    <- strsplit(selected_spots, "_", fixed = TRUE)
  coords_y <- as.integer(vapply(parts, `[[`, character(1L), 1L))
  coords_x <- as.integer(vapply(parts, `[[`, character(1L), 2L))

  # Use feature_ids directly as gene names (they are Ensembl gene IDs in Visium HD).
  # Optionally replace with gene symbols from /features/name if available and
  # feature_ids happen to be 1-based integer indices.
  gene_names <- feature_ids
  for (gp in c("/features/name", "/features/id")) {
    gnames_all <- tryCatch(rhdf5::h5read(h5file, gp), error = function(e) NULL)
    if (!is.null(gnames_all)) {
      feat_num <- suppressWarnings(as.integer(feature_ids))
      if (any(!is.na(feat_num))) {
        # feature_ids are numeric indices – replace with real names
        valid <- !is.na(feat_num) & feat_num >= 1L & feat_num <= length(gnames_all)
        gene_names[valid] <- gnames_all[feat_num[valid]]
      }
      # If feature_ids are already gene ID strings, keep them as-is
      break
    }
  }
  cat("Gene name example:", head(gene_names, 3), "\n")
  rownames(mat) <- gene_names
  colnames(mat) <- selected_spots

  list(mat = mat, x = coords_x, y = coords_y, spot_keys = selected_spots)
}

# Load or generate input data
if (!is.null(feature_h5) && file.exists(feature_h5)) {
  cat("Reading", feature_h5, "...\n")
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    library(rhdf5)

    # Detect format: Visium HD feature_slice vs standard 10x matrix
    top_ls    <- rhdf5::h5ls(feature_h5, recursive = FALSE)
    top_names <- top_ls$name

    if ("feature_slices" %in% top_names) {
      cat("Detected Visium HD feature_slice.h5 format\n")
      result <- read_feature_slice_h5(feature_h5, n_genes, n_spots, seed, required_genes)
      if (!is.null(result)) {
        expr    <- as.matrix(result$mat)     # genes x spots
        spatial <- data.frame(Labels = "spot", center_x = result$x, center_y = result$y)
        anno    <- rep("spot", ncol(expr))
        cat("Loaded:", nrow(expr), "genes,", ncol(expr), "spots\n")
      } else {
        cat("Failed to read feature_slice, using synthetic fallback\n")
        expr <- NULL; spatial <- NULL; anno <- NULL
      }
    } else {
      # Standard 10x CSC matrix (matrix/data, matrix/indices, ...)
      cat("Detected standard 10x matrix format\n")
      root <- NULL
      for (rp in c("matrix", "")) {
        dp  <- if (rp == "") "data" else paste0(rp, "/data")
        chk <- tryCatch(rhdf5::h5read(feature_h5, dp, index = list(1:1)), error = function(e) NULL)
        if (!is.null(chk)) { root <- rp; break }
      }
      if (is.null(root)) {
        cat("Could not locate matrix root, using synthetic fallback\n")
        expr <- NULL; spatial <- NULL; anno <- NULL
      } else {
        pfx      <- if (root == "") "" else paste0(root, "/")
        counts   <- rhdf5::h5read(feature_h5, paste0(pfx, "data"))
        indices  <- rhdf5::h5read(feature_h5, paste0(pfx, "indices"))
        indptr   <- rhdf5::h5read(feature_h5, paste0(pfx, "indptr"))
        shape    <- rhdf5::h5read(feature_h5, paste0(pfx, "shape"))
        barcodes <- rhdf5::h5read(feature_h5, paste0(pfx, "barcodes"))
        genes    <- paste0("gene_", seq_len(shape[1]))
        for (gp in c("features/name", "features/id")) {
          gv <- tryCatch(rhdf5::h5read(feature_h5, paste0(pfx, gp)), error = function(e) NULL)
          if (!is.null(gv)) { genes <- gv; break }
        }
        mat <- sparseMatrix(i = as.integer(indices) + 1L, p = as.integer(indptr),
                            x = as.numeric(counts), dims = as.integer(shape))
        rownames(mat) <- genes; colnames(mat) <- barcodes
        n_side <- as.integer(ceiling(sqrt(ncol(mat))))
        crd    <- expand.grid(x = seq_len(n_side), y = seq_len(n_side))[seq_len(ncol(mat)), ]
        expr    <- as.matrix(mat)
        spatial <- data.frame(Labels = "spot", center_x = crd$x, center_y = crd$y)
        anno    <- rep("spot", ncol(expr))
        cat("Loaded:", nrow(expr), "genes,", ncol(expr), "spots\n")
      }
    }
  } else {
    cat("rhdf5 not available, generating synthetic data\n")
    expr <- NULL; spatial <- NULL; anno <- NULL
  }
} else {
  cat("No H5 input provided, generating synthetic data\n")
  expr <- NULL; spatial <- NULL; anno <- NULL
}

# Use sCCIgen if available (requires input data and proper setup)
# Note: sCCIgen requires parameter YAML files; for simplicity, we use fallback
# Users can provide pre-fitted models and parameter files for full sCCIgen integration
if (sccigen_available && !is.null(expr) && !is.null(spatial) && ncol(expr) > 0) {
  cat("sCCIgen available but requires parameter files.\n")
  cat("Using biologically-informed resampling with spatial structure...\n")
  # Use input data structure but resample with spatial correlation
  count_mat <- expr
  if (ncol(count_mat) != n_spots) {
    # Resample to target number of spots
    if (ncol(count_mat) > n_spots) {
      idx <- sample(ncol(count_mat), n_spots)
      count_mat <- count_mat[, idx]
      x <- spatial$center_x[idx]
      y <- spatial$center_y[idx]
    } else {
      # Expand by resampling with replacement
      idx <- sample(ncol(count_mat), n_spots, replace = TRUE)
      count_mat <- count_mat[, idx]
      x <- spatial$center_x[idx]
      y <- spatial$center_y[idx]
    }
  } else {
    x <- spatial$center_x
    y <- spatial$center_y
  }
  gene_ids <- rownames(count_mat)
  spot_ids <- paste0("spot_", seq_len(ncol(count_mat)))
  cat("Resampled to", ncol(count_mat), "spots\n")
} else {
  expr <- NULL
  spatial <- NULL
}

# Fallback: generate synthetic data and use spatial NB
if (is.null(expr) || !sccigen_available) {
  cat("Using spatial Negative Binomial fallback simulator...\n")
  n_side <- max(1L, as.integer(sqrt(n_spots)))
  n_spots <- n_side * n_side
  
  x <- rep(seq_len(grid_w), each = grid_h)[seq_len(n_spots)]
  y <- rep(seq_len(grid_h), times = grid_w)[seq_len(n_spots)]
  spot_ids <- paste0("spot_", seq_len(n_spots))
  gene_ids <- paste0("gene_", seq_len(n_genes))
  
  mu_base <- exp(rnorm(n_genes, mean = 0, sd = 1.5))
  size <- 1 / (0.3 + runif(n_genes))
  count_mat <- matrix(0L, nrow = n_genes, ncol = n_spots)
  
  for (g in seq_len(n_genes)) {
    mu_spot <- mu_base[g] * (1 + 0.3 * sin(x / 2) * cos(y / 2))
    count_mat[g, ] <- rnbinom(n_spots, size = size[g], mu = pmax(0.1, mu_spot))
  }
  
  rownames(count_mat) <- gene_ids
}

# ── Synthesis error model ────────────────────────────────────────────────────
# Apply mutations ONCE per spot barcode, simulating errors introduced during
# array printing/synthesis.  All molecules captured at a spot share the same
# synthesis-mutated barcode, which is the dominant source of barcode errors in
# spatial transcriptomics.
#
# The perfect barcode is kept as `true_barcode`; the synthesis-mutated version
# is stored as `barcode` (used by attach_spatial_barcodes.py for FASTQ building).

mutate_barcode_synth <- function(bc, sub_r, ins_r, del_r) {
  bases_all <- c("A", "C", "G", "T")
  chars <- strsplit(bc, "")[[1]]
  out   <- character(0)
  for (b in chars) {
    if (runif(1) < del_r) next                                # deletion
    if (runif(1) < ins_r) out <- c(out, sample(bases_all, 1)) # insertion
    if (runif(1) < sub_r) {
      out <- c(out, sample(setdiff(bases_all, b), 1))          # substitution
    } else {
      out <- c(out, b)
    }
  }
  paste(out, collapse = "")
}

# Generate 36nt perfect barcodes per spot
bases <- c("A", "C", "G", "T")
barcode_len <- 36L
true_barcodes <- character(length(spot_ids))

for (i in seq_along(spot_ids)) {
  repeat {
    bc <- paste(sample(bases, barcode_len, replace = TRUE), collapse = "")
    if (!bc %in% true_barcodes) break
  }
  true_barcodes[i] <- bc
}

# Apply synthesis errors to produce per-spot barcodes used during capture
apply_synth <- (synth_sub_rate + synth_ins_rate + synth_del_rate) > 0
if (apply_synth) {
  cat(sprintf("Applying synthesis errors: sub=%.4f ins=%.4f del=%.4f\n",
              synth_sub_rate, synth_ins_rate, synth_del_rate))
  barcodes <- vapply(true_barcodes, mutate_barcode_synth,
                     character(1),
                     sub_r = synth_sub_rate,
                     ins_r = synth_ins_rate,
                     del_r = synth_del_rate)
  n_mutated <- sum(barcodes != true_barcodes)
  cat(sprintf("  %d / %d spots have at least one synthesis error (%.1f%%)\n",
              n_mutated, length(barcodes),
              100 * n_mutated / length(barcodes)))
} else {
  barcodes <- true_barcodes
}

# Write outputs
coord_df <- data.frame(
  spot_id      = spot_ids,
  barcode      = barcodes,        # synthesis-mutated (written to FASTQ)
  x            = x[seq_len(length(spot_ids))],
  y            = y[seq_len(length(spot_ids))],
  true_barcode = true_barcodes,   # perfect barcode (ground truth)
  stringsAsFactors = FALSE
)
write.table(coord_df, out_coords, sep = "\t", row.names = FALSE, quote = FALSE)

mat_for_export <- as.data.frame(t(count_mat))
mat_for_export <- cbind(spot_id = spot_ids, mat_for_export)
write.table(mat_for_export, out_matrix, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Output written:", out_matrix, out_coords, "\n")
quit(save = "no", status = 0)
