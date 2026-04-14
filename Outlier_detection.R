#' RNA-seq Outlier Detection using Isolation Forest
#'
#' Functions to simulate, normalize, and detect outlier samples in RNA-seq data.
#'
#' @author
#' @date 2025

# Load required packages
library(DESeq2)
library(isotree)
library(ggplot2)
library(dplyr)
library(patchwork)

# ----------------------------------------------------------------------
# 1. Simulate RNA-seq count data with outliers
# ----------------------------------------------------------------------

#' Simulate RNA-seq count matrix with contaminated samples
#'
#' @param n_samples Integer, total number of samples.
#' @param n_genes Integer, number of genes.
#' @param n_outliers Integer, number of outlier samples to introduce.
#' @param dispersion Numeric, negative binomial dispersion (default 0.1).
#' @param outlier_strength Numeric vector length 2, range of log2 fold increase for affected genes.
#' @param frac_affected Numeric, proportion of genes altered in each outlier sample.
#' @return Integer matrix of counts (genes x samples).
#' @export
simulate_rnaseq_counts <- function(n_samples = 100,
                                   n_genes = 2000,
                                   n_outliers = 10,
                                   dispersion = 0.1,
                                   outlier_strength = c(3, 5),
                                   frac_affected = 0.2) {
  set.seed(42) # for reproducibility

  # Baseline log2 expression per gene
  gene_means <- runif(n_genes, 3, 12)
  log_expr <- matrix(rnorm(n_samples * n_genes, mean = gene_means, sd = 1),
                     nrow = n_genes, ncol = n_samples)
  colnames(log_expr) <- paste0("Sample", 1:n_samples)
  rownames(log_expr) <- paste0("Gene", 1:n_genes)

  # Introduce outliers
  outlier_idx <- sample(1:n_samples, n_outliers)
  n_affected <- round(frac_affected * n_genes)

  for (i in outlier_idx) {
    affected_genes <- sample(1:n_genes, n_affected)
    shift <- runif(n_affected, outlier_strength[1], outlier_strength[2])
    log_expr[affected_genes, i] <- log_expr[affected_genes, i] + shift
  }

  # Convert to negative binomial counts
  counts <- matrix(NA, nrow = n_genes, ncol = n_samples)
  for (i in 1:n_genes) {
    for (j in 1:n_samples) {
      mu <- 2^log_expr[i, j]  # inverse log2
      counts[i, j] <- rnbinom(1, mu = mu, size = 1/dispersion)
    }
  }
  colnames(counts) <- colnames(log_expr)
  rownames(counts) <- rownames(log_expr)

  # Attach outlier metadata as attribute
  attr(counts, "outlier_indices") <- outlier_idx
  return(counts)
}

# ----------------------------------------------------------------------
# 2. Normalization (VST)
# ----------------------------------------------------------------------

#' Apply variance stabilizing transformation to count matrix
#'
#' @param counts Integer matrix (genes x samples).
#' @return Matrix of normalized expression (genes x samples).
normalize_vst <- function(counts) {
  coldata <- data.frame(condition = factor(rep("A", ncol(counts))))
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ 1)
  vsd <- vst(dds, blind = TRUE)
  assay(vsd)
}

# ----------------------------------------------------------------------
# 3. Isolation Forest outlier detection
# ----------------------------------------------------------------------

#' Detect outlier samples using Isolation Forest
#'
#' @param norm_expr Matrix of normalized expression (genes x samples).
#' @param ntrees Number of trees in the forest.
#' @param sample_size Subsample size for each tree.
#' @param threshold_quantile Quantile to use as outlier cutoff (default 0.95).
#' @return List with scores, outlier flags, and threshold.
detect_outliers_isoforest <- function(norm_expr,
                                      ntrees = 500,
                                      sample_size = 256,
                                      threshold_quantile = 0.95) {
  # Transpose: samples as rows, genes as columns
  X <- t(norm_expr)

  # Train isolation forest
  forest <- isolation_forest(X,
                             ntrees = ntrees,
                             sample_size = min(sample_size, nrow(X)),
                             ndim = ncol(X))

  # Compute scores (higher = more anomalous)
  scores <- predict(forest, X, type = "score")
  outlier_prob <- 2^(-scores)   # convert to probability

  threshold <- quantile(outlier_prob, threshold_quantile)
  outlier_flags <- outlier_prob > threshold

  list(scores = outlier_prob,
       flags = outlier_flags,
       threshold = threshold,
       forest = forest)
}

# ----------------------------------------------------------------------
# 4. Traditional PCA + Mahalanobis (for comparison)
# ----------------------------------------------------------------------

#' Detect outliers using PCA and Mahalanobis distance
#'
#' @param norm_expr Matrix of normalized expression (genes x samples).
#' @param n_pcs Number of principal components to retain.
#' @param alpha Significance level for chi-squared cutoff (default 0.01).
#' @return List with distances, flags, and cutoff.
detect_outliers_mahalanobis <- function(norm_expr, n_pcs = 10, alpha = 0.01) {
  X <- t(norm_expr)
  pca <- prcomp(X, center = TRUE, scale. = FALSE)
  pc_scores <- pca$x[, 1:min(n_pcs, ncol(pca$x))]

  center <- colMeans(pc_scores)
  cov_mat <- cov(pc_scores)
  mahal_dist <- mahalanobis(pc_scores, center, cov_mat)

  df <- ncol(pc_scores)
  cutoff <- qchisq(1 - alpha, df = df)
  flags <- mahal_dist > cutoff

  list(distance = mahal_dist,
       flags = flags,
       cutoff = cutoff,
       pcs = pc_scores)
}

# ----------------------------------------------------------------------
# 5. Visualisation and reporting
# ----------------------------------------------------------------------

#' Plot comparison of Isolation Forest and PCA-based outlier scores
#'
#' @param iso_result Result from detect_outliers_isoforest.
#' @param mahal_result Result from detect_outliers_mahalanobis.
#' @param true_outliers Optional vector of true outlier indices (for simulation).
#' @return ggplot object (combined with patchwork).
plot_outlier_comparison <- function(iso_result, mahal_result, true_outliers = NULL) {
  n <- length(iso_result$scores)
  plot_df <- data.frame(
    Sample = 1:n,
    IsoScore = iso_result$scores,
    MahalDist = mahal_result$distance,
    TrueOutlier = if (!is.null(true_outliers)) 1:n %in% true_outliers else FALSE
  )

  p1 <- ggplot(plot_df, aes(x = Sample, y = IsoScore, color = TrueOutlier)) +
    geom_point(size = 2) +
    geom_hline(yintercept = iso_result$threshold, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("black", "orange"), labels = c("Normal", "True Outlier")) +
    labs(title = "Isolation Forest Outlier Probabilities",
         y = "Outlier Probability", x = "Sample Index") +
    theme_minimal()

  p2 <- ggplot(plot_df, aes(x = Sample, y = MahalDist, color = TrueOutlier)) +
    geom_point(size = 2) +
    geom_hline(yintercept = mahal_result$cutoff, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("black", "orange"), labels = c("Normal", "True Outlier")) +
    labs(title = "PCA + Mahalanobis Distance",
         y = "Mahalanobis Distance", x = "Sample Index") +
    theme_minimal()

  combined <- p1 / p2 + plot_annotation(title = "Outlier Detection in RNA-seq Samples")
  print(combined)
  invisible(plot_df)
}

# ----------------------------------------------------------------------
# 6. Main workflow wrapper
# ----------------------------------------------------------------------

#' Complete outlier detection pipeline
#'
#' @param counts Integer matrix (genes x samples).
#' @param true_outliers Optional vector of true outlier indices (for simulation).
#' @return List containing both detection results and comparison plot.
run_outlier_pipeline <- function(counts, true_outliers = NULL) {
  cat("Normalizing data with VST...\n")
  norm_expr <- normalize_vst(counts)

  cat("Running Isolation Forest...\n")
  iso_res <- detect_outliers_isoforest(norm_expr)

  cat("Running PCA + Mahalanobis...\n")
  mahal_res <- detect_outliers_mahalanobis(norm_expr)

  cat("Generating comparison plot...\n")
  plot_df <- plot_outlier_comparison(iso_res, mahal_res, true_outliers)

  # Print confusion matrices if true outliers are known
  if (!is.null(true_outliers)) {
    cat("\n--- Isolation Forest Performance ---\n")
    print(table(Predicted = iso_res$flags, Actual = 1:ncol(counts) %in% true_outliers))
    cat("\n--- PCA + Mahalanobis Performance ---\n")
    print(table(Predicted = mahal_res$flags, Actual = 1:ncol(counts) %in% true_outliers))
  }

  list(isolation_forest = iso_res,
       mahalanobis = mahal_res,
       plot_data = plot_df)
}
