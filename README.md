# RNA-seq Outlier Detection using Isolation Forest

An R implementation of **unsupervised anomaly detection** for RNA-seq samples.  
This tool flags potentially mislabeled or contaminated samples using **Isolation Forest**, a robust, high‑dimensional outlier detection algorithm.

## Features

- No dimensionality reduction required (works directly on thousands of genes).
- Automatic outlier scoring (probability scale 0‑1).
- Built‑in simulation of contaminated RNA-seq data.
- Comparison with traditional PCA + Mahalanobis distance.

## Installation

Clone the repository and run the R script. Required packages:

```r
install.packages(c("DESeq2", "isotree", "ggplot2", "dplyr", "patchwork"))
