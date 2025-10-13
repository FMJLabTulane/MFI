# KANN-MFI: Reproducible Mitochondrial Fitness Prediction

A production-ready implementation of Kolmogorov-Arnold Networks (KAN) for predicting mitochondrial fitness indices from single-cell RNA-seq data, with complete reproducibility guarantees.

## Features

✅ **Complete Reproducibility**
- SHA-256 checksums of input data
- Full session info capture (R version, packages, torch version)
- Seed management for deterministic results
- Training provenance tracking

✅ **Data Validation**
- Automatic verification that input data hasn't changed
- Feature availability checking
- Cell overlap validation
- Warnings for missing genes or data modifications

✅ **Easy Model Persistence**
- Save/load complete model bundles
- All metadata preserved
- Human-readable provenance files (JSON)
- Automatic README generation

✅ **Professional Workflow**
- Progress bars and informative messages
- Diagnostic plots (predicted vs. observed, residuals)
- Permutation importance with batching
- Print methods for easy inspection

## Installation

```r
# Install required packages
install.packages(c("Seurat", "dplyr", "digest", "jsonlite", "ggplot2"))

# Install torch (if not already installed)
install.packages("torch")
torch::install_torch()

# Source the script
source("kann_from_seurat.R")
```

## Quick Start

```r
library(Seurat)
source("kann_from_seurat.R")

# 1) Train model on your Seurat object
fit <- train_kann_on_seurat(
  seu = combined,
  target_col = "MQI_v1",    # Your phenotype column
  assay = "RNA",
  slot = "data",             # log-normalized counts
  hidden = 48,
  epochs = 300,
  seed = 42
)

# 2) Predict and write back to Seurat
combined <- predict_kann_to_seurat(
  combined, 
  fit, 
  out_col = "KANN_MFI_pred",
  validate = TRUE  # checks data integrity
)

# 3) Save complete model bundle (everything needed to reproduce)
save_kann_fit(fit, dir = "models/kann_mfi_v1")

# 4) Evaluate performance
idx <- which(is.finite(combined$MQI_v1) & is.finite(combined$KANN_MFI_pred))
r2 <- cor(combined$MQI_v1[idx], combined$KANN_MFI_pred[idx])^2
rmse <- sqrt(mean((combined$MQI_v1[idx] - combined$KANN_MFI_pred[idx])^2))
cat(sprintf("R² = %.3f, RMSE = %.4f (n = %d cells)\n", r2, rmse, length(idx)))

# 5) Permutation importance
imp <- perm_importance_seurat(combined, fit, repeats = 3)
head(imp, 20)
```

## Reproducibility Workflow

### First Run: Train and Save

```r
# Set seed and train
fit <- train_kann_on_seurat(
  seu = seu_object,
  target_col = "MQI_v1",
  seed = 42,           # Reproducible initialization
  verbose = TRUE
)

# Save complete bundle
save_kann_fit(fit, dir = "kann_model_20250115")

# This saves:
# - model_state.pt      (PyTorch weights)
# - scaler.rds          (standardization params)
# - config.rds          (hyperparameters)
# - provenance.json     (complete metadata)
# - genes_used.csv      (feature list)
# - train_history.csv   (loss over epochs)
# - README.md           (auto-generated documentation)
```

### Later: Load and Reproduce

```r
# Load model (reconstructs everything exactly)
fit <- load_kann_fit("kann_model_20250115")

# Validate against current Seurat object
validate_kann_fit(seu_object, fit, verbose = TRUE)

# Predict (with automatic validation)
seu_object <- predict_kann_to_seurat(seu_object, fit, validate = TRUE)
```

### What Gets Tracked?

The provenance record includes:
- **Data fingerprint**: SHA-256 hash of exact expression matrix used
- **Feature info**: Which genes, their order, summary statistics
- **Target stats**: Mean, SD, range, n_valid samples
- **Model config**: All hyperparameters (hidden units, bins, epochs, LR)
- **Seeds**: R seed and torch seed for perfect reproduction
- **Session info**: R version, package versions, torch version, platform
- **Timestamps**: When trained, when loaded

## Advanced Usage

### Use Specific Gene Sets

```r
# Define your regulators/markers
mito_genes <- c("MT-ND1", "MT-CO1", "MT-ATP6", "PPARGC1A", "TFAM", ...)

fit <- train_kann_on_seurat(
  seu = seu_object,
  target_col = "MQI_v1",
  features = mito_genes,  # Use only these genes
  min_features = 10
)
```

### Train/Test Split

```r
# Split cells
set.seed(42)
train_cells <- sample(colnames(seu), size = floor(0.8 * ncol(seu)))
test_cells <- setdiff(colnames(seu), train_cells)

# Train on subset
fit <- train_kann_on_seurat(
  seu = seu_object,
  target_col = "MQI_v1",
  cells_filter = train_cells  # Only use these cells
)

# Predict on test set
seu_test <- subset(seu_object, cells = test_cells)
seu_test <- predict_kann_to_seurat(seu_test, fit)

# Evaluate on held-out data
test_idx <- which(is.finite(seu_test$MQI_v1) & is.finite(seu_test$KANN_pred))
test_r2 <- cor(seu_test$MQI_v1[test_idx], seu_test$KANN_pred[test_idx])^2
```

### Batch Training for Large Datasets

```r
fit <- train_kann_on_seurat(
  seu = large_seu,
  target_col = "MQI_v1",
  batch_size = 512,    # Mini-batch training
  epochs = 500,
  lr = 1e-3
)
```

### Diagnostic Plots

```r
plots <- plot_kann_diagnostics(seu_object, fit, out_col = "KANN_MFI_pred")

# View plots
plots$pred_vs_obs      # Predicted vs. Observed scatter
plots$residuals        # Residual plot
plots$residual_hist    # Residual distribution

# Access underlying data
head(plots$stats)
```

### Permutation Importance

```r
# Compute ΔMSE for each feature
imp <- perm_importance_seurat(
  seu_object, 
  fit, 
  repeats = 5,          # More repeats = more stable
  batch_rows = 16384,   # Adjust for memory
  verbose = TRUE
)

# Top 20 most important genes
head(imp, 20)

# Save importance
write.csv(
  data.frame(gene = names(imp), delta_mse = as.numeric(imp)),
  "importance_results.csv",
  row.names = FALSE
)
```

## Validation and QC

### Validate Before Predicting

```r
# Check if your Seurat object is compatible with the saved model
validation <- validate_kann_fit(seu_object, fit, verbose = TRUE)

# Programmatic checks
if (!validation$valid) {
  warning("Validation failed! Check results:")
  print(validation$checks)
}
```

### Inspect Fit Object

```r
# Print summary
print(fit)

# Access components
fit$model              # PyTorch model
fit$scaler             # Standardization parameters
fit$used_genes         # Genes in training
fit$used_cells         # Cells in training
fit$provenance         # Complete metadata
fit$train_history      # Loss over epochs

# Check data hash
fit$provenance$data_source$data_hash
```

## Example: Complete Analysis Pipeline

```r
# ============================================================================
# Complete Reproducible KANN-MFI Analysis
# ============================================================================

library(Seurat)
library(ggplot2)
source("kann_from_seurat.R")

# Load data
combined <- readRDS("seurat_object.rds")
DefaultAssay(combined) <- "RNA"

# Check target
table(is.finite(combined$MQI_v1))

# ---- 1. Train model --------------------------------------------------------
fit <- train_kann_on_seurat(
  seu = combined,
  target_col = "MQI_v1",
  assay = "RNA",
  slot = "data",
  hidden = 48,
  n_bins = 10,
  epochs = 300,
  lr = 2e-3,
  seed = 42,
  verbose = TRUE
)

# ---- 2. Predict ------------------------------------------------------------
combined <- predict_kann_to_seurat(
  combined, 
  fit, 
  out_col = "KANN_MFI",
  resid_col = "KANN_resid",
  validate = TRUE
)

# ---- 3. Evaluate -----------------------------------------------------------
idx <- which(is.finite(combined$MQI_v1) & is.finite(combined$KANN_MFI))
obs <- combined$MQI_v1[idx]
pred <- combined$KANN_MFI[idx]

r2 <- cor(obs, pred)^2
rmse <- sqrt(mean((obs - pred)^2))
mae <- mean(abs(obs - pred))

cat(sprintf("\n=== Performance Metrics ===\n"))
cat(sprintf("n = %d cells\n", length(idx)))
cat(sprintf("R² = %.4f\n", r2))
cat(sprintf("RMSE = %.4f\n", rmse))
cat(sprintf("MAE = %.4f\n", mae))

# ---- 4. Diagnostics --------------------------------------------------------
plots <- plot_kann_diagnostics(combined, fit, out_col = "KANN_MFI")
ggsave("diagnostics_pred_vs_obs.pdf", plots$pred_vs_obs, width = 6, height = 5)
ggsave("diagnostics_residuals.pdf", plots$residuals, width = 6, height = 5)

# ---- 5. Feature importance -------------------------------------------------
imp <- perm_importance_seurat(combined, fit, repeats = 5, verbose = TRUE)
write.csv(
  data.frame(gene = names(imp), delta_mse = as.numeric(imp)),
  "importance_delta_mse.csv",
  row.names = FALSE
)

# Plot top 20
imp_df <- data.frame(
  gene = names(imp)[1:20],
  importance = as.numeric(imp)[1:20]
)
imp_df$gene <- factor(imp_df$gene, levels = rev(imp_df$gene))

p_imp <- ggplot(imp_df, aes(x = importance, y = gene)) +
  geom_col(fill = "steelblue") +
  labs(title = "Top 20 Features by ΔMSE",
       x = "Importance (ΔMSE)", y = NULL) +
  theme_minimal()
ggsave("importance_top20.pdf", p_imp, width = 7, height = 6)

# ---- 6. Save everything ----------------------------------------------------
save_kann_fit(fit, dir = "kann_mfi_final", overwrite = TRUE)
saveRDS(combined, "seurat_with_predictions.rds")

cat("\n✓ Analysis complete! Check kann_mfi_final/ for model artifacts.\n")
```

## Troubleshooting

### "Missing features" warning
```r
# Some genes not in your Seurat object
# Check which are missing:
fit$missing_genes

# Either:
# 1) Use genes that are present in your data
# 2) Add missing genes to Seurat (e.g., from raw counts)
```

### "Data hash mismatch" warning
```r
# Your expression matrix has changed since training
# This could mean:
# - Different normalization
# - Subset of cells
# - Updated Seurat object

# To proceed anyway:
seu <- predict_kann_to_seurat(seu, fit, validate = FALSE)
```

### Memory issues with large datasets
```r
# Use smaller batches
fit <- train_kann_on_seurat(
  seu, target_col = "MQI_v1",
  batch_size = 256  # Reduce if needed
)

# For importance on big data
imp <- perm_importance_seurat(
  seu, fit,
  batch_rows = 8192  # Reduce if OOM
)
```

## Citation

If you use this code in your research, please cite:

```
@software{kann_mfi_2025,
  author = {<Your Name>},
  title = {KANN-MFI: Kolmogorov-Arnold Networks for Mitochondrial Fitness},
  year = {2025},
  url = {https://github.com/<your-username>/kann-mfi}
}
```

## License

MIT License - see LICENSE file

## References

- Kolmogorov-Arnold representation theorem
- Liu et al. (2024) "KAN: Kolmogorov-Arnold Networks" arXiv:2404.19756
- Seurat: Hao et al. (2021) Cell
- torch for R: Falbel & Luraschi (2023)

---

**Questions?** Open an issue on GitHub or contact <your.email@institution.edu>
