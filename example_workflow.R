# ==============================================================================
# Example Workflow: KANN-MFI Reproducible Analysis
# ==============================================================================
# This script demonstrates a complete analysis from training to validation
# with full reproducibility guarantees.

library(Seurat)
library(ggplot2)
source("kann_from_seurat.R")

# ==============================================================================
# 1. PREPARE DATA
# ==============================================================================

# Load your Seurat object
combined <- readRDS("path/to/your/seurat_object.rds")
DefaultAssay(combined) <- "RNA"

# Check that target exists and has values
stopifnot("MQI_v1" %in% colnames(combined@meta.data))
cat(sprintf("Target column 'MQI_v1': %d finite values out of %d cells\n",
           sum(is.finite(combined$MQI_v1)), ncol(combined)))

# Optional: Define specific gene set (or use NULL for HVGs)
# my_genes <- c("MT-ND1", "MT-CO1", "PPARGC1A", "TFAM", ...)
my_genes <- NULL

# ==============================================================================
# 2. TRAIN MODEL (with full reproducibility)
# ==============================================================================

cat("\n=== TRAINING KANN MODEL ===\n")

fit <- train_kann_on_seurat(
  seu = combined,
  target_col = "MQI_v1",
  assay = "RNA",
  slot = "data",           # log-normalized expression
  features = my_genes,     # NULL = use HVGs
  hidden = 48,             # hidden layer size
  n_bins = 10,             # spline knots per feature
  epochs = 300,
  lr = 2e-3,
  batch_size = NULL,       # NULL = full batch
  seed = 42,               # for reproducibility
  verbose = TRUE
)

# Inspect the fit
print(fit)

# ==============================================================================
# 3. PREDICT AND EVALUATE
# ==============================================================================

cat("\n=== MAKING PREDICTIONS ===\n")

combined <- predict_kann_to_seurat(
  seu = combined,
  fit = fit,
  out_col = "KANN_MFI_pred",
  resid_col = "KANN_MFI_resid",
  validate = TRUE,         # check data integrity
  verbose = TRUE
)

# Compute performance metrics
idx <- which(is.finite(combined$MQI_v1) & is.finite(combined$KANN_MFI_pred))
obs <- combined$MQI_v1[idx]
pred <- combined$KANN_MFI_pred[idx]

r2 <- cor(obs, pred)^2
rmse <- sqrt(mean((obs - pred)^2))
mae <- mean(abs(obs - pred))

cat("\n=== PERFORMANCE METRICS ===\n")
cat(sprintf("n = %d cells with predictions\n", length(idx)))
cat(sprintf("R² = %.4f\n", r2))
cat(sprintf("RMSE = %.4f\n", rmse))
cat(sprintf("MAE = %.4f\n\n", mae))

# ==============================================================================
# 4. DIAGNOSTIC PLOTS
# ==============================================================================

cat("=== GENERATING DIAGNOSTIC PLOTS ===\n")

# Create plots
plots <- plot_kann_diagnostics(combined, fit, out_col = "KANN_MFI_pred")

# Save to PDF
dir.create("figures", showWarnings = FALSE)
ggsave("figures/pred_vs_obs.pdf", plots$pred_vs_obs, width = 6, height = 5)
ggsave("figures/residuals.pdf", plots$residuals, width = 6, height = 5)
ggsave("figures/residual_hist.pdf", plots$residual_hist, width = 6, height = 4)

cat("✓ Plots saved to figures/\n")

# ==============================================================================
# 5. FEATURE IMPORTANCE
# ==============================================================================

cat("\n=== COMPUTING FEATURE IMPORTANCE ===\n")

imp <- perm_importance_seurat(
  seu = combined,
  fit = fit,
  repeats = 5,           # increase for more stable estimates
  seed = 73,
  verbose = TRUE
)

# Save full results
dir.create("results", showWarnings = FALSE)
write.csv(
  data.frame(
    gene = names(imp),
    delta_mse = as.numeric(imp),
    rank = seq_along(imp)
  ),
  "results/importance_full.csv",
  row.names = FALSE
)

# Show top 20
cat("\nTop 20 most important features:\n")
print(head(imp, 20))

# Plot top 20
imp_df <- data.frame(
  gene = names(imp)[1:20],
  importance = as.numeric(imp)[1:20]
)
imp_df$gene <- factor(imp_df$gene, levels = rev(imp_df$gene))

p_imp <- ggplot(imp_df, aes(x = importance, y = gene)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(
    title = "Top 20 Features by Permutation Importance",
    x = "Importance (ΔMSE)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("figures/importance_top20.pdf", p_imp, width = 7, height = 6)
cat("✓ Importance saved to results/\n")

# ==============================================================================
# 6. SAVE MODEL BUNDLE (for reproducibility)
# ==============================================================================

cat("\n=== SAVING MODEL ARTIFACTS ===\n")

# Save complete model bundle with all metadata
save_kann_fit(
  fit,
  dir = "kann_model_v1",
  overwrite = TRUE
)

# Save updated Seurat object with predictions
saveRDS(combined, "seurat_with_kann_predictions.rds")

cat("✓ Model bundle saved to kann_model_v1/\n")
cat("✓ Seurat object saved to seurat_with_kann_predictions.rds\n")

# ==============================================================================
# 7. DEMONSTRATE RELOAD AND VALIDATION
# ==============================================================================

cat("\n=== DEMONSTRATING MODEL RELOAD ===\n")

# Clear and reload (simulating a new session)
rm(fit)

# Load model from disk
fit_reloaded <- load_kann_fit("kann_model_v1")

# Validate against Seurat object
cat("\nValidating reloaded model:\n")
validation <- validate_kann_fit(combined, fit_reloaded, verbose = TRUE)

# Predict again (should give identical results)
combined_check <- predict_kann_to_seurat(
  combined,
  fit_reloaded,
  out_col = "KANN_MFI_check",
  validate = TRUE,
  verbose = FALSE
)

# Verify predictions are identical
identical_preds <- all.equal(
  combined$KANN_MFI_pred,
  combined_check$KANN_MFI_check,
  tolerance = 1e-6
)

cat(sprintf("\nPredictions identical after reload: %s\n",
           ifelse(isTRUE(identical_preds), "✓ YES", "✗ NO")))

# ==============================================================================
# 8. SUMMARY AND NEXT STEPS
# ==============================================================================

cat("\n" %||% paste0(rep("=", 70), collapse = "") %||% "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste0(rep("=", 70), collapse = "") %||% "\n")
cat(sprintf("Model: %d genes → %d hidden → 1 output\n",
           length(fit_reloaded$used_genes), fit_reloaded$config$hidden))
cat(sprintf("Performance: R²=%.4f, RMSE=%.4f (n=%d)\n", r2, rmse, length(idx)))
cat(sprintf("Top gene: %s (ΔMSE=%.6f)\n", names(imp)[1], imp[1]))
cat("\nArtifacts saved:\n")
cat("  • kann_model_v1/              (complete model bundle)\n")
cat("  • seurat_with_kann_predictions.rds\n")
cat("  • figures/*.pdf               (diagnostic plots)\n")
cat("  • results/importance_full.csv\n")
cat("\nTo reproduce in a new session:\n")
cat("  fit <- load_kann_fit('kann_model_v1')\n")
cat("  seu <- predict_kann_to_seurat(seu, fit)\n")
cat(paste0(rep("=", 70), collapse = "") %||% "\n\n")

# ==============================================================================
# OPTIONAL: Additional analyses
# ==============================================================================

# # Visualize predictions on UMAP
# if ("UMAP_1" %in% colnames(combined@reductions$umap@cell.embeddings)) {
#   p_umap <- FeaturePlot(combined, features = c("MQI_v1", "KANN_MFI_pred"),
#                        blend = FALSE, ncol = 2)
#   ggsave("figures/umap_predictions.pdf", p_umap, width = 12, height = 5)
# }
# 
# # Compare predictions across cell types
# if ("cell_type" %in% colnames(combined@meta.data)) {
#   library(dplyr)
#   by_type <- combined@meta.data %>%
#     group_by(cell_type) %>%
#     summarise(
#       n = n(),
#       mqi_mean = mean(MQI_v1, na.rm = TRUE),
#       pred_mean = mean(KANN_MFI_pred, na.rm = TRUE),
#       cor = cor(MQI_v1, KANN_MFI_pred, use = "complete.obs")
#     ) %>%
#     arrange(desc(cor))
#   write.csv(by_type, "results/performance_by_celltype.csv", row.names = FALSE)
# }
# 
# # Train/test split for honest evaluation
# set.seed(999)
# train_cells <- sample(colnames(combined), size = floor(0.8 * ncol(combined)))
# test_cells <- setdiff(colnames(combined), train_cells)
# 
# fit_honest <- train_kann_on_seurat(
#   seu = combined,
#   target_col = "MQI_v1",
#   cells_filter = train_cells,
#   seed = 42
# )
# 
# seu_test <- subset(combined, cells = test_cells)
# seu_test <- predict_kann_to_seurat(seu_test, fit_honest)
# 
# test_idx <- which(is.finite(seu_test$MQI_v1) & is.finite(seu_test$KANN_pred))
# test_r2 <- cor(seu_test$MQI_v1[test_idx], seu_test$KANN_pred[test_idx])^2
# cat(sprintf("Held-out test R² = %.4f\n", test_r2))
