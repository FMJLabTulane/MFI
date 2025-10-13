# ==============================================================================
# KANN-MFI Example Workflow
# ==============================================================================

library(Seurat)
library(qs)
library(ggplot2)

# Load functions
source("kann_mfi.R")

# ==============================================================================
# 1. LOAD YOUR DATA
# ==============================================================================

# Load Seurat object (edit path to your data)
combined <- qread("path/to/your/seurat.qs")
# OR: combined <- readRDS("path/to/your/seurat.rds")

DefaultAssay(combined) <- "RNA"

# Check what columns you have
head(combined@meta.data)

# ==============================================================================
# 2. CALCULATE MQI
# ==============================================================================

# This requires: Mitophagy, OXPHOS1, percent.mt
# Optional: TCA1, Glyco1
# UPRmtScore is calculated automatically

combined <- calculate_mqi(combined, target_col = "MQI_v1")

# Check the distribution
summary(combined$MQI_v1)
hist(combined$MQI_v1, breaks = 50, main = "MQI Distribution")

# ==============================================================================
# 3. TRAIN KANN MODEL
# ==============================================================================

# Train on all cells with MQI
fit <- train_kann(
  seu = combined,
  target_col = "MQI_v1",
  features = NULL,      # NULL = use mitochondrial regulators
  hidden = 48,          # hidden layer size
  n_bins = 10,          # spline complexity
  epochs = 300,
  lr = 0.002,
  seed = 42,
  verbose = TRUE
)

# ==============================================================================
# 4. PREDICT
# ==============================================================================

combined <- predict_kann(combined, fit, out_col = "MQI_pred")

# Check performance
idx <- which(is.finite(combined$MQI_v1) & is.finite(combined$MQI_pred))
r2 <- cor(combined$MQI_v1[idx], combined$MQI_pred[idx])^2
rmse <- sqrt(mean((combined$MQI_v1[idx] - combined$MQI_pred[idx])^2))

cat(sprintf("\nPerformance:\n"))
cat(sprintf("  n = %d cells\n", length(idx)))
cat(sprintf("  R² = %.4f\n", r2))
cat(sprintf("  RMSE = %.4f\n\n", rmse))

# ==============================================================================
# 5. PLOTS
# ==============================================================================

# Predicted vs Observed
p1 <- plot_predictions(combined, obs_col = "MQI_v1", pred_col = "MQI_pred")
ggsave("pred_vs_obs.pdf", p1, width = 6, height = 5)

# On UMAP (if you have UMAP)
if ("umap" %in% names(combined@reductions)) {
  p2 <- FeaturePlot(combined, features = c("MQI_v1", "MQI_pred"), ncol = 2)
  ggsave("umap_mqi.pdf", p2, width = 12, height = 5)
}

# ==============================================================================
# 6. PERMUTATION IMPORTANCE
# ==============================================================================

# This takes time! Start with repeats=3
imp <- perm_importance(
  seu = combined,
  fit = fit,
  repeats = 3,      # increase to 5+ for stability
  seed = 42,
  verbose = TRUE
)

# Top 20 genes
head(imp, 20)

# Plot
p3 <- plot_importance(imp, top_n = 30)
ggsave("importance.pdf", p3, width = 8, height = 7)

# Save to CSV
write.csv(
  data.frame(gene = names(imp), delta_mse = as.numeric(imp)),
  "importance_all_genes.csv",
  row.names = FALSE
)

# ==============================================================================
# 7. SAVE EVERYTHING
# ==============================================================================

# Save model
save_kann(fit, dir = "kann_model_v1", overwrite = TRUE)

# Save Seurat with predictions
qsave(combined, "seurat_with_predictions.qs")

# ==============================================================================
# 8. LATER: RELOAD AND PREDICT ON NEW DATA
# ==============================================================================

# Load saved model
# fit <- load_kann("kann_model_v1")

# Apply to new Seurat object
# new_seu <- predict_kann(new_seu, fit)

cat("\n✓ Done! Check the PDF files and CSV.\n")
