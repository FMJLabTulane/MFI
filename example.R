# ==============================================================================
# KANN-MFI Complete Example
# ==============================================================================
# For someone with just a basic Seurat object!
# You only need normalized counts - everything else is calculated automatically
# ==============================================================================

library(Seurat)
library(ggplot2)

# Load the functions
source("kann_mfi.R")

# ==============================================================================
# 1. LOAD YOUR DATA
# ==============================================================================

# Load your Seurat object (edit this path!)
seu <- readRDS("path/to/your/seurat.rds")
# OR if using .qs format:
# library(qs)
# seu <- qread("path/to/your/seurat.qs")

# Make sure you're using RNA assay
DefaultAssay(seu) <- "RNA"

# Check that you have normalized data
# If you don't, run: seu <- NormalizeData(seu)

# ==============================================================================
# 2. CALCULATE MQI (Everything is automatic!)
# ==============================================================================

# This will:
# - Calculate Mitophagy score (from 27 mitophagy genes)
# - Calculate OXPHOS score (from Complex I-V genes)
# - Calculate TCA score (from TCA cycle genes)
# - Calculate Glycolysis score (from glycolysis genes)
# - Calculate UPRmt score (from 15 UPRmt genes)
# - Calculate percent.mt (mitochondrial fraction)
# - Combine everything into MQI

seu <- calculate_mqi(seu, target_col = "MQI_v1")

# Look at the results
summary(seu$MQI_v1)
hist(seu$MQI_v1, breaks = 50, main = "MQI Distribution", col = "lightblue")

# Check individual scores if you want
head(seu@meta.data[, c("Mitophagy", "OXPHOS1", "TCA1", "Glyco1", 
                       "UPRmtScore", "percent.mt", "MQI_v1")])

# ==============================================================================
# 3. TRAIN KANN MODEL
# ==============================================================================

# This trains a neural network to predict MQI from gene expression
# It automatically uses ~200 mitochondrial regulator genes

fit <- train_kann(
  seu = seu,
  target_col = "MQI_v1",
  hidden = 48,        # neural network size
  epochs = 300,       # training iterations
  lr = 0.002,         # learning rate
  seed = 42           # for reproducibility
)

# ==============================================================================
# 4. PREDICT MQI FOR ALL CELLS
# ==============================================================================

# Now use the trained model to predict MQI for all cells
seu <- predict_kann(seu, fit, out_col = "MQI_pred")

# Check how well it worked
idx <- which(is.finite(seu$MQI_v1) & is.finite(seu$MQI_pred))
r2 <- cor(seu$MQI_v1[idx], seu$MQI_pred[idx])^2
rmse <- sqrt(mean((seu$MQI_v1[idx] - seu$MQI_pred[idx])^2))

cat("\n═══════════════════════════════════════\n")
cat("PERFORMANCE:\n")
cat("═══════════════════════════════════════\n")
cat(sprintf("  n cells:  %d\n", length(idx)))
cat(sprintf("  R²:       %.4f\n", r2))
cat(sprintf("  RMSE:     %.4f\n", rmse))
cat("═══════════════════════════════════════\n\n")

# ==============================================================================
# 5. MAKE PLOTS
# ==============================================================================

# Predicted vs Observed
p1 <- plot_predictions(seu, obs_col = "MQI_v1", pred_col = "MQI_pred")
ggsave("01_pred_vs_obs.pdf", p1, width = 6, height = 5)
print(p1)

# Show MQI on UMAP (if you have UMAP)
if ("umap" %in% names(seu@reductions)) {
  p2 <- FeaturePlot(seu, features = c("MQI_v1", "MQI_pred"), ncol = 2) +
    plot_annotation(title = "MQI: Observed vs Predicted")
  ggsave("02_umap_mqi.pdf", p2, width = 12, height = 5)
  print(p2)
}

# Distribution comparison
p3 <- ggplot(seu@meta.data, aes(x = MQI_v1)) +
  geom_density(fill = "blue", alpha = 0.3) +
  geom_density(aes(x = MQI_pred), fill = "red", alpha = 0.3) +
  labs(title = "MQI Distribution: Observed (blue) vs Predicted (red)",
       x = "MQI", y = "Density") +
  theme_classic()
ggsave("03_mqi_distributions.pdf", p3, width = 7, height = 4)
print(p3)

# ==============================================================================
# 6. FIND IMPORTANT GENES (What drives MQI?)
# ==============================================================================

# This takes a while! It tests each gene by shuffling it
# Start with repeats=3 (faster), increase to 5+ for publication

cat("\nCalculating feature importance...\n")
cat("(This will take several minutes)\n\n")

imp <- perm_importance(
  seu = seu,
  fit = fit,
  repeats = 3,      # increase to 5 for more stable results
  seed = 42
)

# Show top 20 genes
cat("\nTop 20 most important genes:\n")
print(head(imp, 20))

# Plot top 30
p4 <- plot_importance(imp, top_n = 30, title = "Top 30 Genes by Importance")
ggsave("04_importance_top30.pdf", p4, width = 8, height = 7)
print(p4)

# Save full importance list to CSV
write.csv(
  data.frame(gene = names(imp), importance = as.numeric(imp)),
  "importance_all_genes.csv",
  row.names = FALSE
)

# ==============================================================================
# 7. SAVE EVERYTHING
# ==============================================================================

# Save the trained model (so you can reuse it later)
save_kann(fit, dir = "kann_model_v1", overwrite = TRUE)

# Save your Seurat object with all the new scores and predictions
saveRDS(seu, "seurat_with_mqi_predictions.rds")
# OR: library(qs); qsave(seu, "seurat_with_mqi_predictions.qs")

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("✓ ANALYSIS COMPLETE!\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("\nFiles created:\n")
cat("  • 01_pred_vs_obs.pdf\n")
cat("  • 02_umap_mqi.pdf\n")
cat("  • 03_mqi_distributions.pdf\n")
cat("  • 04_importance_top30.pdf\n")
cat("  • importance_all_genes.csv\n")
cat("  • kann_model_v1/  (saved model)\n")
cat("  • seurat_with_mqi_predictions.rds\n")
cat("\n")
cat("Your Seurat object now has these new columns:\n")
cat("  • Mitophagy, OXPHOS1, TCA1, Glyco1, UPRmtScore, percent.mt\n")
cat("  • MQI_v1 (calculated MQI)\n")
cat("  • MQI_pred (KANN prediction)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# ==============================================================================
# 8. LATER: RELOAD THE MODEL
# ==============================================================================

# In a future R session, you can reload and reuse the model:

# library(Seurat)
# source("kann_mfi.R")
# 
# # Load your saved model
# fit <- load_kann("kann_model_v1")
# 
# # Apply to a new Seurat object
# new_seu <- readRDS("new_data.rds")
# new_seu <- calculate_mqi(new_seu)  # calculate MQI first
# new_seu <- predict_kann(new_seu, fit)  # then predict

# ==============================================================================
# BONUS: Check which genes were actually used
# ==============================================================================

cat("\nGenes used in the model:\n")
cat(sprintf("  Total: %d genes\n", length(fit$scaler$features)))
cat("\n  First 20:\n")
print(head(fit$scaler$features, 20))

# See which mitochondrial gene sets had good coverage
gene_sets <- get_gene_sets()
coverage <- lapply(gene_sets, function(genes) {
  present <- sum(genes %in% rownames(seu))
  total <- length(genes)
  pct <- 100 * present / total
  c(present = present, total = total, pct = pct)
})

cat("\nGene set coverage in your data:\n")
for (name in names(coverage)) {
  cat(sprintf("  %s: %d/%d (%.0f%%)\n",
             name,
             coverage[[name]]["present"],
             coverage[[name]]["total"],
             coverage[[name]]["pct"]))
}
