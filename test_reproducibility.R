# ==============================================================================
# Reproducibility Test Suite for KANN-MFI
# ==============================================================================
# Run this script to verify that your implementation produces identical results
# across multiple runs with the same seed.

source("kann_from_seurat.R")

cat("=============================================================\n")
cat("KANN-MFI REPRODUCIBILITY TEST SUITE\n")
cat("=============================================================\n\n")

# ==============================================================================
# Test 1: Synthetic Data Reproducibility
# ==============================================================================

cat("Test 1: Synthetic data reproducibility\n")
cat("---------------------------------------\n")

# Create synthetic Seurat object
set.seed(12345)
n_cells <- 1000
n_genes <- 100

# Generate expression data
counts <- matrix(
  rpois(n_cells * n_genes, lambda = 10),
  nrow = n_genes,
  ncol = n_cells,
  dimnames = list(
    paste0("Gene_", 1:n_genes),
    paste0("Cell_", 1:n_cells)
  )
)

# Create Seurat object
seu_test <- CreateSeuratObject(counts = counts)
seu_test <- NormalizeData(seu_test, verbose = FALSE)

# Create synthetic target (linear combination of genes with noise)
coef_true <- rnorm(20)
gene_idx <- sample(n_genes, 20)
X_subset <- t(as.matrix(GetAssayData(seu_test, slot = "data")[gene_idx, ]))
seu_test$MQI_synthetic <- as.numeric(X_subset %*% coef_true + rnorm(n_cells, sd = 0.5))

# Train model twice with same seed
fit1 <- train_kann_on_seurat(
  seu = seu_test,
  target_col = "MQI_synthetic",
  hidden = 24,
  n_bins = 5,
  epochs = 50,
  seed = 42,
  verbose = FALSE
)

fit2 <- train_kann_on_seurat(
  seu = seu_test,
  target_col = "MQI_synthetic",
  hidden = 24,
  n_bins = 5,
  epochs = 50,
  seed = 42,
  verbose = FALSE
)

# Predict with both
seu_test <- predict_kann_to_seurat(seu_test, fit1, out_col = "pred1", verbose = FALSE)
seu_test <- predict_kann_to_seurat(seu_test, fit2, out_col = "pred2", verbose = FALSE)

# Compare predictions
max_diff <- max(abs(seu_test$pred1 - seu_test$pred2), na.rm = TRUE)
identical_check <- max_diff < 1e-6

cat(sprintf("  Max difference: %.2e\n", max_diff))
cat(sprintf("  Result: %s\n\n", ifelse(identical_check, "✓ PASS", "✗ FAIL")))

if (!identical_check) {
  cat("  ERROR: Predictions differ! Check random seed handling.\n\n")
}

# ==============================================================================
# Test 2: Save/Load Consistency
# ==============================================================================

cat("Test 2: Save/load consistency\n")
cat("------------------------------\n")

# Save model
temp_dir <- tempfile("kann_test_")
save_kann_fit(fit1, dir = temp_dir, overwrite = TRUE)
cat(sprintf("  Saved to: %s\n", temp_dir))

# Load model
fit_loaded <- load_kann_fit(temp_dir)

# Predict with loaded model
seu_test <- predict_kann_to_seurat(seu_test, fit_loaded, out_col = "pred_loaded", verbose = FALSE)

# Compare
max_diff_load <- max(abs(seu_test$pred1 - seu_test$pred_loaded), na.rm = TRUE)
load_check <- max_diff_load < 1e-6

cat(sprintf("  Max difference: %.2e\n", max_diff_load))
cat(sprintf("  Result: %s\n\n", ifelse(load_check, "✓ PASS", "✗ FAIL")))

# Clean up
unlink(temp_dir, recursive = TRUE)

# ==============================================================================
# Test 3: Data Hash Validation
# ==============================================================================

cat("Test 3: Data hash validation\n")
cat("-----------------------------\n")

# Compute hash
hash1 <- compute_data_hash(
  seu_test, "RNA", "data",
  fit1$used_genes, fit1$used_cells
)

# Verify hash matches provenance
hash_match <- hash1 == fit1$provenance$data_source$data_hash
cat(sprintf("  Hash match: %s\n", ifelse(hash_match, "✓ PASS", "✗ FAIL")))

# Modify data slightly
seu_test_mod <- seu_test
seu_test_mod@assays$RNA@data[1, 1] <- seu_test_mod@assays$RNA@data[1, 1] + 0.001

# Compute hash of modified data
hash2 <- compute_data_hash(
  seu_test_mod, "RNA", "data",
  fit1$used_genes, fit1$used_cells
)

# Verify hashes differ
hash_differ <- hash1 != hash2
cat(sprintf("  Hash detects change: %s\n\n", ifelse(hash_differ, "✓ PASS", "✗ FAIL")))

# ==============================================================================
# Test 4: Scaler Consistency
# ==============================================================================

cat("Test 4: Scaler consistency\n")
cat("--------------------------\n")

# Extract training data
X_train <- t(as.matrix(GetAssayData(seu_test, slot = "data")[fit1$used_genes, fit1$used_cells]))

# Apply scaler
Xz1 <- scaler_apply(X_train, fit1$scaler)

# Refit scaler independently
sc_test <- scaler_fit(X_train)
Xz2 <- scaler_apply(X_train, sc_test)

# Compare
max_diff_scaler <- max(abs(Xz1 - Xz2))
scaler_check <- max_diff_scaler < 1e-10

cat(sprintf("  Max difference: %.2e\n", max_diff_scaler))
cat(sprintf("  Result: %s\n\n", ifelse(scaler_check, "✓ PASS", "✗ FAIL")))

# ==============================================================================
# Test 5: Validation Function
# ==============================================================================

cat("Test 5: Validation function\n")
cat("---------------------------\n")

# Should pass on original object
val1 <- validate_kann_fit(seu_test, fit1, verbose = FALSE)
cat(sprintf("  Original object: %s\n", ifelse(val1$valid, "✓ PASS", "✗ FAIL")))

# Create incompatible object (remove features)
seu_incompatible <- seu_test
seu_incompatible@assays$RNA@data <- seu_incompatible@assays$RNA@data[-1:-50, ]

# Should detect incompatibility
val2 <- validate_kann_fit(seu_incompatible, fit1, verbose = FALSE)
detects_incompatible <- !val2$valid
cat(sprintf("  Detects incompatible: %s\n\n", 
           ifelse(detects_incompatible, "✓ PASS", "✗ FAIL")))

# ==============================================================================
# Test 6: Permutation Importance Reproducibility
# ==============================================================================

cat("Test 6: Importance reproducibility\n")
cat("-----------------------------------\n")

# Run importance twice with same seed
imp1 <- perm_importance_seurat(seu_test, fit1, repeats = 2, seed = 99, verbose = FALSE)
imp2 <- perm_importance_seurat(seu_test, fit1, repeats = 2, seed = 99, verbose = FALSE)

# Compare
max_diff_imp <- max(abs(imp1 - imp2))
imp_check <- max_diff_imp < 1e-6

cat(sprintf("  Max difference: %.2e\n", max_diff_imp))
cat(sprintf("  Result: %s\n\n", ifelse(imp_check, "✓ PASS", "✗ FAIL")))

# ==============================================================================
# Summary
# ==============================================================================

all_tests <- c(
  identical_check,
  load_check,
  hash_match && hash_differ,
  scaler_check,
  val1$valid && detects_incompatible,
  imp_check
)

cat("=============================================================\n")
cat("SUMMARY\n")
cat("=============================================================\n")
cat(sprintf("Tests passed: %d / %d\n", sum(all_tests), length(all_tests)))

if (all(all_tests)) {
  cat("\n✓ ALL TESTS PASSED\n")
  cat("Your installation is working correctly and produces\n")
  cat("reproducible results.\n")
} else {
  cat("\n✗ SOME TESTS FAILED\n")
  cat("Please review the output above and check:\n")
  cat("  - torch installation\n")
  cat("  - Random seed handling\n")
  cat("  - File I/O operations\n")
}
cat("=============================================================\n")

# Return test results invisibly
invisible(list(
  tests = c(
    "Reproducibility" = identical_check,
    "Save/Load" = load_check,
    "Hash validation" = hash_match && hash_differ,
    "Scaler" = scaler_check,
    "Validation" = val1$valid && detects_incompatible,
    "Importance" = imp_check
  ),
  all_passed = all(all_tests)
))
