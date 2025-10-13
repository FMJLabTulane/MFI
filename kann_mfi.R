# ==============================================================================
# KANN-MFI: Kolmogorov-Arnold Networks for Mitochondrial Fitness Index
# ==============================================================================
# Complete pipeline: Auto-calculate all scores → KANN training → Importance
# Author: <Your Name> | MIT License | 2025
# 
# BEGINNER-FRIENDLY: Just give it a Seurat object, it handles the rest!
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(torch)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# GENE SETS FOR ALL METABOLIC PATHWAYS
# ==============================================================================

get_gene_sets <- function() {
  list(
    # Mitophagy genes
    mitophagy = c(
      "PINK1", "PRKN", "SQSTM1", "OPTN", "CALCOCO2", "TAX1BP1", "NBR1",
      "BNIP3", "BNIP3L", "FUNDC1", "PHB2", "FKBP8", "BCL2L13",
      "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "GABARAP", "GABARAPL1", "GABARAPL2",
      "ULK1", "ULK2", "ATG5", "ATG7", "ATG12", "TOMM20", "TOMM7"
    ),
    
    # OXPHOS genes (Complexes I-V)
    oxphos = c(
      # Complex I
      "NDUFB8", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2",
      # Complex II
      "SDHA", "SDHB", "SDHC", "SDHD",
      # Complex III
      "UQCRC1", "UQCRC2", "UQCRFS1", "CYC1",
      # Complex IV
      "COX4I1", "COX5A", "COX5B", "COX6C", "COX7A2", "COX7C",
      # Complex V
      "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5F1", "ATP5G1", "ATP5H", "ATP5O"
    ),
    
    # TCA cycle genes
    tca = c(
      "CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G",
      "OGDH", "SUCLA2", "SUCLG1", "SUCLG2", "SDHA", "FH", "MDH1", "MDH2"
    ),
    
    # Glycolysis genes
    glycolysis = c(
      "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP",
      "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2",
      "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "PKLR", "LDHA", "LDHB"
    ),
    
    # UPRmt genes
    uprmt = c(
      "LONP1", "CLPP", "CLPX", "HTRA2", "HSPD1", "HSPE1", "HSPA9",
      "DNAJA3", "YME1L1", "SPG7", "AFG3L2", "LRPPRC", "ATF5", "ATF4", "DDIT3"
    )
  )
}

# ==============================================================================
# AUTO-CALCULATE ALL METABOLIC SCORES
# ==============================================================================

#' Calculate a module score safely (handles missing genes)
calc_score_safe <- function(seu, genes, score_name, verbose = TRUE) {
  present <- genes[genes %in% rownames(seu)]
  
  if (length(present) == 0) {
    if (verbose) message(sprintf("  ⚠ %s: NO genes found! Setting to NA", score_name))
    seu[[score_name]] <- NA_real_
    return(seu)
  }
  
  if (length(present) < length(genes)) {
    missing_pct <- 100 * (length(genes) - length(present)) / length(genes)
    if (verbose) message(sprintf("  ⚠ %s: %d/%d genes found (%.0f%% missing)",
                                score_name, length(present), length(genes), missing_pct))
  } else {
    if (verbose) message(sprintf("  ✓ %s: %d/%d genes found",
                                score_name, length(present), length(genes)))
  }
  
  seu <- AddModuleScore(
    seu,
    features = list(present),
    name = paste0(score_name, "_temp"),
    nbin = 24,
    ctrl = 100,
    seed = 42
  )
  
  seu[[score_name]] <- seu[[paste0(score_name, "_temp1")]]
  seu[[paste0(score_name, "_temp1")]] <- NULL
  
  seu
}

#' Calculate percent.mt if missing
calc_percent_mt <- function(seu, verbose = TRUE) {
  if ("percent.mt" %in% colnames(seu@meta.data)) {
    if (verbose) message("  ✓ percent.mt: already present")
    return(seu)
  }
  
  # Try to find MT genes (works for human/mouse)
  mt_genes <- grep("^MT-|^mt-", rownames(seu), value = TRUE, ignore.case = TRUE)
  
  if (length(mt_genes) == 0) {
    if (verbose) message("  ⚠ percent.mt: no MT genes found! Setting to NA")
    seu$percent.mt <- NA_real_
    return(seu)
  }
  
  seu <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-", col.name = "percent.mt")
  if (verbose) message(sprintf("  ✓ percent.mt: calculated from %d genes", length(mt_genes)))
  
  seu
}

#' Calculate ALL required metabolic scores
#' 
#' This is the main function that checks what's missing and calculates it
calculate_metabolic_scores <- function(seu, force_recalc = FALSE, verbose = TRUE) {
  if (verbose) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  CALCULATING METABOLIC SCORES\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  gene_sets <- get_gene_sets()
  
  scores_to_calc <- list(
    list(name = "Mitophagy", genes = gene_sets$mitophagy),
    list(name = "OXPHOS1", genes = gene_sets$oxphos),
    list(name = "TCA1", genes = gene_sets$tca),
    list(name = "Glyco1", genes = gene_sets$glycolysis),
    list(name = "UPRmtScore", genes = gene_sets$uprmt)
  )
  
  for (score_info in scores_to_calc) {
    score_name <- score_info$name
    
    # Skip if already present and not forcing recalc
    if (score_name %in% colnames(seu@meta.data) && !force_recalc) {
      if (verbose) message(sprintf("  ✓ %s: already present (use force_recalc=TRUE to recalculate)", score_name))
      next
    }
    
    seu <- calc_score_safe(seu, score_info$genes, score_name, verbose = verbose)
  }
  
  # percent.mt
  seu <- calc_percent_mt(seu, verbose = verbose)
  
  if (verbose) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  ✓ ALL SCORES CALCULATED\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  seu
}

# ==============================================================================
# MQI CALCULATION
# ==============================================================================

#' Calculate Mitochondrial Quality Index
compute_MQI <- function(
  md,
  weights = list(mitophagy = 1.0, oxphos = 1.0, uprmt = 0.7, 
                 tca_pen = 0.2, mtfrac = 0.5, imbalance = 0.5),
  flips = list(mitophagy = FALSE, oxphos = FALSE, uprmt = FALSE, 
               tca = FALSE, glyco = FALSE, mtfrac = TRUE)
) {
  z <- function(x) if (all(!is.finite(x))) rep(NA_real_, length(x)) else as.numeric(scale(x))
  
  a <- list(
    mitophagy = md$Mitophagy,
    oxphos = md$OXPHOS1,
    tca = md$TCA1 %||% NA_real_,
    glyco = md$Glyco1 %||% NA_real_,
    mtfrac = md$percent.mt %||% NA_real_,
    uprmt = md$UPRmtScore %||% NA_real_
  )
  
  for (nm in names(a)) {
    if (isTRUE(flips[[nm]])) a[[nm]] <- -a[[nm]]
    a[[nm]] <- z(a[[nm]])
  }
  
  imb <- if (all(!is.finite(a$oxphos)) || all(!is.finite(a$glyco))) {
    rep(0, nrow(md))
  } else z(abs(a$oxphos - a$glyco))
  
  good <- weights$mitophagy * a$mitophagy + weights$oxphos * a$oxphos + weights$uprmt * a$uprmt
  pen <- 0
  if (any(is.finite(a$tca))) pen <- pen + weights$tca_pen * pmax(0, a$tca)
  if (any(is.finite(a$mtfrac))) pen <- pen + weights$mtfrac * z(a$mtfrac)
  if (any(is.finite(imb))) pen <- pen + weights$imbalance * imb
  
  as.numeric(good - pen)
}

#' Calculate MQI for Seurat object (with auto-calculation of missing scores)
calculate_mqi <- function(seu, target_col = "MQI_v1", force_recalc = FALSE, verbose = TRUE) {
  # Auto-calculate any missing metabolic scores
  seu <- calculate_metabolic_scores(seu, force_recalc = force_recalc, verbose = verbose)
  
  # Check what we have
  required <- c("Mitophagy", "OXPHOS1", "percent.mt")
  optional <- c("TCA1", "Glyco1", "UPRmtScore")
  
  has_required <- required %in% colnames(seu@meta.data)
  has_optional <- optional %in% colnames(seu@meta.data)
  
  if (verbose) {
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  MQI CALCULATION\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("\nRequired scores:\n")
    for (i in seq_along(required)) {
      cat(sprintf("  %s %s\n", ifelse(has_required[i], "✓", "✗"), required[i]))
    }
    cat("\nOptional scores (improve MQI):\n")
    for (i in seq_along(optional)) {
      cat(sprintf("  %s %s\n", ifelse(has_optional[i], "✓", "✗"), optional[i]))
    }
    cat("\n")
  }
  
  if (!all(has_required)) {
    stop("Missing required scores! Something went wrong in score calculation.")
  }
  
  # Calculate MQI
  seu@meta.data[[target_col]] <- compute_MQI(seu@meta.data)
  
  # Remove cells with NA MQI
  n_na <- sum(is.na(seu@meta.data[[target_col]]))
  if (n_na > 0 && verbose) {
    message(sprintf("  ⚠ Removing %d cells with NA MQI", n_na))
  }
  
  if (verbose) {
    finite_vals <- seu@meta.data[[target_col]][is.finite(seu@meta.data[[target_col]])]
    cat(sprintf("✓ MQI calculated for %d cells\n", length(finite_vals)))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(finite_vals), max(finite_vals)))
    cat(sprintf("  Mean: %.3f, SD: %.3f\n", mean(finite_vals), sd(finite_vals)))
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  seu
}

# ==============================================================================
# MITOCHONDRIAL REGULATOR GENE SETS
# ==============================================================================

get_regulators <- function() {
  unique(c(
    # Mitophagy core
    "PINK1","PRKN","TOMM7","TOMM20","TOMM22","SQSTM1","OPTN","CALCOCO2","TAX1BP1","NBR1","TBK1","TBKBP1","ULK1","ULK2","RB1CC1","ATG13",
    # Receptors
    "BNIP3","BNIP3L","FUNDC1","PHB2","FKBP8","BCL2L13","NLRX1","NIPSNAP1","NIPSNAP2","AMBRA1",
    # E3/DUB
    "MARCHF5","MUL1","RNF185","HUWE1","USP30","USP35","USP15",
    # Dynamics
    "DNM1L","MFF","MIEF1","MIEF2","FIS1","GDAP1","OPA1","MFN1","MFN2","OMA1","YME1L1","IMMT",
    # Biogenesis
    "PPARGC1A","PPARGC1B","TFAM","TFB2M","TFB1M","POLG","POLG2","POLRMT","NRF1","NFE2L2","ESRRA","ESRRB","ESRRG","PPARA","PPARD","PPARG","CREB1",
    # UPRmt
    "HSPD1","HSPE1","HSPA9","CLPP","LONP1","HTRA2","SPG7","AFG3L2","DNAJA3","LRPPRC","ATF4","DDIT3","ATF5","ATF2",
    # Antioxidant
    "SOD2","PRDX3","PRDX5","GPX1","GPX4","TXN2","GLRX2","CAT","HMOX1","NQO1","KEAP1",
    # Autophagy
    "BECN1","PIK3C3","PIK3R4","WIPI1","WIPI2","ATG5","ATG7","ATG3","ATG12","ATG16L1","UVRAG","TECPR1","GABARAP","GABARAPL1","GABARAPL2","MAP1LC3A","MAP1LC3B","MAP1LC3C",
    # ERphagy
    "FAM134A","FAM134B","FAM134C","RTN3","ATL3","SEC62","CCPG1","TEX264","CALCOCO1","VAPA","VAPB","UFM1","UFL1","DDRGK1","UFSP2","CDK5RAP3",
    # Signaling
    "PRKAA1","PRKAA2","PRKAB1","PRKAB2","PRKAG1","PRKAG2","PRKAG3","HIF1A","EPAS1","ARNT","ESR1","ESR2","AR","FOXO1","FOXO3","TFEB","TFE3","MITF","RELA","NFKB1","STAT3","JUN","FOS"
  ))
}

# ==============================================================================
# KANN ARCHITECTURE
# ==============================================================================

KANSpline1D <- nn_module(
  "KANSpline1D",
  initialize = function(n_bins = 10, xmin = -3, xmax = 3) {
    self$n_bins <- n_bins
    self$register_buffer("knots", torch_linspace(xmin, xmax, n_bins))
    self$coef <- nn_parameter(torch_zeros(n_bins))
  },
  forward = function(x) {
    x_exp <- x$unsqueeze(2)
    d <- (x_exp - self$knots$unsqueeze(1))$abs()
    w <- (self$knots[2] - self$knots[1])
    phi <- (1 - d / w)$clamp(min = 0)
    (phi * self$coef)$sum(dim = 2)
  }
)

KANLayer <- nn_module(
  "KANLayer",
  initialize = function(in_features, out_features, n_bins = 10, xmin = -3, xmax = 3) {
    self$in_features <- in_features
    self$out_features <- out_features
    self$splines <- nn_module_list(lapply(1:(out_features * in_features), function(i)
      KANSpline1D(n_bins = n_bins, xmin = xmin, xmax = xmax)
    ))
    self$linear <- nn_linear(in_features, out_features)
    self$bias <- nn_parameter(torch_zeros(out_features))
  },
  forward = function(x) {
    N <- x$size(1)
    outs <- vector("list", self$out_features)
    idx <- 1
    for (j in 1:self$out_features) {
      acc <- torch_zeros(N)
      for (i in 1:self$in_features) {
        acc <- acc + self$splines[[idx]](x[, i])
        idx <- idx + 1
      }
      outs[[j]] <- acc
    }
    S <- torch_stack(outs, dim = 2)
    S + self$linear(x) + self$bias$unsqueeze(1)
  }
)

KANNet <- nn_module(
  "KANNet",
  initialize = function(in_features, hidden = 48, out_features = 1, n_bins = 10) {
    self$layer1 <- KANLayer(in_features, hidden, n_bins = n_bins)
    self$act1 <- nn_relu()
    self$layer2 <- KANLayer(hidden, out_features, n_bins = n_bins)
  },
  forward = function(x) self$layer2(self$act1(self$layer1(x)))
)

# ==============================================================================
# TRAINING
# ==============================================================================

#' Train KANN model
train_kann <- function(
  seu,
  target_col = "MQI_v1",
  features = NULL,
  hidden = 48,
  n_bins = 10,
  epochs = 300,
  lr = 0.002,
  seed = 42,
  verbose = TRUE
) {
  set.seed(seed)
  torch_manual_seed(as.integer(seed))
  
  if (verbose) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  TRAINING KANN MODEL\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  # Get features
  if (is.null(features)) {
    features <- get_regulators()
    features <- intersect(features, rownames(seu))
    if (verbose) message(sprintf("Using %d mitochondrial regulators", length(features)))
    
    if (length(features) < 12) {
      if (verbose) message("Too few regulators, falling back to highly variable genes...")
      if (length(VariableFeatures(seu)) == 0) {
        seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      }
      features <- VariableFeatures(seu)
      if (verbose) message(sprintf("Using %d variable features instead", length(features)))
    }
  }
  
  # Build X, y
  if (!target_col %in% colnames(seu@meta.data)) {
    stop(sprintf("Target column '%s' not found! Run calculate_mqi() first.", target_col))
  }
  
  cells <- colnames(seu)[is.finite(seu@meta.data[[target_col]])]
  if (length(cells) == 0) {
    stop("No cells with finite MQI values!")
  }
  
  expr <- GetAssayData(seu, slot = "data")
  X <- t(as.matrix(expr[features, cells]))
  y <- seu@meta.data[cells, target_col]
  
  ok <- apply(X, 1, function(r) all(is.finite(r)))
  X <- X[ok, ]
  y <- y[ok]
  
  if (verbose) message(sprintf("Training data: %d cells × %d features", nrow(X), ncol(X)))
  
  # Standardize
  mu <- colMeans(X)
  sd <- apply(X, 2, sd)
  sd[sd == 0 | !is.finite(sd)] <- 1
  X <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
  
  # Train
  X_t <- torch_tensor(X, dtype = torch_float())
  y_t <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
  
  model <- KANNet(ncol(X), hidden = hidden, out_features = 1, n_bins = n_bins)
  opt <- optim_adam(model$parameters, lr = lr)
  loss_fn <- nn_mse_loss()
  
  if (verbose) message(sprintf("\nArchitecture: %d → %d → 1", ncol(X), hidden))
  if (verbose) message(sprintf("Training for %d epochs...\n", epochs))
  
  for (e in 1:epochs) {
    model$train()
    opt$zero_grad()
    pred <- model(X_t)
    loss <- loss_fn(pred, y_t)
    loss$backward()
    opt$step()
    
    if (verbose && (e %% 50 == 0 || e == 1)) {
      cat(sprintf("  Epoch %3d/%d: MSE = %.6f\n", e, epochs, loss$item()))
    }
  }
  
  if (verbose) {
    cat("\n✓ Training complete!\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  list(
    model = model,
    scaler = list(mu = mu, sd = sd, features = features),
    cells_used = names(y),
    config = list(hidden = hidden, n_bins = n_bins, epochs = epochs, lr = lr, seed = seed),
    target_col = target_col
  )
}

#' Predict with KANN
predict_kann <- function(seu, fit, out_col = "MQI_pred", verbose = TRUE) {
  if (verbose) message("Predicting MQI...")
  
  cells <- colnames(seu)
  expr <- GetAssayData(seu, slot = "data")
  
  missing <- setdiff(fit$scaler$features, rownames(expr))
  if (length(missing) > 0) {
    warning(sprintf("%d/%d features missing", length(missing), length(fit$scaler$features)))
  }
  
  features <- intersect(fit$scaler$features, rownames(expr))
  X <- t(as.matrix(expr[features, cells]))
  X <- sweep(sweep(X, 2, fit$scaler$mu[features], "-"), 2, fit$scaler$sd[features], "/")
  
  fit$model$eval()
  pred <- as.numeric(fit$model(torch_tensor(X, dtype = torch_float())))
  
  seu[[out_col]] <- NA_real_
  seu[[out_col]][cells] <- pred
  
  if (verbose) message(sprintf("✓ Predictions saved to '%s'", out_col))
  
  seu
}

# ==============================================================================
# PERMUTATION IMPORTANCE
# ==============================================================================

predict_batches <- function(model, X, batch = 32768) {
  n <- nrow(X)
  out <- numeric(n)
  i <- 1
  while (i <= n) {
    j <- min(i + batch - 1, n)
    out[i:j] <- as.numeric(model(torch_tensor(X[i:j, , drop = FALSE], dtype = torch_float())))
    i <- j + 1
  }
  out
}

#' Permutation importance
perm_importance <- function(
  seu,
  fit,
  repeats = 3,
  seed = 42,
  batch = 32768,
  verbose = TRUE
) {
  set.seed(seed)
  torch_manual_seed(as.integer(seed))
  
  if (verbose) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  PERMUTATION IMPORTANCE\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  # Build X
  cells <- fit$cells_used[fit$cells_used %in% colnames(seu)]
  expr <- GetAssayData(seu, slot = "data")
  X <- t(as.matrix(expr[fit$scaler$features, cells]))
  X <- sweep(sweep(X, 2, fit$scaler$mu, "-"), 2, fit$scaler$sd, "/")
  
  y <- seu@meta.data[cells, fit$target_col]
  
  if (verbose) message(sprintf("Computing importance for %d features (%d repeats)...", ncol(X), repeats))
  
  fit$model$eval()
  base_pred <- predict_batches(fit$model, X, batch)
  base_mse <- mean((base_pred - y)^2)
  
  p <- ncol(X)
  imp <- numeric(p)
  names(imp) <- colnames(X)
  
  if (verbose) pb <- txtProgressBar(0, p, style = 3)
  
  for (j in 1:p) {
    mses <- numeric(repeats)
    for (r in 1:repeats) {
      Xp <- X
      Xp[, j] <- X[sample(nrow(X)), j]
      pred <- predict_batches(fit$model, Xp, batch)
      mses[r] <- mean((pred - y)^2)
    }
    imp[j] <- mean(mses) - base_mse
    if (verbose) setTxtProgressBar(pb, j)
  }
  
  if (verbose) {
    close(pb)
    cat("\n\n✓ Importance calculation complete!\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  sort(imp, decreasing = TRUE)
}

# ==============================================================================
# SAVE/LOAD
# ==============================================================================

save_kann <- function(fit, dir, overwrite = FALSE) {
  if (dir.exists(dir) && !overwrite) stop("Directory exists. Use overwrite=TRUE")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  torch_save(fit$model$state_dict(), file.path(dir, "model.pt"))
  qsave(fit$scaler, file.path(dir, "scaler.qs"))
  qsave(fit$config, file.path(dir, "config.qs"))
  qsave(fit$target_col, file.path(dir, "target_col.qs"))
  
  writeLines(c(
    "KANN-MFI Model",
    sprintf("Features: %d", length(fit$scaler$features)),
    sprintf("Hidden: %d", fit$config$hidden),
    sprintf("Target: %s", fit$target_col),
    sprintf("Saved: %s", Sys.time())
  ), file.path(dir, "README.txt"))
  
  message(sprintf("✓ Model saved to: %s/", dir))
}

load_kann <- function(dir) {
  scaler <- qread(file.path(dir, "scaler.qs"))
  config <- qread(file.path(dir, "config.qs"))
  target_col <- qread(file.path(dir, "target_col.qs"))
  
  model <- KANNet(
    length(scaler$features),
    hidden = config$hidden,
    out_features = 1,
    n_bins = config$n_bins
  )
  
  state <- torch_load(file.path(dir, "model.pt"))
  model$load_state_dict(state)
  model$eval()
  
  message(sprintf("✓ Model loaded from: %s/", dir))
  
  list(model = model, scaler = scaler, config = config, target_col = target_col)
}

# ==============================================================================
# PLOTTING
# ==============================================================================

plot_predictions <- function(seu, obs_col = "MQI_v1", pred_col = "MQI_pred") {
  df <- seu@meta.data[, c(obs_col, pred_col)]
  df <- df[complete.cases(df), ]
  
  r2 <- cor(df[[obs_col]], df[[pred_col]])^2
  rmse <- sqrt(mean((df[[obs_col]] - df[[pred_col]])^2))
  
  ggplot(df, aes_string(x = obs_col, y = pred_col)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
            label = sprintf("R² = %.3f\nRMSE = %.3f\nn = %d", r2, rmse, nrow(df))) +
    labs(title = "KANN Predictions", x = "Observed MQI", y = "Predicted MQI") +
    theme_classic(base_size = 12)
}

plot_importance <- function(imp, top_n = 30, title = "Permutation Importance") {
  df <- data.frame(
    gene = names(head(imp, top_n)),
    importance = as.numeric(head(imp, top_n))
  )
  df$gene <- factor(df$gene, levels = rev(df$gene))
  
  ggplot(df, aes(x = gene, y = importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = title, x = NULL, y = "ΔMSE (increase in prediction error)") +
    theme_classic(base_size = 12)
}

cat("\n✓ KANN-MFI functions loaded!\n\n")
cat("Quick start:\n")
cat("  1. seu <- calculate_mqi(seu)\n")
cat("  2. fit <- train_kann(seu, epochs=300)\n")
cat("  3. seu <- predict_kann(seu, fit)\n")
cat("  4. plot_predictions(seu)\n\n")
