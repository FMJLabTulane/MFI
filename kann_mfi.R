# ==============================================================================
# KANN-MFI: Kolmogorov-Arnold Networks for Mitochondrial Fitness Index
# ==============================================================================
# Complete pipeline: MQI calculation → KANN training → Permutation importance
# Author: <Your Name> | MIT License | 2025
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
# 1) MQI CALCULATION
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

#' Add UPRmt score to Seurat
add_uprmt <- function(seu) {
  if ("UPRmtScore" %in% colnames(seu@meta.data)) return(seu)
  
  genes <- c("LONP1","CLPP","CLPX","HTRA2","HSPD1","HSPE1","HSPA9",
            "DNAJA3","YME1L1","SPG7","AFG3L2","LRPPRC","ATF5","ATF4","DDIT3")
  present <- genes[genes %in% rownames(seu)]
  
  if (length(present) == 0) {
    seu$UPRmtScore <- NA_real_
    return(seu)
  }
  
  seu <- AddModuleScore(seu, features = list(present), name = "UPRmt", nbin = 24, ctrl = 100, seed = 42)
  seu$UPRmtScore <- seu$UPRmt1
  seu$UPRmt1 <- NULL
  seu
}

#' Calculate MQI for Seurat object
calculate_mqi <- function(seu, target_col = "MQI_v1") {
  seu <- add_uprmt(seu)
  required <- c("Mitophagy", "OXPHOS1", "percent.mt")
  missing <- setdiff(required, colnames(seu@meta.data))
  if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = ", "))
  
  seu@meta.data[[target_col]] <- compute_MQI(seu@meta.data)
  message(sprintf("✓ MQI calculated: range [%.2f, %.2f], mean %.2f",
                 min(seu[[target_col]], na.rm = TRUE),
                 max(seu[[target_col]], na.rm = TRUE),
                 mean(seu[[target_col]], na.rm = TRUE)))
  seu
}

# ==============================================================================
# 2) REGULATOR GENE SETS
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
# 3) KANN ARCHITECTURE
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
# 4) TRAINING
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
  
  # Get features
  if (is.null(features)) {
    features <- get_regulators()
    features <- intersect(features, rownames(seu))
    if (length(features) < 12) {
      if (length(VariableFeatures(seu)) == 0) {
        seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      }
      features <- VariableFeatures(seu)
    }
  }
  
  # Build X, y
  cells <- colnames(seu)[is.finite(seu@meta.data[[target_col]])]
  expr <- GetAssayData(seu, slot = "data")
  X <- t(as.matrix(expr[features, cells]))
  y <- seu@meta.data[cells, target_col]
  
  ok <- apply(X, 1, function(r) all(is.finite(r)))
  X <- X[ok, ]
  y <- y[ok]
  
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
  
  if (verbose) message(sprintf("Training: %d cells × %d features → %d hidden → 1", nrow(X), ncol(X), hidden))
  
  for (e in 1:epochs) {
    model$train()
    opt$zero_grad()
    pred <- model(X_t)
    loss <- loss_fn(pred, y_t)
    loss$backward()
    opt$step()
    
    if (verbose && (e %% 50 == 0 || e == 1)) {
      cat(sprintf("Epoch %3d: MSE=%.6f\n", e, loss$item()))
    }
  }
  
  list(
    model = model,
    scaler = list(mu = mu, sd = sd, features = features),
    cells_used = cells,
    config = list(hidden = hidden, n_bins = n_bins, epochs = epochs, lr = lr, seed = seed)
  )
}

#' Predict with KANN
predict_kann <- function(seu, fit, out_col = "MQI_pred") {
  cells <- colnames(seu)
  expr <- GetAssayData(seu, slot = "data")
  
  missing <- setdiff(fit$scaler$features, rownames(expr))
  if (length(missing) > 0) {
    warning(sprintf("%d features missing", length(missing)))
  }
  
  features <- intersect(fit$scaler$features, rownames(expr))
  X <- t(as.matrix(expr[features, cells]))
  X <- sweep(sweep(X, 2, fit$scaler$mu[features], "-"), 2, fit$scaler$sd[features], "/")
  
  fit$model$eval()
  pred <- as.numeric(fit$model(torch_tensor(X, dtype = torch_float())))
  
  seu[[out_col]] <- NA_real_
  seu[[out_col]][cells] <- pred
  seu
}

# ==============================================================================
# 5) PERMUTATION IMPORTANCE
# ==============================================================================

#' Batched prediction helper
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
  
  # Build X using same cells/features as training
  cells <- fit$cells_used[fit$cells_used %in% colnames(seu)]
  expr <- GetAssayData(seu, slot = "data")
  X <- t(as.matrix(expr[fit$scaler$features, cells]))
  X <- sweep(sweep(X, 2, fit$scaler$mu, "-"), 2, fit$scaler$sd, "/")
  
  y <- seu@meta.data[cells, "MQI_v1"]
  
  fit$model$eval()
  base_pred <- predict_batches(fit$model, X, batch)
  base_mse <- mean((base_pred - y)^2)
  
  p <- ncol(X)
  imp <- numeric(p)
  names(imp) <- colnames(X)
  
  if (verbose) {
    message(sprintf("Computing importance: %d features, %d repeats", p, repeats))
    pb <- txtProgressBar(0, p, style = 3)
  }
  
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
  
  if (verbose) close(pb)
  
  sort(imp, decreasing = TRUE)
}

# ==============================================================================
# 6) SAVE/LOAD
# ==============================================================================

#' Save model bundle
save_kann <- function(fit, dir, overwrite = FALSE) {
  if (dir.exists(dir) && !overwrite) stop("Directory exists. Use overwrite=TRUE")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  torch_save(fit$model$state_dict(), file.path(dir, "model.pt"))
  qsave(fit$scaler, file.path(dir, "scaler.qs"))
  qsave(fit$config, file.path(dir, "config.qs"))
  
  writeLines(c(
    sprintf("KANN-MFI Model"),
    sprintf("Features: %d", length(fit$scaler$features)),
    sprintf("Hidden: %d", fit$config$hidden),
    sprintf("Trained: %s", Sys.time())
  ), file.path(dir, "README.txt"))
  
  message("✓ Saved to: ", dir)
}

#' Load model bundle
load_kann <- function(dir) {
  scaler <- qread(file.path(dir, "scaler.qs"))
  config <- qread(file.path(dir, "config.qs"))
  
  model <- KANNet(
    length(scaler$features),
    hidden = config$hidden,
    out_features = 1,
    n_bins = config$n_bins
  )
  
  state <- torch_load(file.path(dir, "model.pt"))
  model$load_state_dict(state)
  model$eval()
  
  list(model = model, scaler = scaler, config = config)
}

# ==============================================================================
# 7) PLOTTING
# ==============================================================================

#' Quick diagnostic plot
plot_predictions <- function(seu, obs_col = "MQI_v1", pred_col = "MQI_pred") {
  df <- seu@meta.data[, c(obs_col, pred_col)]
  df <- df[complete.cases(df), ]
  
  r2 <- cor(df[[obs_col]], df[[pred_col]])^2
  rmse <- sqrt(mean((df[[obs_col]] - df[[pred_col]])^2))
  
  ggplot(df, aes_string(x = obs_col, y = pred_col)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
            label = sprintf("R² = %.3f\nRMSE = %.3f", r2, rmse)) +
    labs(title = "KANN Predictions", x = "Observed MQI", y = "Predicted MQI") +
    theme_classic()
}

#' Importance barplot
plot_importance <- function(imp, top_n = 30, title = "Permutation Importance") {
  df <- data.frame(
    gene = names(head(imp, top_n)),
    importance = as.numeric(head(imp, top_n))
  )
  df$gene <- factor(df$gene, levels = rev(df$gene))
  
  ggplot(df, aes(x = gene, y = importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = title, x = NULL, y = "ΔMSE") +
    theme_classic()
}
