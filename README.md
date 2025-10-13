# KANN-MFI

Kolmogorov-Arnold Networks for Mitochondrial Fitness Index prediction from single-cell RNA-seq.

<p align="center">
  <img src="logo.png" width="200"/>
</p>

## What it does

1. **Calculates MQI** from metabolic axes (Mitophagy, OXPHOS, UPRmt, etc.)
2. **Trains a KANN** to predict MQI from mitochondrial regulator genes
3. **Permutation importance** to find key regulators
4. **Saves everything** for reproducibility

## Quick Start

```r
# Install
install.packages(c("Seurat", "qs", "dplyr", "torch", "ggplot2"))
torch::install_torch()

# Load
source("kann_mfi.R")

# Your Seurat object needs these columns:
# - Mitophagy, OXPHOS1, percent.mt (required)
# - TCA1, Glyco1 (optional)
# UPRmtScore will be calculated automatically

# 1. Calculate MQI
seu <- calculate_mqi(seu, target_col = "MQI_v1")

# 2. Train KANN
fit <- train_kann(
  seu,
  target_col = "MQI_v1",
  hidden = 48,
  epochs = 300,
  seed = 42
)

# 3. Predict
seu <- predict_kann(seu, fit, out_col = "MQI_pred")

# 4. Evaluate
plot_predictions(seu, obs_col = "MQI_v1", pred_col = "MQI_pred")

# 5. Feature importance
imp <- perm_importance(seu, fit, repeats = 3)
plot_importance(imp, top_n = 30)

# 6. Save for later
save_kann(fit, dir = "my_model")
```

## Load saved model

```r
fit <- load_kann("my_model")
seu <- predict_kann(seu, fit)
```

## File Structure

```
kann-mfi/
├── kann_mfi.R          # All functions (source this!)
├── example.R           # Complete working example
├── README.md           # This file
├── LICENSE             # MIT
└── logo.png            # Hex sticker
```

## What you need in your Seurat object

**Required columns:**
- `Mitophagy` - mitophagy score
- `OXPHOS1` - oxidative phosphorylation score  
- `percent.mt` - mitochondrial fraction

**Optional (improves MQI):**
- `TCA1` - TCA cycle score
- `Glyco1` - glycolysis score

**Auto-calculated:**
- `UPRmtScore` - from 15 UPRmt genes

## How MQI is calculated

```r
MQI = (Mitophagy + OXPHOS + 0.7×UPRmt) - penalties
```

Penalties include high mt%, OXPHOS/glycolysis imbalance, excess TCA.

## Citation

```bibtex
@software{kann_mfi_2025,
  author = {<Your Name>},
  title = {KANN-MFI: Predicting Mitochondrial Fitness with Kolmogorov-Arnold Networks},
  year = {2025},
  url = {https://github.com/<your-username>/kann-mfi}
}
```

## License

MIT - see LICENSE file

## Questions?

Open an issue on GitHub
