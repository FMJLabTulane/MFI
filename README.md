# KANN-MFI

**Predict Mitochondrial Fitness from single-cell RNA-seq** using Kolmogorov-Arnold Neural Networks

<p align="center">
  <img src="logo_final.png" width="200"/>
</p>

## What it does

**You give it:** A Seurat object with normalized counts  
**It gives you:** Mitochondrial Function Index (MFI) for every cell

### The pipeline:
1. **Auto-calculates metabolic scores** (Mitophagy, OXPHOS, TCA, Glycolysis, UPRmt, mt%)
2. **Computes MFI** from those scores
3. **Trains a neural network** to predict MQI from ~200 mitochondrial genes
4. **Finds important genes** that drive mitochondrial quality

## Quick Start (Complete Noob)

### Install (one time)
```r
install.packages(c("Seurat", "qs", "dplyr", "torch", "ggplot2"))
torch::install_torch()  # Downloads PyTorch (~500MB)
```

### Run (every time)
```r
library(Seurat)
source("kann_mfi.R")

# Load your data
seu <- readRDS("my_seurat_object.rds")

# Step 1: Calculate MQI (auto-calculates everything!)
seu <- calculate_mqi(seu)

# Step 2: Train neural network
fit <- train_kann(seu, epochs = 300)

# Step 3: Predict for all cells
seu <- predict_kann(seu, fit)

# Step 4: Plot
plot_predictions(seu)

# Step 5: Find important genes
imp <- perm_importance(seu, fit)
plot_importance(imp)

# Step 6: Save everything
save_kann(fit, dir = "my_model")
saveRDS(seu, "seurat_with_mqi.rds")
```

**That's it!** Check out `example.R` for a complete annotated walkthrough.

---

## Requirements

### What you need:
- A Seurat object with **normalized counts** (run `NormalizeData()` if you haven't)
- That's it!

### What gets calculated automatically:
- Mitophagy score (27 genes)
- OXPHOS score (Complex I-V genes)
- TCA cycle score
- Glycolysis score
- UPRmt score (15 genes)
- Mitochondrial % (from MT- genes)
- MQI (combines all of the above)

**All you need is raw RNA counts!** 🎉

---

## What is MFI?

**Mitochondrial Function Index** = a single number that summarizes mitochondrial health

```
MQI = (Mitophagy + OXPHOS + UPRmt) - penalties
```

**Higher MFI = healthier mitochondria**

Penalties for:
- High mitochondrial DNA %
- OXPHOS/Glycolysis imbalance  
- Excess TCA activity

---

## 🔬 Output

### New columns in your Seurat object:
- `Mitophagy` - mitophagy activity
- `OXPHOS1` - oxidative phosphorylation
- `TCA1` - TCA cycle activity
- `Glyco1` - glycolysis
- `UPRmtScore` - mitochondrial unfolded protein response
- `percent.mt` - mitochondrial fraction
- **`MQI_v1`** - calculated mitochondrial quality index ⭐
- **`MQI_pred`** - KANN prediction of MQI ⭐

### Files created:
- `pred_vs_obs.pdf` - How well the model works
- `importance.pdf` - Top genes that determine MQI
- `importance_all_genes.csv` - Complete gene rankings
- `my_model/` - Saved model (reusable!)

---

## Examples

### Basic usage
```r
source("kann_mfi.R")
seu <- readRDS("data.rds")

# Everything in 4 lines:
seu <- calculate_mqi(seu)
fit <- train_kann(seu, epochs = 300)
seu <- predict_kann(seu, fit)
plot_predictions(seu)
```

### Reuse saved model
```r
source("kann_mfi.R")

# Load model
fit <- load_kann("my_model")

# Apply to new data
new_seu <- readRDS("new_data.rds")
new_seu <- calculate_mqi(new_seu)
new_seu <- predict_kann(new_seu, fit)
```

### Compare cell types
```r
# After running the pipeline
library(ggplot2)

# MQI by cell type
ggplot(seu@meta.data, aes(x = cell_type, y = MQI_v1)) +
  geom_boxplot() +
  theme_classic()

# Top genes for each cell type
for (ct in unique(seu$cell_type)) {
  seu_sub <- subset(seu, cell_type == ct)
  imp <- perm_importance(seu_sub, fit, repeats = 3)
  print(head(imp, 10))
}
```

---

## 📊 Understanding the Results

### Good performance:
- **R² > 0.7** - Model captures most MQI variation
- **RMSE < 0.5** - Predictions are close to observed

### Importance scores (ΔMSE):
- **> 0.01** - Gene strongly affects MQI
- **0.001-0.01** - Moderate effect
- **< 0.001** - Weak effect

### Top genes usually include:
- Mitophagy receptors (BNIP3, FUNDC1, PINK1, PRKN)
- OXPHOS genes (Complex I, Complex V)
- Quality control (LONP1, CLPP, HSPD1)

---

## Files in this repo

```
kann-mfi/
├── kann_mfi.R          # 🔴 Main file - source this!
├── example.R           # Complete working example
├── README.md           # This file
├── LICENSE             # MIT License
└── logo.png            # Hex sticker
```

**Just download `kann_mfi.R` and you're good to go!**

---

## ❓ FAQ

**Q: I'm so lost what do I do?  
A: Create an issue on the issue tab :)

**Q: I only have a count matrix, not a Seurat object**  
A: 
```r
seu <- CreateSeuratObject(counts = your_matrix)
seu <- NormalizeData(seu)
# Now use the pipeline!
```

**Q: Some genes are missing from my dataset**  
A: No problem! The pipeline uses whatever genes are available and tells you the coverage %.

**Q: Can I use my own gene sets?**  
A: Yes! Pass your own genes to `train_kann()`:
```r
my_genes <- c("PINK1", "PRKN", "MFN1", "MFN2", ...)
fit <- train_kann(seu, features = my_genes)
```

**Q: How long does training take?**  
A: ~2-5 minutes for 10,000 cells on a laptop. Permutation importance takes longer (10-30 min).

**Q: Can I use this on mouse data?**  
A: Yes! Gene names should be in the same format as your Seurat object (e.g., `Ppargc1a` for mouse).

**Q: What if I already have metabolic scores?**  
A: The pipeline checks for existing scores and skips them. Use `force_recalc=TRUE` to recalculate:
```r
seu <- calculate_mqi(seu, force_recalc = TRUE)
```

---

## Citation

```bibtex
@software{kann_mfi_2025,
  author = {<Fahd Qadir aka Dragonmasterx87>},
  title = {KANN-MFI: Kolmogorov-Arnold Networks for Mitochondrial Fitness Prediction},
  year = {2025},
  url = {https://github.com/FMJLabTulane/kann-mfi}
}
```

---

## 📄 License

MIT License - Free to use, modify, and share!

---

## 🆘 Help

- **Issues?** Open an issue on GitHub
- **Questions?** Check `example.R` for detailed walkthrough
- **Email:** mqadir@tulane.edu

---

## 🙏 Acknowledgments

Based on:
- Kolmogorov-Arnold representation theorem
- Liu et al. (2024) "KAN: Kolmogorov-Arnold Networks" arXiv:2404.19756
- Seurat: Hao et al. (2021) Cell
