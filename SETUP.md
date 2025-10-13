# Setup Guide

## 📥 Step-by-Step Installation

### 1. Clone or Download Repository

**Option A: Using Git**
```bash
git clone https://github.com/YOUR-USERNAME/kann-mfi.git
cd kann-mfi
```

**Option B: Download ZIP**
1. Click the green "Code" button on GitHub
2. Select "Download ZIP"
3. Extract to your desired location
4. Open the folder in your terminal/RStudio

### 2. Install R Dependencies

Open R or RStudio and run:

```r
# Install CRAN packages
install.packages(c("Seurat", "dplyr", "digest", "jsonlite", "ggplot2"))

# Install torch (special installation)
install.packages("torch")
torch::install_torch()  # Downloads PyTorch libraries (~500MB)
```

**Verify torch installation:**
```r
library(torch)
torch_is_installed()  # Should return TRUE
```

### 3. Load the KANN Functions

```r
# Set working directory to the repository folder
setwd("path/to/kann-mfi")

# Load all functions
source("kann_from_seurat.R")
```

### 4. Test Installation

Run the reproducibility test suite:

```r
source("test_reproducibility.R")
```

You should see:
```
✓ ALL TESTS PASSED
Your installation is working correctly and produces
reproducible results.
```

---

## 🚀 Quick Start

### Example 1: Basic Training

```r
library(Seurat)
source("kann_from_seurat.R")

# Load your Seurat object
seu <- readRDS("your_seurat_object.rds")

# Train model
fit <- train_kann_on_seurat(
  seu = seu,
  target_col = "MQI_v1",  # Your target column name
  seed = 42
)

# Predict
seu <- predict_kann_to_seurat(seu, fit)

# Check performance
idx <- which(is.finite(seu$MQI_v1) & is.finite(seu$KANN_pred))
r2 <- cor(seu$MQI_v1[idx], seu$KANN_pred[idx])^2
cat(sprintf("R² = %.3f\n", r2))
```

### Example 2: Complete Workflow

```r
# Run the complete example
source("example_workflow.R")
```

Edit `example_workflow.R` to point to your data file.

---

## 📁 Repository Structure

```
kann-mfi/
├── kann_from_seurat.R       # Main implementation (source this)
├── example_workflow.R        # Complete analysis example
├── test_reproducibility.R   # Validation suite
├── README.md                 # Documentation
├── SETUP.md                  # This file
├── LICENSE                   # MIT License
└── .gitignore               # Git ignore rules
```

---

## 🔧 Troubleshooting

### Torch Installation Issues

**Error: "torch is not installed"**
```r
# Reinstall torch
remove.packages("torch")
install.packages("torch")
torch::install_torch()
```

**On Windows:** May need [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

**On Mac with Apple Silicon (M1/M2):**
```r
# Use CPU version (GPU not needed for this)
torch::install_torch(type = "cpu")
```

**On Linux:**
```r
# If you have CUDA GPU
torch::install_torch(type = "cuda")

# Otherwise
torch::install_torch(type = "cpu")
```

### Seurat Issues

**Error: "Can't find Seurat"**
```r
# Install from CRAN
install.packages("Seurat")

# Or latest version
remotes::install_github("satijalab/seurat", ref = "master")
```

### Memory Issues

For large datasets (>50k cells):
```r
fit <- train_kann_on_seurat(
  seu = seu,
  target_col = "MQI_v1",
  batch_size = 256,  # Reduce batch size
  epochs = 200       # May need fewer epochs
)
```

---

## 🎯 What Each File Does

| File | Purpose |
|------|---------|
| `kann_from_seurat.R` | Core implementation - source this to load all functions |
| `example_workflow.R` | Complete example from training to evaluation |
| `test_reproducibility.R` | Verifies installation works correctly |
| `README.md` | Full documentation with examples |
| `SETUP.md` | This file - installation guide |

---

## ✅ Next Steps

Once setup is complete:

1. **Read the README**: Full documentation at `README.md`
2. **Run example**: Edit and run `example_workflow.R`
3. **Adapt to your data**: Replace file paths and column names
4. **Check provenance**: Every model saves complete metadata

---

## 🆘 Getting Help

- **Issues**: Open an issue on GitHub
- **Questions**: Check README.md troubleshooting section
- **Email**: your.email@institution.edu

---

## 📊 System Requirements

- **R**: Version 4.0 or higher
- **RAM**: 8GB minimum (16GB+ for large datasets)
- **Disk**: ~2GB for torch + dependencies
- **OS**: Windows, macOS, or Linux
