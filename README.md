# Insite: A New Cancer Progression Model
### *From synthetic tumors to real data and back*

This repository contains the core **Insite** R package, the data analysis pipelines, and the scripts required to reproduce all results and figures presented in our publication.

---

## 🌐 Project Ecosystem

To ensure both user-friendliness and scientific reproducibility, the project is divided into three components:

* **Insite (This Repository):** The core R package, processed data, and scripts for analysis and figure generation.
* **[Insite_Interface](LASCIA_QUI_IL_LINK_A_INSITE_INTERFACE):** A Docker-based Shiny application providing a graphical user interface for interactive simulations.
* **[Zenodo Dataset](LASCIA_QUI_IL_LINK_A_ZENODO):** A permanent archive containing the complete raw simulation outputs used in the study.

---

## 📂 Repository Structure

### 1. Data (`/Data`)
This folder contains the processed data required to generate the paper's figures immediately.
* `best_fit/`: Results from model fitting on real clinical cohorts (AML, Breast, Lung, etc.) - **Main Text**.
* `Clonal_Architectures/`: Data and configurations for Experiment 3 - **Main Text**.
* `Experiment1/` & `Experiment2/`: Data regarding synthetic experiments - **Supplementary Material**.

### 2. Analysis & Figures (`/scripts`)
* `PaperAnalysis/`: Scripts that bridge the gap between raw simulation outputs and processed data. These scripts perform the statistical processing that justifies our findings.
* `PaperFigures/`: R scripts dedicated to generating the specific figures (e.g., `Fig3.R`, `Fig5.R`). 
    * *Note:* These scripts are designed for interactive use (e.g., within RStudio). They will render plots directly to your active graphics device rather than saving files to disk.

### 3. CLI Simulation Wrappers (`/scripts` root)
The `.R` files located directly in the `/scripts` folder (e.g., `run_simulation.R`) are command-line wrappers for the package functions. They use `optparse` and can be executed from a terminal with flags.

---

## 🚀 Getting Started

### Option 1: Docker (Recommended)
To ensure a consistent environment with all dependencies pre-installed:
```bash
docker pull tuo-username/insite:latest
docker run -it tuo-username/insite R
