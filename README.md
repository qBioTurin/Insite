# A New Cancer Progression Model: From synthetic tumors to real data and back

This repository contains the core **Insite** R package, the data analysis pipelines, and the scripts required to reproduce all results and figures presented in our publication.

---

## 🌐 Project Ecosystem

To ensure both user-friendliness and scientific reproducibility, the project is divided into three components:

* **Insite (This Repository):** The core R package, processed data, and scripts for analysis and figure generation.
* **[Insite_Interface](https://github.com/qBioTurin/Insite_Interface):** A Docker-based NextJS application providing a graphical user interface for interactive simulations.
* **[Zenodo Dataset](LASCIA_QUI_IL_LINK_A_ZENODO):** A permanent archive containing the complete raw simulation outputs used in the study.

---

## 📂 Repository Structure

### 1. Data (`/Data`)
This folder contains the processed data required to generate the paper's figures immediately.
* `best_fit/`: Results from model fitting on real clinical cohorts (AML, Breast, Lung, etc.) - **Main Text**.
* `Clonal_Architectures/`: Data and configurations for Experiment 3 - **Main Text**.
* `Experiment1/` & `Experiment2/`: Data regarding synthetic experiments - **Supplementary Material**.

### 2. Analysis & Figures (`/scripts`)
* `PaperAnalysis/`: Scripts that bridge the gap between raw simulation outputs contained in the Zenodo repo and processed data.
* `PaperFigures/`: R scripts dedicated to generating the specific figures (e.g., `Fig3.R`, `Fig5.R`). 
    * *Note:* These scripts are designed for interactive use (e.g., within RStudio). They will render plots directly to your active graphics device rather than saving files to disk.

### 3. CLI Simulation Wrappers (`/scripts` root)
The `.R` files located directly in the `/scripts` folder (e.g., `run_simulation.R`) are command-line wrappers for the package functions. They use `optparse` and can be executed from a terminal with flags.

---

## 🚀 Getting Started

### Option 1: Docker (Recommended)
To ensure a consistent environment with all dependencies pre-installed:
```bash
docker pull qbioturin/insite:latest
docker run -it qbioturin/insite R
```

### Option 2: Local R Installation
```R
# install.packages("devtools")
devtools::install_github("qBioTurin/Insite")
```

## 🛠 Reproducibility Guide
### Level 1: Regenerating Figures (From elaborated data)
To view the figures using the processed data already included in this repository:
1. Open the project in RStudio.
2. Run any script located in scripts/PaperFigures/.
3. The plot will be rendered in the "Plots" pane.

### Level 2: Full Pipeline (From Raw Data)
To verify the analysis by processing the raw simulation outputs yourself:
1. Download the raw dataset from Zenodo.
2. Create a folder named OutputSim in the root directory of this repository.
3. Place the downloaded raw files into OutputSim.
4. Execute the scripts in scripts/PaperAnalysis/. These will read from OutputSim and update the processed files in the Data/ folder.

## 🧪 Running New Simulations (CLI)

Beyond the R package functions, we provide four command-line wrappers in the `/scripts` directory. These allow you to run the full pipeline—from simulation to analysis—directly from the terminal using `Rscript`.

### 1. `run_simulation.R`
Executes a new cancer progression simulation based on a JSON configuration. It saves the state of the tumor at each timestep as a timestamped `.RData` file.

| Option | Description | Default |
| :--- | :--- | :--- |
| `--params` | Path to the `.json` file containing model parameters. | `params.json` |
| `--dir` | Directory where simulation timesteps will be stored. | `raw` |
| `--seed` | Random seed for the simulation. | `Sys.time()` |
| `--Nexp` | Unique simulation ID/number. | `1` |

### 2. `draw_plot.R`
Generates visual representations of the simulated tumor's evolutionary history, saving them as PDFs. It produces both **Muller plots** (clonal dynamics over time) and **Phylogenetic trees**.

| Option | Description | Default |
| :--- | :--- | :--- |
| `--sim_dir` | Path to the folder containing simulation `.RData` files. | `raw/sim1` |
| `--params` | Path to the elaborated `Parameters.RData` file. | `raw/Parameters.RData` |
| `--path_out` | Folder where the output PDF plots will be saved. | `raw` |
| `--depth` | Prevalence threshold (retrieves clones reaching at least $1/10^{depth}$). | `3` |
| `--relative` | Logical. Use relative (`TRUE`) or absolute (`FALSE`) abundance. | `FALSE` |

### 3. `sequence_tumor.R`
Performs the synthetic sequencing procedure described in the paper. It generates a VCF-like output saved as a `.txt` file.

| Option | Description | Default |
| :--- | :--- | :--- |
| `--sim_dir` | Folder containing simulation outputs. | `raw/sim1` |
| `--params` | Path to elaborated `Parameters.RData`. | `raw/Parameters.RData` |
| `--seq_day` | Day of sequencing (use `Inf` for the final state). | `Inf` |
| `--nregions` | Number of regions (1 for Bulk, >1 for Multi-regional). | `1` |
| `--ncells` | Cells per region (use `Inf` for the whole mass). | `1000` |
| `--neighborhood`| Path to pre-computed neighborhood file to save time. | `NULL` |
| `--repeat` | Number of sequencing repetitions to perform. | `1` |
| `--dens` | Path to DP probability density for coverage. | `Data/dens.RData` |
| `--seed` | Fixed seed for reproducible sequencing. | `NULL` |

### 4. `derive_nD_indices.R`
Calculates the **Clonal Nesting ($n$)** and **Clonal Diversity ($D$)** indices as defined by *Noble et al. (2022)*. Results are exported to a `.txt` file.

| Option | Description | Default |
| :--- | :--- | :--- |
| `--sim_dir` | Folder containing simulation outputs. | `raw/sim1` |
| `--path_out` | Folder where the index results will be stored. | `raw` |
| `--seq_day` | Simulation day to analyze. | `Inf` |

---

## 📄 Citation
If you use Insite in your research, please cite:
> *Daniela Volpatto et al. (2026). A new cancer progression model: from synthetic tumors to real data and back.*

---
