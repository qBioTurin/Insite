# A New Cancer Progression Model: From synthetic tumors to real data and back

This repository contains the core **Insite** R package, the data analysis pipelines, and the scripts required to reproduce all results and figures presented in our publication.

---

## 🌐 Project Ecosystem

To ensure both user-friendliness and scientific reproducibility, the project is divided into three components:

* **Insite (This Repository):** The core R package, processed data, and scripts for analysis and figure generation.
* **[Insite_Interface](https://github.com/qBioTurin/Insite_Interface):** A Docker-based NextJS application providing a graphical user interface for interactive simulations.
* **[Zenodo Dataset](10.5281/zenodo.19821839):** A permanent archive containing the complete raw simulation outputs used in the study.

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

Clone the Github repository via

```
git clone https://github.com/qBioTurin/Insite.git
```

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
To recreate the paper figures using the processed data included in this repository:
1. Open the project in RStudio and ensure the working directory is set to the repository root.
2. Run any script in `scripts/PaperFigures/`. The plots will appear in the *Plots* pane.

### Level 2: Full Pipeline (From Raw Data)
To verify the analysis by processing the raw simulation outputs yourself:
1. Download the raw dataset from Zenodo (a compressed folder named OutputSim).
2. Unzip it and place the `OutputSim` folder in the root directory of this repository. Do not rename it, otherwise the scripts will not work.
3. Run the scripts in `scripts/PaperAnalysis/`. These will read from `OutputSim` and update the processed files in the `Data/` folder.

## 🧪 Running New Simulations

The R package allows you to run new simulations and analyze the results.

In addition to the package functions, four command-line wrappers are provided in the `scripts/` directory. These let you run the full pipeline—from simulation to analysis—directly from the terminal using `Rscript`.
For testing and as a reference for how parameters should be passed to the simulator, an example `params.json` file is included in the repository root. By default, the command-line wrappers use this file as input and create a `raw/` directory where simulation outputs and all subsequent processing steps are stored. Both the input file and output directory can be changed.

In detail:

### 1. `run_simulation.R`
Executes a new cancer simulation based on parameters specified in a `.json` file. At each timestep, the tumor state is saved as a timestamped `.RData` file in a folder named `simNexp`, where `Nexp` is provided by the user. This folder is created in the specified output directory.
The design supports running multiple simulations and facilitates parallelization over the `Nexp` argument.
A seed can be specified to reproduce results. The effective seed is computed as `seed + Nexp`, allowing multiple runs to be reproducible with minimal information. Even if no seed is provided, the script records the seed used for the simulation in a `seed.txt` file.


| Option | Description | Default |
| :--- | :--- | :--- |
| `--params` | Path to the `.json` file containing model parameters. | `params.json` |
| `--dir` | Output directory where simulation results (e.g. `simNexp/`) will be stored. | `raw` |
| `--seed` | Base random seed. The effective seed is computed as seed + Nexp. If not provided, a seed is generated and saved to `seed.txt`. | `Sys.time()` |
| `--Nexp` | Simulation index. Used to name the output folder (`simNexp`) and to derive the effective seed. | `1` |

Along with the simulation `.RData` files, the script will produce an auxiliary file `Parameters.RData` and place it in the directory indicated (`raw` by default).

### 2. `draw_plot.R`
Generates visual representations of the simulated tumor's evolutionary history, saving them as PDFs images. It produces both **Muller plots** (clonal dynamics over time) and **Phylogenetic trees**.

| Option | Description | Default |
| :--- | :--- | :--- |
| `--sim_dir` | Path to the folder containing simulation `.RData` files. | `raw/sim1` |
| `--params` | Path to the auxiliary `Parameters.RData` file. | `raw/Parameters.RData` |
| `--path_out` | Folder where the output PDF plots will be saved. | `raw` |
| `--depth` | Prevalence threshold (plots retrieves clones reaching at least $1/10^{depth}$). | `3` |
| `--relative` | Logical. Use relative (`TRUE`) or absolute (`FALSE`) abundance for the mullerplot visualization. | `FALSE` |

### 3. `sequence_tumor.R`
Performs the synthetic sequencing procedure described in the paper at a given time point (`seq_day`).

This procedure exploits the neighborhood structure of clones based on their parental relationships: this is the most computationally expensive step. To avoid recomputation, once executed for a given timestep the script saves the result as a `Clones_ordered_seq_day.RData` file. This file can be reused in subsequent runs by providing its path via the `--neighborhood` argument, otherwise optional.

An additional optional argument, `--dens`, allows overriding the default probability distribution used to sample per-base coverage (DP). The function expects an object of class `Density`, which can be stored in a `.RData` file and passed to the script. A default density derived from TCGA pancancer data is included in the package. The script used to generate it is available in `scripts/PaperAnalysis/` as `density_coverage_TCGA.R`.

The function generates one or more VCF-like outputs (as specified by `--repeat`), each saved as a `.txt` file.

| Option           | Description                                                        | Default                |
| :--------------- | :----------------------------------------------------------------- | :--------------------- |
| `--sim_dir`      | Directory containing simulation outputs.                           | `raw/sim1`             |
| `--params`       | Path to auxiliary `Parameters.RData` file.                         | `raw/Parameters.RData` |
| `--seq_day`      | Time point for sequencing (use `Inf` for the final state).         | `Inf`                  |
| `--nregions`     | Number of sampled regions (1 = bulk, >1 = multi-regional).         | `1`                    |
| `--ncells`       | Number of cells per region (use `Inf` for the full tumor mass).    | `1000`                 |
| `--neighborhood` | Path to a precomputed neighborhood file to skip recomputation.     | `NULL`                 |
| `--repeat`       | Number of sequencing replicates (one VCF-like file per replicate). | `1`                    |
| `--dens`         | Path to the DP coverage density (`.RData` file).                   | `Data/dens.RData`      |
| `--seed`         | Seed for reproducible sequencing.                                  | `NULL`                 |

### 4. `derive_nD_indices.R`
Calculates the **Clonal Nesting ($n$)** and **Clonal Diversity ($D$)** indices as defined by *Noble et al. (2022)*. Results are exported to a `.txt` file.

| Option | Description | Default |
| :--- | :--- | :--- |
| `--sim_dir` | Folder containing simulation outputs. | `raw/sim1` |
| `--path_out` | Folder where the index results will be stored. | `raw` |
| `--seq_day` | Simulation day to analyze (use `Inf` for the final state). | `Inf` |

---

## 📄 Citation
If you use Insite in your research, please cite:
> *Daniela Volpatto et al. (2026). A new cancer progression model: from synthetic tumors to real data and back.*

---
