## Setup

Install the required R libraries using `renv`.

### Option 1 — Automatic (recommended)

Run the following in your terminal:

```
bash
git clone git@github.com:qBioTurin/Insite.git
cd Insite

R -e "install.packages('renv', repos = 'https://cloud.r-project.org'); \
options(repos = c(CRAN = 'https://cloud.r-project.org')); \
renv::restore(prompt = FALSE)"
```

### Option 2 — Manual (RStudio)
1. Clone the repository using your preferred method
2. Open the project in RStudio (or create an .Rproj file in the folder)
3. Run the following in the R console:
```install.packages('renv', repos = 'https://cloud.r-project.org')
options(repos = c(CRAN = 'https://cloud.r-project.org'))
renv::restore(prompt = FALSE)```