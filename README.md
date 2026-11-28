# Y1 Temporal Proteomics and Phosphoproteomics

Code and analysis workflows accompanying the study:

*A comprehensive temporal proteome and phosphoproteome dataset of Y1 adrenocortical cells under FGF2 stimulation and SOS inhibition.*

The repository contains the preprocessing, batch-correction, dimensionality-reduction, and differential-analysis workflows. The original proteomics data are available from the jPOST Repository and are not included in this repository.

## Data availability

Download the datasets from jPOST:

- **Proteome:** [JPST004922](https://repository.jpostdb.org/entry/JPST004922) / ProteomeXchange **PXD084835**
- **Phosphoproteome:** [JPST004842](https://repository.jpostdb.org/entry/JPST004842) / ProteomeXchange **PXD084379**

After downloading and preparing the corresponding input files, place them at the following paths, using the filenames expected by the preprocessing notebook:

```text
Combat_CoVariates/
├── Normal_Proteome/
│   └── input_files/
│       └── result.pg_matrix.tsv
└── Phosphoproteome/
    └── input_files/
        ├── Phospho (STY)Sites_PIF_0.5_Rearranged.txt
        └── batch_info_wReferences.csv
```

The `batch_info_wReferences.csv` file is provided in this repository and should be retained in the indicated location. The `model.xlsx` file in `outer_functions/DimRed_data/` is also required for the PCA and UMAP notebook.

## Workflow

Run the notebooks and scripts in the following order.

### 1. Preprocessing and batch correction

Run `Combat_CoVariates/01_preprocessing_normalizing_bpca_combat.ipynb`.

This notebook preprocesses the phosphoproteome and normal proteome datasets. The phosphoproteome workflow includes filtering, log2 transformation, BPCA imputation, ComBat batch correction, and mean-centering. The normal proteome workflow includes preprocessing, BPCA imputation, and mean-centering.

The notebook writes its processed files to the `output_files/` directories under `Phosphoproteome/` and `Normal_Proteome/`.

The BPCA implementation used in this workflow is based on Kitamoto et al. (2023): https://github.com/Idesan/bpca. The BPCA-related folders under `outer_functions/` are retained for this workflow.

### 2. PCA and UMAP

After preprocessing, run `Combat_CoVariates/02_PCA_UMAP.ipynb`.

This notebook performs PCA and UMAP analyses on the processed proteome and phosphoproteome data. It uses `outer_functions/DimRed_data/model.xlsx` for sample metadata and grouping.

### 3. Differential analysis with limma

Run the R scripts separately:

- `Limma_NormalProteome/limma.R`
- `Limma_PhosphoProteome/limma.R`

The scripts expect the processed input tables at these paths:

- Normal proteome: `Limma_NormalProteome/Data/03_NormalProteome_result.pg_matrix_BPCA_Centered.tsv`
- Phosphoproteome: `Limma_PhosphoProteome/Data/05_Phospho (STY)Sites_PIF_0.5_DownstreamReady.tsv`

These tables are generated during preprocessing. Copy the corresponding generated files from the preprocessing `output_files/` directories into the `Data/` directories using the filenames listed above before running the R scripts. The `Data/` directories do not need to contain these generated tables in the published repository.

The R scripts save their results under the corresponding `Outputs/` directories.

## Dependencies

### Python

Install the packages listed in `Combat_CoVariates/requirements.txt` from the `Combat_CoVariates/` directory:

```bash
pip install -r requirements.txt
```

The preprocessing notebook also imports `pycombat` using `from combat.pycombat import pycombat`

### R

The limma scripts require `limma` and `openxlsx`. Install them in R with:

```r
install.packages("BiocManager")
BiocManager::install("limma")
install.packages("openxlsx")
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24629325.svg)](https://doi.org/10.5281/zenodo.24629325)
