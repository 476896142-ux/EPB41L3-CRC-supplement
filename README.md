# EPB41L3 colorectal cancer multi-omics supplementary materials

This repository accompanies a study of EPB41L3 and related candidate genes in colorectal cancer using bulk transcriptomics, eQTL/pQTL Mendelian randomization, network analysis, machine learning, immune deconvolution, single-cell and spatial transcriptomics, and exploratory survival analyses.

It contains the analysis scripts archived for the manuscript, tabular results and figure source data, and the 13 manuscript figures. It is a **publication companion repository**, not a turnkey rerun package: the archived scripts retain portions of the original project-specific workflow and require external source data, software environments, reference panels, and access credentials described below.

## Contents

- `code/`: analysis scripts grouped in the manuscript workflow order. A credential embedded in the archived pQTL script has been removed; set `IEUGWASR_TOKEN` locally when that API is required.
- `results/`: derived tables, source data, QC records, model outputs, and statistical results underlying the manuscript.
- `figures/`: the 13 PNG figures extracted from the manuscript and SHA-256 matched to the manuscript copy during the original material audit.
- `data/reference/`: small reference files needed by the GSEA/PPI components.
- `metadata/`: figure mapping, source-file checksums, excluded-data register, and publication checklist.

## Data access and reuse

This repository deliberately does **not** redistribute large raw/processed matrices, R serialized objects, third-party software/signatures, private access tokens, or files whose redistribution terms were not verified. See [metadata/DATA_SOURCES.md](metadata/DATA_SOURCES.md) and [metadata/EXCLUDED_CONTENT.md](metadata/EXCLUDED_CONTENT.md).

Public datasets should be retrieved from their original repositories, rather than copied from this repository. The scripts and the associated result tables identify the public cohorts used, including GEO series GSE196006, GSE146009, GSE110224, GSE44076, GSE20916, GSE132257, GSE132465, and GSE146771, plus TCGA data where indicated.

## Reproducing components

1. Clone the repository and inspect the relevant script and result folder.
2. Retrieve the upstream public data at the source and place it in a local project layout appropriate for the selected script.
3. Install the R/Python dependencies required by that script. Existing software-version tables are retained in relevant `results/` folders where available.
4. For OpenGWAS access, set `IEUGWASR_TOKEN` or `OPENGWAS_JWT` in your local environment. Do not place credentials in source files.
5. Run only after verifying input paths, genome build, reference panels, and software versions. The archived code has not been rerun as part of this repository assembly.

## Important interpretation note

The pQTL manuscript candidate set used an unadjusted IVW P < 0.05 intersected with differential-expression results. The file `results/02_pQTL_MR/MR_IVW_FDR05_DEG_intersection.csv` is a separate full-scan FDR-filtered result and must not be substituted for the manuscript candidate set.

## Release procedure

Before public release, complete [metadata/RELEASE_CHECKLIST.md](metadata/RELEASE_CHECKLIST.md), choose an explicit code license, create a tagged GitHub release, and archive that exact release with Zenodo to obtain a DOI. Replace the placeholders in `CITATION.cff` and the manuscript availability statements with the final repository URL and DOI.
