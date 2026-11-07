# Data sources and access notes

## Public cohorts named in the archived workflow

| Analysis component | Source datasets/services | Repository handling |
|---|---|---|
| Bulk differential expression | GEO: GSE196006, GSE146009, GSE110224 | Retrieve from GEO; derived DEG tables are provided. |
| Machine learning and validation | GEO: GSE44076, GSE20916 | Retrieve from GEO; processed result tables and annotation mappings are provided. |
| Single-cell localization and CellChat | GEO: GSE132257, GSE132465 | Retrieve deposited data and author annotations from GEO; derived summary tables are provided. |
| Independent single-cell localization | GEO: GSE146771 | Retrieve from GEO; derived summary tables are provided. |
| Spatial expression and RCTD | Source data and spatial images are not redistributed here | Use the original deposition identified in the analysis records; derived QC and weight tables are provided. |
| eQTL/pQTL MR | GTEx-derived eQTL inputs, pQTL inputs, OpenGWAS/VCF outcome resources, LD reference panel | Do not redistribute large association data, VCFs, reference panels, or API credentials. |
| Bulk microenvironment and survival | TCGA and GEO cohorts | Retrieve from original providers; derived association and model tables are provided. |

## Third-party methods and resources

The repository does not include CIBERSORT source code/signature resources or MCP-counter program files. Obtain them from their official sources and comply with their respective license or access terms. STRING-derived interaction tables are retained only as the small analysis input used for the documented network step; cite STRING when reusing them. KEGG pathway files are retained as small historical analysis inputs; users should verify current KEGG reuse terms before redistribution or new analyses.

## Repository data classes

`results/` contains derived outputs and figure source data. It does not represent a replacement for original public datasets. Any sample identifiers, clinical variables, or expression summaries should be used under the terms of the original data repository.
