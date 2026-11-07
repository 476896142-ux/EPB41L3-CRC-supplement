# Content intentionally excluded from the public release

| Excluded class | Reason | How to obtain or handle it |
|---|---|---|
| Large raw and processed expression matrices, single-cell objects, CellChat objects, and spatial objects | GitHub is not the source archive; redistribution and repository-size suitability were not verified | Retrieve the public source data; rebuild intermediate objects locally. |
| `eqtl_result.RData` (approximately 278 MB) and other serialized R objects | Large binary artifact; provenance and redistribution terms must be checked before release | Regenerate from documented upstream inputs, or archive separately in Zenodo only after provenance review. |
| CIBERSORT and LM22 resources | Third-party distribution and license conditions | Obtain directly from the official CIBERSORT channel. |
| MCP-counter source/reference files | Third-party package/license handling | Install or obtain from the official MCP-counter distribution. |
| OpenGWAS JWT/API token | Credential | Create and export an individual token locally; never commit it. |
| Manuscript DOCX and duplicate PDF figures | Not needed for source-code/data provision; journal submission system retains the manuscript files | Upload through the journal workflow as required. |

The archival source folder is preserved unchanged outside this repository. Exclusion here is a publication-safety decision, not deletion of the underlying material.
