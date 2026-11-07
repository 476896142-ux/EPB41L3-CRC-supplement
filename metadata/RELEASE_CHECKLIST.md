# Public-release checklist

- [ ] Confirm the manuscript title, author list, corresponding author, and final repository name.
- [ ] Review every tracked file for personal data, access tokens, emails, local credentials, and unpublished data.
- [ ] Confirm that every public data source and third-party resource is redistributed lawfully; otherwise retain only its accession/link and acquisition script.
- [ ] Check that the stated eQTL/pQTL instrument selection, outcome source, genome build, and LD reference panel match the final manuscript.
- [ ] Verify each figure and supplementary table against the submitted manuscript version.
- [ ] Run at least the documented reproducible components in a clean environment, or state clearly which steps remain archival rather than rerunnable.
- [ ] Choose and add an explicit code license only after confirming author/institutional ownership.
- [ ] Replace the placeholders in `CITATION.cff` with the final title, authors, repository URL, and version.
- [ ] Create a GitHub release tag (for example, `v1.0.0`) corresponding exactly to the submitted or accepted manuscript.
- [ ] Archive that tag in Zenodo and insert the returned DOI into the manuscript and `CITATION.cff`.
