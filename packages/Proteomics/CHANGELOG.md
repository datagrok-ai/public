# Proteomics changelog

## 0.1.0

Initial release — mass spectrometry-based proteomics analysis, carrying a study from
search-engine output to biological interpretation inside Datagrok.

- **Import** — MaxQuant (`proteinGroups.txt`), Spectronaut Report (PG-level, including
  streaming aggregation of long-form precursor exports), Spectronaut Candidates
  (pre-computed differential expression), FragPipe (`combined_protein.tsv`), and a generic
  CSV/TSV column-mapping importer. Contaminant/decoy filtering, log2 transformation, and
  semantic-type tagging are applied on import.
- **Pipeline** — experiment annotation (group assignment), normalization
  (median / quantile / VSN), missing-value imputation (kNN / MinProb / mean / median /
  zero), and differential expression cascading DEqMS → limma → client-side Welch's t-test.
- **Visualization** — volcano, clustered expression heatmap, sample-level PCA, group-mean
  correlation, QC dashboard (MA / CV / missingness / intensity trend), and an SPC dashboard
  with Nelson-rule instrument monitoring.
- **Biological interpretation** — g:Profiler over-representation analysis (GO, KEGG,
  Reactome, WikiPathways across 9 organisms) with dot/bar charts cross-linked to the
  volcano, plus a UniProt context panel for per-protein metadata.
- **Sharing** — publish read-only, access-verified, versioned review snapshots
  (table + volcano + enrichment) to a reviewer group.
- R scripts (limma, DEqMS, VSN) run server-side with client-side fallbacks throughout, so
  the package is fully functional without a configured R environment.
