# Proteomics changelog

## 1.3.0 (2026-07-13)

- **Rank–Abundance (dynamic-range) viewer** — visualize the abundance dynamic range across the
  experiment.
- **Abundance Correlation viewer** — the separate sample-correlation views are consolidated into
  one viewer.
- **Share analyses by project** — publish a read-only analysis snapshot to a reviewer team under a
  controlled, admin-maintained project vocabulary. The Share dialog pre-selects a default reviewer
  team, and administrators add/edit/remove project names from package settings.
- **Enrichment ↔ volcano linked selection** — selecting one or more enrichment terms highlights the
  union of their member proteins on the volcano.
- **Numerical robustness** — imputation and quantile normalization now guard against NaN/Infinity on
  edge-case inputs (single-value columns, degenerate random draws).
- **Package icon** plus expanded documentation (publishing design rationale, the administrator role,
  and projects/teams/package settings).

## 1.2.0

- **Shared enrichment charts now survive reopen** — a published/shared analysis's enrichment
  dot and bar charts no longer lose their value column (the "Column negLog10FDR does not
  exist" error) when the project is reopened, and the Up/Down term tables no longer leak
  into the project as extra tables. The charts are now self-contained on the enrichment
  frame and split by direction via a per-viewer filter.
- **Volcano opens automatically on Candidates import** — importing a Spectronaut Candidates
  file now lands you on the volcano immediately, matching the Report → Differential
  Expression path (Candidates arrive with the result already computed, so the pipeline's DE
  step is skipped).
- **Personas & capabilities documentation** — a reference describing what the proteomics
  analyst vs. the biology scientist can do and how that boundary is enforced.

## 0.2.0

- **Organism-aware subcellular-location coloring** — the volcano's "color by subcellular
  location" now fetches UniProt annotations for the experiment's organism (persisted as a
  `proteomics.organism` tag from the Enrichment dialog), so non-human data is no longer
  mis-colored by human entries. The taxonomy filter applies only to the reviewed-by-gene
  fallback; accession lookups stay unfiltered, and all 9 supported organisms are covered.
- **Internal consistency** — significance thresholds, the direction palette, and column-name
  constants are now single-sourced instead of duplicated across files; the differential-
  expression FDR cutoff is wired through to the BH correction.
- **README** — illustrated with the main analysis view, an authored pipeline diagram, and an
  enrichment→volcano cross-link walkthrough.

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
