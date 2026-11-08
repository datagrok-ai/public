# EDA analyze — what the features do not claim

## Resolved

- public#4473 (MCL on number columns only never finished: the dialog passed no column): the MCL editor builds its
  per-column parameters for number columns too (libraries/ml); `clustering.feature` runs MCL on two number columns.

## Still open

- Reduce Dimensionality greys OK with no column picked by a CSS class only (EDA `src/package.ts`, `reduceDimensionality`
  `validate`), with no `disabled` or `aria-disabled`; claim when fixed: `OK button in "Reduce Dimensionality" dialog
  should be disabled` before a column is picked.

## Won't do

- DBSCAN labels, the t-test numbers, the embeddings: EDA package tests ('DBSCAN', ttest-tests, 'Dimensionality
  reduction: UMAP / T-SNE'); MCL clusters: Peptides 'Model: MCL clusters'.
- Option permutations of one run (t-SNE, Postprocessing None, Plot Embeddings off; Student, Full report off; a t-test
  on overlapping groups): one rule each, the dialog's OK is claimed once.
- The dialogs' defaults and captions, and the reason a disabled t-test Run gives in its tooltip: no test gap.
- Use WebGPU (Reduce Dimensionality, MCL): headless Chromium has no WebGPU adapter.
- core#416, core#418 (Control Comparisons tooltip jump, dropdown flash): cosmetic and closed, no state to read.
