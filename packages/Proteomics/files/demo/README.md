# Demo Datasets

Sample search-engine outputs that exercise every import path and the full analysis
pipeline (annotate → normalize → impute → differential expression → enrichment →
visualize), plus the fixtures the package's tests pin against. They are browsable in the
platform under **Files | App Data | Proteomics | demo**.

Most carry an engineered or known-ground-truth signal — spike-in standards or species
mixed at fixed ratios between conditions — so the differential-expression and volcano
results are checkable against the truth, not just plausible.

## Which dataset should I use?

| Dataset | Vendor / format | Acquisition | Best for |
|---------|-----------------|-------------|----------|
| `cptac-spike-in.txt` | MaxQuant `proteinGroups.txt` | DDA | The full pipeline on a real benchmark with known UPS spike-in truth |
| `proteinGroups.txt` | MaxQuant `proteinGroups.txt` | DDA | A fast, tiny MaxQuant import / UI smoke test |
| `fragpipe-smoke-test.tsv` | FragPipe `combined_protein.tsv` | DDA | The FragPipe path + UniProt panel (real accessions) with an engineered volcano |
| `spectronaut-hye-mix.tsv` | Spectronaut report (long) | DIA | Long-format DIA import and the three-species (HYE) ground-truth story |
| `spectronaut-hye-demo.tsv` | Spectronaut report (long, slim) | DIA | **The registered `Proteomics Demo`** — the full pipeline's richest endpoint, run end to end |
| `enrichment-demo.csv` | Engineered human DE result | — | **The registered `Proteomics Enrichment Demo`** — a single-organism pathway-enrichment story |
| `spectronaut-hye-candidates.tsv` | Spectronaut Candidates | DIA | The pre-computed-DE shortcut (straight to volcano/enrichment) |
| `spectronaut-hye-precursor*.{tsv,json}` | Spectronaut precursor / aggregated | DIA | Streaming-aggregation import and its golden-file equivalence tests |

The two HYE files (`spectronaut-hye-mix.tsv` and `spectronaut-hye-candidates.tsv`) tell the
same biological story from two angles — per-sample quantities vs. pre-computed
fold-change — so they pair well for comparing the package's own DE against a vendor-shaped
result.

## proteinGroups.txt

| Property | Value |
|----------|-------|
| **Format** | MaxQuant `proteinGroups.txt` (wide, tab-separated) |
| **Source** | Synthetic — generated for development and testing |
| **Organism** | Human |
| **Proteins** | 27 |
| **Samples** | 6 (Sample1–Sample6) |
| **Quantification** | LFQ intensity |
| **Columns** | 17 (Protein IDs, Gene names, Protein names, filter flags, LFQ intensities, coverage, peptide counts) |
| **Size** | 3.7 KB |

Minimal MaxQuant output suitable for quick import testing and UI development. Includes
`Potential contaminant`, `Reverse`, and `Only identified by site` filter columns.

## cptac-spike-in.txt

| Property | Value |
|----------|-------|
| **Format** | MaxQuant `proteinGroups.txt` (wide, tab-separated) |
| **Source** | CPTAC spike-in study — processed through MaxQuant |
| **Organism** | Human (with UPS spike-in standards) |
| **Proteins** | 1,569 |
| **Samples** | 6 (two conditions: 6A and 6B, 3 replicates each) |
| **Quantification** | LFQ intensity and raw Intensity |
| **Columns** | 80 (full MaxQuant output including peptide counts, sequence coverage, identification types, MS/MS counts) |
| **Size** | 4.6 MB |

Full-scale MaxQuant output from the CPTAC spike-in benchmark. Two conditions (6A vs 6B) with
known differentially abundant UPS proteins spiked into a constant yeast background. Ideal for
testing the complete pipeline: import, normalization, imputation, and differential expression.

## fragpipe-smoke-test.tsv

| Property | Value |
|----------|-------|
| **Format** | FragPipe `combined_protein.tsv` (wide, tab-separated) |
| **Source** | Synthetic — generated to exercise the FragPipe parser |
| **Organism** | Human (plus one bovine contaminant) |
| **Proteins** | 10 (8 targets + 1 `contam_` + 1 `rev_` decoy) |
| **Samples** | 4 (Ctrl_1, Ctrl_2, Treat_1, Treat_2 — two conditions × two replicates) |
| **Quantification** | MaxLFQ Intensity, Intensity, Razor Intensity (all three per sample) |
| **Columns** | 19 (Protein, Protein ID, Entry Name, Gene, Length, Organism, Description + 3 × 4 sample intensity columns) |
| **Size** | 2.1 KB |

Minimal FragPipe output for smoke-testing import and the full analysis pipeline. Real UniProt
accessions (TP53, BRCA1, EGFR, AKT1, MTOR, F2, PRKCD, HRAS) so the UniProt context panel
resolves. Includes a `contam_BOVIN_CASEIN` row and a `rev_sp|P12345|...` decoy row to exercise
the parser's contaminant/decoy filter. Differential abundance is engineered: TP53/AKT1/F2 up
in Treat, BRCA1/PRKCD down, EGFR/MTOR/HRAS unchanged — gives a predictable volcano-plot pattern.

## spectronaut-hye-mix.tsv

| Property | Value |
|----------|-------|
| **Format** | Spectronaut report (long, tab-separated) |
| **Source** | [SpectroPipeR](https://github.com/stemicha/SpectroPipeR) — HYE species mix benchmark (`inst/extdata/SN_test_HYE_mix_file.tsv`) |
| **Organism** | Human / Yeast / *E. coli* three-species mix |
| **Proteins** | 93 protein groups |
| **Samples** | 8 runs (two conditions: HYE mix A and HYE mix B, 4 replicates each) |
| **Quantification** | PG.Quantity (protein group), FG.MS2Quantity (fragment), PG.IBAQ |
| **Columns** | 35 (R.FileName, R.Condition, R.Replicate, PG.ProteinGroups, PG.Organisms, EG.Qvalue, EG.ModifiedPeptide, FG.Charge, and more) |
| **Instrument** | timsTOF Pro2 |
| **Size** | 5.0 MB |

Native Spectronaut DIA-MS long-format report from a three-species spike-in experiment. Each row
is a peptide precursor measurement across runs. Human, yeast, and *E. coli* proteins are mixed at
known ratios between conditions A and B, providing ground-truth differential abundance. Useful for
testing long-format import and DIA-specific workflows.

## spectronaut-hye-demo.tsv

| Property | Value |
|----------|-------|
| **Format** | Spectronaut report (long, tab-separated) — slimmed to the 7 columns the parser reads |
| **Source** | **Derived** from `spectronaut-hye-mix.tsv`: kept `R.FileName, R.Condition, R.Replicate, PG.ProteinGroups, PG.Organisms, PG.IBAQ, EG.Qvalue`, and scaled `PG.IBAQ` ×1000 into realistic intensity magnitudes. |
| **Organism** | Human / Yeast / *E. coli* three-species mix |
| **Proteins** | 93 protein groups |
| **Samples** | 8 runs (HYE mix A and HYE mix B, 4 replicates each) |
| **Quantification** | PG.IBAQ (scaled to realistic magnitude) |
| **Columns** | 7 |
| **Signal** | client-side Welch t-test + BH FDR yields ~45 significant (E. coli 27 + yeast 17 + human 1), log2FC up to ≈ 2.5 |
| **Size** | 1.1 MB |

Backs the package's single registered **`Proteomics Demo`** (`meta.demoPath: Bioinformatics |
Proteomics Differential Expression`). The demo imports this file and runs the whole pipeline —
import → impute → differential expression → volcano + heatmap — landing on the richest analysis
endpoint with no R environment and no manual steps. It also opens a **QC dashboard** on a second
tab and pre-selects the top-significant protein so its **UniProt** entry is already showing in the
context panel. Enrichment is intentionally excluded: g:Profiler enrichment assumes one organism,
and this is a three-species mix, so it would only map the human subset and mislead.

Why a separate file rather than reusing `spectronaut-hye-mix.tsv`: that fixture's `PG.IBAQ`
values were scaled *down* for size (max ≈ 34 k), which lands them in the ambiguous band where the
parser's log2-vs-raw heuristic (`detectLog2Status`) wrongly concludes the data is already
log2-transformed and skips the transform — producing meaningless linear "fold changes" (± thousands).
Scaling every intensity by a constant leaves fold-changes and significance **unchanged** (log2FC is
a difference of logs) but restores realistic magnitudes so the parser log2-transforms correctly.
The result is the textbook HYE outcome: the two spiked species (E. coli, yeast) separate on the
volcano while the constant human background stays flat. Regenerate with:

```
awk -F'\t' 'BEGIN{OFS="\t"; split("R.FileName,R.Condition,R.Replicate,PG.ProteinGroups,PG.Organisms,PG.IBAQ,EG.Qvalue", w, ",")}
NR==1{for(i=1;i<=NF;i++)h[$i]=i; s=""; for(k=1;k<=length(w);k++)s=s (k>1?OFS:"") w[k]; print s; next}
{o=""; for(k=1;k<=length(w);k++){c=w[k]; v=(c in h)?$(h[c]):""; if(c=="PG.IBAQ"&&v!=""&&v!="NA"&&v!="NaN")v=v*1000; o=o (k>1?OFS:"") v} print o}' \
  spectronaut-hye-mix.tsv > spectronaut-hye-demo.tsv
```

## enrichment-demo.csv

| Property | Value |
|----------|-------|
| **Format** | Pre-computed DE result (CSV): `Protein ID, Gene, log2FC, p-value, adj.p-value, significant` |
| **Source** | **Engineered for demonstration.** Real HGNC gene symbols + UniProt accessions; fold-changes and p-values are synthetic. |
| **Organism** | Human (single organism — enrichment requires one) |
| **Proteins** | 119 (30 up, 25 down, 64 non-significant background) |
| **Signal** | UP = cell cycle / mitosis; DOWN = oxidative phosphorylation; background = a diverse, unrelated set |
| **Size** | 4.5 KB |

Backs the registered **`Proteomics Enrichment Demo`** (`meta.demoPath: Proteomics | Enrichment
Analysis`). The demo loads this as a pre-computed DE result, opens the volcano, runs g:Profiler
over-representation analysis (GO / KEGG / Reactome / WikiPathways), and docks the enrichment
charts cross-linked to the volcano — selecting a term highlights its proteins.

Why engineered rather than a real dataset: enrichment assumes a **single organism**, which rules
out the three-species HYE files, and it needs a **coherent** regulated set so the returned terms
are recognizable. The two regulated sets are drawn from textbook, strongly-annotated pathways
(cell cycle up, oxidative phosphorylation down), so g:Profiler returns crisp terms — e.g. KEGG
"Cell cycle" and "Oxidative phosphorylation" at FDR < 1e-15. Gene symbols are real (g:Profiler
maps them); accessions are real (the UniProt panel resolves); only the statistics are synthetic.
Regenerate with `bash tools/generate-enrichment-fixture.sh` (the gene lists live in that script).

## spectronaut-hye-candidates.tsv

| Property | Value |
|----------|-------|
| **Format** | Spectronaut Candidates report (wide, tab-separated, pre-computed DE) |
| **Source** | **Derived** from `spectronaut-hye-mix.tsv` via this package's own Welch's t-test + BH FDR. Regenerate with `node tools/generate-spectronaut-candidates-fixture.mjs`. |
| **Comparison** | HYE mix B / HYE mix A (single comparison) |
| **Proteins** | 93 protein groups (30 Human + 30 Yeast + 30 *E. coli* + 3 other) |
| **Columns** | 18 (Valid, Comparison, Condition Numerator/Denominator, # of Ratios, Group, AVG Group Quantity Num/Den, AVG Log2 Ratio, Absolute AVG Log2 Ratio, % Change, Ratio, Pvalue, Qvalue, ProteinGroups, Genes, UniProtIds, Organisms) |
| **Signal** | log2FC range ≈ [−2.46, +3.66]; 43 significant at \|log2FC\| ≥ 1 ∧ q ≤ 0.05 (26 up, 17 down) |
| **Size** | 20 KB |

Exercises the `Proteomics | Import | Spectronaut Candidates...` path with values that line up
with the source PG file, so the two HYE fixtures tell the same biological story from different
angles. **Not a verbatim Spectronaut export** — we don't have Spectronaut here, so the script
does its own DE math and lays it out in Spectronaut's Candidates column shape. Schema matches
real Spectronaut Candidates output closely enough to round-trip through the importer; numerical
values are ours, not Biognosys's.

## spectronaut-hye-precursor.tsv

| Property | Value |
|----------|-------|
| **Format** | Spectronaut report (long, tab-separated, precursor/fragment-level) |
| **Source** | Synthetic — generated to exercise the streaming Spectronaut import path. Regenerate with `node tools/generate-spectronaut-precursor-fixture.mjs`. |
| **Organism** | Human / Yeast / *E. coli* three-species mix (synthetic) |
| **Proteins** | 41 protein groups (39 passing + 1 `CON__` + 1 `REV__`) |
| **Samples** | 6 runs (two conditions: CondA / CondB, 3 replicates each) |
| **Quantification** | PG.Quantity (constant per protein × condition × replicate) |
| **Columns** | 11 (R.FileName, R.Condition, R.Replicate, PG.ProteinGroups, PG.Organisms, PG.Quantity, EG.Qvalue, PEP.StrippedSequence, EG.ModifiedPeptide, FG.Charge, FG.Id) |
| **Size** | ~44 KB |

Small synthetic precursor-level Spectronaut report: ≥2 distinct precursor/fragment
rows per protein × sample. Carries the D-01 precursor signature columns
(`EG.ModifiedPeptide` / `FG.Charge` / `PEP.StrippedSequence`) so the header-sniff
routes it to the streaming aggregator, uses `CondA`/`CondB`, and carries **no**
`PG.Genes`/`PG.ProteinAccessions` columns (the real `spectronaut-hye-mix.tsv`
lacks them too). Deliberately includes one `CON__` row, one `REV__` row, one
protein mixing a passing precursor with a >0.01 sub-threshold one, one
non-numeric (`Profiled`) q-value protein, and one empty-string q-value protein so
every Spectronaut filter branch is exercised by a single fixture.

## spectronaut-hye-precursor-golden.tsv

| Property | Value |
|----------|-------|
| **Format** | Aggregated Spectronaut report (one row per protein × condition × replicate, tab-separated) |
| **Source** | **Derived** — verbatim duckdb output of `tools/spectronaut-aggregate.sh` over `spectronaut-hye-precursor.tsv`. Regenerate with `tools/spectronaut-aggregate.sh files/demo/spectronaut-hye-precursor.tsv files/demo/spectronaut-hye-precursor-golden.tsv`. |
| **Proteins** | 39 (the `CON__`/`REV__` decoys are dropped by the WHERE filter) |
| **Samples** | 6 (CondA / CondB × 3 replicates) |
| **Rows** | 234 data rows + header |
| **Columns** | 6 (PG.ProteinGroups, R.Condition, R.Replicate, R.FileName, PG.Quantity, EG.Qvalue) |
| **Size** | ~12 KB |

The D-04 equivalence oracle: the real duckdb aggregation of the fixture, not a
hand-derived approximation. The Plan-03 golden test compares the streaming
TypeScript aggregator's output against this file.

## spectronaut-hye-precursor-golden.json

| Property | Value |
|----------|-------|
| **Format** | JSON map `{ "<protein><condition>_<replicate>": { "quantity": number, "qvalue": number } }` |
| **Source** | **Transcribed verbatim** from `spectronaut-hye-precursor-golden.tsv` — no re-aggregation. Regenerate with `node tools/derive-precursor-golden-sidecar.mjs`. |
| **Entries** | 234 (one per surviving protein × condition × replicate) |
| **Size** | ~9 KB |

Read by `grok test` as the in-test oracle when committed-file reads are
unavailable in the test runner. It is **derived FROM** the golden `.tsv` and
never re-aggregated, so the equivalence test stays pinned to real duckdb output;
re-running the deriver against an unchanged golden produces a byte-identical JSON.

### Flip caveat & regeneration order

The committed `tools/spectronaut-aggregate.sql` carries a **reference-file-only**
DMD↔WT `R.Condition` correction (a one-off fix for the mislabeled
`2026-05-13 BP DMD WT.tsv` reference file). It is a structural **no-op** on this
`CondA`/`CondB` fixture, which is exactly why the same script doubles as both the
documented manual fallback **and** the test oracle without perturbing the golden.
The committed SQL also **drops** the two unused `any_value` carry-along SELECT
terms for the gene-symbol and protein-accession columns present in the `/tmp`
original — a documented divergence, because the package's Spectronaut data does
not carry those columns and keeping them would Binder-Error duckdb over this
fixture.

The fixture → golden `.tsv` → golden `.json` chain **must be regenerated in
order**, or the artifacts silently drift:

1. `node tools/generate-spectronaut-precursor-fixture.mjs`
2. `tools/spectronaut-aggregate.sh files/demo/spectronaut-hye-precursor.tsv files/demo/spectronaut-hye-precursor-golden.tsv`
3. `node tools/derive-precursor-golden-sidecar.mjs`

## License

- **proteinGroups.txt** — Synthetic data, no restrictions.
- **cptac-spike-in.txt** — Public benchmark data from the Clinical Proteomic Tumor Analysis Consortium.
- **fragpipe-smoke-test.tsv** — Synthetic data, no restrictions.
- **spectronaut-hye-mix.tsv** — From SpectroPipeR, licensed under MIT.
- **spectronaut-hye-candidates.tsv** — Derived from `spectronaut-hye-mix.tsv` (MIT). DE values computed by this package; no restrictions.
- **spectronaut-hye-precursor.tsv** — Synthetic data, no restrictions.
- **spectronaut-hye-precursor-golden.tsv** — Derived (duckdb aggregation) from `spectronaut-hye-precursor.tsv` (synthetic); no restrictions.
- **spectronaut-hye-precursor-golden.json** — Transcribed from `spectronaut-hye-precursor-golden.tsv` (synthetic); no restrictions.
