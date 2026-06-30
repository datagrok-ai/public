---
status: diagnosed
phase: 12-spectronaut-input-coverage
source: [12-01-SUMMARY.md, 12-02-SUMMARY.md, 12-03-SUMMARY.md]
started: 2026-05-15T17:40:53Z
updated: 2026-05-15T18:30:09Z
---

## Current Test

[testing complete]

## Tests

### 1. Streaming precursor import via menu (committed fixture)
expected: Open Proteomics | Import | Spectronaut Report on files/demo/spectronaut-hye-precursor.tsv. Imports with no error, brief progress indicator, wide DataFrame ~39 proteins × 6 sample columns, CON__/REV__ dropped, tags proteomics.source=spectronaut + proteomics.preNormalized=true + proteomics.groups (2 conditions CondA/CondB).
result: pass
note: "User confirmed 39 rows × 6 sample columns. User observed info message 'Spectronaut import: skipped 30 malformed line(s)'. Investigation: all 492 fixture lines are well-formed (11/11 fields); the 30 'skipped' are correctly-filtered CON__/REV__ decoys + q>0.01 rows mislabeled as 'malformed' by the shared skipped counter in handleFields. Data integrity OK (matches duckdb golden). Mislabeled-message defect logged as a gap — see Gaps."

### 2. Pre-existing spectronaut-hye-mix.tsv imports to an equivalent result (regression)
expected: Importing the canonical existing demo files/demo/spectronaut-hye-mix.tsv via Proteomics | Import | Spectronaut Report produces the same protein × sample result as before Phase 12 (no silent numeric regression). [Test-premise correction: this file is precursor/fragment-level — has PEP.StrippedSequence/EG.ModifiedPeptide/FG.Charge/FG.Id — so sniffIsPrecursor correctly routes it to the NEW streaming path, not the text path as originally written.]
result: pass
note: "User: '93 rows, 8 samples, 18 columns' + info 'skipped 159 malformed line(s)'. Verified: 93 distinct non-decoy proteins, 2 conditions x 4 replicates = 8 samples (18 cols = id + organism + 8 raw + 8 log2) — correct shape. Numeric-equivalence proof: pivotSpectronaut (text) and parseSpectronautStream apply IDENTICAL q<=0.01 + CON__/REV__ filters; they differ only first-encountered (text) vs max (stream) quantity per group. PG.IBAQ verified group-constant across ALL 741 (protein x cond x rep) groups (0 with >1 distinct value) -> first == max -> path switch is a numeric no-op for this file. No data regression. The 159 'malformed' == EXACTLY the 159 EG.Qvalue>0.01 rows in the file -> corroborates the Test-1 mislabel gap with hard evidence (correctly q-filtered precursors, not malformed)."

### 3. 2.6 GB reference file E2E + full pipeline
expected: Import ~/Downloads/2026-05-13 BP DMD WT.tsv (2.6 GB precursor report) via Proteomics | Import | Spectronaut Report. No V8 string-length / OOM error; TaskBar progress advances monotonically; no Chrome "Page Unresponsive" dialog; tab switching works mid-import; result ≈ 8,328 proteins × 24 samples with correct tags; Annotate → Normalize → Impute → DE → Volcano runs end-to-end. (Already executed and PASSED in 12-HUMAN-UAT.md — confirm it still reflects reality.)
result: pass
note: "User: '8328 rows, 50 columns, 24 samples'; context panel confirms proteomics.source=spectronaut, proteomics.preNormalized=true, proteomics.groups auto-populated (DMD/WT). Pipeline Normalize/Impute/DE/Volcano all ran E2E. Headline streaming-import + full-pipeline acceptance MET; 12-HUMAN-UAT.md still reflects reality. Two pre-existing UX issues surfaced during this run (NOT Phase-12 regressions — experiment-setup.ts / normalization.ts untouched by Phase 12) — logged as out-of-scope gaps (obs-A, obs-B)."

### 4. duckdb fallback committed and discoverable
expected: tools/spectronaut-aggregate.sql and tools/spectronaut-aggregate.sh are committed in the package git tree; files/demo/README.md documents them with a regen one-liner and the R.Condition flip caveat; a failed Spectronaut import surfaces a warning pointing at tools/spectronaut-aggregate.sh + files/demo/README.md as the manual fallback.
result: pass
note: "Independently verified: git ls-files shows tools/spectronaut-aggregate.{sh,sql} tracked; README.md:123,165 has the regen one-liner; README.md:150-151 documents the reference-file-only DMD<->WT R.Condition flip caveat; package.ts:144-149 emits grok.shell.error then grok.shell.warning citing tools/spectronaut-aggregate.sh + files/demo/README.md on Spectronaut import failure. All R5 sub-claims hold."

## Summary

total: 4
passed: 4
issues: 0
pending: 0
skipped: 0
blocked: 0
gaps_open: 3
gaps_in_phase_12_scope: 1
gaps_out_of_scope_preexisting: 2

## Gaps

- truth: "Spectronaut precursor import reports only genuinely malformed lines as 'malformed'; correctly-filtered decoy/low-confidence rows are not surfaced as errors"
  status: failed
  reason: "parseSpectronautStream.handleFields (src/parsers/spectronaut-parser.ts:307-344) funnels four distinct return-false reasons — too-few-fields, CON__/REV__ decoy filter, q-value>threshold filter, both-casts-null — into a single `skipped` counter surfaced via grok.shell.info:426 as 'skipped N malformed line(s)'. These are correctly-filtered decoys + q>0.01 precursors, NOT malformed; aggregation output is correct. CORROBORATED across two independent files with exact filter-count matches: (Test 1) committed fixture, all 492 lines well-formed (11/11 fields), reported '30 malformed'; (Test 2) spectronaut-hye-mix.tsv reported '159 malformed' == EXACTLY its 159 EG.Qvalue>0.01 rows. On the 2.6 GB reference file this mislabels ~millions of by-design-filtered precursor rows as 'malformed', falsely signaling a corrupt import on a provably-correct one — most damaging on the exact headline use case Phase 12 delivers."
  severity: major
  test: 1
  also_observed_test: 2
  fix_hint: "Split the single `skipped` counter into (a) genuinely malformed (f.length < expectedFields, !protein, or qty===null && q===null) and (b) by-design filtered (decoy/contaminant, q>threshold). The user-facing grok.shell.info should report only (a) as 'malformed/unparseable'; (b) is normal filtering and should either be silent or phrased as 'filtered N low-confidence/decoy rows'. handleFields should return an enum/discriminant instead of a bare boolean. pivotSpectronaut (text path) applies the same filters silently — matching that silence for category (b) keeps the two paths' UX consistent."
  root_cause: "parseSpectronautStream collapses four semantically distinct handleFields() return-false outcomes (genuinely-malformed: f.length<expectedFields, !protein, qty===null&&q===null; by-design-filtered: CON__/REV__ decoy, q>threshold) into one `skipped` counter (src/parsers/spectronaut-parser.ts:392-393,417-418), then reports it as 'skipped N malformed line(s)' (line 426). Categories 3-4 are intentional duckdb-parity filtering, not malformed. Empirically proven: Test-1 fixture has 0 structurally-malformed lines (492/492 = 11 fields) yet reported 30; Test-2 mix file reported 159 == EXACTLY its 159 EG.Qvalue>0.01 rows. Aggregation output is correct; defect is purely the user-facing message/categorization."
  artifacts:
    - path: "src/parsers/spectronaut-parser.ts"
      issue: "handleFields (307-344) returns a bare boolean; caller loops (392-393, 417-418) funnel every false into one `skipped`; progress msg (402-405) and completion msg (425-426) both label the lot 'malformed'. pivotSpectronaut (47-80) applies the same q-filter/decoy-drop SILENTLY — the UX precedent to match for the by-design category."
  missing:
    - "handleFields returns a discriminated outcome ('kept' | 'malformed' | 'filtered') instead of bare boolean"
    - "Two separate counters (malformed, filtered); aggregation/output logic unchanged"
    - "Completion + in-progress messages report only genuine-malformed as 'malformed'; filtered rows silent (match text path) or neutrally phrased, never 'malformed'"
    - "Regression test: a fixture of only decoy + q>threshold rows yields malformed == 0; existing Spectronaut tests stay green"
  debug_session: ".planning/debug/spectronaut-malformed-mislabel.md"

- id: obs-A
  truth: "Annotate Experiment dialog reflects already-populated proteomics.groups instead of presenting empty Control/Treatment defaults"
  status: backlogged
  routed_to: "ROADMAP.md Backlog Phase 999.1 (.planning/phases/999.1-annotate-dialog-prefill-from-groups/) — DO NOT include in Phase-12 gap-closure"
  reason: "On the 2.6 GB import the Spectronaut parser auto-populated proteomics.groups (DMD/WT, visible in context panel), but showAnnotationDialog (src/analysis/experiment-setup.ts:35-38) hard-codes value:'Control'/'Treatment' and value:[] and never calls the existing getGroups(df) (line 19) to pre-fill. DATA-LOSS RISK: clicking OK as-presented (empty column selections) calls setGroups and overwrites the correct auto-inferred DMD/WT groups with empty groups, silently corrupting the experimental design feeding DE. User raised this directly during Test 3."
  severity: major
  test: 3
  scope: "OUT-OF-PHASE-12 — experiment-setup.ts was not modified by Phase 12 (Phase 12 = streaming input path only; SPEC explicitly scopes downstream pipeline as unchanged). Pre-existing UX/data-loss bug surfaced by Phase-12 UAT. Recommend routing to backlog (gsd-capture), NOT folding into Phase-12 gap-closure. Phase 12's auto-group deliverable (R3) makes this interaction newly visible, so it is worth tracking prominently."
  fix_hint: "In showAnnotationDialog, call getGroups(df); if non-null, seed group1Name/group1Cols/group2Name/group2Cols from it instead of the Control/Treatment/[] defaults."
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""

- id: obs-B
  truth: "Pre-normalized banner wording matches the certainty of the proteomics.preNormalized=true gate it is shown behind"
  status: backlogged
  routed_to: "ROADMAP.md Backlog Phase 999.2 (.planning/phases/999.2-prenormalized-banner-wording/) — DO NOT include in Phase-12 gap-closure"
  reason: "normalization.ts:189 only displays the banner when df.getTag('proteomics.preNormalized') === 'true' (a definite assertion the Spectronaut parser sets unconditionally), yet the copy (normalization.ts:192) hedges: 'This data may be pre-normalized (Spectronaut). Additional normalization may distort results.' Since the banner only renders when pre-normalization is asserted, 'This data IS pre-normalized (Spectronaut)...' is accurate. User raised this directly during Test 3. (Trailing 'may distort' is fine — distortion is genuinely possible, not certain.)"
  severity: cosmetic
  test: 3
  scope: "OUT-OF-PHASE-12 — normalization.ts was not modified by Phase 12 (downstream pipeline explicitly out of scope per SPEC). Pre-existing wording issue surfaced by Phase-12 UAT. Recommend routing to backlog (gsd-capture), NOT Phase-12 gap-closure."
  fix_hint: "normalization.ts:192 — change 'This data may be pre-normalized (Spectronaut).' to 'This data is pre-normalized (Spectronaut).' Keep the second sentence ('Additional normalization may distort results.') as-is."
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""
