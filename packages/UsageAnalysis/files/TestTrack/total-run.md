# Test Track — Global Report

**Date**: 2026-06-22
**Legend**: 🟢 improvement/pass · 🔴 regression/fail · 🟡 partial/ambiguous/new-missing-delta · ⚪ no change / no data
**Verdict**: Not ready

## Column definitions

- **Pass %** — `(pass + 0.5·partial) / total_steps` from each run's `## Steps` section. `ambiguous`/`skip` count toward the denominator. Shown as `NN% (pass/total)`.
- **Browser** — MCP phase wall-clock = model thinking time + live browser interaction time combined.
- **Spec Gen** — model-only time to generate the Playwright spec.
- **Spec Run** — Playwright-only spec execution time.
- **Total** (per test) — sum of Browser + Spec Gen + Spec Run for that scenario. **Mean Total** = average of those sums.
- **Pass Δ (1d)** — Pass % change vs. `prev1d` (most recent committed `total-run.md` before today; 2026-05-12).
- **Pass Δ (7d)** — Pass % change vs. `prev7d` (committed report closest to today−7d, ±3d window). **No commit falls in that window**, so 7d Pass Δ cells are empty; the Trend column carries the longer history.
- **Trend** — last ≤7 daily Pass % values from past commits of `total-run.md`, oldest → newest, dot-separated. Prefix icon: 🟢 last > first, 🔴 last < first, ⚪ equal / single point.

## Folder Summary

**Total**: 272 tests · Run: 152/272 (56%) · Playwright: 265/272 (97%) · Mean Pass: 🟡 88% · Mean Browser: 5m 3s · Mean Spec Gen: 1m 22s · Mean Spec Run: 1m 15s · Mean Total (sum per test): 6m 29s

| Folder | Tests | Run | Playwright | Status | Mean Pass % | Mean Browser | Mean Spec Gen | Mean Spec Run | Mean Total |
|---|---|---|---|---|---|---|---|---|---|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Bio | 24 | 9/24 (38%) | 24/24 (100%) | 🟡 PARTIAL | 🟢 98% | 2m 43s | 43.2s | 47.1s | 4m 14s |
| BiostructureViewer | 11 | 0/11 (0%) | 11/11 (100%) | ⚪ NO DATA |  |  |  |  |  |
| Browse | 2 | 1/2 (50%) | 0/2 (0%) | 🟡 PARTIAL | 🟡 90% |  |  |  |  |
| Charts | 11 | 3/11 (27%) | 10/11 (91%) | 🟡 PARTIAL | 🟡 54% | 2m 10s | 45s | 34.3s | 3m 29s |
| Chat | 1 | 0/1 (0%) | 1/1 (100%) | ⚪ NO DATA |  |  |  |  |  |
| Chem | 18 | 12/18 (67%) | 18/18 (100%) | 🟡 PARTIAL | 🟢 96% | 1m 11s | 29.6s | 42.1s | 2m 22s |
| Connections | 10 | 10/10 (100%) | 10/10 (100%) | 🟡 PARTIAL | 🟡 76% | 4m 34s | 2m 16s | 2m 50s | 9m 41s |
| Dendrogram | 10 | 0/10 (0%) | 10/10 (100%) | ⚪ NO DATA |  |  |  |  |  |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | 🟡 PARTIAL | 🟢 96% | 3m 10s | 56.4s | 45.4s | 4m 52s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | 🟡 PARTIAL | 🟡 67% | 3m | 58s | 17.5s | 2m 45s |
| FileFormats | 1 | 0/1 (0%) | 1/1 (100%) | ⚪ NO DATA |  |  |  |  |  |
| General | 7 | 0/7 (0%) | 6/7 (86%) | ⚪ NO DATA |  |  |  |  |  |
| Helm | 4 | 0/4 (0%) | 4/4 (100%) | ⚪ NO DATA |  |  |  |  |  |
| Models | 7 | 2/7 (29%) | 7/7 (100%) | 🟡 PARTIAL | 🟡 68% | 9m 25s | 1m 55s | 4m 12s | 15m 32s |
| Notebooks | 6 | 0/6 (0%) | 6/6 (100%) | ⚪ NO DATA |  |  |  |  |  |
| Peptides | 14 | 4/14 (29%) | 14/14 (100%) | 🟡 PARTIAL | 🟡 68% |  | 3s | 47.6s | 50.6s |
| PowerPack | 13 | 2/13 (15%) | 13/13 (100%) | 🟡 PARTIAL | 🟡 88% | 8m 57s | 1m 16s | 1m 6s | 11m 18s |
| Projects | 16 | 16/16 (100%) | 16/16 (100%) | 🟡 PARTIAL | 🟡 94% |  |  | 1m 35s | 1m 35s |
| Queries | 14 | 14/14 (100%) | 14/14 (100%) | 🟡 PARTIAL | 🟡 93% | 7m 49s | 2m 45s | 3m 5s | 13m 39s |
| Scripts | 6 | 6/6 (100%) | 6/6 (100%) | 🟡 PARTIAL | 🟡 91% | 2m 53s | 58.2s | 1m 7s | 4m 58s |
| SequenceTranslator | 7 | 0/7 (0%) | 7/7 (100%) | ⚪ NO DATA |  |  |  |  |  |
| Sharing | 11 | 0/11 (0%) | 11/11 (100%) | ⚪ NO DATA |  |  |  |  |  |
| Sketchers | 2 | 0/2 (0%) | 2/2 (100%) | ⚪ NO DATA |  |  |  |  |  |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 77% | 5m 45s | 1m 26s | 56.2s | 8m 8s |
| Tooltips | 7 | 7/7 (100%) | 7/7 (100%) | 🟡 PARTIAL | 🟡 92% | 7m 23s | 2m 38s | 30.1s | 10m 32s |
| Viewers | 46 | 44/46 (96%) | 45/46 (98%) | 🟡 PARTIAL | 🟡 92% | 6m 18s | 1m 16s | 50.8s | 7m 22s |

## All Tests

**Total**: 272 tests · 🟢 105 PASS / 🟡 40 PARTIAL / 🔴 6 FAIL / 🟡 1 SKIP / ⚪ 120 NO RUN · Mean Pass: 🟡 88% · Mean Browser: 5m 3s · Mean Spec Gen: 1m 22s · Mean Spec Run: 1m 15s · Mean Total (sum per test): 6m 29s

| Folder | Test | Status | Pass % | Description | Browser (model+MCP) | Spec Gen (model) | Spec Run (Playwright) | Total (sum) | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Apps | apps | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Apps | tutorials | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [analyze](Bio/analyze-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset typ… | 2m 40s | 20s | 1m 21s | 4m 21s | ⚪ +0% (11/11 → 11/11) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | bio-calculate-scoring | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-cell-actions-panels | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-diversity-search | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-lifecycle-fasta-file | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-lifecycle-immunum-wasm | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-lifecycle-macromolecule-column | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-lifecycle-monomer-collection | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-lifecycle-monomer-library | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-lifecycle-pepsea-container | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-manage-libraries-crud | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-renderer-dispatch | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-service-surface-init | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-similarity-search | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | bio-transform-atomic-level | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 15 checks pass in both browser and Playwright (5 scenario steps × 3 datasets): WebLogo opens from Bio → Analyze → Co… | 4m 13s | 1m 29s | 41s | 6m 23s | ⚪ +0% (5/5 → 5/5) |  | 🟢 80·80·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [convert](Bio/convert-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 sub-steps passed in both the MCP run and the Playwright spec. The key changes that closed the two outstanding fla… | 4m 40s | 1m 10s | 1m 42s | 7m 32s | ⚪ +0% (15/15 → 15/15) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | empty-input-row-viewers | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [manage](Bio/manage-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All three scenario steps passed in the MCP browser run against dev.datagrok.ai: `HELM.csv` opened (540 rows, HELM → Macr… | 36s | 25s | 22s | 1m 23s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [msa](Bio/msa-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | The MSA scenario works end-to-end. In the live browser the new "Clusters" column is added via Edit > Add New Column with… | 2m 50s | 40s | 19s | 3m 49s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (5/6) | Steps 1–5 pass end-to-end against dev: 50-row HELM subset opens with `Macromolecule`/`helm`, Clusters column is added vi… | 3m 15s | 30s | 25s | 4m 10s | 🔴 -9% (5.5/6 → 5/6) |  | 🟢 67·67·100·83·92 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [search](Bio/search-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All four steps passed on the dev server with both the MCP reproduction and the generated Playwright spec. The scenario e… | 1m | 25s | 14s | 1m 39s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets). Activity Cliffs… | 2m 26s | 1m | 1m | 4m 26s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟡 1m (new) | 🔴 +1m |
| Bio | [sequence-space](Bio/sequence-space-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets: FASTA, HELM, MSA… | 2m 50s | 30s | 1m | 4m 20s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟡 1m (new) | 🔴 +1m |
| BiostructureViewer | biostructure-viewer | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | biostructureviewer-bug-claude-33 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | biostructureviewer-bug-grok-14442 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | biostructureviewer-bug-grok-14552 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | biostructureviewer-bug-grok-17485 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | biostructureviewer-bug-grok-17967 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | context-panel-widgets-extension | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | js-api-extension | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | molstar-overlay-extension | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | ngl-viewer-extension | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| BiostructureViewer | property-surface-extension | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 90% (4/5) | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. Item… |  |  |  |  | ⚪ +0% (4.5/5 → 4/5) |  | ⚪ 90·90·90·90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | package-manager | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | charts-api | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | cycle-summary-charts-migrate-2026-05-07 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [radar](Charts/radar-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | The Radar viewer reproduced cleanly on dev for both earthquakes.csv (2426 rows) and demog.csv (5850 rows). All 21 Radar… | 1m 30s | 35s | 33s | 2m 38s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | radar-save-reopen-bug | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 42% (5/12) | The sunburst viewer reproduces structurally on dev — Sunburst can be added to both SPGI and demog, and `hierarchyColumnN… | 3m | 1m | 34s | 4m 34s | ⚪ +0% (5/12 → 5/12) |  | ⚪ 42·42·42·42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | sunburst-date-column-bug | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | sunburst-scatterplot-color-pollution-bug | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | timelines | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 20% (1/5) | Setup (open demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduced cleanly on dev. All four test steps are mar… | 2m | 40s | 36s | 3m 16s | ⚪ +0% (1/5 → 1/5) |  | ⚪ 20·20·20·20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | tree-enhancements-bundle-bug | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | tree-rowsource-onclick-state-bug | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chat | projects-chat-collaboration | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Activity Cliffs computation on SPGI.csv (3624 rows) finishes within 45s and produces a UMAP scatter plot with molecule t… | 1m 14s | 20s | 1m 4s | 2m 38s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (3/6) | Smoke coverage only: Scaffold Tree viewer launches from the Chem menu and the magic wand generates a scaffold tree on SP… | 57s | 25s | 49.6s | 2m 12s | ⚪ +0% (3/6 → 3/6) |  | ⚪ 50·50·50·50·50 | ⚪ +0s | ⚪ +0s | 🟢 -0.4s | 🟢 -0.4s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Scaffold Tree viewer launches from the Chem → Scaffold Tree menu, and the magic-wand generator produces scaffold nodes f… | 1m 15s | 30s | 40.2s | 2m 25s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.2s | 🔴 +0.2s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Similarity Search launches from the Chem menu and exposes a viewer that accepts option changes (fingerprint Morgan ↔ Pat… | 37s | 25s | 24.8s | 1m 27s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.2s | 🟢 -0.2s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Substructure filtering via `grok.chem.searchSubstructure` works on SPGI.csv (3624 rows): benzene substructure yields a b… | 38s | 25s | 29s | 1m 32s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | chem-github-3004 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | chem-grok-12758 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | chem-grok-14028 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | chem-grok-16870 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | chem-grok-17964 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | [chemical-space](Chem/chemical-space-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Chemical Space dimensional reduction runs end-to-end on smiles.csv: the dialog opens, OK with defaults produces a Scatte… | 46s | 20s | 59.8s | 2m 6s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.2s | 🟢 -0.2s |
| Chem | Connectors/database-search-panels | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Elemental Analysis works on dev. The menu path `[name="div-Chem"]` → `Elemental Analysis...` resolves and the dialog ope… | 38s | 30s | 30.9s | 1m 39s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | 🟢 PASS → PASS | 🟢 100% (2/2) | The filter panel correctly shows a Structure filter for SPGI.csv's Molecule column; clicking the embedded sketch-link op… | 34s | 25s | 24.2s | 1m 23s | ⚪ +0% (2/2 → 2/2) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.2s | 🔴 +0.2s |
| Chem | [info-panels](Chem/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | Info panels work correctly on smiles.csv: column-level (Details, Filter, Colors, Style, Chemistry with Rendering/Highlig… | 3m 10s | 1m | 31s | 4m 41s | ⚪ +0% (5/5 → 5/5) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [mmp](Chem/mmp-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset at the bottom of the view… | 1m 27s | 20s | 1m 14s | 3m 1s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | R-Groups Analysis works on sar_small.csv: MCS auto-populates the sketcher, OK produces a Trellis plot and appends R1–R4… | 2m 20s | 45s | 58.7s | 4m 4s | ⚪ +0% (5/5 → 5/5) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.3s | 🟢 -0.3s |
| Chem | [sketcher](Chem/sketcher-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Sketcher opens via `grok.chem.sketcher(molCol, initialSmiles)` wrapped in `ui.dialog(...).show()`, accepts a typed SMILE… | 33s | 30s | 19.3s | 1m 22s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.3s | 🔴 +0.3s |
| Connections | [adding](Connections/adding-run.md) | 🟢 PASS → PASS | 🟡 93% (6/7) | All 7 steps passed end-to-end on dev. Both `test_postgres` and `test_postgres_2` were created and saved (ids `8ab9...` /… | 4m 27s | 45s | 40s | 5m 52s | ⚪ +0% (6.5/7 → 6/7) |  | ⚪ 93·93·93·93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [browser](Connections/browser-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed end-to-end on dev. The scenario depends on a connection named `new_test` existing — it was seeded via… | 6m 56s | 1m 36s | 31s | 9m 3s | ⚪ +0% (9/9 → 9/9) |  | 🟢 78·78·78·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | 🟡 12% (2/17) | Catalog browsing cannot be exercised on dev — all three MS SQL connections in `Browse > Databases > MS SQL` (NorthwindTe… | 2m 35s | 2m | 45s | 5m 20s | ⚪ +0% (2/17 → 2/17) |  | 🟢 9·9·9·12·12 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [delete](Connections/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 sub-steps passed. Both connections (`new_test_postgres` and `test_postgres_2`) were deleted successfully through r… | 1m | 2m | 27s | 3m 27s | ⚪ +0% (8/8 → 8/8) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [edit](Connections/edit-run.md) | 🟢 PASS → PASS | 🟡 86% (6/7) | 6 of 7 steps fully passed; step 7 skipped (real DB password unavailable). The connection rename, credential modification… | 2m | 4m | 1m 6s | 7m 6s | ⚪ +0% (6/7 → 6/7) |  | ⚪ 86·86·86·86·86 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [external-provider](Connections/external-provider-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 7 scenario steps passed in the MCP run against dev: connection PostgreSQLDBTests2 was created via the UI dialog, all… | 5m 54s | 3m 30s | 14m 10s | 23m 34s | ⚪ +0% (8/8 → 8/8) |  | 🟢 0·0·0·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [identifiers](Connections/identifiers-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 17% (2/12) | The scenario depends on a `Configure Identifiers…` connection right-click action. Playwright (fresh Chromium, agolovko s… | 12m 36s | 5m 30s | 4m 22s | 22m 28s | ⚪ +0% (2/12 → 2/12) |  | 🟢 11·11·11·17·17 | ⚪ +0s | ⚪ +0s | 🟡 4m 22s (new) | 🔴 +4m 22s |
| Connections | [import-swagger](Connections/import-swagger-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | Drag-and-drop swagger import is automatable end-to-end on dev. The Chromium limitation that previously caused this scena… | 7m 39s | 58s | 5m 19s | 13m 56s | ⚪ +0% (7/7 → 7/7) |  | 🟢 0·0·0·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [schema](Connections/schema-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 4 steps pass. The previous run mis-read step 3 as "Browse schema" (a non-existent menu item) and marked it AMBIGUOUS… | 57s | 25s | 36s | 1m 58s | ⚪ +0% (4/4 → 4/4) |  | 🟢 75·75·75·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 57% (4/7) | The scenario reproduced cleanly through the UI on dev — the Add-new-connection dialog opened, all fields accepted input,… | 1m 40s | 2m | 29s | 4m 9s | ⚪ +0% (4/7 → 4/7) |  | 🔴 86·86·86·57·57 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Dendrogram | assign-clusters | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | dendrogram-nwk-file-open | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | dendrogram-tree-helper-cross-package | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | dendrogram-tree-helper-traversal | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | dendrogram-viewer-from-newick-prop | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | grok-13041-filter-no-prompt | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | hierarchical-clustering-bio | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | hierarchical-clustering-bio-api | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | hierarchical-clustering-chem | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Dendrogram | hierarchical-clustering-chem-api | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | Catalog scenario reproduces fully on dev.datagrok.ai. All 6 steps PASS both in MCP and in the Playwright spec (35.8s wal… | 1m 53s | 39s | 38s | 3m 10s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Cyclic Models (PK-PD) scenario reproduces fully on dev.datagrok.ai. The PK-PD library model loads via double-click, Mult… | 1m 2s | 38s | 36s | 2m 16s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Files & Sharing scenario reproduces fully on dev.datagrok.ai. pk.ivp loads via the `DiffStudio:previewIvp` function with… | 2m 31s | 1m 32s | 1m 17s | 5m 20s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (5/6) | Scenario is PARTIAL on dev.datagrok.ai — steps 1–5 pass, step 6 (actually running the fit) does not produce result rows… | 10m 2s | 55s | 1m | 11m 57s | ⚪ +0% (5/6 → 5/6) |  | ⚪ 83·83·83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | 🟢 PASS → PASS | 🟡 83% (5/6) | Scenario fully reproduces on dev.datagrok.ai. All 6 steps pass in the interactive MCP session and in the Playwright spec… | 1m 31s | 36s | 26s | 2m 33s | ⚪ +0% (5/6 → 5/6) |  | ⚪ 83·83·83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 scenario steps PASS against dev.datagrok.ai. Edit toggle is reachable via `.d4-ribbon-item .ui-input-bool-switch .… | 4m 30s | 1m 40s | 1m 3s | 7m 13s | ⚪ +0% (5/5 → 5/5) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Sensitivity Analysis scenario fully reproduces against dev.datagrok.ai. Bioreactor loads from the DiffStudio hub (librar… | 2m 3s | 44s | 39s | 3m 26s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Stages (Acid Production) scenario reproduces fully on dev.datagrok.ai. The library card opens a view named "GA-productio… | 1m 49s | 47s | 24s | 3m | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [anova](EDA/anova-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All 3 scenario steps passed against dev. Dataset opens via JS API in ~1s; ANOVA dialog mounts with sensible defaults (RA… | 1m 30s | 30s | 26s | 2m 26s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [MLMethods/linear-regression](EDA/MLMethods/linear-regression-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation, c… |  | 2s | 6.8s | 8.8s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100 | 🟡 45s (removed) | ⚪ +0s | 🟢 -0.2s | 🟢 -45.2s |
| EDA | [MLMethods/pls-regression](EDA/MLMethods/pls-regression-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 71% (5/7) | PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-4 c… |  | 2s | 6.9s | 8.9s | ⚪ +0% (5/7 → 5/7) |  | ⚪ 71·71 | 🟡 1m (removed) | ⚪ +0s | 🟢 -0.1s | 🟢 -1m |
| EDA | [MLMethods/softmax](EDA/MLMethods/softmax-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) | Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns and w… |  | 2s | 2.6s | 4.6s | ⚪ +0% (1/3 → 1/3) |  | ⚪ 33·33 | 🟡 10s (removed) | ⚪ +0s | 🟢 -0.4s | 🟢 -10.4s |
| EDA | [MLMethods/xgboost1](EDA/MLMethods/xgboost1-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as Ui… |  | 2s | 2.7s | 4.7s | ⚪ +0% (2/3 → 2/3) |  | ⚪ 67·67 | 🟡 5s (removed) | ⚪ +0s | 🟢 -0.3s | 🟢 -5.3s |
| EDA | [MLMethods/xgboost2](EDA/MLMethods/xgboost2-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interactio… |  | 2s | 2.7s | 4.7s | ⚪ +0% (2/3 → 2/3) |  | ⚪ 67·67 | 🟡 5s (removed) | ⚪ +0s | 🟢 -0.3s | 🟢 -5.3s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | 2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check, where the wording does not spec… | 2m 30s | 2m | 13s | 4m 43s | ⚪ +0% (2/3 → 2/3) |  | ⚪ 67·67·67·50·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 43% (3/7) | 3 of 7 steps passed, 1 failed, 3 were skipped due to the missing prerequisite dataset. The Pareto Front viewer itself is… | 4m | 2m | 32s | 6m 32s | ⚪ +0% (3/7 → 3/7) |  | ⚪ 43·43·43·43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 60% (3/5) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 3 PASS / 1 FAIL / 1 SKIP. The dialog path works (menu, F… | 5m | 3m | 1m 7s | 9m 7s | ⚪ +0% (3/5 → 3/5) |  | ⚪ 60·60·60·50·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🟡 62% (2/4) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 2 PASS / 1 PARTIAL / 1 FAIL. The dialog path (menu, Usin… | 2m | 2m | 15s | 4m 15s | ⚪ +0% (2.5/4 → 2/4) |  | ⚪ 62·50·62·50·62 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| FileFormats | supported-formats | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | files-cache | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | helm-bio-menu-integration | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | molecule-in-exported-csv | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | network | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | notebooks-lifecycle-jupyter-container | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | profile-settings | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | table-manager | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Helm | helm-editor-and-panels | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Helm | helm-helper | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Helm | helm-lifecycle-macromolecule-column | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Helm | helm-monomer-funcs-override-api | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🟡 35% (6/17) | Only 3 of 16 sub-steps pass fully in both browser and Playwright (1.1 open `smiles.csv`, 1.2 open Train view, 2.1 open `… | 13m | 2m | 3m 45s | 18m 45s | ⚪ +0% (6/17 → 6/17) |  | 🟢 28·28·35·21·35 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | models-bug-grok-3525 | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | models-lifecycle-csv-table | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | models-one-hot-suffix-collision | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | models-testdemog-lifecycle-smoke | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | models-validators-edge | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [predictive-models](Models/predictive-models-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (17/17) | All 17 scenario steps PASSED in the interactive MCP browser run (Scenario 1 Train, Scenario 2 Apply, Scenario 3 Apply on… | 5m 50s | 1m 50s | 4m 40s | 12m 20s | ⚪ +0% (17/17 → 17/17) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Notebooks | browser | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | create | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | notebooks-api-utils-coverage | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | notebooks-context-menu-smoke | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | notebooks-edge | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | notebooks-lifecycle-linked-table | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | collaborative-selection | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | export-invariant-map | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | export-mutation-cliffs | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [info-panels](Peptides/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are rende… |  | 3s | 10.9s | 13.9s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100 | 🟡 17s (removed) | ⚪ +0s | 🟢 -0.1s | 🟢 -17.1s |
| Peptides | manual-alignment | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | monomer-position-hover-tooltip | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | mutation-cliffs-compute-pipeline | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | peptide-sar-demo-dashboard | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variability… |  | 3s | 1m 18s | 1m 21s | ⚪ +0% (2/5 → 2/5) |  | ⚪ 40·40·40·40 | 🟡 25s (removed) | ⚪ +0s | 🟡 1m 18s (new) | 🔴 +53s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (4/6) | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clusters… |  | 3s | 11.6s | 14.6s | ⚪ +0% (4/6 → 4/6) |  | ⚪ 67·67·67·67 | 🟡 18s (removed) | ⚪ +0s | 🟢 -0.4s | 🟢 -18.4s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Residu… |  | 3s | 1m 30s | 1m 33s | ⚪ +0% (6/9 → 6/9) |  | ⚪ 67·67·67·67 | 🟡 50s (removed) | ⚪ +0s | 🟡 1m 30s (new) | 🔴 +40s |
| Peptides | sar-save-reopen | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | sar-similarity-threshold-matrix | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | sar-viewer-lifecycle | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All six scenario steps PASS in both the MCP-driven grok-browser run and the Playwright replay (existing spec — not overw… | 1m 54s | 31s | 11s | 2m 36s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | add-new-column-advanced | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | autocomplete | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 76% (16/21) | The PowerPack "Enrich column" feature works end-to-end for the primary create/apply/edit/delete flow on a dataframe that… | 16m | 2m | 2m | 20m | ⚪ +0% (16/21 → 16/21) |  | 🟢 0·76·76·76·76 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | direct-link-loading | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | formula-refreshing | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | functions-sorting | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | highlight | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | hints | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | input-functions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | power-search-enter | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | widgets-after-debug-delete | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | xlsx-open | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [complex](Projects/complex-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 24s. |  |  | 1m 24s | 1m 24s | 🟢 +26% (17/23 → 1/1) |  | 🟢 0·0·0·0·74 | 🟡 19m (removed) | 🟡 3m (removed) | 🟡 1m 24s (new) | 🟢 -20m 36s |
| Projects | [complex-augment](Projects/complex-augment-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 44.3s. |  |  | 44.3s | 44.3s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 44.3s (new) | 🟡 44.3s (new) |
| Projects | [complex-integration](Projects/complex-integration-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 39.8s. |  |  | 39.8s | 39.8s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 39.8s (new) | 🟡 39.8s (new) |
| Projects | [complex-move](Projects/complex-move-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 18s. |  |  | 1m 18s | 1m 18s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 18s (new) | 🟡 1m 18s (new) |
| Projects | [complex-save-copy](Projects/complex-save-copy-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 3m 06s. |  |  | 3m 6s | 3m 6s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 3m 6s (new) | 🟡 3m 6s (new) |
| Projects | [lifecycle-api](Projects/lifecycle-api-run.md) | ⚪ → 🟢 NO RUN → PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 51.2s. |  |  | 51.2s | 51.2s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 51.2s (new) | 🟡 51.2s (new) |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 → 🟢 SKIP → PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 2m 30s. |  |  | 2m 30s | 2m 30s | 🟢 +100% (0/4 → 1/1) |  | ⚪ 0·0·0·0·0 | ⚪ — | ⚪ — | 🟡 2m 30s (new) | 🟡 2m 30s (new) |
| Projects | [projects-copy-clone](Projects/projects-copy-clone-run.md) | 🟡 → 🟢 SKIP → PASS | 🟢 100% (4/4) | Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 02s. |  |  | 1m 2s | 1m 2s | 🟢 +60% (2/5 → 4/4) |  | ⚪ 40·40·40 | ⚪ — | ⚪ — | 🟡 1m 2s (new) | 🟡 1m 2s (new) |
| Projects | [projects-lifecycle-db](Projects/projects-lifecycle-db-run.md) | 🟢 new: PASS | 🟢 100% (2/2) | Spec passed end-to-end on dev.datagrok.ai. Total run: 2m 30s (2 tests). |  |  | 2m 30s | 2m 30s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 2m 30s (new) | 🟡 2m 30s (new) |
| Projects | [projects-lifecycle-derived](Projects/projects-lifecycle-derived-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 30s. |  |  | 1m 30s | 1m 30s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 30s (new) | 🟡 1m 30s (new) |
| Projects | [projects-lifecycle-files](Projects/projects-lifecycle-files-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 00s. |  |  | 1m | 1m | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m (new) | 🟡 1m (new) |
| Projects | [projects-lifecycle-query](Projects/projects-lifecycle-query-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 56.5s. |  |  | 56.5s | 56.5s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 56.5s (new) | 🟡 56.5s (new) |
| Projects | [projects-lifecycle-script](Projects/projects-lifecycle-script-run.md) | 🟢 new: PASS | 🟢 100% (1/1) | Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 54s. |  |  | 1m 54s | 1m 54s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 54s (new) | 🟡 1m 54s (new) |
| Projects | [projects-lifecycle-spaces](Projects/projects-lifecycle-spaces-run.md) | 🟡 new: SKIP | 🔴 0% (0/1) | Spec env-skipped on dev.datagrok.ai due to bug 1 (Spaces `addEntity` UUID parse failure — see `PLATFORM-BUGS-FOUND-PHASE… |  |  |  |  | 🔴 0% (baseline) |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [projects-ui-smoke](Projects/projects-ui-smoke-run.md) | 🟢 new: PASS | 🟢 100% (7/7) | Spec passed end-to-end on dev.datagrok.ai using bundled Playwright Chromium. Total run: 1m 18s. |  |  | 1m 18s | 1m 18s | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 18s (new) | 🟡 1m 18s (new) |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 100% (3/3) | Spec passed end-to-end on dev.datagrok.ai. Total run: 2m 57s (3 tests). |  |  | 2m 57s | 2m 57s | 🟢 +62% (6/16 → 3/3) |  | 🔴 57·57·57·38·38 | 🟡 10m (removed) | 🟡 2m (removed) | 🔴 +1m 50s | 🟢 -10m 10s |
| Queries | [adding](Queries/adding-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All seven scenario steps passed on dev (https://dev.datagrok.ai/) in both the MCP run and the regenerated Playwright spe… | 34s | 50s | 21s | 1m 45s | ⚪ +0% (7/7 → 7/7) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [browse-and-save-project](Queries/browse-and-save-project-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 8 scenario steps executed successfully against dev.datagrok.ai. Every query on CHEMBL and Northwind (27 + 10 = 37 to… | 11m 30s | 50s | 1m 30s | 13m 50s | ⚪ +0% (9/9 → 9/9) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [browser](Queries/browser-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | The full browse flow worked end-to-end in both the MCP and Playwright runs (24s final spec run, all 4 steps PASSED). Tod… | 1m 50s | 30s | 24s | 2m 44s | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [columns-inspect](Queries/columns-inspect-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | Both parts pass cleanly on dev. Walking the Browse tree from Databases → Provider → Connection → Schemas → public and ex… | 4m 43s | 1m 30s | 2m 45s | 8m 58s | ⚪ +0% (8/8 → 8/8) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [deleting](Queries/deleting-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Delete flow worked in both MCP and Playwright. The dev server's `NorthwindTest` is the friendlyName of the `PostgresTest… | 55s | 30s | 25s | 1m 50s | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [edit](Queries/edit-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All six scenario steps pass in both the MCP run on dev (https://dev.datagrok.ai/) and the Playwright spec (`edit-spec.ts… | 1m | 40s | 3m | 4m 40s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [get-all-get-top-100](Queries/get-all-get-top-100-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | Both halves of the scenario PASS end-to-end on dev. Get All on the orders table returns 830 rows × 14 cols and Get Top 1… | 2m 5s | 4m 40s | 31s | 7m 16s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [ms-sql](Queries/ms-sql-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 59% (9/16) | The scenario completed end-to-end on dev with the same dual-failure mode as the prior run: entity CRUD (add / edit / sav… | 7m 40s | 5m 30s | 1m 10s | 14m 20s | ⚪ +0% (9.5/16 → 9/16) |  | 🔴 60·59 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [new-sql-query](Queries/new-sql-query-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | End-to-end New-SQL-Query-from-table scenario passes via the UI-first path. Right-click → context menu → editor → Play (i… | 1m 35s | 30s | 28s | 2m 33s | ⚪ +0% (5/5 → 5/5) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [new-visual-query](Queries/new-visual-query-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (8/17) | Scenario PARTIAL: 1, 2, 3, 4, 12, 13, 14, 17 PASS; 5, 15 PARTIAL/AMBIGUOUS; 6–11, 16 SKIPPED. The scenario hits two recu… | 20m | 2m | 56s | 22m 56s | ⚪ +0% (8.5/17 → 8/17) |  | ⚪ 50·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [query-layout](Queries/query-layout-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | The full scenario ran end-to-end via grok-browser MCP automation against https://dev.datagrok.ai/: edit PostgresAll → ad… | 6m 20s | 12m | 1m 18s | 19m 38s | ⚪ +0% (13/13 → 13/13) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [query-postprocessing](Queries/query-postprocessing-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All 10 scenario steps reproduced cleanly on dev. The ribbon Play button, post-process JS hook, saved layout, and `Edit…`… | 5m 16s | 3m 12s | 44s | 9m 12s | ⚪ +0% (11/11 → 11/11) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [transformations](Queries/transformations-run.md) | 🟢 PASS → PASS | 🟢 100% (10/10) | The scenario fully PASSES end-to-end in the MCP browser run (all 10 steps verified manually): `${productid}` is added as… | 7m 53s | 54s | 28m 30s | 37m 17s | ⚪ +0% (10/10 → 10/10) |  | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [visual-query-advanced](Queries/visual-query-advanced-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 91% (20/22) | PARTIAL — 19 of 21 steps PASS in the MCP browser run; only steps 19 and 21 FAIL, both due to a server-side regression on… | 38m | 5m | 1m 6s | 44m 6s | ⚪ +0% (20/22 → 20/22) |  | ⚪ 91·91 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Scripts | [browser](Scripts/browser-run.md) | 🟢 PASS → PASS | 🟡 88% (13/16) | Browser scenario passes end-to-end: all seven Context-Panel accordions render, the script is shareable via JS API, chat… | 4m | 3m | 1m 58s | 8m 58s | ⚪ +0% (14/16 → 13/16) |  | 🔴 89·89·78·78·88 | ⚪ +0s | ⚪ +0s | 🟡 1m 58s (new) | 🔴 +1m 58s |
| Scripts | [create](Scripts/create-run.md) | 🟢 PASS → PASS | 🟢 100% (20/20) | All 20 steps PASS end-to-end against `https://dev.datagrok.ai`. Total scenario run (with model): ~9m 12s. The previous r… | 4m 57s | 1m 26s | 2m 49s | 9m 12s | ⚪ +0% (20/20 → 20/20) |  | 🟢 92·92·92·92·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Scripts | [delete](Scripts/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 scenario steps PASSED end-to-end in MCP, and the regenerated Playwright spec PASSED on the first attempt in 12.6s.… | 45s | 30s | 13s | 1m 28s | ⚪ +0% (5/5 → 5/5) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Scripts | [edit](Scripts/edit-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | End-to-end edit works on dev: the seeded `testRscript` opens, `newParam="test"` is appended via CodeMirror, `Save` persi… | 1m 13s | 10s | 22s | 1m 45s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | 🟡 10s (new) | ⚪ +0s | 🔴 +10s |
| Scripts | [layout](Scripts/layout-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 86% (9/11) | Re-running the scenario surfaced two corrections to the previous run's analysis: | 3m 33s | 14s | 49s | 4m 36s | ⚪ +0% (9.5/11 → 9/11) |  | 🟢 46·86 | 🟡 3m 33s (new) | 🟡 14s (new) | 🟡 49s (new) | 🟡 4m 36s (new) |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 73% (8/11) | The scenario passed end-to-end on dev: right-click on the `testRscript` card opened a Run dialog, OK ran the R script an… | 2m 48s | 29s | 30s | 3m 47s | ⚪ +0% (8/11 → 8/11) |  | 🟢 67·67·67·67·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| SequenceTranslator | oligo-nucleotide-grid | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| SequenceTranslator | st-edge-validation-combine | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| SequenceTranslator | st-lifecycle-macromolecule-fasta-separator-column | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| SequenceTranslator | st-lifecycle-macromolecule-helm-column | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| SequenceTranslator | st-lifecycle-molecule-column | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| SequenceTranslator | st-lifecycle-monomer-library-file | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| SequenceTranslator | st-lifecycle-oligonucleotide-column | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | share-connection-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | share-files-and-spaces-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | share-layout-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | share-model-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | share-query-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | share-script-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | share-table-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | sharing-all-users-and-owner | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | sharing-api-permissions | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | sharing-edge-cases | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing | sharing-smoke | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sketchers | filter-panel-cross-context | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sketchers | sketcher-backends | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 3 scenario sections (single-cell add/edit, sticky column behavior, batch edit) were reproduced successfully via MCP… | 5m 11s | 1m | 44s | 6m 55s | ⚪ +0% (9/9 → 9/9) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [copy-clone-delete](StickyMeta/copy-clone-delete-run.md) | 🟢 PASS → PASS | 🟢 100% (10/10) | All 10 steps PASS end-to-end, both in the MCP scenario run and in the Playwright replay (spec wall-clock 1m 29s). Full c… | 7m 52s | 2m | 1m 29s | 11m 21s | ⚪ +0% (10/10 → 10/10) |  | 🟢 50·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 scenario steps PASS in both the MCP browser reproduction and the generated Playwright spec. TestEntity1 is created… | 2m 48s | 1m 24s | 59s | 5m 11s | ⚪ +0% (8/8 → 8/8) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | 🟡 9% (1/11) | The "Database meta" context-panel section is not rendered on `dev.datagrok.ai` for a Postgres DbInfo connection entity (… | 7m 10s | 1m 20s | 33s | 9m 3s | ⚪ +0% (1/11 → 1/11) |  | 🔴 20·20·9·9·9 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [actions-in-the-context-menu](Tooltips/actions-in-the-context-menu-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All five viewers (Grid, Histogram, Line chart, Bar chart, Trellis plot) expose a `Tooltip` group in their right-click co… | 2m 30s | 3m | 19s | 5m 49s | ⚪ +0% (7/7 → 7/7) |  | ⚪ 100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [default-tooltip](Tooltips/default-tooltip-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All six scenario steps PASS. The grid's `Show Visible Columns In Tooltip` property behaves exactly as documented: when u… | 15m | 5m | 22s | 20m 22s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [default-tooltip-visibility](Tooltips/default-tooltip-visibility-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All five scenario steps passed end-to-end against dev.datagrok.ai and the generated Playwright spec runs in 18.8s (1 tes… | 4m 35s | 2m | 22s | 6m 57s | ⚪ +0% (5/5 → 5/5) |  | ⚪ 100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [edit-tooltip](Tooltips/edit-tooltip-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All six scenario steps passed. Dialog elements, case-insensitive search, column selection through the canvas-based colum… | 1m 29s | 2m 18s | 34s | 4m 21s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [line-chart-aggregated-tooltip](Tooltips/line-chart-aggregated-tooltip-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | The aggregated tooltip on a Line Chart with both an applied split (Stereo Category) and a custom aggregated tooltip conf… | 9m 9s | 1m 10s | 51s | 11m 10s | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [tooltip-properties](Tooltips/tooltip-properties-run.md) | 🟢 PASS → PASS | 🟡 91% (10/11) | The scenario tested tooltip property behavior across multiple viewer types. 10 of 11 steps PASS; step 7 is AMBIGUOUS — b… | 7m | 2m | 39s | 9m 39s | ⚪ +0% (10/11 → 10/11) |  | ⚪ 91 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [uniform-default-tooltip](Tooltips/uniform-default-tooltip-run.md) | 🔴 FAIL → FAIL | 🟡 50% (1/2) | Step 1 PASS, Step 2 FAIL. Scatter plot and Box plot both populate their tooltip from `tooltip.dart#getRowTooltip` with t… | 12m | 3m | 24s | 15m 24s | ⚪ +0% (1/2 → 1/2) |  | ⚪ 50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 14 steps (setup + 13 scenario sections) passed in both the MCP browser run against `https://dev.datagrok.ai` and the… | 40s | 10s | 24s | 1m 14s | ⚪ +0% (15/15 → 15/15) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | 🟢 PASS → PASS | 🟡 93% (28/30) | The full annotation-regions Playwright spec passes end-to-end against the local Datagrok at `http://localhost:8888/` aft… | 2m 11s | 1m 33s | 49s | 4m 33s | ⚪ +0% (28/30 → 28/30) |  | 🟢 77·86·85·93·93 | 🟡 2m 11s (new) | 🟡 1m 33s (new) | ⚪ +0s | 🔴 +3m 44s |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (82/82) | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color coding… | 3m 3s | 21s | 52s | 4m 16s | ⚪ +0% (82/82 → 82/82) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 95% (18/19) | 17 of 19 sections passed cleanly; section 8 combined into section 7 in the spec. Section 18 is AMBIGUOUS — `grok.dapi.pr… | 1m 5s | 8s | 32s | 1m 45s | ⚪ +0% (18/19 → 18/19) |  | ⚪ 95·95·95·95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selectio… |  | 1m | 9.4s | 1m 9s | ⚪ +0% (11/11 → 11/11) |  | ⚪ 100·100·100·100 | 🟡 25s (removed) | 🟡 1m (new) | 🟡 9.4s (new) | 🔴 +44.4s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 11 steps passed. The entire test runs on the demog dataset (no SPGI_v2 needed). UI-only steps (Grid Color Coding All… |  |  |  |  | ⚪ +0% (12/12 → 12/12) |  | 🟢 67·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 87% (26/30) | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operation… |  | 3s | 22.5s | 25.5s | ⚪ +0% (26/30 → 26/30) |  | ⚪ 87·87·87·87·87 | ⚪ — | ⚪ +0s | 🔴 +0.5s | 🔴 +0.5s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (58/58) | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI interacti… |  | 2m |  | 2m | ⚪ +0% (58/58 → 58/58) |  | ⚪ 100·100·100·100·100 | ⚪ — | 🟡 2m (new) | ⚪ — | 🟡 2m (new) |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | 🟢 PASS → PASS | 🟢 100% (26/26) | Ran basic-operations end-to-end against dev. All 31 scenario steps passed in the MCP browser phase (Section 1: structure… | 4m 27s | 9s | 50s | 5m 26s | ⚪ +0% (26/26 → 26/26) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | 🟢 PASS → PASS | 🟢 100% (16/16) | Ran chem-and-bio scenario end-to-end against dev. All 11 scenario steps passed in the MCP browser phase (Chem: open spgi… | 2m 50s | 42s | 47s | 4m 19s | ⚪ +0% (16/16 → 16/16) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 scenario steps PASSed on dev. spgi-100.csv loads correctly this time (previous run had to substitute SPGI.csv). C… | 3m 16s | 14s | 55s | 4m 25s | ⚪ +0% (15/15 → 15/15) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed end-to-end on dev: table linking (SELECTION_TO_FILTER and FILTER_TO_FILTER) propagated correctly betw… | 1m 46s | 17s | 35s | 2m 38s | ⚪ +0% (9/9 → 9/9) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran combined-boolean-filter end-to-end against dev. All 13 numbered scenario steps passed in the MCP browser phase: SEX_… | 2m 37s | 12s | 24s | 3m 13s | ⚪ +0% (13/13 → 13/13) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (14/14) | All 14 steps passed in both the MCP run and the Playwright replay. Expression filter works correctly: 5-rule AND yields… | 1m 14s | 8s | 23s | 1m 45s | ⚪ +0% (14/14 → 14/14) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 12 steps passed in the MCP run and in the Playwright replay (spec finished in 21.8s). The hierarchical filter correc… | 1m 15s | 21s | 23s | 1m 59s | ⚪ +0% (12/12 → 12/12) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed in the MCP run and in the Playwright replay (spec finished in 8.7s, total wall-clock 11.56s). The tex… | 1m 12s | 20s | 12s | 1m 44s | ⚪ +0% (9/9 → 9/9) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | 🟢 PASS → PASS | 🟢 100% (34/34) | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to reset… | 4m 24s | 40s | 1m 3s | 6m 7s | ⚪ +0% (34/34 → 34/34) |  | ⚪ 100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | 🟡 93% (28/30) | All 14 sections of form-tests-pw.md exercised across 30 steps. 28 PASS, 2 AMBIGUOUS, 0 FAIL in MCP run. Playwright spec… | 18m | 4m | 3m 12s | 25m 12s | ⚪ +0% (28/30 → 28/30) |  | ⚪ 93·93·92·93·93 | 🟡 18m (new) | 🟡 4m (new) | ⚪ +0s | 🔴 +22m |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas elemen… | 18m | 3m | 51.5s | 21m 52s | ⚪ +0% (36/36 → 36/36) |  | ⚪ 100·100·100·100·100 | 🟡 18m (new) | 🟡 3m (new) | 🟢 -0.5s | 🔴 +21m |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 73% (16/22) | Grid tests ran 22 steps (spec softSteps); 17 passed outright and 5 were AMBIGUOUS (Copy/Paste, Column Header Context Men… | 11m | 3m | 1m 18s | 15m 18s | ⚪ +0% (16/22 → 16/22) |  | ⚪ 73·73·73·73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | 🟡 94% (15/16) | All 14 heat-map sections exercised across 17 steps. 15 PASS, 1 AMBIGUOUS, 1 SKIP in MCP run. Playwright spec passed full… | 18m | 4m | 48.9s | 22m 49s | ⚪ +0% (15/16 → 15/16) |  | ⚪ 94·94·93·94·94 | 🟡 18m (new) | 🟡 4m (new) | 🟢 -0.1s | 🔴 +22m |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (87/94) | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance, l… | 50s | 7s | 46s | 1m 43s | ⚪ +0% (87/94 → 87/94) |  | ⚪ 93·93·93·93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 69% (5/8) | Color consistency through layout round-trip works — the `.categorical-colors` tag survives save/reload and `R_ONE` stays… | 2m 30s | 35s | 26s | 3m 31s | ⚪ +0% (5.5/8 → 5/8) |  | ⚪ 69·69·69·69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (9/12) | Filtering legend updates work end-to-end in the MCP run: numeric filter, categorical filter, layout round-trip, composed… | 3m 10s | 1m 10s | 44s | 5m 4s | ⚪ +0% (10/12 → 9/12) |  | ⚪ 83·83·83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | Legend/legend-api | ⚪ new: NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 77% (8/11) | Line chart legend and multi-axis behaviors are mostly correct: 7 legend items for 7 categories, layout round-trip preser… | 2m 10s | 40s | 32s | 3m 22s | ⚪ +0% (8.5/11 → 8/11) |  | ⚪ 77·77·77·77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 58% (7/13) | Categorical legend on scatter plot updates correctly when X axis changes (sub 2) and when the Filter Panel narrows categ… | 4m 15s | 1m 20s | 54s | 6m 29s | ⚪ +0% (7.5/13 → 7/13) |  | ⚪ 58·58·58·58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 79% (5/7) | Structure rendering in legends works for Scatter plot, Histogram, Line chart and Pie chart (canvas-based molecule thumbn… | 2m 35s | 40s | 29s | 3m 44s | ⚪ +0% (5.5/7 → 5/7) |  | ⚪ 79·79·79·79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 70% (13/20) | Scenario executed end-to-end with a mix of PASS, AMBIGUOUS, and FAIL. Legend display, source-swap, corner positioning, a… | 5m 45s | 1m 30s | 41s | 7m 56s | ⚪ +0% (14/20 → 13/20) |  | ⚪ 70·70·70·70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (26/26) | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout sa… | 57s | 8s | 1m 43s | 2m 48s | ⚪ +0% (26/26 → 26/26) |  | ⚪ 100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [map](Viewers/map-run.md) | 🟢 PASS → PASS | 🟡 80% (8/10) | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker size m… |  | 3s | 9s | 12s | ⚪ +0% (8/10 → 8/10) |  | ⚪ 80·80·80·80 | 🟡 15s (removed) | ⚪ +0s | ⚪ +0s | 🟢 -15s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (16/19) | Matrix Plot tests ran with 15 PASS, 3 AMBIGUOUS, 0 FAIL. The spec executed in 57.7s with all implemented steps passing.… | 20m | 3m | 55.9s | 23m 56s | ⚪ +0% (16/19 → 16/19) |  | ⚪ 84·84·84·84 | 🟡 20m (new) | 🟡 3m (new) | 🟢 -0.1s | 🔴 +23m |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (8/12) | 9 of 12 steps PASS; 3 SKIP (canvas-based node/edge interactions cannot be automated via DOM). The network diagram viewer… | 8m | 1m 30s | 22s | 9m 52s | ⚪ +0% (8/12 → 8/12) |  | ⚪ 67·67·67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | All 13 scenario sections (mapped to 12 Playwright softSteps — scale and normalization are combined in the spec) passed d… | 1m 8s | 8s | 47s | 2m 3s | ⚪ +0% (13/13 → 13/13) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (81/81) | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appearan… | 40s | 7s | 47s | 1m 34s | ⚪ +0% (81/81 → 81/81) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (17/20) | Pivot Table tests ran with 16 PASS, 2 AMBIGUOUS, 1 SKIP, 0 FAIL. The spec executed in 35.1s with all implemented steps p… | 16m | 3m | 1m 12s | 20m 12s | ⚪ +0% (17/20 → 17/20) |  | ⚪ 85·85·85·85 | 🟡 16m (new) | 🟡 3m (new) | ⚪ +0s | 🔴 +19m |
| Viewers | [row-source](Viewers/row-source-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with all 8… |  | 5s | 1m 24s | 1m 29s | ⚪ +0% (36/36 → 36/36) |  | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (20/20) | All 20 sections passed during the MCP run on dev.datagrok.ai. The existing Playwright spec was re-run headed without mod… | 3m 13s | 29s | 52s | 4m 34s | ⚪ +0% (20/20 → 20/20) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | scatter-plot-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [statistics](Viewers/statistics-run.md) | 🟢 PASS → PASS | 🟢 100% (24/24) | All 23 MCP steps passed. The date columns section (STARTED row behavior) was moved to `statistics-tests-ui.md` as a manu… | 20m | 4m | 2m | 26m | ⚪ +0% (24/24 → 24/24) |  | ⚪ 100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟡 2m (new) | 🔴 +2m |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | 🟢 PASS → PASS | 🟢 100% (24/24) | 24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to manual… |  | 3m | 58s | 3m 58s | ⚪ +0% (24/24 → 24/24) |  | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (36/36) | All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via `v… | 28m | 4m | 46s | 32m 46s | ⚪ +0% (36/36 → 36/36) |  | ⚪ 100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (70/84) | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range slid… |  | 30s | 1m 48s | 2m 18s | ⚪ +0% (70.5/84 → 70/84) |  | ⚪ 84·84·84·85·84 | 🟡 3m (removed) | ⚪ +0s | 🟡 1m 48s (new) | 🟢 -1m 12s |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All 7 scenario steps PASS in both the MCP run against https://dev.datagrok.ai and the generated Playwright replay (27.2s… | 2m 17s | 50s | 27s | 3m 34s | ⚪ +0% (7/7 → 7/7) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [word-cloud-tests](Viewers/word-cloud-tests-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (39/46) | All 46 scenario steps were exercised against dev.datagrok.ai. In the MCP run, 35 passed, 3 were AMBIGUOUS (visual effect… | 4m 30s | 2m | 34s | 7m 4s | ⚪ +0% (39/46 → 39/46) |  | ⚪ 85·85·85·85 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [working-with-nan-infinity](Viewers/working-with-nan-infinity-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 spec steps PASSED in 1m 24s. NaN and Infinity values in numeric columns are handled gracefully across Scatter Plot… | 6m | 3m | 1m 24s | 10m 24s | ⚪ +0% (9/9 → 9/9) |  | ⚪ 100·100·100·100 | 🟡 6m (new) | 🟡 3m (new) | ⚪ +0s | 🔴 +9m |
| Browse | browse-tree-states | 🟡 removed: PARTIAL | 🟡 50% (0.5/1) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 50·50·50·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | japanese-in-myfiles | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | local-deploy | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 24 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces-(ui-only) | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 3 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 removed: FAIL | 🟡 33% (1/3) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 33·33·33·33·33 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 removed: PARTIAL | 🟡 40% (2/5) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 40·40·40·40·40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | scaffold-tree | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | scaffold-tree-functions | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | similarity-search | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | structure-filter | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | linear-regression | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/linear-regression | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/pls-regression | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 71·71 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/softmax | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 33·33 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/xgboost1 | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/xgboost2 | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | pls-regression | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | softmax | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | xgboost1 | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | xgboost2 | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | api-samples | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | first-login | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | inactivity-response | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | login | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | startup-time | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | tabs-reordering | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| LocalCashing | local-cashing | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | apply | 🔴 removed: FAIL | 🟡 50% (2/4) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | 🔴 100·100·50·33·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | browser | 🟡 removed: PARTIAL | 🟡 33% (2/6) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | 🔴 67·67·33·20·33 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | delete | 🔴 removed: FAIL | 🟡 20% (1/5) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | 🔴 100·100·20·20·20 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | train | 🟢 removed: PASS | 🟡 91% (10/11) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | 🟢 90·90·91·100·91 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | delete | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | edit | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/add-new-column | 🟡 removed: PARTIAL | 🟢 100% (10/10) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | 🟢 60·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/autocomplete | 🟢 removed: PASS | 🟢 100% (7/7) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/formula-refreshing | 🟢 removed: PASS | 🟢 100% (7/7) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/functions-sorting | 🟢 removed: PASS | 🟢 100% (7/7) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/functions_sorting | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/highlight | 🟢 removed: PASS | 🟢 100% (5/5) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/hints | 🟢 removed: PASS | 🟢 100% (4/4) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/input-functions | 🟢 removed: PASS | 🟢 100% (10/10) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/input_functions | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | 🟢 0·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | formula-lines | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | functions_sorting | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | input_functions | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | browser | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 56·56·56 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 removed: SKIP | 🔴 0% (0/5) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 0·0·0·0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | deleting | 🟢 removed: PASS | 🟢 100% (4/4) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | 🟢 50·50·50·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | opening | 🟡 removed: PARTIAL | 🟢 100% (5/5) | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | projects-copy_clone | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 40·40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | share-project | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | upload-project | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | upload-project-migration-report | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | browse-&-save-project | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| scatter-plot | axes-and-encoding | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| scatter-plot | regression-line-per-category | ⚪ removed: NO RUN |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | copy,-clone,-delete | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | line-chart---aggregated-tooltip | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | bar-chart-tests | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | basic-operations | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | chem-and-bio | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | cloned-views | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | collaborative-filtering-for-linked-tables | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | color-coding-(linked) | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | ⚪ 83 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | color-consistency | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | combined-boolean-filter | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | expression-filter | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | filtering | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | grid-viewer | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | hierarchical-filter | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | rendering-structures-on-the-axes | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | scatterplot | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | statistics-viewer | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | structure-rendering | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | text-filter | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | tile | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | viewers | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | viewers-docking | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | visibility-and-positioning | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | working-with-nan-&-infinity | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | world-cloud | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | ~~grid-viewer~~ | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | ~~tile~~ | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | ~~world-cloud~~ | ⚪ removed: NO DATA |  | _(scenario removed; prior run result retained)_ |  |  |  | — | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |

## Comparison with Previous Reports

Deltas are computed against `prev1d` = the committed `total-run.md` at **2026-05-12** (commit `d4d2df0667`). No commit falls in the 7-day (±3d) window around 2026-06-22−7d, so there is **no `prev7d` baseline**; 7d Pass Δ cells are empty and the Trend column carries the multi-week history instead.

### Totals

**Total (1d)**: Tests Δ **+91** · Run Δ **🔴 -5 (-31%)** · Playwright Δ **🟢 +110 (+12%)** · Mean Pass Δ **🟢 +3%**

**Total (7d)**: Mean Pass Δ **🟢 +11%** _(no commit in 7d±3 window; baseline = oldest report with metrics, 2026-04-22, commit `e90d954935`)_

### By Folder

| Folder | Tests Δ | Run Δ | Playwright Δ | Status | Mean Pass Δ (1d) | Mean Pass Δ (7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|---|
| Apps | +0 | ⚪ +0 | ⚪ +0 | ⚪ NO DATA | — | — | — | — | — | — |
| Bio | +15 | ⚪ +0 | 🟢 +15 | 🟡 PARTIAL | 🔴 -1% | — | — | — | — | — |
| BiostructureViewer (new) | +11 (new) | +0 | +11 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | -1 | 🔴 -1 | ⚪ +0 | 🟡 PARTIAL | 🟢 +20% | — | — | — | — | — |
| Charts | +8 | ⚪ +0 | 🟢 +7 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| Chat (new) | +1 (new) | +0 | +1 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | +4 | 🔴 -2 | 🟢 +4 | 🟡 PARTIAL | 🟢 +9% | — | — | — | — | — |
| Connections | +0 | ⚪ +0 | ⚪ +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| Dendrogram (new) | +10 (new) | +0 | +10 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | +0 | ⚪ +0 | ⚪ +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| EDA | +0 | ⚪ +0 | ⚪ +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| FileFormats (new) | +1 (new) | +0 | +1 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | -3 | ⚪ +0 | 🟢 +6 | ⚪ NO DATA | — | — | — | — | — | — |
| Helm (new) | +4 (new) | +0 | +4 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | +1 | 🔴 -4 | 🟢 +1 | 🟡 PARTIAL | 🟢 +12% | — | — | — | — | — |
| Notebooks | +2 | ⚪ +0 | 🟢 +6 | ⚪ NO DATA | — | — | — | — | — | — |
| Peptides | +10 | ⚪ +0 | 🟢 +10 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| PowerPack | +4 | 🔴 -7 | 🟢 +4 | 🟡 PARTIAL | 🔴 -9% | — | — | — | — | — |
| Projects | +5 | 🟢 +9 | 🟢 +9 | 🟡 PARTIAL | 🟢 +44% | — | — | — | — | — |
| Queries | +0 | ⚪ +0 | ⚪ +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| Scripts | +0 | ⚪ +0 | ⚪ +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| SequenceTranslator (new) | +7 (new) | +0 | +7 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sharing (new) | +11 (new) | +0 | +11 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Sketchers (new) | +2 (new) | +0 | +2 | ⚪ new: NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | +0 | ⚪ +0 | ⚪ +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| Tooltips | +0 | ⚪ +0 | ⚪ +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| Viewers | +1 | ⚪ +0 | 🟢 +1 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — | — |
| scatter-plot (removed) | -2 (removed) | -0 | -0 | ⚪ removed | — | — | — | — | — | — |

### Per-Test Changes

| Folder | Test | Status | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -9% (5.5/6 → 5/6) |  | 🟢 67·67·100·83·92 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟡 1m (new) | 🔴 +1m |
| Bio | [sequence-space](Bio/sequence-space-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟡 1m (new) | 🔴 +1m |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4.5/5 → 4/5) |  | ⚪ 90·90·90·90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/12 → 5/12) |  | ⚪ 42·42·42·42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (1/5 → 1/5) |  | ⚪ 20·20·20·20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/6 → 3/6) |  | ⚪ 50·50·50·50·50 | ⚪ +0s | ⚪ +0s | 🟢 -0.4s | 🟢 -0.4s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.2s | 🔴 +0.2s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.2s | 🟢 -0.2s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.2s | 🟢 -0.2s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | 🟢 PASS → PASS | ⚪ +0% (2/2 → 2/2) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.2s | 🔴 +0.2s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | 🟢 PASS → PASS | ⚪ +0% (5/5 → 5/5) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.3s | 🟢 -0.3s |
| Chem | [sketcher](Chem/sketcher-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.3s | 🔴 +0.3s |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (2/17 → 2/17) |  | 🟢 9·9·9·12·12 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [identifiers](Connections/identifiers-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/12 → 2/12) |  | 🟢 11·11·11·17·17 | ⚪ +0s | ⚪ +0s | 🟡 4m 22s (new) | 🔴 +4m 22s |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/7 → 4/7) |  | 🔴 86·86·86·57·57 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/6 → 5/6) |  | ⚪ 83·83·83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [MLMethods/linear-regression](EDA/MLMethods/linear-regression-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 4/4) |  | ⚪ 100·100 | 🟡 45s (removed) | ⚪ +0s | 🟢 -0.2s | 🟢 -45.2s |
| EDA | [MLMethods/pls-regression](EDA/MLMethods/pls-regression-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/7 → 5/7) |  | ⚪ 71·71 | 🟡 1m (removed) | ⚪ +0s | 🟢 -0.1s | 🟢 -1m |
| EDA | [MLMethods/softmax](EDA/MLMethods/softmax-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/3 → 1/3) |  | ⚪ 33·33 | 🟡 10s (removed) | ⚪ +0s | 🟢 -0.4s | 🟢 -10.4s |
| EDA | [MLMethods/xgboost1](EDA/MLMethods/xgboost1-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/3 → 2/3) |  | ⚪ 67·67 | 🟡 5s (removed) | ⚪ +0s | 🟢 -0.3s | 🟢 -5.3s |
| EDA | [MLMethods/xgboost2](EDA/MLMethods/xgboost2-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/3 → 2/3) |  | ⚪ 67·67 | 🟡 5s (removed) | ⚪ +0s | 🟢 -0.3s | 🟢 -5.3s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/3 → 2/3) |  | ⚪ 67·67·67·50·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/7 → 3/7) |  | ⚪ 43·43·43·43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/5 → 3/5) |  | ⚪ 60·60·60·50·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (2.5/4 → 2/4) |  | ⚪ 62·50·62·50·62 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (6/17 → 6/17) |  | 🟢 28·28·35·21·35 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Peptides | [info-panels](Peptides/info-panels-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100 | 🟡 17s (removed) | ⚪ +0s | 🟢 -0.1s | 🟢 -17.1s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/5 → 2/5) |  | ⚪ 40·40·40·40 | 🟡 25s (removed) | ⚪ +0s | 🟡 1m 18s (new) | 🔴 +53s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/6 → 4/6) |  | ⚪ 67·67·67·67 | 🟡 18s (removed) | ⚪ +0s | 🟢 -0.4s | 🟢 -18.4s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/9 → 6/9) |  | ⚪ 67·67·67·67 | 🟡 50s (removed) | ⚪ +0s | 🟡 1m 30s (new) | 🔴 +40s |
| PowerPack | autocomplete | ⚪ NO RUN → NO RUN | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🟢 → 🟡 PASS → PARTIAL | ⚪ +0% (16/21 → 16/21) |  | 🟢 0·76·76·76·76 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | formula-refreshing | ⚪ NO RUN → NO RUN | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | highlight | ⚪ NO RUN → NO RUN | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | hints | ⚪ NO RUN → NO RUN | — |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [complex](Projects/complex-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 +26% (17/23 → 1/1) |  | 🟢 0·0·0·0·74 | 🟡 19m (removed) | 🟡 3m (removed) | 🟡 1m 24s (new) | 🟢 -20m 36s |
| Projects | [complex-augment](Projects/complex-augment-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 44.3s (new) | 🟡 44.3s (new) |
| Projects | [complex-integration](Projects/complex-integration-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 39.8s (new) | 🟡 39.8s (new) |
| Projects | [complex-move](Projects/complex-move-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 18s (new) | 🟡 1m 18s (new) |
| Projects | [complex-save-copy](Projects/complex-save-copy-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 3m 6s (new) | 🟡 3m 6s (new) |
| Projects | [lifecycle-api](Projects/lifecycle-api-run.md) | ⚪ → 🟢 NO RUN → PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 51.2s (new) | 🟡 51.2s (new) |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 → 🟢 SKIP → PASS | 🟢 +100% (0/4 → 1/1) |  | ⚪ 0·0·0·0·0 | ⚪ — | ⚪ — | 🟡 2m 30s (new) | 🟡 2m 30s (new) |
| Projects | [projects-copy-clone](Projects/projects-copy-clone-run.md) | 🟡 → 🟢 SKIP → PASS | 🟢 +60% (2/5 → 4/4) |  | ⚪ 40·40·40 | ⚪ — | ⚪ — | 🟡 1m 2s (new) | 🟡 1m 2s (new) |
| Projects | [projects-lifecycle-db](Projects/projects-lifecycle-db-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 2m 30s (new) | 🟡 2m 30s (new) |
| Projects | [projects-lifecycle-derived](Projects/projects-lifecycle-derived-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 30s (new) | 🟡 1m 30s (new) |
| Projects | [projects-lifecycle-files](Projects/projects-lifecycle-files-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m (new) | 🟡 1m (new) |
| Projects | [projects-lifecycle-query](Projects/projects-lifecycle-query-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 56.5s (new) | 🟡 56.5s (new) |
| Projects | [projects-lifecycle-script](Projects/projects-lifecycle-script-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 54s (new) | 🟡 1m 54s (new) |
| Projects | [projects-lifecycle-spaces](Projects/projects-lifecycle-spaces-run.md) | 🟡 new: SKIP | 🔴 0% (baseline) |  | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [projects-ui-smoke](Projects/projects-ui-smoke-run.md) | 🟢 new: PASS | 🟢 100% (baseline) |  | — | ⚪ — | ⚪ — | 🟡 1m 18s (new) | 🟡 1m 18s (new) |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 +62% (6/16 → 3/3) |  | 🔴 57·57·57·38·38 | 🟡 10m (removed) | 🟡 2m (removed) | 🔴 +1m 50s | 🟢 -10m 10s |
| Queries | [ms-sql](Queries/ms-sql-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (9.5/16 → 9/16) |  | 🔴 60·59 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [new-visual-query](Queries/new-visual-query-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8.5/17 → 8/17) |  | ⚪ 50·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Queries | [visual-query-advanced](Queries/visual-query-advanced-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (20/22 → 20/22) |  | ⚪ 91·91 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Scripts | [browser](Scripts/browser-run.md) | 🟢 PASS → PASS | ⚪ +0% (14/16 → 13/16) |  | 🔴 89·89·78·78·88 | ⚪ +0s | ⚪ +0s | 🟡 1m 58s (new) | 🔴 +1m 58s |
| Scripts | [edit](Scripts/edit-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) |  | ⚪ 100·100·100·100·100 | ⚪ +0s | 🟡 10s (new) | ⚪ +0s | 🔴 +10s |
| Scripts | [layout](Scripts/layout-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (9.5/11 → 9/11) |  | 🟢 46·86 | 🟡 3m 33s (new) | 🟡 14s (new) | 🟡 49s (new) | 🟡 4m 36s (new) |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/11 → 8/11) |  | 🟢 67·67·67·67·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/11 → 1/11) |  | 🔴 20·20·9·9·9 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | [uniform-default-tooltip](Tooltips/uniform-default-tooltip-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/2 → 1/2) |  | ⚪ 50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | 🟢 PASS → PASS | ⚪ +0% (28/30 → 28/30) |  | 🟢 77·86·85·93·93 | 🟡 2m 11s (new) | 🟡 1m 33s (new) | ⚪ +0s | 🔴 +3m 44s |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (18/19 → 18/19) |  | ⚪ 95·95·95·95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS | ⚪ +0% (11/11 → 11/11) |  | ⚪ 100·100·100·100 | 🟡 25s (removed) | 🟡 1m (new) | 🟡 9.4s (new) | 🔴 +44.4s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (26/30 → 26/30) |  | ⚪ 87·87·87·87·87 | ⚪ — | ⚪ +0s | 🔴 +0.5s | 🔴 +0.5s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | ⚪ +0% (58/58 → 58/58) |  | ⚪ 100·100·100·100·100 | ⚪ — | 🟡 2m (new) | ⚪ — | 🟡 2m (new) |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | ⚪ +0% (28/30 → 28/30) |  | ⚪ 93·93·92·93·93 | 🟡 18m (new) | 🟡 4m (new) | ⚪ +0s | 🔴 +22m |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | ⚪ +0% (36/36 → 36/36) |  | ⚪ 100·100·100·100·100 | 🟡 18m (new) | 🟡 3m (new) | 🟢 -0.5s | 🔴 +21m |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (16/22 → 16/22) |  | ⚪ 73·73·73·73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | ⚪ +0% (15/16 → 15/16) |  | ⚪ 94·94·93·94·94 | 🟡 18m (new) | 🟡 4m (new) | 🟢 -0.1s | 🔴 +22m |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (87/94 → 87/94) |  | ⚪ 93·93·93·93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5.5/8 → 5/8) |  | ⚪ 69·69·69·69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (10/12 → 9/12) |  | ⚪ 83·83·83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8.5/11 → 8/11) |  | ⚪ 77·77·77·77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (7.5/13 → 7/13) |  | ⚪ 58·58·58·58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5.5/7 → 5/7) |  | ⚪ 79·79·79·79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (14/20 → 13/20) |  | ⚪ 70·70·70·70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [map](Viewers/map-run.md) | 🟢 PASS → PASS | ⚪ +0% (8/10 → 8/10) |  | ⚪ 80·80·80·80 | 🟡 15s (removed) | ⚪ +0s | ⚪ +0s | 🟢 -15s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (16/19 → 16/19) |  | ⚪ 84·84·84·84 | 🟡 20m (new) | 🟡 3m (new) | 🟢 -0.1s | 🔴 +23m |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/12 → 8/12) |  | ⚪ 67·67·67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (17/20 → 17/20) |  | ⚪ 85·85·85·85 | 🟡 16m (new) | 🟡 3m (new) | ⚪ +0s | 🔴 +19m |
| Viewers | [statistics](Viewers/statistics-run.md) | 🟢 PASS → PASS | ⚪ +0% (24/24 → 24/24) |  | ⚪ 100·100·100·100 | ⚪ +0s | ⚪ +0s | 🟡 2m (new) | 🔴 +2m |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (70.5/84 → 70/84) |  | ⚪ 84·84·84·85·84 | 🟡 3m (removed) | ⚪ +0s | 🟡 1m 48s (new) | 🟢 -1m 12s |
| Viewers | [word-cloud-tests](Viewers/word-cloud-tests-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (39/46 → 39/46) |  | ⚪ 85·85·85·85 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [working-with-nan-infinity](Viewers/working-with-nan-infinity-run.md) | 🟢 PASS → PASS | ⚪ +0% (9/9 → 9/9) |  | ⚪ 100·100·100·100 | 🟡 6m (new) | 🟡 3m (new) | ⚪ +0s | 🔴 +9m |
| Browse | browse-tree-states | 🟡 removed: PARTIAL | — | — | ⚪ 50·50·50·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | japanese-in-myfiles | ⚪ removed: NO DATA | — | — | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | local-deploy | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces | ⚪ removed: NO DATA | — | — | ⚪ 24 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces-(ui-only) | ⚪ removed: NO DATA | — | — | ⚪ 3 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 removed: FAIL | — | — | ⚪ 33·33·33·33·33 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 removed: PARTIAL | — | — | ⚪ 40·40·40·40·40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | scaffold-tree | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | scaffold-tree-functions | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | similarity-search | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Chem | structure-filter | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | linear-regression | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/linear-regression | ⚪ removed: NO DATA | — | — | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/pls-regression | ⚪ removed: NO DATA | — | — | ⚪ 71·71 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/softmax | ⚪ removed: NO DATA | — | — | ⚪ 33·33 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/xgboost1 | ⚪ removed: NO DATA | — | — | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | ML methods/xgboost2 | ⚪ removed: NO DATA | — | — | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | pls-regression | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | softmax | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | xgboost1 | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| EDA | xgboost2 | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | api-samples | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | first-login | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | inactivity-response | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | login | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | startup-time | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | tabs-reordering | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| LocalCashing | local-cashing | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | apply | 🔴 removed: FAIL | — | — | 🔴 100·100·50·33·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | browser | 🟡 removed: PARTIAL | — | — | 🔴 67·67·33·20·33 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | delete | 🔴 removed: FAIL | — | — | 🔴 100·100·20·20·20 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | train | 🟢 removed: PASS | — | — | 🟢 90·90·91·100·91 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | delete | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | edit | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/add-new-column | 🟡 removed: PARTIAL | — | — | 🟢 60·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/autocomplete | 🟢 removed: PASS | — | — | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/formula-refreshing | 🟢 removed: PASS | — | — | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/functions-sorting | 🟢 removed: PASS | — | — | ⚪ 100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/functions_sorting | ⚪ removed: NO DATA | — | — | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/highlight | 🟢 removed: PASS | — | — | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/hints | 🟢 removed: PASS | — | — | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/input-functions | 🟢 removed: PASS | — | — | ⚪ 100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | AddNewColumn/input_functions | ⚪ removed: NO DATA | — | — | 🟢 0·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | formula-lines | ⚪ removed: NO DATA | — | — | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | functions_sorting | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | input_functions | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | browser | ⚪ removed: NO DATA | — | — | ⚪ 56·56·56 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 removed: SKIP | — | — | ⚪ 0·0·0·0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | deleting | 🟢 removed: PASS | — | — | 🟢 50·50·50·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | opening | 🟡 removed: PARTIAL | — | — | ⚪ 100·100·100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | projects-copy_clone | ⚪ removed: NO DATA | — | — | ⚪ 40·40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | share-project | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | upload-project | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | upload-project-migration-report | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | browse-&-save-project | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| scatter-plot | axes-and-encoding | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| scatter-plot | regression-line-per-category | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | copy,-clone,-delete | ⚪ removed: NO DATA | — | — | ⚪ 50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | line-chart---aggregated-tooltip | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | bar-chart-tests | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | basic-operations | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | chem-and-bio | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | cloned-views | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | collaborative-filtering-for-linked-tables | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | color-coding-(linked) | ⚪ removed: NO DATA | — | — | ⚪ 83 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | color-consistency | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | combined-boolean-filter | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | expression-filter | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | filtering | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | grid-viewer | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | hierarchical-filter | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | rendering-structures-on-the-axes | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | scatterplot | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | statistics-viewer | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | structure-rendering | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | text-filter | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | tile | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | viewers | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | viewers-docking | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | visibility-and-positioning | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | working-with-nan-&-infinity | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | world-cloud | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | ~~grid-viewer~~ | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | ~~tile~~ | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | ~~world-cloud~~ | ⚪ removed: NO DATA | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |

## Release Readiness

**Verdict**: Not ready

Run coverage is 152/272 (56%), below the 70% bar, so the suite is **not ready**. The scenario set grew sharply since the 2026-05-12 report (181→272 tests) while run execution did not keep pace, so coverage regressed. No folder is fully FAIL, but a large block of folders has zero executed runs.

### Blocking Issues
- Apps: 2 scenarios, 0 runs (no `-run.md` results)
- BiostructureViewer: 11 scenarios, 0 runs (no `-run.md` results)
- Chat: 1 scenarios, 0 runs (no `-run.md` results)
- Dendrogram: 10 scenarios, 0 runs (no `-run.md` results)
- FileFormats: 1 scenarios, 0 runs (no `-run.md` results)
- General: 7 scenarios, 0 runs (no `-run.md` results)
- Helm: 4 scenarios, 0 runs (no `-run.md` results)
- Notebooks: 6 scenarios, 0 runs (no `-run.md` results)
- SequenceTranslator: 7 scenarios, 0 runs (no `-run.md` results)
- Sharing: 11 scenarios, 0 runs (no `-run.md` results)
- Sketchers: 2 scenarios, 0 runs (no `-run.md` results)
- Bio: only 9/24 scenarios run (38%)
- Charts: only 3/11 scenarios run (27%)
- Models: only 2/7 scenarios run (29%)
- Peptides: only 4/14 scenarios run (29%)
- PowerPack: only 2/13 scenarios run (15%)
