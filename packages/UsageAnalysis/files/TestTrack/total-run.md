# Test Track — Global Report

**Date**: 2026-05-05
**Legend**: 🟢 improvement/pass · 🔴 regression/fail · 🟡 partial/ambiguous/new-missing-delta · ⚪ no change / no data
**Verdict**: Conditionally ready

## Column definitions

- **Pass %** — `(pass + 0.5·partial) / total_steps` from each run's `## Steps` section.
- **Browser** — MCP phase wall-clock = model thinking time + live browser interaction time.
- **Spec Gen** — model-only time to generate the Playwright spec.
- **Spec Run** — Playwright-only spec execution time.
- **Total** (per test) — sum of Browser + Spec Gen + Spec Run for that scenario. **Mean Total** = average of those sums.
- **Pass Δ (1d)** — Pass % change vs. `prev1d` (the most recent committed `total-run.md` strictly before today).
- **Pass Δ (7d)** — Pass % change vs. `prev7d` (the committed `total-run.md` closest to today − 7 days, ±3-day window). Empty when no commit falls in that window.
- **Trend** — last ≤7 daily Pass % values, oldest → newest, dot-separated. Prefix icon: 🟢 last > first, 🔴 last < first, ⚪ equal or only one point.

## Folder Summary

**Total**: 181 tests · Run: 150/181 (83%) · Playwright: 144/181 (80%) · Mean Pass: 🟡 83% · Mean Browser: 4m 7s · Mean Spec Gen: 1m 12s · Mean Spec Run: 1m 7s · Mean Total (sum per test): 6m

| Folder | Tests | Run | Playwright | Status | Mean Pass % | Mean Browser | Mean Spec Gen | Mean Spec Run | Mean Total |
|---|---|---|---|---|---|---|---|---|---|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Bio | 9 | 9/9 (100%) | 9/9 (100%) | 🟡 PARTIAL | 🟢 98% | 2m 43s | 43s | 43s | 4m |
| Browse | 3 | 2/3 (67%) | 0/3 (0%) | 🟡 PARTIAL |  |  |  |  |  |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | 🟡 PARTIAL | 🟡 54% | 2m 10s | 45s | 34s | 3m 29s |
| Chem | 14 | 14/14 (100%) | 14/14 (100%) | 🟡 PARTIAL | 🟡 87% | 1m 37s | 32s | 40s | 2m 46s |
| Connections | 10 | 10/10 (100%) | 10/10 (100%) | 🟡 PARTIAL | 🟡 76% | 4m 34s | 2m 16s | 2m 40s | 9m 15s |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | 🟡 PARTIAL | 🟢 96% | 3m 10s | 56s | 45s | 4m 52s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | 🟡 PARTIAL | 🟡 59% | 1m 42s | 58s | 17s | 2m 58s |
| General | 10 | 0/10 (0%) | 0/10 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | 🟡 PARTIAL | 🟡 49% | 8m 20s | 2m 3s | 1m 45s | 12m 8s |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL |  | 28s | 3s | 11s | 36s |
| PowerPack | 9 | 9/9 (100%) | 9/9 (100%) | 🟡 PARTIAL | 🟢 97% | 5m 16s | 1m 9s | 30s | 6m 54s |
| Projects | 11 | 7/11 (64%) | 3/11 (27%) | 🟡 PARTIAL | 🟡 40% | 6m 38s | 2m | 1m 7s | 9m 12s |
| Queries | 14 | 14/14 (100%) | 14/14 (100%) | 🟡 PARTIAL | 🟡 93% | 7m 49s | 2m 45s | 3m 5s | 13m 39s |
| Scripts | 6 | 6/6 (100%) | 6/6 (100%) | 🟡 PARTIAL | 🟡 80% | 3m 16s | 34s | 32s | 4m 22s |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 77% | 5m 45s | 1m 26s | 56s | 8m 8s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Viewers | 45 | 44/45 (98%) | 44/45 (98%) | 🟡 PARTIAL | 🟡 92% | 4m 8s | 52s | 49s | 4m 37s |
| scatter-plot | 2 | 0/2 (0%) | 0/2 (0%) | ⚪ NO DATA |  |  |  |  |  |

## All Tests

**Total**: 181 tests · 🟢 89 PASS / 🟡 46 PARTIAL / 🔴 8 FAIL / 🟡 4 SKIP / ⚪ 34 NO RUN · Mean Pass: 🟡 83% · Mean Browser: 4m 7s · Mean Spec Gen: 1m 12s · Mean Spec Run: 1m 7s · Mean Total (sum per test): 6m

| Folder | Test | Status | Pass % | Description | Browser (model+MCP) | Spec Gen (model) | Spec Run (Playwright) | Total (sum) | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Apps | apps | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Apps | tutorials | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [analyze](Bio/analyze-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset typ | 2m 40s | 20s | 1m 21s | 4m 21s | ⚪ +0% (11/11 → 11/11) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 15 checks pass in both browser and Playwright (5 scenario steps × 3 datasets): WebLogo opens from Bio → Analyze → Co | 4m 13s | 1m 29s | 41s | 6m 23s | ⚪ +0% (5/5 → 4/4) | ⚪ — (no 7d baseline) | 🟢 80·80·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [convert](Bio/convert-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 sub-steps passed in both the MCP run and the Playwright spec. The key changes that closed the two outstanding fla | 4m 40s | 1m 10s | 1m 42s | 7m 32s | ⚪ +0% (15/15 → 15/15) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [manage](Bio/manage-run.md) | 🟢 PASS → PASS | 🟢 100% (2/2) | All three scenario steps passed in the MCP browser run against dev.datagrok.ai: `HELM.csv` opened (540 rows, HELM → Macr | 36s | 25s | 22s | 1m 23s | ⚪ +0% (3/3 → 2/2) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [msa](Bio/msa-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | The MSA scenario works end-to-end. In the live browser the new "Clusters" column is added via Edit > Add New Column with | 2m 50s | 40s | 19s | 3m 49s | ⚪ +0% (6/6 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (5/6) | Steps 1–5 pass end-to-end against dev: 50-row HELM subset opens with `Macromolecule`/`helm`, Clusters column is added vi | 3m 15s | 30s | 25s | 4m 10s | 🔴 -17% (5/5 → 5/6) | ⚪ — (no 7d baseline) | 🟢 67·67·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [search](Bio/search-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All four steps passed on the dev server with both the MCP reproduction and the generated Playwright spec. The scenario e | 1m | 25s | 14s | 1m 39s | ⚪ +0% (4/4 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets). Activity Cliffs | 2m 26s | 1m |  | 3m 26s | ⚪ +0% (6/6 → 6/6) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Bio | [sequence-space](Bio/sequence-space-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets: FASTA, HELM, MSA | 2m 50s | 30s |  | 3m 20s | ⚪ +0% (6/6 → 6/6) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL |  | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. Item |  |  |  |  | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 90·90·90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | 🟡 PARTIAL → PARTIAL |  | 1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a single session |  |  |  |  | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 50·50·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | package-manager | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [radar](Charts/radar-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | The Radar viewer reproduced cleanly on dev for both earthquakes.csv (2426 rows) and demog.csv (5850 rows). All 21 Radar | 1m 30s | 35s | 33s | 2m 38s | ⚪ +0% (3/3 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 42% (5/12) | The sunburst viewer reproduces structurally on dev — Sunburst can be added to both SPGI and demog, and `hierarchyColumnN | 3m | 1m | 34s | 4m 34s | ⚪ 0% (5/12 → 5/12) | ⚪ — (no 7d baseline) | ⚪ 42·42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 20% (1/5) | Setup (open demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduced cleanly on dev. All four test steps are mar | 2m | 40s | 36s | 3m 16s | ⚪ +0% (1/5 → 1/5) | ⚪ — (no 7d baseline) | ⚪ 20·20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (3/6) | Smoke coverage only: Scaffold Tree viewer launches from the Chem menu and the magic wand generates a scaffold tree on SP | 57s | 25s | 50s | 2m 12s | ⚪ +0% (3/6 → 3/6) | ⚪ — (no 7d baseline) | ⚪ 50·50·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Scaffold Tree viewer launches from the Chem → Scaffold Tree menu, and the magic-wand generator produces scaffold nodes f | 1m 15s | 30s | 40s | 2m 25s | ⚪ +0% (3/3 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Similarity Search launches from the Chem menu and exposes a viewer that accepts option changes (fingerprint Morgan ↔ Pat | 37s | 25s | 25s | 1m 27s | ⚪ +0% (3/3 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Substructure filtering via `grok.chem.searchSubstructure` works on SPGI.csv (3624 rows): benzene substructure yields a b | 38s | 25s | 29s | 1m 32s | ⚪ +0% (4/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Activity Cliffs computation on SPGI.csv (3624 rows) finishes within 45s and produces a UMAP scatter plot with molecule t | 1m 14s | 20s | 1m 4s | 2m 38s | ⚪ +0% (4/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) | Calculate Descriptors cannot be exercised on `dev` right now. The Chem top menu fails to open its popup — both through D | 8m | 1m |  | 9m | ⚪ +0% (1/3 → 1/3) | ⚪ — (no 7d baseline) | ⚪ 33·33·33 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Chemical Space dimensional reduction runs end-to-end on smiles.csv: the dialog opens, OK with defaults produces a Scatte | 46s | 20s | 1m | 2m 6s | ⚪ +0% (3/3 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | ChemProp scenario is partially automated: the spec confirms mol1K.sdf opens and the Train Model view is reachable from t | 26s | 30s | 19s | 1m 15s | ⚪ +0% (2/5 → 2/5) | ⚪ — (no 7d baseline) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Elemental Analysis works on dev. The menu path `[name="div-Chem"]` → `Elemental Analysis...` resolves and the dialog ope | 38s | 30s | 31s | 1m 39s | ⚪ +0% (3/3 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | 🟢 PASS → PASS | 🟢 100% (2/2) | The filter panel correctly shows a Structure filter for SPGI.csv's Molecule column; clicking the embedded sketch-link op | 34s | 25s | 24s | 1m 23s | ⚪ +0% (2/2 → 2/2) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [info-panels](Chem/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | Info panels work correctly on smiles.csv: column-level (Details, Filter, Colors, Style, Chemistry with Rendering/Highlig | 3m 10s | 1m | 31s | 4m 41s | ⚪ +0% (5/5 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [mmp](Chem/mmp-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset at the bottom of the view | 1m 27s | 20s | 1m 14s | 3m 1s | ⚪ +0% (3/3 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | R-Groups Analysis works on sar_small.csv: MCS auto-populates the sketcher, OK produces a Trellis plot and appends R1–R4 | 2m 20s | 45s | 59s | 4m 4s | ⚪ +0% (5/5 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [sketcher](Chem/sketcher-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Sketcher opens via `grok.chem.sketcher(molCol, initialSmiles)` wrapped in `ui.dialog(...).show()`, accepts a typed SMILE | 33s | 30s | 19s | 1m 22s | ⚪ +0% (3/3 → 3/3) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [adding](Connections/adding-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟡 93% (6.5/7) | All 7 steps passed end-to-end on dev. Both `test_postgres` and `test_postgres_2` were created and saved (ids `8ab9...` / | 4m 27s | 45s | 40s | 5m 52s | ⚪ 0% (6/7 → 6.5/7) | ⚪ — (no 7d baseline) | ⚪ 93·93·93 | 🟡 4m 27s (new) | 🟡 45s (new) | 🟡 40s (new) | 🟡 5m 52s (new) |
| Connections | [browser](Connections/browser-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 100% (9/9) | All 9 steps passed end-to-end on dev. The scenario depends on a connection named `new_test` existing — it was seeded via | 6m 56s | 1m 36s | 31s | 9m 3s | 🟢 +22% (7/9 → 9/9) | ⚪ — (no 7d baseline) | ⚪ 78·78·78 | 🟡 6m 56s (new) | 🟡 1m 36s (new) | 🟡 31s (new) | 🟡 9m 3s (new) |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | 🟡 12% (2/17) | Catalog browsing cannot be exercised on dev — all three MS SQL connections in `Browse > Databases > MS SQL` (NorthwindTe | 2m 35s | 2m | 45s | 5m 20s | 🟢 +3% (1/11 → 2/17) | ⚪ — (no 7d baseline) | ⚪ 9·9·9 | 🟡 2m 35s (new) | 🟡 2m (new) | 🟡 45s (new) | 🟡 5m 20s (new) |
| Connections | [delete](Connections/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 sub-steps passed. Both connections (`new_test_postgres` and `test_postgres_2`) were deleted successfully through r | 1m | 2m | 27s | 3m 27s | ⚪ +0% (8/8 → 8/8) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | 🟡 1m (new) | 🟡 2m (new) | 🟡 27s (new) | 🟡 3m 27s (new) |
| Connections | [edit](Connections/edit-run.md) | 🟢 PASS → PASS | 🟡 86% (6/7) | 6 of 7 steps fully passed; step 7 skipped (real DB password unavailable). The connection rename, credential modification | 2m | 4m | 1m 6s | 7m 6s | ⚪ 0% (6/7 → 6/7) | ⚪ — (no 7d baseline) | ⚪ 86·86·86 | 🟡 2m (new) | 🟡 4m (new) | 🟡 1m 6s (new) | 🟡 7m 6s (new) |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 → 🟢 FAIL → PASS | 🟢 100% (8/8) | All 7 scenario steps passed in the MCP run against dev: connection PostgreSQLDBTests2 was created via the UI dialog, all | 5m 54s | 3m 30s | 14m 10s | 23m 34s | 🟢 +100% (0/7 → 8/8) | ⚪ — (no 7d baseline) | ⚪ 0·0·0 | 🟡 5m 54s (new) | 🟡 3m 30s (new) | 🟡 14m 10s (new) | 🟡 23m 34s (new) |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 → 🟡 FAIL → PARTIAL | 🟡 17% (2/12) | The scenario depends on a `Configure Identifiers…` connection right-click action. Playwright (fresh Chromium, agolovko s | 12m 36s | 5m 30s |  | 18m 6s | 🟢 +6% (1/9 → 2/12) | ⚪ — (no 7d baseline) | ⚪ 11·11·11 | 🟡 12m 36s (new) | 🟡 5m 30s (new) | ⚪ — | 🟡 18m 6s (new) |
| Connections | [import-swagger](Connections/import-swagger-run.md) | ⚪ removed: FAIL | 🟢 100% (7/7) | Drag-and-drop swagger import is automatable end-to-end on dev. The Chromium limitation that previously caused this scena | 7m 39s | 58s | 5m 19s | 13m 56s | 🟢 +100% (0/7 → 7/7) | ⚪ — (no 7d baseline) | ⚪ 0·0·0 | 🟡 7m 39s (new) | 🟡 58s (new) | 🟡 5m 19s (new) | 🟡 13m 56s (new) |
| Connections | [schema](Connections/schema-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 100% (4/4) | All 4 steps pass. The previous run mis-read step 3 as "Browse schema" (a non-existent menu item) and marked it AMBIGUOUS | 57s | 25s | 36s | 1m 58s | 🟢 +25% (3/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 75·75·75 | 🟡 57s (new) | 🟡 25s (new) | 🟡 36s (new) | 🟡 1m 58s (new) |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 57% (4/7) | The scenario reproduced cleanly through the UI on dev — the Add-new-connection dialog opened, all fields accepted input, | 1m 40s | 2m | 29s | 4m 9s | 🔴 -29% (6/7 → 4/7) | ⚪ — (no 7d baseline) | ⚪ 86·86·86 | 🟡 1m 40s (new) | 🟡 2m (new) | 🟡 29s (new) | 🟡 4m 9s (new) |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | Catalog scenario reproduces fully on dev.datagrok.ai. All 6 steps PASS both in MCP and in the Playwright spec (35.8s wal | 1m 53s | 39s | 38s | 3m 10s | ⚪ +0% (6/6 → 6/6) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Cyclic Models (PK-PD) scenario reproduces fully on dev.datagrok.ai. The PK-PD library model loads via double-click, Mult | 1m 2s | 38s | 36s | 2m 16s | ⚪ +0% (4/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Files & Sharing scenario reproduces fully on dev.datagrok.ai. pk.ivp loads via the `DiffStudio:previewIvp` function with | 2m 31s | 1m 32s | 1m 17s | 5m 20s | ⚪ +0% (4/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (5/6) | Scenario is PARTIAL on dev.datagrok.ai — steps 1–5 pass, step 6 (actually running the fit) does not produce result rows | 10m 2s | 55s | 1m | 11m 57s | ⚪ +0% (5/6 → 5/6) | ⚪ — (no 7d baseline) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | 🟢 PASS → PASS | 🟡 83% (5/6) | Scenario fully reproduces on dev.datagrok.ai. All 6 steps pass in the interactive MCP session and in the Playwright spec | 1m 31s | 36s | 26s | 2m 33s | ⚪ +0% (5/6 → 5/6) | ⚪ — (no 7d baseline) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 scenario steps PASS against dev.datagrok.ai. Edit toggle is reachable via `.d4-ribbon-item .ui-input-bool-switch . | 4m 30s | 1m 40s | 1m 3s | 7m 13s | ⚪ +0% (5/5 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Sensitivity Analysis scenario fully reproduces against dev.datagrok.ai. Bioreactor loads from the DiffStudio hub (librar | 2m 3s | 44s | 39s | 3m 26s | ⚪ +0% (4/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Stages (Acid Production) scenario reproduces fully on dev.datagrok.ai. The library card opens a view named "GA-productio | 1m 49s | 47s | 24s | 3m | ⚪ +0% (4/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [MLMethods/linear-regression](EDA/MLMethods/linear-regression-run.md) | 🟢 PASS → PASS |  | Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation, c | 45s | 2s | 7s | 54s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [MLMethods/pls-regression](EDA/MLMethods/pls-regression-run.md) | 🟡 PARTIAL → PARTIAL |  | PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-4 c | 1m | 2s | 7s | 1m 9s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 71 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [MLMethods/softmax](EDA/MLMethods/softmax-run.md) | 🔴 FAIL → FAIL |  | Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns and w | 10s | 2s | 3s | 15s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 33 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [MLMethods/xgboost1](EDA/MLMethods/xgboost1-run.md) | 🟡 PARTIAL → PARTIAL |  | XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as Ui | 5s | 2s | 3s | 10s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [MLMethods/xgboost2](EDA/MLMethods/xgboost2-run.md) | 🟡 PARTIAL → PARTIAL |  | XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interactio | 5s | 2s | 3s | 10s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [anova](EDA/anova-run.md) | 🟢 PASS → PASS | 🟢 100% (2/2) | All 3 scenario steps passed against dev. Dataset opens via JS API in ~1s; ANOVA dialog mounts with sensible defaults (RA | 1m 30s | 30s | 26s | 2m 26s | ⚪ +0% (3/3 → 2/2) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (1/2) | 2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check, where the wording does not spec | 2m 30s | 2m | 13s | 4m 43s | 🔴 -17% (2/3 → 1/2) | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 43% (3/7) | 3 of 7 steps passed, 1 failed, 3 were skipped due to the missing prerequisite dataset. The Pareto Front viewer itself is | 4m | 2m | 32s | 6m 32s | ⚪ 0% (3/7 → 3/7) | ⚪ — (no 7d baseline) | ⚪ 43·43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (2/4) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 3 PASS / 1 FAIL / 1 SKIP. The dialog path works (menu, F | 5m | 3m | 1m 7s | 9m 7s | 🔴 -10% (3/5 → 2/4) | ⚪ — (no 7d baseline) | ⚪ 60·60·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🟡 50% (1.5/3) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 2 PASS / 1 PARTIAL / 1 FAIL. The dialog path (menu, Usin | 2m | 2m | 15s | 4m 15s | 🔴 -12% (2/4 → 1.5/3) | ⚪ — (no 7d baseline) | ⚪ 62·50·62 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| General | files-cache | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | first-login | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | inactivity-response | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | login | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | molecule-in-exported-csv | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | network | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | profile-settings | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | startup-time | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | table-manager | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | tabs-reordering | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [apply](Models/apply-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) | Steps 1 and 2 passed — demog.csv opened and the "Apply predictive model" dialog opened correctly via the top menu. Step | 2m 12s | 30s | 20s | 3m 2s | 🔴 -17% (2/4 → 1/3) | ⚪ — (no 7d baseline) | 🔴 100·100·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 20% (1/5) | 2 of 6 steps pass fully in both browser and Playwright (step 1: navigation; step 4: Filter Templates). Steps 2, 3, 5, an | 4m 45s | 2m 10s | 51s | 7m 46s | 🔴 -13% (2/6 → 1/5) | ⚪ — (no 7d baseline) | 🔴 67·67·33 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🟡 21% (3/14) | Only 3 of 16 sub-steps pass fully in both browser and Playwright (1.1 open `smiles.csv`, 1.2 open Train view, 2.1 open ` | 13m | 2m | 3m 45s | 18m 45s | 🔴 -14% (6/17 → 3/14) | ⚪ — (no 7d baseline) | 🟢 28·28·35 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [delete](Models/delete-run.md) | 🔴 FAIL → FAIL | 🟡 20% (1/5) | The delete UI itself could not be exercised on dev because the prerequisite predictive model from train.md/browser.md do | 18m | 4m | 24s | 22m 24s | ⚪ +0% (1/5 → 1/5) | ⚪ — (no 7d baseline) | 🔴 100·100·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [predictive-models](Models/predictive-models-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (16/16) | All 17 scenario steps PASSED in the interactive MCP browser run (Scenario 1 Train, Scenario 2 Apply, Scenario 3 Apply on | 5m 50s | 1m 50s | 4m 40s | 12m 20s | ⚪ +0% (17/17 → 16/16) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [train](Models/train-run.md) | 🟢 PASS → PASS | 🟢 100% (10/10) | All 10 scenario steps passed in the interactive MCP browser run. Both models — classification (TestDemog, Predict probab | 6m 10s | 1m 50s | 28s | 8m 28s | 🟢 +9% (10/11 → 10/10) | ⚪ — (no 7d baseline) | 🟢 90·90·91 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Notebooks | browser | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | create | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | delete | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | edit | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [info-panels](Peptides/info-panels-run.md) | 🟢 PASS → PASS |  | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are rende | 17s | 3s | 11s | 31s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL |  | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variability | 25s | 3s |  | 28s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | 🟡 1m 18s (removed) | 🟢 -1m 18s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL |  | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clusters | 18s | 3s | 12s | 33s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL |  | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Residu | 50s | 3s |  | 53s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | 🟡 1m 30s (removed) | 🟢 -1m 30s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (10/10) | All 10 scenario steps reproduce successfully in the MCP run; 9/10 pass in the Playwright replay. The one FAILED Playwrig | 7m 45s | 2m | 32s | 10m 17s | ⚪ +0% (10/10 → 10/10) | ⚪ — (no 7d baseline) | 🟢 60·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All six autocomplete behaviours PASS in both the MCP run and the Playwright replay. `.cm-tooltip-autocomplete` appears o | 1m 45s | 20s | 8s | 2m 13s | ⚪ +0% (7/7 → 7/7) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All seven sub-steps pass in both the MCP run and the Playwright replay. Dependency propagation across calculated columns | 2m 15s | 1m | 28s | 3m 43s | ⚪ +0% (7/7 → 7/7) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/functions-sorting](PowerPack/AddNewColumn/functions-sorting-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All five scenario steps PASS in both the MCP run and the Playwright replay (17s, 1 test, 0 failures). The previous run's | 8m 20s | 2m | 17s | 10m 37s | ⚪ +0% (7/7 → 6/6) | ⚪ — (no 7d baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps PASS in both the MCP run and the Playwright replay. `${AGE}`, `$[AGE]` and the autocomplete-inserted `${HEIG | 4m 43s | 1m 13s | 25s | 6m 21s | ⚪ +0% (5/5 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All four steps pass in both MCP and Playwright. The CodeMirror formula editor shows a `.cm-tooltip-hover` on hover with | 1m 10s | 15s | 8s | 1m 33s | ⚪ +0% (4/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/input-functions](PowerPack/AddNewColumn/input-functions-run.md) | 🟢 PASS → PASS | 🟢 100% (10/10) | All 10 scenario steps pass end-to-end — both in interactive MCP driving and in the Playwright replay (18.8s). A single i | 3m 30s | 1m | 19s | 4m 49s | ⚪ +0% (10/10 → 10/10) | ⚪ — (no 7d baseline) | ⚪ 100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All six scenario steps PASS in both the MCP-driven grok-browser run and the Playwright replay (existing spec — not overw | 1m 54s | 31s | 11s | 2m 36s | ⚪ +0% (6/6 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | ⚪ removed: PARTIAL | 🟡 76% (16/21) | The PowerPack "Enrich column" feature works end-to-end for the primary create/apply/edit/delete flow on a dataframe that | 16m | 2m | 2m | 20m | ⚪ +0% (16/21 → 16/21) | ⚪ — (no 7d baseline) | 🟢 0·76·76 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Projects | [complex](Projects/complex-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/13) | All 13 steps skipped. This is the most complex scenario requiring tables from 7+ different sources, drag-and-drop, entit |  |  |  |  | ⚪ +0% (0/13 → 0/13), broken | ⚪ — (no 7d baseline) | ⚪ 0·0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/5) | All 5 steps skipped. This scenario requires running a custom JavaScript script with Data Sync enabled, then modifying fi |  |  |  |  | ⚪ +0% (0/5 → 0/5), broken | ⚪ — (no 7d baseline) | ⚪ 0·0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 100% (4/4) | All four scenario steps PASS in the live MCP run on `dev.datagrok.ai`: two test projects were created via API, located i | 3m 17s | 2m |  | 5m 17s | 🟢 +50% (2/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 50·50·50 | 🟡 3m 17s (new) | 🟡 2m (new) | ⚪ — | 🟡 5m 17s (new) |
| Projects | lifecycle-api | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [opening](Projects/opening-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (5/5) | All 5 steps passed. Projects from the Uploading step are accessible in Browse > Dashboards. Context Panel correctly show |  |  |  |  | ⚪ +0% (5/5 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/4) | All steps skipped. This scenario depends on Projects copy_clone.md (order 5) which was not fully executed. The Link/Clon |  |  |  |  | ⚪ +0% (0/4 → 0/4), broken | ⚪ — (no 7d baseline) | ⚪ 0·0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [projects-copy-clone](Projects/projects-copy-clone-run.md) | 🟡 SKIP → SKIP | 🟡 40% (2/5) | 2 of 5 steps passed, 3 skipped. Project preview and opening work. Copy/clone/link operations were not tested because the |  |  |  |  | ⚪ +0% (2/5 → 2/5) | ⚪ — (no 7d baseline) | ⚪ 40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | share-project | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | upload-project | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | upload-project-migration-report | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 38% (6/16) | Total scenario run (with model): ~14m. Of the 9 cases, only Case 1 was runnable on dev — and only after substituting `Sy | 10m | 2m | 1m 7s | 13m 7s | 🔴 -20% (8/14 → 6/16) | ⚪ — (no 7d baseline) | ⚪ 57·57·57 | 🟡 10m (new) | 🟡 2m (new) | 🟡 1m 7s (new) | 🟡 13m 7s (new) |
| Queries | [adding](Queries/adding-run.md) | 🟢 new: PASS | 🟢 100% (7/7) | All seven scenario steps passed on dev (https://dev.datagrok.ai/) in both the MCP run and the regenerated Playwright spe | 34s | 50s | 21s | 1m 45s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 34s (new) | 🟡 50s (new) | 🟡 21s (new) | 🟡 1m 45s (new) |
| Queries | [browse-and-save-project](Queries/browse-and-save-project-run.md) | 🟢 new: PASS | 🟢 100% (9/9) | All 8 scenario steps executed successfully against dev.datagrok.ai. Every query on CHEMBL and Northwind (27 + 10 = 37 to | 11m 30s | 50s | 1m 30s | 13m 50s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 11m 30s (new) | 🟡 50s (new) | 🟡 1m 30s (new) | 🟡 13m 50s (new) |
| Queries | [browser](Queries/browser-run.md) | 🟢 new: PASS | 🟢 100% (4/4) | The full browse flow worked end-to-end in both the MCP and Playwright runs (24s final spec run, all 4 steps PASSED). Tod | 1m 50s | 30s | 24s | 2m 44s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 1m 50s (new) | 🟡 30s (new) | 🟡 24s (new) | 🟡 2m 44s (new) |
| Queries | [columns-inspect](Queries/columns-inspect-run.md) | 🟢 new: PASS | 🟢 100% (4/4) | Both parts pass cleanly on dev. Walking the Browse tree from Databases → Provider → Connection → Schemas → public and ex | 4m 43s | 1m 30s | 2m 45s | 8m 58s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 4m 43s (new) | 🟡 1m 30s (new) | 🟡 2m 45s (new) | 🟡 8m 58s (new) |
| Queries | [deleting](Queries/deleting-run.md) | 🟢 new: PASS | 🟢 100% (3/3) | Delete flow worked in both MCP and Playwright. The dev server's `NorthwindTest` is the friendlyName of the `PostgresTest | 55s | 30s | 25s | 1m 50s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 55s (new) | 🟡 30s (new) | 🟡 25s (new) | 🟡 1m 50s (new) |
| Queries | [edit](Queries/edit-run.md) | 🟢 new: PASS | 🟢 100% (6/6) | All six scenario steps pass in both the MCP run on dev (https://dev.datagrok.ai/) and the Playwright spec (`edit-spec.ts | 1m | 40s | 3m | 4m 40s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 1m (new) | 🟡 40s (new) | 🟡 3m (new) | 🟡 4m 40s (new) |
| Queries | [get-all-get-top-100](Queries/get-all-get-top-100-run.md) | 🟢 new: PASS | 🟢 100% (3/3) | Both halves of the scenario PASS end-to-end on dev. Get All on the orders table returns 830 rows × 14 cols and Get Top 1 | 2m 5s | 4m 40s | 31s | 7m 16s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 2m 5s (new) | 🟡 4m 40s (new) | 🟡 31s (new) | 🟡 7m 16s (new) |
| Queries | [ms-sql](Queries/ms-sql-run.md) | 🟡 new: PARTIAL | 🟡 60% (3/5) | The scenario completed end-to-end on dev with the same dual-failure mode as the prior run: entity CRUD (add / edit / sav | 7m 40s | 5m 30s | 1m 10s | 14m 20s | 🟡 60% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 7m 40s (new) | 🟡 5m 30s (new) | 🟡 1m 10s (new) | 🟡 14m 20s (new) |
| Queries | [new-sql-query](Queries/new-sql-query-run.md) | 🟢 new: PASS | 🟢 100% (5/5) | End-to-end New-SQL-Query-from-table scenario passes via the UI-first path. Right-click → context menu → editor → Play (i | 1m 35s | 30s | 28s | 2m 33s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 1m 35s (new) | 🟡 30s (new) | 🟡 28s (new) | 🟡 2m 33s (new) |
| Queries | [new-visual-query](Queries/new-visual-query-run.md) | 🟡 new: PARTIAL | 🟡 50% (8.5/17) | Scenario PARTIAL: 1, 2, 3, 4, 12, 13, 14, 17 PASS; 5, 15 PARTIAL/AMBIGUOUS; 6–11, 16 SKIPPED. The scenario hits two recu | 20m | 2m | 56s | 22m 56s | 🟡 50% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 20m (new) | 🟡 2m (new) | 🟡 56s (new) | 🟡 22m 56s (new) |
| Queries | [query-layout](Queries/query-layout-run.md) | 🟢 new: PASS | 🟢 100% (10/10) | The full scenario ran end-to-end via grok-browser MCP automation against https://dev.datagrok.ai/: edit PostgresAll → ad | 6m 20s | 12m | 1m 18s | 19m 38s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 6m 20s (new) | 🟡 12m (new) | 🟡 1m 18s (new) | 🟡 19m 38s (new) |
| Queries | [query-postprocessing](Queries/query-postprocessing-run.md) | 🟢 new: PASS | 🟢 100% (11/11) | All 10 scenario steps reproduced cleanly on dev. The ribbon Play button, post-process JS hook, saved layout, and `Edit…` | 5m 16s | 3m 12s | 44s | 9m 12s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 5m 16s (new) | 🟡 3m 12s (new) | 🟡 44s (new) | 🟡 9m 12s (new) |
| Queries | [transformations](Queries/transformations-run.md) | 🟢 new: PASS | 🟢 100% (10/10) | The scenario fully PASSES end-to-end in the MCP browser run (all 10 steps verified manually): `${productid}` is added as | 7m 53s | 54s | 28m 30s | 37m 17s | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 7m 53s (new) | 🟡 54s (new) | 🟡 28m 30s (new) | 🟡 37m 17s (new) |
| Queries | [visual-query-advanced](Queries/visual-query-advanced-run.md) | 🟡 new: PARTIAL | 🟡 91% (20/22) | PARTIAL — 19 of 21 steps PASS in the MCP browser run; only steps 19 and 21 FAIL, both due to a server-side regression on | 38m | 5m | 1m 6s | 44m 6s | 🟡 91% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 38m (new) | 🟡 5m (new) | 🟡 1m 6s (new) | 🟡 44m 6s (new) |
| Scripts | [browser](Scripts/browser-run.md) | 🟢 PASS → PASS | 🟡 78% (7/9) | Browser scenario passed in the MCP session — all accordions (Details, Script, Run, Activity, Sharing, Chats, Dev) render | 1m 17s | 25s | 16s | 1m 58s | ⚪ 0% (7/9 → 7/9) | ⚪ — (no 7d baseline) | 🔴 89·89·78 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 92% (12/13) | The Create scenario mostly works end-to-end on dev: R script is created, signature editor lets the name be set, a parame | 8m 39s | 50s | 1m 14s | 10m 43s | ⚪ +0% (10/12 → 12/13) | ⚪ — (no 7d baseline) | ⚪ 92·92·92 | 🟡 8m 39s (new) | 🟡 50s (new) | 🟡 1m 14s (new) | 🟡 10m 43s (new) |
| Scripts | [delete](Scripts/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | Delete passed end-to-end in the MCP session: confirmation dialog was correctly shown ("Are you sure? Delete script 'test | 57s | 25s | 22s | 1m 44s | ⚪ +0% (5/5 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | 🟡 57s (new) | 🟡 25s (new) | 🟡 22s (new) | 🟡 1m 44s (new) |
| Scripts | [edit](Scripts/edit-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | Edit works end-to-end in the MCP run: the script reopens, a new line is added via CodeMirror, Save persists, and re-open | 1m 25s | 30s | 21s | 2m 16s | ⚪ +0% (6/6 → 6/6) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | 🟡 1m 25s (new) | 🟡 30s (new) | 🟡 21s (new) | 🟡 2m 16s (new) |
| Scripts | [layout](Scripts/layout-run.md) | 🟡 new: PARTIAL | 🟡 46% (5.5/12) | Script creation and the Layout tab's "Run script" flow work — the preview dataframe is built from `System:DemoFiles/chem | 5m 9s | 35s | 40s | 6m 24s | 🟡 46% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 5m 9s (new) | 🟡 35s (new) | 🟡 40s (new) | 🟡 6m 24s (new) |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Run scenario passed end-to-end in the MCP session: right-click → Run opened the dialog with `cars` pre-selected, clickin | 2m 10s | 40s | 19s | 3m 9s | ⚪ 0% (6/9 → 6/9) | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | 🟡 2m 10s (new) | 🟡 40s (new) | 🟡 19s (new) | 🟡 3m 9s (new) |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 3 scenario sections (single-cell add/edit, sticky column behavior, batch edit) were reproduced successfully via MCP | 5m 11s | 1m | 44s | 6m 55s | ⚪ +0% (9/9 → 9/9) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [copy-clone-delete](StickyMeta/copy-clone-delete-run.md) | 🟢 PASS → PASS | 🟢 100% (10/10) | All 10 steps PASS end-to-end, both in the MCP scenario run and in the Playwright replay (spec wall-clock 1m 29s). Full c | 7m 52s | 2m | 1m 29s | 11m 21s | ⚪ +0% (10/10 → 10/10) | ⚪ — (no 7d baseline) | 🟢 50·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 scenario steps PASS in both the MCP browser reproduction and the generated Playwright spec. TestEntity1 is created | 2m 48s | 1m 24s | 59s | 5m 11s | ⚪ +0% (8/8 → 8/8) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | 🟡 9% (1/11) | The "Database meta" context-panel section is not rendered on `dev.datagrok.ai` for a Postgres DbInfo connection entity ( | 7m 10s | 1m 20s | 33s | 9m 3s | ⚪ +0% (1/11 → 1/11) | ⚪ — (no 7d baseline) | 🔴 20·20·9 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | actions-in-the-context-menu | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip-visibility | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | edit-tooltip | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | line-chart-aggregated-tooltip | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | tooltip-properties | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | uniform-default-tooltip | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 14 steps (setup + 13 scenario sections) passed in both the MCP browser run against `https://dev.datagrok.ai` and the | 40s | 10s | 24s | 1m 14s | ⚪ +0% (15/15 → 15/15) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran basic-operations end-to-end against dev. All 31 scenario steps passed in the MCP browser phase (Section 1: structure | 4m 27s | 9s | 50s | 5m 26s | ⚪ +0% (13/13 → 13/13) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | Ran chem-and-bio scenario end-to-end against dev. All 11 scenario steps passed in the MCP browser phase (Chem: open spgi | 2m 50s | 42s | 47s | 4m 19s | ⚪ +0% (11/11 → 11/11) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 scenario steps PASSed on dev. spgi-100.csv loads correctly this time (previous run had to substitute SPGI.csv). C | 3m 16s | 14s | 55s | 4m 25s | ⚪ +0% (15/15 → 15/15) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed end-to-end on dev: table linking (SELECTION_TO_FILTER and FILTER_TO_FILTER) propagated correctly betw | 1m 46s | 17s | 35s | 2m 38s | ⚪ +0% (9/9 → 9/9) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran combined-boolean-filter end-to-end against dev. All 13 numbered scenario steps passed in the MCP browser phase: SEX_ | 2m 37s | 12s | 24s | 3m 13s | ⚪ +0% (13/13 → 13/13) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (14/14) | All 14 steps passed in both the MCP run and the Playwright replay. Expression filter works correctly: 5-rule AND yields | 1m 14s | 8s | 23s | 1m 45s | ⚪ +0% (14/14 → 14/14) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 12 steps passed in the MCP run and in the Playwright replay (spec finished in 21.8s). The hierarchical filter correc | 1m 15s | 21s | 23s | 1m 59s | ⚪ +0% (12/12 → 12/12) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed in the MCP run and in the Playwright replay (spec finished in 8.7s, total wall-clock 11.56s). The tex | 1m 12s | 20s | 12s | 1m 44s | ⚪ +0% (9/9 → 9/9) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | 🟢 PASS → PASS |  | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to reset | 4m 24s | 40s | 1m 3s | 6m 7s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 69% (5.5/8) | Color consistency through layout round-trip works — the `.categorical-colors` tag survives save/reload and `R_ONE` stays | 2m 30s | 35s | 26s | 3m 31s | ⚪ 0% (5/8 → 5.5/8) | ⚪ — (no 7d baseline) | ⚪ 69·69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (10/12) | Filtering legend updates work end-to-end in the MCP run: numeric filter, categorical filter, layout round-trip, composed | 3m 10s | 1m 10s | 44s | 5m 4s | ⚪ +0% (9/12 → 10/12) | ⚪ — (no 7d baseline) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 77% (8.5/11) | Line chart legend and multi-axis behaviors are mostly correct: 7 legend items for 7 categories, layout round-trip preser | 2m 10s | 40s | 32s | 3m 22s | ⚪ +0% (8/11 → 8.5/11) | ⚪ — (no 7d baseline) | ⚪ 77·77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 58% (7.5/13) | Categorical legend on scatter plot updates correctly when X axis changes (sub 2) and when the Filter Panel narrows categ | 4m 15s | 1m 20s | 54s | 6m 29s | ⚪ 0% (7/13 → 7.5/13) | ⚪ — (no 7d baseline) | ⚪ 58·58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 79% (5.5/7) | Structure rendering in legends works for Scatter plot, Histogram, Line chart and Pie chart (canvas-based molecule thumbn | 2m 35s | 40s | 29s | 3m 44s | ⚪ 0% (5/7 → 5.5/7) | ⚪ — (no 7d baseline) | ⚪ 79·79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 70% (14/20) | Scenario executed end-to-end with a mix of PASS, AMBIGUOUS, and FAIL. Legend display, source-swap, corner positioning, a | 5m 45s | 1m 30s | 41s | 7m 56s | ⚪ +0% (13/20 → 14/20) | ⚪ — (no 7d baseline) | ⚪ 70·70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | ⚪ removed: PASS | 🟡 93% (28/30) | The full annotation-regions Playwright spec passes end-to-end against the local Datagrok at `http://localhost:8888/` aft |  |  | 49s | 49s | 🟢 +8% (11/13 → 28/30) | ⚪ — (no 7d baseline) | 🟢 77·86·85 | ⚪ — | ⚪ — | 🟡 49s (new) | 🟡 49s (new) |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (82/82) | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color coding | 3m 3s | 21s | 52s | 4m 16s | ⚪ +0% (82/82 → 82/82) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 95% (18/19) | 17 of 19 sections passed cleanly; section 8 combined into section 7 in the spec. Section 18 is AMBIGUOUS — `grok.dapi.pr | 1m 5s | 8s | 32s | 1m 45s | ⚪ 0% (18/19 → 18/19) | ⚪ — (no 7d baseline) | ⚪ 95·95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS |  | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selectio | 25s |  |  | 25s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ — | ⚪ — | ⚪ +0s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | 🟢 PASS → PASS |  | All 11 steps passed. The entire test runs on the demog dataset (no SPGI_v2 needed). UI-only steps (Grid Color Coding All |  |  |  |  | ⚪ removed | ⚪ — (no 7d baseline) | 🟢 67·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 87% (26/30) | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operation |  | 3s | 22s | 26s | ⚪ 0% (26/30 → 26/30) | ⚪ — (no 7d baseline) | ⚪ 87·87·87 | ⚪ — | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (58/58) | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI interacti |  |  |  |  | ⚪ +0% (58/58 → 58/58) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | 🟡 93% (28/30) | All 14 sections of form-tests-pw.md exercised across 30 steps. 28 PASS, 2 AMBIGUOUS, 0 FAIL in MCP run. Playwright spec |  |  | 3m 12s | 3m 12s | 🟢 +1% (24/26 → 28/30) | ⚪ — (no 7d baseline) | 🔴 93·93·92 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas elemen |  |  | 52s | 52s | ⚪ +0% (30/30 → 36/36) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 73% (16/22) | Grid tests ran 22 steps (spec softSteps); 17 passed outright and 5 were AMBIGUOUS (Copy/Paste, Column Header Context Men | 11m | 3m | 1m 18s | 15m 18s | ⚪ 0% (16/22 → 16/22) | ⚪ — (no 7d baseline) | ⚪ 73·73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | 🟡 94% (15/16) | All 14 heat-map sections exercised across 17 steps. 15 PASS, 1 AMBIGUOUS, 1 SKIP in MCP run. Playwright spec passed full |  |  | 49s | 49s | 🟢 +1% (13/14 → 15/16) | ⚪ — (no 7d baseline) | 🔴 94·94·93 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (87/94) | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance, l | 50s | 7s | 46s | 1m 43s | ⚪ 0% (87/94 → 87/94) | ⚪ — (no 7d baseline) | ⚪ 93·93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | 🟢 PASS → PASS |  | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout sa | 57s | 8s | 1m 43s | 2m 48s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [map](Viewers/map-run.md) | 🟢 PASS → PASS |  | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker size m | 15s | 3s | 9s | 27s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 80·80·80 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL |  | Matrix Plot tests ran with 15 PASS, 3 AMBIGUOUS, 0 FAIL. The spec executed in 57.7s with all implemented steps passing. |  |  | 56s | 56s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 84·84·84 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (8/12) | 9 of 12 steps PASS; 3 SKIP (canvas-based node/edge interactions cannot be automated via DOM). The network diagram viewer | 8m | 1m 30s | 22s | 9m 52s | ⚪ 0% (8/12 → 8/12) | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | All 13 scenario sections (mapped to 12 Playwright softSteps — scale and normalization are combined in the spec) passed d | 1m 8s | 8s | 47s | 2m 3s | ⚪ +0% (13/13 → 13/13) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (81/81) | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appearan | 40s | 7s | 47s | 1m 34s | ⚪ +0% (81/81 → 81/81) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL |  | Pivot Table tests ran with 16 PASS, 2 AMBIGUOUS, 1 SKIP, 0 FAIL. The spec executed in 35.1s with all implemented steps p |  |  | 1m 12s | 1m 12s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 85·85·85 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| Viewers | [row-source](Viewers/row-source-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with all 8 |  | 5s | 1m 24s | 1m 29s | ⚪ +0% (36/36 → 36/36) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ — | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (20/20) | All 20 sections passed during the MCP run on dev.datagrok.ai. The existing Playwright spec was re-run headed without mod | 3m 13s | 29s | 52s | 4m 34s | ⚪ +0% (20/20 → 20/20) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | scatter-plot-tests | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [statistics](Viewers/statistics-run.md) | 🟢 PASS → PASS | 🟢 100% (23/23) | All 23 MCP steps passed. The date columns section (STARTED row behavior) was moved to `statistics-tests-ui.md` as a manu | 20m | 4m |  | 24m | ⚪ +0% (24/24 → 23/23) | ⚪ — (no 7d baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | 🟢 PASS → PASS | 🟢 100% (24/24) | 24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to manual |  | 3m | 58s | 3m 58s | ⚪ +0% (24/24 → 24/24) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ — | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | 🟡 PARTIAL → PARTIAL |  | All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via `v | 28m | 4m | 46s | 32m 46s | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (71/84) | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range slid | 3m | 30s |  | 3m 30s | 🟢 +1% (70/84 → 71/84) | ⚪ — (no 7d baseline) | ⚪ 84·84·84 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All 7 scenario steps PASS in both the MCP run against https://dev.datagrok.ai and the generated Playwright replay (27.2s | 2m 17s | 50s | 27s | 3m 34s | ⚪ +0% (7/7 → 7/7) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [word-cloud-tests](Viewers/word-cloud-tests-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (39/46) | All 46 scenario steps were exercised against dev.datagrok.ai. In the MCP run, 35 passed, 3 were AMBIGUOUS (visual effect | 4m 30s | 2m | 34s | 7m 4s | ⚪ 0% (39/46 → 39/46) | ⚪ — (no 7d baseline) | ⚪ 85·85 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [working-with-nan-infinity](Viewers/working-with-nan-infinity-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 spec steps PASSED in 1m 24s. NaN and Infinity values in numeric columns are handled gracefully across Scatter Plot |  |  | 1m 24s | 1m 24s | ⚪ +0% (9/9 → 9/9) | ⚪ — (no 7d baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| scatter-plot | axes-and-encoding | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| scatter-plot | regression-line-per-category | ⚪ NO RUN |  |  |  |  |  |  | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |

## Comparison with Previous Reports

Deltas are computed against `prev1d` (2026-04-24, 185 test entries) and `prev7d` (`null` — no commit fell in the 2026-04-25 → 2026-05-01 window). 7d-only metrics fall back to `prev7d_metrics` = the latest commit ≥4 days old that has Pass %/timing data (2026-04-24, the same as `prev1d` here).

### Totals

**Total (1d)**: Tests Δ **-4** · Run Δ **🟢 +7** · Mean Pass Δ **🟢 +3%** · Browser Δ **🔴 +36s** · Spec Gen Δ **🔴 +18s** · Spec Run Δ **🔴 +22s** · Total Δ **🔴 +1m 13s**

**Total (7d)**: same baseline as 1d (2026-04-24, no commit in the strict 7d±3 window). _(baseline 2026-04-24, ±3-day window vs target 2026-04-28 yielded no commit; using `prev7d_metrics` fallback.)_

### By Folder

| Folder | Tests Δ | Run Δ | Playwright Δ | Status | Mean Pass Δ (1d) | Mean Pass Δ (7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|---|
| Apps | +0 | ⚪ +0 (+0%) | — | ⚪ NO DATA → NO DATA | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | 🔴 -2% | 🔴 -2% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Browse | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | 🟡 (removed) | 🟡 (removed) | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | ⚪ 0% | ⚪ 0% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | ⚪ +0% | ⚪ +0% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ -0.03s |
| Connections | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | 🟢 +23% | 🟢 +23% (same baseline) | 🟡 4m 34s (new) | 🟡 2m 16s (new) | 🟡 2m 40s (new) | 🟡 9m 15s (new) |
| DiffStudio | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | ⚪ +0% | ⚪ +0% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | -5 | 🔴 -5 (-50%) | — | 🟡 PARTIAL → PARTIAL | 🔴 -8% | 🔴 -8% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ -0.01s |
| General | +0 | ⚪ +0 (+0%) | — | ⚪ NO DATA → NO DATA | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | 🔴 -6% | 🔴 -6% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Notebooks | +0 | ⚪ +0 (+0%) | — | ⚪ NO DATA → NO DATA | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | 🟡 (removed) | 🟡 (removed) | ⚪ +0s | ⚪ +0s | 🟢 -36s | 🟢 -42s |
| PowerPack | -1 | 🔴 -1 (-11%) | — | 🟡 PARTIAL → PARTIAL | ⚪ +0% | ⚪ +0% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ -0.02s |
| Projects | +2 | 🔴 -2 (-18%) | — | 🟡 PARTIAL → PARTIAL | 🟢 +2% | 🟢 +2% (same baseline) | 🟡 6m 38s (new) | 🟡 2m (new) | 🟡 1m 7s (new) | 🟡 9m 12s (new) |
| Queries | -1 | 🟢 +14 (+100%) | — | ⚪ → 🟡 NO DATA → PARTIAL | 🟡 93% (new) | 🟡 93% (new) | 🟡 7m 49s (new) | 🟡 2m 45s (new) | 🟡 3m 5s (new) | 🟡 13m 39s (new) |
| Scripts | +0 | 🟢 +1 (+17%) | — | 🟡 PARTIAL → PARTIAL | 🔴 -7% | 🔴 -7% (same baseline) | 🔴 +1m 59s | 🔴 +9s | 🔴 +16s | 🔴 +2m 24s |
| StickyMeta | +0 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | ⚪ +0% | ⚪ +0% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | +0 | ⚪ +0 (+0%) | — | ⚪ NO DATA → NO DATA | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | -1 | ⚪ +0 (+0%) | — | 🟡 PARTIAL → PARTIAL | ⚪ 0% | ⚪ 0% (same baseline) | ⚪ +0s | ⚪ +0s | ⚪ +0.01s | 🟢 -6s |
| scatter-plot | +2 (new) | ⚪ +0 (+0%) | — | 🟡 new: NO DATA | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |

### Per-Test Changes

| Folder | Test | Status | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|
| Apps | apps | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Apps | tutorials | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -17% (5/5 → 5/6) | ⚪ — (no 7d baseline) | 🟢 67·67·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Browse | package-manager | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (5/12 → 5/12) | ⚪ — (no 7d baseline) | ⚪ 42·42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (1/5 → 1/5) | ⚪ — (no 7d baseline) | ⚪ 20·20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/6 → 3/6) | ⚪ — (no 7d baseline) | ⚪ 50·50·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/3 → 1/3) | ⚪ — (no 7d baseline) | ⚪ 33·33·33 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/5 → 2/5) | ⚪ — (no 7d baseline) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [adding](Connections/adding-run.md) | 🟡 → 🟢 PARTIAL → PASS | ⚪ 0% (6/7 → 6.5/7) | ⚪ — (no 7d baseline) | ⚪ 93·93·93 | 🟡 4m 27s (new) | 🟡 45s (new) | 🟡 40s (new) | 🟡 5m 52s (new) |
| Connections | [browser](Connections/browser-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 +22% (7/9 → 9/9) | ⚪ — (no 7d baseline) | ⚪ 78·78·78 | 🟡 6m 56s (new) | 🟡 1m 36s (new) | 🟡 31s (new) | 🟡 9m 3s (new) |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | 🟢 +3% (1/11 → 2/17) | ⚪ — (no 7d baseline) | ⚪ 9·9·9 | 🟡 2m 35s (new) | 🟡 2m (new) | 🟡 45s (new) | 🟡 5m 20s (new) |
| Connections | [delete](Connections/delete-run.md) | 🟢 PASS → PASS | ⚪ +0% (8/8 → 8/8) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | 🟡 1m (new) | 🟡 2m (new) | 🟡 27s (new) | 🟡 3m 27s (new) |
| Connections | [edit](Connections/edit-run.md) | 🟢 PASS → PASS | ⚪ 0% (6/7 → 6/7) | ⚪ — (no 7d baseline) | ⚪ 86·86·86 | 🟡 2m (new) | 🟡 4m (new) | 🟡 1m 6s (new) | 🟡 7m 6s (new) |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 → 🟢 FAIL → PASS | 🟢 +100% (0/7 → 8/8) | ⚪ — (no 7d baseline) | ⚪ 0·0·0 | 🟡 5m 54s (new) | 🟡 3m 30s (new) | 🟡 14m 10s (new) | 🟡 23m 34s (new) |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 → 🟡 FAIL → PARTIAL | 🟢 +6% (1/9 → 2/12) | ⚪ — (no 7d baseline) | ⚪ 11·11·11 | 🟡 12m 36s (new) | 🟡 5m 30s (new) | ⚪ — | 🟡 18m 6s (new) |
| Connections | [import-swagger](Connections/import-swagger-run.md) | ⚪ removed: FAIL | 🟢 +100% (0/7 → 7/7) | ⚪ — (no 7d baseline) | ⚪ 0·0·0 | 🟡 7m 39s (new) | 🟡 58s (new) | 🟡 5m 19s (new) | 🟡 13m 56s (new) |
| Connections | [schema](Connections/schema-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 +25% (3/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 75·75·75 | 🟡 57s (new) | 🟡 25s (new) | 🟡 36s (new) | 🟡 1m 58s (new) |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -29% (6/7 → 4/7) | ⚪ — (no 7d baseline) | ⚪ 86·86·86 | 🟡 1m 40s (new) | 🟡 2m (new) | 🟡 29s (new) | 🟡 4m 9s (new) |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/6 → 5/6) | ⚪ — (no 7d baseline) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -17% (2/3 → 1/2) | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (3/7 → 3/7) | ⚪ — (no 7d baseline) | ⚪ 43·43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -10% (3/5 → 2/4) | ⚪ — (no 7d baseline) | ⚪ 60·60·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🔴 -12% (2/4 → 1.5/3) | ⚪ — (no 7d baseline) | ⚪ 62·50·62 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| General | files-cache | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | first-login | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | inactivity-response | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | login | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | molecule-in-exported-csv | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | network | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | profile-settings | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | startup-time | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | table-manager | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | tabs-reordering | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [apply](Models/apply-run.md) | 🔴 FAIL → FAIL | 🔴 -17% (2/4 → 1/3) | ⚪ — (no 7d baseline) | 🔴 100·100·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -13% (2/6 → 1/5) | ⚪ — (no 7d baseline) | 🔴 67·67·33 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🔴 -14% (6/17 → 3/14) | ⚪ — (no 7d baseline) | 🟢 28·28·35 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [delete](Models/delete-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/5 → 1/5) | ⚪ — (no 7d baseline) | 🔴 100·100·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [train](Models/train-run.md) | 🟢 PASS → PASS | 🟢 +9% (10/11 → 10/10) | ⚪ — (no 7d baseline) | 🟢 90·90·91 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Notebooks | browser | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | create | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | delete | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | edit | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | 🟡 1m 18s (removed) | 🟢 -1m 18s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ removed | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | 🟡 1m 30s (removed) | 🟢 -1m 30s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | ⚪ removed: PARTIAL | ⚪ +0% (16/21 → 16/21) | ⚪ — (no 7d baseline) | 🟢 0·76·76 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 +50% (2/4 → 4/4) | ⚪ — (no 7d baseline) | ⚪ 50·50·50 | 🟡 3m 17s (new) | 🟡 2m (new) | ⚪ — | 🟡 5m 17s (new) |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -20% (8/14 → 6/16) | ⚪ — (no 7d baseline) | ⚪ 57·57·57 | 🟡 10m (new) | 🟡 2m (new) | 🟡 1m 7s (new) | 🟡 13m 7s (new) |
| Queries | [adding](Queries/adding-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 34s (new) | 🟡 50s (new) | 🟡 21s (new) | 🟡 1m 45s (new) |
| Queries | [browse-and-save-project](Queries/browse-and-save-project-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 11m 30s (new) | 🟡 50s (new) | 🟡 1m 30s (new) | 🟡 13m 50s (new) |
| Queries | [browser](Queries/browser-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 1m 50s (new) | 🟡 30s (new) | 🟡 24s (new) | 🟡 2m 44s (new) |
| Queries | [columns-inspect](Queries/columns-inspect-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 4m 43s (new) | 🟡 1m 30s (new) | 🟡 2m 45s (new) | 🟡 8m 58s (new) |
| Queries | [deleting](Queries/deleting-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 55s (new) | 🟡 30s (new) | 🟡 25s (new) | 🟡 1m 50s (new) |
| Queries | [edit](Queries/edit-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 1m (new) | 🟡 40s (new) | 🟡 3m (new) | 🟡 4m 40s (new) |
| Queries | [get-all-get-top-100](Queries/get-all-get-top-100-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 2m 5s (new) | 🟡 4m 40s (new) | 🟡 31s (new) | 🟡 7m 16s (new) |
| Queries | [ms-sql](Queries/ms-sql-run.md) | 🟡 new: PARTIAL | 🟡 60% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 7m 40s (new) | 🟡 5m 30s (new) | 🟡 1m 10s (new) | 🟡 14m 20s (new) |
| Queries | [new-sql-query](Queries/new-sql-query-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 1m 35s (new) | 🟡 30s (new) | 🟡 28s (new) | 🟡 2m 33s (new) |
| Queries | [new-visual-query](Queries/new-visual-query-run.md) | 🟡 new: PARTIAL | 🟡 50% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 20m (new) | 🟡 2m (new) | 🟡 56s (new) | 🟡 22m 56s (new) |
| Queries | [query-layout](Queries/query-layout-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 6m 20s (new) | 🟡 12m (new) | 🟡 1m 18s (new) | 🟡 19m 38s (new) |
| Queries | [query-postprocessing](Queries/query-postprocessing-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 5m 16s (new) | 🟡 3m 12s (new) | 🟡 44s (new) | 🟡 9m 12s (new) |
| Queries | [transformations](Queries/transformations-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 7m 53s (new) | 🟡 54s (new) | 🟡 28m 30s (new) | 🟡 37m 17s (new) |
| Queries | [visual-query-advanced](Queries/visual-query-advanced-run.md) | 🟡 new: PARTIAL | 🟡 91% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 38m (new) | 🟡 5m (new) | 🟡 1m 6s (new) | 🟡 44m 6s (new) |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (10/12 → 12/13) | ⚪ — (no 7d baseline) | ⚪ 92·92·92 | 🟡 8m 39s (new) | 🟡 50s (new) | 🟡 1m 14s (new) | 🟡 10m 43s (new) |
| Scripts | [delete](Scripts/delete-run.md) | 🟢 PASS → PASS | ⚪ +0% (5/5 → 5/5) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | 🟡 57s (new) | 🟡 25s (new) | 🟡 22s (new) | 🟡 1m 44s (new) |
| Scripts | [edit](Scripts/edit-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) | ⚪ — (no 7d baseline) | ⚪ 100·100·100 | 🟡 1m 25s (new) | 🟡 30s (new) | 🟡 21s (new) | 🟡 2m 16s (new) |
| Scripts | [layout](Scripts/layout-run.md) | 🟡 new: PARTIAL | 🟡 46% (baseline) | ⚪ — (no 7d baseline) |  | 🟡 5m 9s (new) | 🟡 35s (new) | 🟡 40s (new) | 🟡 6m 24s (new) |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (6/9 → 6/9) | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | 🟡 2m 10s (new) | 🟡 40s (new) | 🟡 19s (new) | 🟡 3m 9s (new) |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/11 → 1/11) | ⚪ — (no 7d baseline) | 🔴 20·20·9 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | actions-in-the-context-menu | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip-visibility | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | edit-tooltip | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | line-chart-aggregated-tooltip | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | tooltip-properties | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | uniform-default-tooltip | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (5/8 → 5.5/8) | ⚪ — (no 7d baseline) | ⚪ 69·69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (9/12 → 10/12) | ⚪ — (no 7d baseline) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/11 → 8.5/11) | ⚪ — (no 7d baseline) | ⚪ 77·77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (7/13 → 7.5/13) | ⚪ — (no 7d baseline) | ⚪ 58·58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (5/7 → 5.5/7) | ⚪ — (no 7d baseline) | ⚪ 79·79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (13/20 → 14/20) | ⚪ — (no 7d baseline) | ⚪ 70·70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | ⚪ removed: PASS | 🟢 +8% (11/13 → 28/30) | ⚪ — (no 7d baseline) | 🟢 77·86·85 | ⚪ — | ⚪ — | 🟡 49s (new) | 🟡 49s (new) |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (18/19 → 18/19) | ⚪ — (no 7d baseline) | ⚪ 95·95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (26/30 → 26/30) | ⚪ — (no 7d baseline) | ⚪ 87·87·87 | ⚪ — | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | 🟢 +1% (24/26 → 28/30) | ⚪ — (no 7d baseline) | 🔴 93·93·92 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (16/22 → 16/22) | ⚪ — (no 7d baseline) | ⚪ 73·73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | 🟢 +1% (13/14 → 15/16) | ⚪ — (no 7d baseline) | 🔴 94·94·93 | ⚪ — | ⚪ — | ⚪ +0s | ⚪ +0s |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (87/94 → 87/94) | ⚪ — (no 7d baseline) | ⚪ 93·93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (8/12 → 8/12) | ⚪ — (no 7d baseline) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | scatter-plot-tests | ⚪ NO RUN | ⚪ — | ⚪ — (no 7d baseline) |  | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 +1% (70/84 → 71/84) | ⚪ — (no 7d baseline) | ⚪ 84·84·84 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Viewers | [word-cloud-tests](Viewers/word-cloud-tests-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ 0% (39/46 → 39/46) | ⚪ — (no 7d baseline) | ⚪ 85·85 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage at 83% with no FAIL folders, but 15 folders are PARTIAL. Review their failing/skipped steps before release.

### PARTIAL folders
- **Bio** — 9/9 runs, mean Pass 🟢 98%
- **Browse** — 2/3 runs, mean Pass 
- **Charts** — 3/3 runs, mean Pass 🟡 54%
- **Chem** — 14/14 runs, mean Pass 🟡 87%
- **Connections** — 10/10 runs, mean Pass 🟡 76%
- **DiffStudio** — 8/8 runs, mean Pass 🟢 96%
- **EDA** — 10/10 runs, mean Pass 🟡 59%
- **Models** — 6/6 runs, mean Pass 🟡 49%
- **Peptides** — 4/4 runs, mean Pass 
- **PowerPack** — 9/9 runs, mean Pass 🟢 97%
- **Projects** — 7/11 runs, mean Pass 🟡 40%
- **Queries** — 14/14 runs, mean Pass 🟡 93%
- **Scripts** — 6/6 runs, mean Pass 🟡 80%
- **StickyMeta** — 4/4 runs, mean Pass 🟡 77%
- **Viewers** — 44/45 runs, mean Pass 🟡 92%