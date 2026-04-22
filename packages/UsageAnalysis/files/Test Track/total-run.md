# Test Track — Global Report

**Date**: 2026-04-23
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

**Total**: 177 tests · Run: 135/177 (76%) · Playwright: 111/177 (63%) · Mean Pass: 🟡 81% · Mean Browser: 3m 57s · Mean Spec Gen: 57.6s · Mean Spec Run: 41.3s · Mean Total (sum per test): 5m 22s

| Folder | Tests | Run | Playwright | Status | Mean Pass % | Mean Browser | Mean Spec Gen | Mean Spec Run | Mean Total |
|---|---|---|---|---|---|---|---|---|---|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Bio | 9 | 9/9 (100%) | 1/9 (11%) | 🟡 PARTIAL | 🟡 94% | 2m 45s | 5s |  | 2m 46s |
| Browse | 3 | 2/3 (67%) | 0/3 (0%) | 🟡 PARTIAL | 🟡 70% |  |  |  |  |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | 🟡 PARTIAL | 🟡 54% | 2m 10s | 45s | 34.3s | 3m 29s |
| Chem | 14 | 14/14 (100%) | 14/14 (100%) | 🟡 PARTIAL | 🟡 87% | 1m 37s | 31.8s | 40.2s | 2m 49s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | 🟡 PARTIAL | 🟡 54% |  |  |  |  |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | 🟡 PARTIAL | 🟢 96% | 3m 10s | 56.4s | 45.4s | 4m 52s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | 🟡 PARTIAL | 🟡 66% | 1m 42s | 58s | 17.5s | 2m 58s |
| General | 10 | 0/10 (0%) | 0/10 (0%) | ⚪ NO DATA |  |  |  |  |  |
| LocalCashing | 0 | 0/0 | 0/0 | ⚪ NO DATA |  |  |  |  |  |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | 🟡 PARTIAL | 🟡 81% |  |  |  |  |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 68% | 27.5s | 3s | 11.2s | 36.1s |
| PowerPack | 9 | 9/9 (100%) | 9/9 (100%) | 🟡 PARTIAL | 🟢 97% | 5m 16s | 1m 9s | 29.8s | 6m 54s |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | 🟡 PARTIAL | 🟡 38% |  |  |  |  |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | 🟡 PARTIAL | 🟡 89% |  |  |  |  |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 68% | 23.8s | 2.2s | 15s | 37.2s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Viewers | 46 | 43/46 (93%) | 43/46 (93%) | 🟡 PARTIAL | 🟡 92% | 6m 11s | 1m 17s | 53s | 7m 52s |

## All Tests

**Total**: 177 tests · 🟢 76 PASS / 🟡 45 PARTIAL / 🔴 9 FAIL / 🟡 1 AMBIGUOUS / 🟡 4 SKIP / ⚪ 42 NO RUN · Mean Pass: 🟡 81% · Mean Browser: 3m 57s · Mean Spec Gen: 57.6s · Mean Spec Run: 41.3s · Mean Total (sum per test): 5m 22s

| Folder | Test | Status | Pass % | Description | Browser (model+MCP) | Spec Gen (model) | Spec Run (Playwright) | Total (sum) | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Apps | apps | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Apps | tutorials | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [analyze](Bio/analyze-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset typ… | 8m | 5s |  | 8m 5s | ⚪ +0% (11/11 → 11/11) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 80% (4/5) | 4 of 5 steps passed. The Composition/WebLogo viewer opens correctly and properties are accessible. The letter-click sele… |  |  |  |  | ⚪ +0% (4/5 → 4/5) | 🟡 80% (baseline) | ⚪ 80·80 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [convert](Bio/convert-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 4 Bio convert/transform functions work correctly on FASTA data. Get Region extracts a subsequence region. PolyTool >… | 2m |  |  | 2m | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ — | ⚪ — | ⚪ +0s |
| Bio | [manage](Bio/manage-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All 3 steps passed. The Manage Monomer Libraries view opens as a full view showing 5 monomer library JSON files (increas… | 30s |  |  | 30s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ — | ⚪ — | ⚪ +0s |
| Bio | [msa](Bio/msa-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. MSA dialog opens with correct fields, Alignment Parameters button adds gap penalty inputs as expecte… |  |  |  |  | ⚪ +0% (6/6 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 AMBIGUOUS → AMBIGUOUS | 🟡 67% (4/6) | The MSA dialog opens correctly for HELM data and shows MAFFT-based method options (mafft --auto, linsi, ginsi, etc.) ins… |  |  |  |  | ⚪ +0% (4/6 → 4/6) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [search](Bio/search-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 4 steps passed. Bio > Search > Subsequence Search opens a filter panel with a Sequence bio substructure filter. Typi… | 30s |  |  | 30s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ — | ⚪ — | ⚪ +0s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. Activity Cliffs works correctly with both default parameters (UMAP/Hamming) and custom parameters (t… |  |  |  |  | ⚪ +0% (6/6 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [sequence-space](Bio/sequence-space-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. Sequence Space works correctly with both default (UMAP/Hamming) and custom (t-SNE/Needlemann-Wunsch)… |  |  |  |  | ⚪ +0% (6/6 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 90% (4/5) | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. Item… |  |  |  |  | ⚪ +0% (4/5 → 4/5) | 🟡 90% (baseline) | ⚪ 90·90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (0/1) | 1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a single session… |  |  |  |  | ⚪ +0% (0/1 → 0/1) | 🟡 50% (baseline) | ⚪ 50·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | package-manager | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | japanese-in-myfiles | 🟢 removed: PASS |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | local-deploy | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces | 🟡 removed: PARTIAL |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces-(ui-only) | 🔴 removed: FAIL |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [radar](Charts/radar-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | The Radar viewer reproduced cleanly on dev for both earthquakes.csv (2426 rows) and demog.csv (5850 rows). All 21 Radar … | 1m 30s | 35s | 33s | 2m 38s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 42% (5/12) | The sunburst viewer reproduces structurally on dev — Sunburst can be added to both SPGI and demog, and `hierarchyColumnN… | 3m | 1m | 34s | 4m 34s | ⚪ +0% (5/12 → 5/12) | 🟡 42% (baseline) | ⚪ 42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 20% (1/5) | Setup (open demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduced cleanly on dev. All four test steps are mar… | 2m | 40s | 36s | 3m 16s | ⚪ +0% (1/5 → 1/5) | 🟡 20% (baseline) | ⚪ 20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (3/6) | Smoke coverage only: Scaffold Tree viewer launches from the Chem menu and the magic wand generates a scaffold tree on SP… | 57s | 25s | 49.6s | 2m 12s | ⚪ +0% (3/6 → 3/6) | 🟡 50% (baseline) | ⚪ 50·50 | 🔴 +17s | ⚪ +0s | 🟢 -0.4s | 🔴 +16.6s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Scaffold Tree viewer launches from the Chem → Scaffold Tree menu, and the magic-wand generator produces scaffold nodes f… | 1m 15s | 30s | 40.2s | 2m 25s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.2s | 🔴 +0.2s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Similarity Search launches from the Chem menu and exposes a viewer that accepts option changes (fingerprint Morgan ↔ Pat… | 37s | 25s | 24.8s | 1m 27s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +12s | ⚪ +0s | 🔴 +1.8s | 🔴 +13.8s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Substructure filtering via `grok.chem.searchSubstructure` works on SPGI.csv (3624 rows): benzene substructure yields a b… | 38s | 25s | 29s | 1m 32s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +8s | ⚪ +0s | 🔴 +5s | 🔴 +13s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Activity Cliffs computation on SPGI.csv (3624 rows) finishes within 45s and produces a UMAP scatter plot with molecule t… | 1m 14s | 20s | 1m 4s | 2m 38s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +44s | ⚪ +0s | 🔴 +4s | 🔴 +48s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) | Calculate Descriptors cannot be exercised on `dev` right now. The Chem top menu fails to open its popup — both through D… | 8m | 1m | 38s | 9m 38s | ⚪ +0% (1/3 → 1/3) | 🟡 33% (baseline) | ⚪ 33·33 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Chemical Space dimensional reduction runs end-to-end on smiles.csv: the dialog opens, OK with defaults produces a Scatte… | 46s | 20s | 59.8s | 2m 6s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +26s | ⚪ +0s | 🔴 +3.8s | 🔴 +29.8s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | ChemProp scenario is partially automated: the spec confirms mol1K.sdf opens and the Train Model view is reachable from t… | 26s | 30s | 19.1s | 1m 15s | ⚪ +0% (2/5 → 2/5) | 🟡 40% (baseline) | ⚪ 40·40 | 🟢 -34s | ⚪ +0s | 🔴 +1.1s | 🟢 -32.9s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Elemental Analysis works on dev. The menu path `[name="div-Chem"]` → `Elemental Analysis...` resolves and the dialog ope… | 38s | 30s | 30.9s | 1m 39s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🟢 -22s | ⚪ +0s | 🔴 +1.9s | 🟢 -20.1s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | 🟢 PASS → PASS | 🟢 100% (2/2) | The filter panel correctly shows a Structure filter for SPGI.csv's Molecule column; clicking the embedded sketch-link op… | 34s | 25s | 24.2s | 1m 23s | ⚪ +0% (2/2 → 2/2) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +4s | ⚪ +0s | 🔴 +3.2s | 🔴 +7.2s |
| Chem | [info-panels](Chem/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | Info panels work correctly on smiles.csv: column-level (Details, Filter, Colors, Style, Chemistry with Rendering/Highlig… | 3m 10s | 1m | 31s | 4m 41s | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [mmp](Chem/mmp-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset at the bottom of the view… | 1m 27s | 20s | 1m 14s | 3m 1s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +57s | ⚪ +0s | 🔴 +2s | 🔴 +59s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | R-Groups Analysis works on sar_small.csv: MCS auto-populates the sketcher, OK produces a Trellis plot and appends R1–R4 … | 2m 20s | 45s | 58.7s | 4m 4s | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.3s | 🟢 -0.3s |
| Chem | [sketcher](Chem/sketcher-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Sketcher opens via `grok.chem.sketcher(molCol, initialSmiles)` wrapped in `ui.dialog(...).show()`, accepts a typed SMILE… | 33s | 30s | 19.3s | 1m 22s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +3s | ⚪ +0s | 🔴 +3.3s | 🔴 +6.3s |
| Connections | [adding](Connections/adding-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (6/7) | 6 of 7 steps fully passed, Step 5 was partial (TEST button works but actual connection test fails without real credentia… |  |  |  |  | ⚪ +0% (6/7 → 6/7) | 🟡 93% (baseline) | ⚪ 93·93 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [browser](Connections/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 78% (7/9) | 7 of 9 steps passed, 2 ambiguous. Search filtering works correctly and the Context Pane shows all expected tabs (Details… |  |  |  |  | ⚪ +0% (7/9 → 7/9) | 🟡 78% (baseline) | ⚪ 78·78 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | 🟡 9% (1/11) | 1 step passed, 1 failed, 15 skipped. The required `NorthwindTest` MS SQL connection is not present on public.datagrok.ai… |  |  |  |  | ⚪ +0% (1/11 → 1/11) | 🟡 9% (baseline) | ⚪ 9·9 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [delete](Connections/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 steps passed. Both connections were deleted successfully. The confirmation dialog uses a red "DELETE" button (not … |  |  |  |  | ⚪ +0% (8/8 → 8/8) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [edit](Connections/edit-run.md) | 🟢 PASS → PASS | 🟡 86% (6/7) | 6 of 7 steps passed (1 skipped due to missing real credentials). The connection rename, credential modification, and err… |  |  |  |  | ⚪ +0% (6/7 → 6/7) | 🟡 86% (baseline) | ⚪ 86·86 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 FAIL → FAIL | 🔴 0% (0/7) | All 7 steps skipped. This scenario requires a specific Postgres connection at db.datagrok.ai:54327 with superuser creden… |  |  |  |  | ⚪ +0% (0/7 → 0/7, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 FAIL → FAIL | 🟡 11% (1/9) | 1 step passed, 1 failed, 7 skipped. This scenario depends on a working Postgres connection to the Northwind database. Th… |  |  |  |  | ⚪ +0% (1/9 → 1/9) | 🟡 11% (baseline) | ⚪ 11·11 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [import-swagger](Connections/import-swagger-run.md) | 🔴 FAIL → FAIL | 🔴 0% (0/7) | All 7 steps skipped. This scenario requires manual interaction: downloading a YAML file to the local machine and drag-dr… |  |  |  |  | ⚪ +0% (0/7 → 0/7, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [schema](Connections/schema-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 75% (3/4) | 3 of 4 steps passed, 1 ambiguous. The "Browse schema" context menu option was not found in the current UI, but the schem… |  |  |  |  | ⚪ +0% (3/4 → 3/4) | 🟡 75% (baseline) | ⚪ 75·75 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 86% (6/7) | 6 of 7 steps passed (1 failed). All UI steps worked correctly. The SPARQL connection was created and deleted successfull… |  |  |  |  | ⚪ +0% (6/7 → 6/7) | 🟡 86% (baseline) | ⚪ 86·86 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | Catalog scenario reproduces fully on dev.datagrok.ai. All 6 steps PASS both in MCP and in the Playwright spec (35.8s wal… | 1m 53s | 39s | 38s | 3m 10s | ⚪ +0% (6/6 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Cyclic Models (PK-PD) scenario reproduces fully on dev.datagrok.ai. The PK-PD library model loads via double-click, Mult… | 1m 2s | 38s | 36s | 2m 16s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Files & Sharing scenario reproduces fully on dev.datagrok.ai. pk.ivp loads via the `DiffStudio:previewIvp` function with… | 2m 31s | 1m 32s | 1m 17s | 5m 20s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (5/6) | Scenario is PARTIAL on dev.datagrok.ai — steps 1–5 pass, step 6 (actually running the fit) does not produce result rows … | 10m 2s | 55s | 1m | 11m 57s | ⚪ +0% (5/6 → 5/6) | 🟡 83% (baseline) | ⚪ 83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | 🟢 PASS → PASS | 🟡 83% (5/6) | Scenario fully reproduces on dev.datagrok.ai. All 6 steps pass in the interactive MCP session and in the Playwright spec… | 1m 31s | 36s | 26s | 2m 33s | ⚪ +0% (5/6 → 5/6) | 🟡 83% (baseline) | ⚪ 83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 scenario steps PASS against dev.datagrok.ai. Edit toggle is reachable via `.d4-ribbon-item .ui-input-bool-switch .… | 4m 30s | 1m 40s | 1m 3s | 7m 13s | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Sensitivity Analysis scenario fully reproduces against dev.datagrok.ai. Bioreactor loads from the DiffStudio hub (librar… | 2m 3s | 44s | 39s | 3m 26s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Stages (Acid Production) scenario reproduces fully on dev.datagrok.ai. The library card opens a view named "GA-productio… | 1m 49s | 47s | 24s | 3m | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [ML methods/linear-regression](EDA/ML methods/linear-regression-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) |  | 45s | 2s | 6.8s | 53.8s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🟢 -0.2s |
| EDA | [ML methods/pls-regression](EDA/ML methods/pls-regression-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 71% (5/7) |  | 1m | 2s | 6.9s | 1m 9s | ⚪ +0% (5/7 → 5/7) | 🟡 71% (baseline) | ⚪ 71·71 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🟢 -0.1s |
| EDA | [ML methods/softmax](EDA/ML methods/softmax-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) |  | 10s | 2s | 2.6s | 14.6s | ⚪ +0% (1/3 → 1/3) | 🟡 33% (baseline) | ⚪ 33·33 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🟢 -0.4s |
| EDA | [ML methods/xgboost1](EDA/ML methods/xgboost1-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) |  | 5s | 2s | 2.7s | 9.7s | ⚪ +0% (2/3 → 2/3) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [ML methods/xgboost2](EDA/ML methods/xgboost2-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) |  | 5s | 2s | 2.7s | 9.7s | ⚪ +0% (2/3 → 2/3) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [anova](EDA/anova-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All 3 scenario steps passed against dev. Dataset opens via JS API in ~1s; ANOVA dialog mounts with sensible defaults (RA… | 1m 30s | 30s | 26s | 2m 26s | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | 2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check, where the wording does not spec… | 2m 30s | 2m | 13s | 4m 43s | ⚪ +0% (2/3 → 2/3) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 43% (3/7) | 3 of 7 steps passed, 1 failed, 3 were skipped due to the missing prerequisite dataset. The Pareto Front viewer itself is… | 4m | 2m | 32s | 6m 32s | ⚪ +0% (3/7 → 3/7) | 🟡 43% (baseline) | ⚪ 43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 60% (3/5) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 3 PASS / 1 FAIL / 1 SKIP. The dialog path works (menu, F… | 5m | 3m | 1m 7s | 9m 7s | ⚪ +0% (3/5 → 3/5) | 🟡 60% (baseline) | ⚪ 60·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🟡 50% (2/5) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 2 PASS / 1 PARTIAL / 1 FAIL. The dialog path (menu, Usin… | 2m | 2m | 15s | 4m 15s | 🔴 -12% (2/4 → 2/5) | 🟡 50% (baseline) | 🔴 62·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| General | files-cache | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | first-login | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | inactivity-response | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | login | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | molecule-in-exported-csv | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | network | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | profile-settings | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | startup-time | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | table-manager | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | tabs-reordering | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | api-samples | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| LocalCashing | local-cashing | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [apply](Models/apply-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 4 steps passed. The Apply Model workflow functions correctly end-to-end. The TestDemog model (trained in Train.md) w… |  |  |  |  | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (4/6) | 4 of 6 steps passed; 2 skipped because only 1 model was available (TestDemog was the only model — the second numeric mod… |  |  |  |  | ⚪ +0% (4/6 → 4/6) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🟡 28% (5/18) | 2 of 17 sub-steps passed, 13 skipped, 1 failed, 1 ambiguous. The scenario fails entirely due to the Chemprop Docker cont… |  |  |  |  | ⚪ +0% (5/18 → 5/18) | 🟡 28% (baseline) | ⚪ 28·28 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [delete](Models/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps passed. The model deletion workflow works correctly end-to-end. The right-click context menu, confirmation d… |  |  |  |  | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [predictive-models](Models/predictive-models-run.md) | 🟢 PASS → PASS | 🟢 100% (21/21) | All 20 sub-steps passed. The full lifecycle (Train → Apply → Apply on new dataset → Delete) for EDA-based predictive mod… |  |  |  |  | ⚪ +0% (21/21 → 21/21) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [train](Models/train-run.md) | 🟢 PASS → PASS | 🟡 90% (9/10) | All 10 steps passed. The Train Model workflow functions correctly on public.datagrok.ai using the built-in EDA engines. … |  |  |  |  | ⚪ +0% (9/10 → 9/10) | 🟡 90% (baseline) | ⚪ 90·90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | browser | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | create | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | delete | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | edit | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [info-panels](Peptides/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are rende… | 17s | 3s | 10.9s | 30.9s | ⚪ +0% (6/6 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variability… | 25s | 3s |  | 28s | ⚪ +0% (2/5 → 2/5) | 🟡 40% (baseline) | ⚪ 40·40 | ⚪ +0s | ⚪ +0s | 🟡 1m 18s (removed) | 🟢 -1m 18s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (4/6) | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clusters… | 18s | 3s | 11.6s | 32.6s | ⚪ +0% (4/6 → 4/6) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | 🟢 -0.4s | 🟢 -0.4s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Residu… | 50s | 3s |  | 53s | ⚪ +0% (6/9 → 6/9) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | 🟡 1m 30s (removed) | 🟢 -1m 30s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (10/10) | All 10 scenario steps reproduce successfully in the MCP run; 9/10 pass in the Playwright replay. The one FAILED Playwrig… | 7m 45s | 2m | 32s | 10m 17s | 🟢 +40% (6/10 → 10/10) | 🟢 100% (baseline) | 🟢 60·100 | 🟡 7m 45s (new) | 🟡 2m (new) | 🟡 32s (new) | 🟡 10m 17s (new) |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All six autocomplete behaviours PASS in both the MCP run and the Playwright replay. `.cm-tooltip-autocomplete` appears o… | 1m 45s | 20s | 8s | 2m 13s | ⚪ +0% (8/8 → 7/7) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 1m 45s (new) | 🟡 20s (new) | 🟡 8s (new) | 🟡 2m 13s (new) |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All seven sub-steps pass in both the MCP run and the Playwright replay. Dependency propagation across calculated columns… | 2m 15s | 1m | 28s | 3m 43s | ⚪ +0% (5/5 → 7/7) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 2m 15s (new) | 🟡 1m (new) | 🟡 28s (new) | 🟡 3m 43s (new) |
| PowerPack | [AddNewColumn/functions-sorting](PowerPack/AddNewColumn/functions-sorting-run.md) | 🟢 new: PASS | 🟢 100% (7/7) | All five scenario steps PASS in both the MCP run and the Playwright replay (17s, 1 test, 0 failures). The previous run's… | 8m 20s | 2m | 17s | 10m 37s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 8m 20s (new) | 🟡 2m (new) | 🟡 17s (new) | 🟡 10m 37s (new) |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps PASS in both the MCP run and the Playwright replay. `${AGE}`, `$[AGE]` and the autocomplete-inserted `${HEIG… | 4m 43s | 1m 13s | 25s | 6m 21s | ⚪ +0% (4/4 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 4m 43s (new) | 🟡 1m 13s (new) | 🟡 25s (new) | 🟡 6m 21s (new) |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All four steps pass in both MCP and Playwright. The CodeMirror formula editor shows a `.cm-tooltip-hover` on hover with … | 1m 10s | 15s | 8s | 1m 33s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 1m 10s (new) | 🟡 15s (new) | 🟡 8s (new) | 🟡 1m 33s (new) |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | 🟡 → 🟢 SKIP → PASS | 🟢 100% (10/10) | All 10 scenario steps pass end-to-end — both in interactive MCP driving and in the Playwright replay (18.8s). A single i… | 3m 30s | 1m | 18.8s | 4m 49s | 🟢 +100% (0/6 → 10/10) | 🟢 100% (baseline) | 🟢 0·100 | 🟡 3m 30s (new) | 🟡 1m (new) | 🟡 18.8s (new) | 🟡 4m 49s (new) |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All six scenario steps PASS in both the MCP-driven grok-browser run and the Playwright replay (existing spec — not overw… | 1m 54s | 31s | 11s | 2m 36s | ⚪ +0% (5/5 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +1m 32s | 🔴 +29s | 🔴 +5.7s | 🔴 +2m 7s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🔴 → 🟡 FAIL → PARTIAL | 🟡 76% (16/21) | The PowerPack "Enrich column" feature works end-to-end for the primary create/apply/edit/delete flow on a dataframe that… | 16m | 2m | 2m | 20m | 🟢 +76% (0/10 → 16/21) | 🟡 76% (baseline) | 🟢 0·76 | 🔴 +8m 50s | 🔴 +45s | 🔴 +1m 17s | 🔴 +10m 52s |
| PowerPack | AddNewColumn/functions_sorting | 🟡 removed: SKIP |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | formula-lines | 🟡 removed: SKIP |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [browser](Projects/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 56% (5/9) | 5 of 9 steps passed, 3 skipped, 1 ambiguous. Browse > Dashboards view works correctly: projects are listed, searchable, … |  |  |  |  | ⚪ +0% (5/9 → 5/9) | 🟡 56% (baseline) | ⚪ 56·56 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [complex](Projects/complex-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/13) | All 13 steps skipped. This is the most complex scenario requiring tables from 7+ different sources, drag-and-drop, entit… |  |  |  |  | ⚪ +0% (0/13 → 0/13, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/5) | All 5 steps skipped. This scenario requires running a custom JavaScript script with Data Sync enabled, then modifying fi… |  |  |  |  | ⚪ +0% (0/5 → 0/5, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (2/4) | 2 of 4 steps passed, 1 skipped, 1 ambiguous. Project deletion works via the API (`grok.dapi.projects.delete()`). The rig… |  |  |  |  | ⚪ +0% (2/4 → 2/4) | 🟡 50% (baseline) | ⚪ 50·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [opening](Projects/opening-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (5/5) | All 5 steps passed. Projects from the Uploading step are accessible in Browse > Dashboards. Context Panel correctly show… |  |  |  |  | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/4) | All steps skipped. This scenario depends on Projects copy_clone.md (order 5) which was not fully executed. The Link/Clon… |  |  |  |  | ⚪ +0% (0/4 → 0/4, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | 🟡 SKIP → SKIP | 🟡 40% (2/5) | 2 of 5 steps passed, 3 skipped. Project preview and opening work. Copy/clone/link operations were not tested because the… |  |  |  |  | ⚪ +0% (2/5 → 2/5) | 🟡 40% (baseline) | ⚪ 40·40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 57% (8/14) | 8 of 14 steps passed, 6 skipped. Core project creation from local tables, file shares, query results, and join results a… |  |  |  |  | ⚪ +0% (8/14 → 8/14) | 🟡 57% (baseline) | ⚪ 57·57 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | adding | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | browse-&-save-project | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | browser | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | columns-inspect | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | deleting | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | edit | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | get-all-get-top-100 | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | ms-sql | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | new-sql-query | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | new-visual-query | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | query-layout | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | query-postprocessing | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | transformations | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | visual-query-advanced | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [browser](Scripts/browser-run.md) | 🟢 PASS → PASS | 🟡 89% (8/9) | The Scripts Browser scenario passed well. The context pane shows all expected accordions (Details, Script, Run, Activity… |  |  |  |  | ⚪ +0% (8/9 → 8/9) | 🟡 89% (baseline) | ⚪ 89·89 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 92% (10/12) | The Create scenario completed successfully overall. The script `testRscript` was created, parameters configured, saved, … |  |  |  |  | ⚪ +0% (10/12 → 10/12) | 🟡 92% (baseline) | ⚪ 92·92 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [delete](Scripts/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps passed. The delete flow works correctly with a confirmation dialog and immediate removal from the scripts li… |  |  |  |  | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [edit](Scripts/edit-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. The Edit scenario works correctly — edits are saved persistently and visible on re-open. |  |  |  |  | ⚪ +0% (6/6 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | layout | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Core run functionality works: the script can be triggered from context menu with a table selection and from the console.… |  |  |  |  | ⚪ +0% (6/9 → 6/9) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All tested steps passed. SPGI.csv opened with TestSchema1 sticky metadata schema pre-configured. The Sticky meta panel i… | 35s | 3s | 17s | 55s | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | 🟢 PASS → PASS | 🟡 50% (2/4) | Steps 1-2 passed: SPGI.csv opened with TestSchema1 sticky metadata schema, and cloning the table preserves the schema an… | 25s | 3s | 21s | 49s | ⚪ +0% (2/4 → 2/4) | 🟡 50% (baseline) | ⚪ 50·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All steps passed. The Sticky Meta Schemas browser at `/meta/schemas` shows 20 schemas including TestSchema1. The "NEW SC… | 10s | 3s | 7s | 20s | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | 🟡 20% (1/5) | Step 1 passed: navigated to Databases > Postgres and found CHEMBL connection. Step 2 failed: the "Database meta" section… | 25s | 0s |  | 25s | ⚪ +0% (1/5 → 1/5) | 🟡 20% (baseline) | ⚪ 20·20 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Tooltips | actions-in-the-context-menu | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip-visibility | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | edit-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | line-chart---aggregated-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | tooltip-properties | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | uniform-default-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 14 steps (setup + 13 scenario sections) passed in both the browser-driven MCP run against https://dev.datagrok.ai an… | 7m | 2m | 20s | 9m 20s | ⚪ +0% (15/15 → 15/15) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -50s | 🟢 -50s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | 🟢 PASS → PASS | 🟢 100% (26/26) | Ran basic-operations end-to-end against dev. All 31 scenario steps passed in the MCP browser phase (Section 1: structure… | 4m 27s | 9s | 50s | 5m 26s | ⚪ +0% (13/13 → 26/26) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | 🟢 PASS → PASS | 🟢 100% (16/16) | Ran chem-and-bio scenario end-to-end against dev. All 11 scenario steps passed in the MCP browser phase (Chem: open spgi… | 2m 50s | 42s | 47s | 4m 19s | ⚪ +0% (11/11 → 16/16) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 scenario steps PASSed on dev. spgi-100.csv loads correctly this time (previous run had to substitute SPGI.csv). C… | 3m 16s | 14s | 55s | 4m 25s | ⚪ +0% (15/15 → 15/15) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed end-to-end on dev: table linking (SELECTION_TO_FILTER and FILTER_TO_FILTER) propagated correctly betw… | 1m 46s | 17s | 35s | 2m 38s | ⚪ +0% (9/9 → 9/9) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran combined-boolean-filter end-to-end against dev. All 13 numbered scenario steps passed in the MCP browser phase: SEX_… | 2m 37s | 12s | 24s | 3m 13s | ⚪ +0% (13/13 → 13/13) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (14/14) | All 14 steps passed in both the MCP run and the Playwright replay. Expression filter works correctly: 5-rule AND yields … | 1m 14s | 8s | 23s | 1m 45s | ⚪ +0% (14/14 → 14/14) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 12 steps passed in the MCP run and in the Playwright replay (spec finished in 21.8s). The hierarchical filter correc… | 1m 15s | 21s | 23s | 1m 59s | ⚪ +0% (12/12 → 12/12) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed in the MCP run and in the Playwright replay (spec finished in 8.7s, total wall-clock 11.56s). The tex… | 1m 12s | 20s | 12s | 1m 44s | ⚪ +0% (9/9 → 9/9) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | 🟢 PASS → PASS | 🟢 100% (34/34) | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to reset… | 4m 24s | 40s | 1m 3s | 6m 7s | ⚪ +0% (34/34 → 34/34) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 69% (5/8) | Color consistency through layout round-trip works — the `.categorical-colors` tag survives save/reload and `R_ONE` stays… | 2m 30s | 35s | 26s | 3m 31s | ⚪ +0% (5/8 → 5/8) | 🟡 69% (baseline) | ⚪ 69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (9/12) | Filtering legend updates work end-to-end in the MCP run: numeric filter, categorical filter, layout round-trip, composed… | 3m 10s | 1m 10s | 44s | 5m 4s | ⚪ +0% (9/12 → 9/12) | 🟡 83% (baseline) | ⚪ 83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 77% (8/11) | Line chart legend and multi-axis behaviors are mostly correct: 7 legend items for 7 categories, layout round-trip preser… | 2m 10s | 40s | 32s | 3m 22s | ⚪ +0% (8/11 → 8/11) | 🟡 77% (baseline) | ⚪ 77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 58% (7/13) | Categorical legend on scatter plot updates correctly when X axis changes (sub 2) and when the Filter Panel narrows categ… | 4m 15s | 1m 20s | 54s | 6m 29s | ⚪ +0% (7/13 → 7/13) | 🟡 58% (baseline) | ⚪ 58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 79% (5/7) | Structure rendering in legends works for Scatter plot, Histogram, Line chart and Pie chart (canvas-based molecule thumbn… | 2m 35s | 40s | 29s | 3m 44s | ⚪ +0% (5/7 → 5/7) | 🟡 79% (baseline) | ⚪ 79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 70% (13/20) | Scenario executed end-to-end with a mix of PASS, AMBIGUOUS, and FAIL. Legend display, source-swap, corner positioning, a… | 5m 45s | 1m 30s | 41s | 7m 56s | ⚪ +0% (13/20 → 13/20) | 🟡 70% (baseline) | ⚪ 70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | 🟢 PASS → PASS | 🟡 86% (12/14) |  |  |  |  |  | 🟢 +9% (10/13 → 12/14) | 🟡 86% (baseline) | 🟢 77·86 | 🟡 7m (removed) | 🟡 1m (removed) | 🟡 17s (removed) | 🟡 8m 17s (removed) |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (82/82) | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color coding… | 3m 3s | 21s | 52s | 4m 16s | ⚪ +0% (82/82 → 82/82) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | bar-chart-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 95% (18/19) | 17 of 19 sections passed cleanly; section 8 combined into section 7 in the spec. Section 18 is AMBIGUOUS — `grok.dapi.pr… | 1m 5s | 8s | 32s | 1m 45s | ⚪ +0% (18/19 → 18/19) | 🟡 95% (baseline) | ⚪ 95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selectio… | 25s | 1m | 9.4s | 1m 34s | ⚪ +0% (11/11 → 11/11) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🔴 +0.4s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 11 steps passed. The entire test runs on the demog dataset (no SPGI_v2 needed). UI-only steps (Grid Color Coding All… |  |  | 29s | 29s | 🟢 +33% (8/12 → 12/12) | 🟢 100% (baseline) | 🟢 67·100 | 🟡 15s (removed) | 🟡 3s (removed) | 🟢 -55s | 🟢 -1m 13s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 87% (26/30) | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operation… | 5m | 3s | 22.5s | 5m 26s | ⚪ +0% (26/30 → 26/30) | 🟡 87% (baseline) | ⚪ 87·87 | ⚪ +0s | ⚪ +0s | 🔴 +0.5s | 🟢 -0.5s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (58/58) | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI interacti… |  | 2m |  | 2m | ⚪ +0% (58/58 → 58/58) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 18m (removed) | ⚪ +0s | ⚪ — | 🟢 -18m |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | 🟡 93% (28/30) | All 14 sections of form-tests-pw.md exercised across 30 steps. 28 PASS, 2 AMBIGUOUS, 0 FAIL in MCP run. Playwright spec … | 18m | 4m | 3m 12s | 25m 12s | ⚪ +0% (28/30 → 28/30) | 🟡 93% (baseline) | ⚪ 93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas elemen… | 18m | 3m | 51.5s | 21m 52s | ⚪ +0% (36/36 → 36/36) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.5s | 🟢 -0.5s |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 73% (16/22) | Grid tests ran 22 steps (spec softSteps); 17 passed outright and 5 were AMBIGUOUS (Copy/Paste, Column Header Context Men… | 11m | 3m | 1m 18s | 15m 18s | ⚪ +0% (16/22 → 16/22) | 🟡 73% (baseline) | ⚪ 73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | 🟡 94% (15/16) | All 14 heat-map sections exercised across 17 steps. 15 PASS, 1 AMBIGUOUS, 1 SKIP in MCP run. Playwright spec passed full… | 18m | 4m | 48.9s | 22m 49s | ⚪ +0% (15/16 → 15/16) | 🟡 94% (baseline) | ⚪ 94·94 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (87/94) | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance, l… | 50s | 7s | 46s | 1m 43s | ⚪ +0% (87/94 → 87/94) | 🟡 93% (baseline) | ⚪ 93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (26/26) | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout sa… | 57s | 8s | 1m 43s | 2m 48s | ⚪ +0% (26/26 → 26/26) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [map](Viewers/map-run.md) | 🟢 PASS → PASS | 🟡 80% (8/10) | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker size m… | 15s | 3s | 9s | 27s | ⚪ +0% (8/10 → 8/10) | 🟡 80% (baseline) | ⚪ 80·80 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (16/19) | Matrix Plot tests ran with 15 PASS, 3 AMBIGUOUS, 0 FAIL. The spec executed in 57.7s with all implemented steps passing. … | 20m | 3m | 55.9s | 23m 56s | ⚪ +0% (16/19 → 16/19) | 🟡 84% (baseline) | ⚪ 84·84 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (8/12) | 9 of 12 steps PASS; 3 SKIP (canvas-based node/edge interactions cannot be automated via DOM). The network diagram viewer… | 8m | 1m 30s | 22s | 9m 52s | ⚪ +0% (8/12 → 8/12) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | All 13 scenario sections (mapped to 12 Playwright softSteps — scale and normalization are combined in the spec) passed d… | 1m 8s | 8s | 47s | 2m 3s | ⚪ +0% (13/13 → 13/13) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (81/81) | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appearan… | 40s | 7s | 47s | 1m 34s | ⚪ +0% (81/81 → 81/81) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (17/20) | Pivot Table tests ran with 16 PASS, 2 AMBIGUOUS, 1 SKIP, 0 FAIL. The spec executed in 35.1s with all implemented steps p… | 16m | 3m | 1m 12s | 20m 12s | ⚪ +0% (17/20 → 17/20) | 🟡 85% (baseline) | ⚪ 85·85 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [row-source](Viewers/row-source-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with all 8… | 4m | 5s | 1m 24s | 5m 29s | ⚪ +0% (36/36 → 36/36) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (20/20) | All 20 sections passed during the MCP run on dev.datagrok.ai. The existing Playwright spec was re-run headed without mod… | 3m 13s | 29s | 52s | 4m 34s | ⚪ +0% (20/20 → 20/20) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | scatter-plot-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [statistics](Viewers/statistics-run.md) | 🟢 new: PASS | 🟢 100% (24/24) | All 23 MCP steps passed. The date columns section (STARTED row behavior) was moved to `statistics-tests-ui.md` as a manu… | 20m | 4m | 2m | 26m | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 20m (new) | 🟡 4m (new) | 🟡 2m (new) | 🟡 26m (new) |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | 🟢 PASS → PASS | 🟢 100% (24/24) | 24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to manual… |  | 3m | 58s | 3m 58s | ⚪ +0% (24/24 → 24/24) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 4m (removed) | ⚪ +0s | ⚪ +0s | 🟢 -4m |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (36/36) | All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via `v… | 28m | 4m | 46s | 32m 46s | ⚪ +0% (36/36 → 36/36) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (70/84) | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range slid… | 3m | 30s | 1m 48s | 5m 18s | ⚪ +0% (70/84 → 70/84) | 🟡 84% (baseline) | ⚪ 84·84 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All 7 MCP scenario steps PASS. The Word Cloud viewer adds via both entry points (Add-Viewer gallery and Toolbox icon), t… | 4m 15s | 1m | 2m 7s | 7m 22s | ⚪ +0% (7/7 → 7/7) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | word-cloud-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [working-with-nan-infinity](Viewers/working-with-nan-infinity-run.md) | 🟢 new: PASS | 🟢 100% (9/9) | All 9 spec steps PASSED in 1m 24s. NaN and Infinity values in numeric columns are handled gracefully across Scatter Plot… | 6m | 3m | 1m 24s | 10m 24s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 6m (new) | 🟡 3m (new) | 🟡 1m 24s (new) | 🟡 10m 24s (new) |
| Viewers | color-coding-(linked) | 🟢 removed: PASS |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | rendering-structures-on-the-axes | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | statistics-viewer | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | viewers-docking | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | working-with-nan-&-infinity | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |

## Comparison with Previous Reports

Deltas are computed against two baselines pulled from git history of `total-run.md`: `prev1d` (latest commit before today) and `prev7d` (commit closest to today − 7 days, ±3-day window). Signed with `+`/`-`; time deltas use the same format as the values. All status and delta cells carry the Legend icons.

### Totals

**Total (1d)**: Tests Δ **-10** · Run Δ **🔴 -3** · Status **🟡 PARTIAL → 🟡 PARTIAL** · Mean Pass Δ **🟢 +4%** · Browser Δ **🔴 +14.9s** · Spec Gen Δ **🔴 +6.3s** · Spec Run Δ **🟢 -0.6s** · Total Δ **🔴 +12s**

**Total (7d)**: Mean Pass Δ **—** · Browser Δ **🟡 3m 57s (new)** · Spec Gen Δ **🟡 57.6s (new)** · Spec Run Δ **🟡 41.3s (new)** · Total Δ **🟡 5m 22s (new)** _(7d-only deltas — count and status deltas live in the 1d row to avoid double-counting.)_

### By Folder

| Folder | Tests Δ | Run Δ | Status | Mean Pass Δ (1d) | Mean Pass Δ (7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|
| Apps | +0 | ⚪ +0 | ⚪ NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Browse | -4 | 🔴 -3 | 🟡 PARTIAL → PARTIAL | 🟢 +17% | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | 🔴 +8.2s | ⚪ +0s | 🔴 +1.8s | 🔴 +10s |
| Connections | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | 🔴 -1% | — | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🟢 -0.1s |
| General | -1 | ⚪ +0 | ⚪ NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| LocalCashing | -1 | ⚪ +0 | ⚪ NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | +0 | ⚪ +0 | ⚪ NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ +0s | ⚪ +0s | 🟢 -36.5s | 🟢 -42.1s |
| PowerPack | -1 | 🔴 -1 | 🟡 PARTIAL → PARTIAL | 🟢 +41% | — | 🔴 +1m 30s | 🔴 +30.3s | 🔴 +5.6s | 🔴 +2m 6s |
| Projects | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | +0 | ⚪ +0 | ⚪ NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | +0 | ⚪ +0 | 🟡 PARTIAL → PARTIAL | ⚪ +0% | — | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Tooltips | +0 | ⚪ +0 | ⚪ NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | -3 | 🟢 +1 | 🟡 PARTIAL → PARTIAL | 🟢 +2% | — | 🔴 +21.6s | 🔴 +10.3s | 🔴 +1.8s | 🔴 +5.3s |

### Per-Test Changes

Lists tests where Pass % (1d or 7d), status, or any timing component changed vs. either baseline — plus persistent-failure rows where Pass % is still < 100% even with a flat Δ. Tests at a clean 100% with no timing change are omitted.

| Folder | Test | Status | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|
| Apps | apps | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Apps | tutorials | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/5 → 4/5) | 🟡 80% (baseline) | ⚪ 80·80 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 AMBIGUOUS → AMBIGUOUS | ⚪ +0% (4/6 → 4/6) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/5 → 4/5) | 🟡 90% (baseline) | ⚪ 90·90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (0/1 → 0/1) | 🟡 50% (baseline) | ⚪ 50·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | package-manager | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | japanese-in-myfiles | 🟢 removed: PASS | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | local-deploy | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces | 🟡 removed: PARTIAL | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | spaces-(ui-only) | 🔴 removed: FAIL | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/12 → 5/12) | 🟡 42% (baseline) | ⚪ 42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (1/5 → 1/5) | 🟡 20% (baseline) | ⚪ 20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/6 → 3/6) | 🟡 50% (baseline) | ⚪ 50·50 | 🔴 +17s | ⚪ +0s | 🟢 -0.4s | 🔴 +16.6s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🔴 +0.2s | 🔴 +0.2s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +12s | ⚪ +0s | 🔴 +1.8s | 🔴 +13.8s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +8s | ⚪ +0s | 🔴 +5s | 🔴 +13s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +44s | ⚪ +0s | 🔴 +4s | 🔴 +48s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/3 → 1/3) | 🟡 33% (baseline) | ⚪ 33·33 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +26s | ⚪ +0s | 🔴 +3.8s | 🔴 +29.8s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/5 → 2/5) | 🟡 40% (baseline) | ⚪ 40·40 | 🟢 -34s | ⚪ +0s | 🔴 +1.1s | 🟢 -32.9s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🟢 -22s | ⚪ +0s | 🔴 +1.9s | 🟢 -20.1s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | 🟢 PASS → PASS | ⚪ +0% (2/2 → 2/2) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +4s | ⚪ +0s | 🔴 +3.2s | 🔴 +7.2s |
| Chem | [mmp](Chem/mmp-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +57s | ⚪ +0s | 🔴 +2s | 🔴 +59s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | 🟢 PASS → PASS | ⚪ +0% (5/5 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.3s | 🟢 -0.3s |
| Chem | [sketcher](Chem/sketcher-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +3s | ⚪ +0s | 🔴 +3.3s | 🔴 +6.3s |
| Connections | [adding](Connections/adding-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/7 → 6/7) | 🟡 93% (baseline) | ⚪ 93·93 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [browser](Connections/browser-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (7/9 → 7/9) | 🟡 78% (baseline) | ⚪ 78·78 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/11 → 1/11) | 🟡 9% (baseline) | ⚪ 9·9 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (0/7 → 0/7, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/9 → 1/9) | 🟡 11% (baseline) | ⚪ 11·11 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [import-swagger](Connections/import-swagger-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (0/7 → 0/7, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [schema](Connections/schema-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/4 → 3/4) | 🟡 75% (baseline) | ⚪ 75·75 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/7 → 6/7) | 🟡 86% (baseline) | ⚪ 86·86 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/6 → 5/6) | 🟡 83% (baseline) | ⚪ 83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [ML methods/linear-regression](EDA/ML methods/linear-regression-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🟢 -0.2s |
| EDA | [ML methods/pls-regression](EDA/ML methods/pls-regression-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/7 → 5/7) | 🟡 71% (baseline) | ⚪ 71·71 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🟢 -0.1s |
| EDA | [ML methods/softmax](EDA/ML methods/softmax-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/3 → 1/3) | 🟡 33% (baseline) | ⚪ 33·33 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🟢 -0.4s |
| EDA | [ML methods/xgboost1](EDA/ML methods/xgboost1-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/3 → 2/3) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [ML methods/xgboost2](EDA/ML methods/xgboost2-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/3 → 2/3) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/3 → 2/3) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/7 → 3/7) | 🟡 43% (baseline) | ⚪ 43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/5 → 3/5) | 🟡 60% (baseline) | ⚪ 60·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🔴 -12% (2/4 → 2/5) | 🟡 50% (baseline) | 🔴 62·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| General | files-cache | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | first-login | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | inactivity-response | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | login | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | molecule-in-exported-csv | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | network | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | profile-settings | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | startup-time | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | table-manager | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | tabs-reordering | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| General | api-samples | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| LocalCashing | local-cashing | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/6 → 4/6) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (5/18 → 5/18) | 🟡 28% (baseline) | ⚪ 28·28 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | browser | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | create | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | delete | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | edit | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [info-panels](Peptides/info-panels-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/5 → 2/5) | 🟡 40% (baseline) | ⚪ 40·40 | ⚪ +0s | ⚪ +0s | 🟡 1m 18s (removed) | 🟢 -1m 18s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/6 → 4/6) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | 🟢 -0.4s | 🟢 -0.4s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/9 → 6/9) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | 🟡 1m 30s (removed) | 🟢 -1m 30s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 +40% (6/10 → 10/10) | 🟢 100% (baseline) | 🟢 60·100 | 🟡 7m 45s (new) | 🟡 2m (new) | 🟡 32s (new) | 🟡 10m 17s (new) |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | 🟢 PASS → PASS | ⚪ +0% (8/8 → 7/7) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 1m 45s (new) | 🟡 20s (new) | 🟡 8s (new) | 🟡 2m 13s (new) |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | 🟢 PASS → PASS | ⚪ +0% (5/5 → 7/7) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 2m 15s (new) | 🟡 1m (new) | 🟡 28s (new) | 🟡 3m 43s (new) |
| PowerPack | [AddNewColumn/functions-sorting](PowerPack/AddNewColumn/functions-sorting-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 8m 20s (new) | 🟡 2m (new) | 🟡 17s (new) | 🟡 10m 37s (new) |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 5/5) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 4m 43s (new) | 🟡 1m 13s (new) | 🟡 25s (new) | 🟡 6m 21s (new) |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 4/4) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 1m 10s (new) | 🟡 15s (new) | 🟡 8s (new) | 🟡 1m 33s (new) |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | 🟡 → 🟢 SKIP → PASS | 🟢 +100% (0/6 → 10/10) | 🟢 100% (baseline) | 🟢 0·100 | 🟡 3m 30s (new) | 🟡 1m (new) | 🟡 18.8s (new) | 🟡 4m 49s (new) |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | 🟢 PASS → PASS | ⚪ +0% (5/5 → 6/6) | 🟢 100% (baseline) | ⚪ 100·100 | 🔴 +1m 32s | 🔴 +29s | 🔴 +5.7s | 🔴 +2m 7s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🔴 → 🟡 FAIL → PARTIAL | 🟢 +76% (0/10 → 16/21) | 🟡 76% (baseline) | 🟢 0·76 | 🔴 +8m 50s | 🔴 +45s | 🔴 +1m 17s | 🔴 +10m 52s |
| PowerPack | AddNewColumn/functions_sorting | 🟡 removed: SKIP | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | formula-lines | 🟡 removed: SKIP | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [browser](Projects/browser-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/9 → 5/9) | 🟡 56% (baseline) | ⚪ 56·56 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [complex](Projects/complex-run.md) | 🟡 SKIP → SKIP | ⚪ +0% (0/13 → 0/13, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 SKIP → SKIP | ⚪ +0% (0/5 → 0/5, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/4 → 2/4) | 🟡 50% (baseline) | ⚪ 50·50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 SKIP → SKIP | ⚪ +0% (0/4 → 0/4, still broken) | 🔴 0% (baseline) | ⚪ 0·0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | 🟡 SKIP → SKIP | ⚪ +0% (2/5 → 2/5) | 🟡 40% (baseline) | ⚪ 40·40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/14 → 8/14) | 🟡 57% (baseline) | ⚪ 57·57 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | adding | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | browse-&-save-project | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | browser | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | columns-inspect | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | deleting | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | edit | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | get-all-get-top-100 | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | ms-sql | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | new-sql-query | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | new-visual-query | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | query-layout | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | query-postprocessing | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | transformations | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | visual-query-advanced | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (10/12 → 10/12) | 🟡 92% (baseline) | ⚪ 92·92 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | layout | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/9 → 6/9) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/5 → 1/5) | 🟡 20% (baseline) | ⚪ 20·20 | ⚪ +0s | ⚪ +0s | ⚪ — | ⚪ +0s |
| Tooltips | actions-in-the-context-menu | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip-visibility | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | edit-tooltip | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | line-chart---aggregated-tooltip | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | tooltip-properties | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | uniform-default-tooltip | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | ⚪ +0% (15/15 → 15/15) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -50s | 🟢 -50s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/8 → 5/8) | 🟡 69% (baseline) | ⚪ 69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (9/12 → 9/12) | 🟡 83% (baseline) | ⚪ 83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/11 → 8/11) | 🟡 77% (baseline) | ⚪ 77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (7/13 → 7/13) | 🟡 58% (baseline) | ⚪ 58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/7 → 5/7) | 🟡 79% (baseline) | ⚪ 79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (13/20 → 13/20) | 🟡 70% (baseline) | ⚪ 70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | 🟢 PASS → PASS | 🟢 +9% (10/13 → 12/14) | 🟡 86% (baseline) | 🟢 77·86 | 🟡 7m (removed) | 🟡 1m (removed) | 🟡 17s (removed) | 🟡 8m 17s (removed) |
| Viewers | bar-chart-tests | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (18/19 → 18/19) | 🟡 95% (baseline) | ⚪ 95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS | ⚪ +0% (11/11 → 11/11) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | 🔴 +0.4s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | 🟢 PASS → PASS | 🟢 +33% (8/12 → 12/12) | 🟢 100% (baseline) | 🟢 67·100 | 🟡 15s (removed) | 🟡 3s (removed) | 🟢 -55s | 🟢 -1m 13s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (26/30 → 26/30) | 🟡 87% (baseline) | ⚪ 87·87 | ⚪ +0s | ⚪ +0s | 🔴 +0.5s | 🟢 -0.5s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | ⚪ +0% (58/58 → 58/58) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 18m (removed) | ⚪ +0s | ⚪ — | 🟢 -18m |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | ⚪ +0% (36/36 → 36/36) | 🟢 100% (baseline) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟢 -0.5s | 🟢 -0.5s |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (16/22 → 16/22) | 🟡 73% (baseline) | ⚪ 73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | ⚪ +0% (15/16 → 15/16) | 🟡 94% (baseline) | ⚪ 94·94 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (87/94 → 87/94) | 🟡 93% (baseline) | ⚪ 93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (16/19 → 16/19) | 🟡 84% (baseline) | ⚪ 84·84 | ⚪ +0s | ⚪ +0s | 🟢 -0.1s | 🟢 -0.1s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/12 → 8/12) | 🟡 67% (baseline) | ⚪ 67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (17/20 → 17/20) | 🟡 85% (baseline) | ⚪ 85·85 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | scatter-plot-tests | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [statistics](Viewers/statistics-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 20m (new) | 🟡 4m (new) | 🟡 2m (new) | 🟡 26m (new) |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | 🟢 PASS → PASS | ⚪ +0% (24/24 → 24/24) | 🟢 100% (baseline) | ⚪ 100·100 | 🟡 4m (removed) | ⚪ +0s | ⚪ +0s | 🟢 -4m |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (70/84 → 70/84) | 🟡 84% (baseline) | ⚪ 84·84 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | word-cloud-tests | ⚪ NO RUN → NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [working-with-nan-infinity](Viewers/working-with-nan-infinity-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 6m (new) | 🟡 3m (new) | 🟡 1m 24s (new) | 🟡 10m 24s (new) |
| Viewers | color-coding-(linked) | 🟢 removed: PASS | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | rendering-structures-on-the-axes | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | statistics-viewer | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | viewers-docking | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | working-with-nan-&-infinity | ⚪ removed: NO RUN | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage is 76% (≥70%) and no folders are fully failing, but 14 folder(s) are PARTIAL:
- Bio, Browse, Charts, Chem, Connections, DiffStudio, EDA, Models, Peptides, PowerPack, Projects, Scripts, StickyMeta, Viewers

### Blocking Issues
- Bio: PARTIAL (Bio/composition-analysis, Bio/pepsea)
- Browse: PARTIAL (Browse/browse, Browse/browse-tree-states)
- Charts: PARTIAL (Charts/sunburst, Charts/tree)
- Chem: PARTIAL (Chem/Advanced/scaffold-tree, Chem/calculate, Chem/chemprop)
- Connections: PARTIAL (Connections/adding, Connections/browser, Connections/catalogs, Connections/external-provider, Connections/identifiers, Connections/import-swagger, …)
- DiffStudio: PARTIAL (DiffStudio/fitting)
- EDA: PARTIAL (EDA/ML methods/pls-regression, EDA/ML methods/softmax, EDA/ML methods/xgboost1, EDA/ML methods/xgboost2, EDA/multivariate-analysis, EDA/pareto-front-viewer, …)
- Models: PARTIAL (Models/browser, Models/chemprop)
- Peptides: PARTIAL (Peptides/peptide-space, Peptides/peptides, Peptides/sar)
- PowerPack: PARTIAL (PowerPack/AddNewColumn/add-new-column, PowerPack/data-enrichment)
- Projects: PARTIAL (Projects/browser, Projects/complex, Projects/custom-creation-scripts, Projects/deleting, Projects/opening, Projects/project-url, …)
- Scripts: PARTIAL (Scripts/create, Scripts/run)
- StickyMeta: PARTIAL (StickyMeta/database-meta)
- Viewers: PARTIAL (Viewers/Legend/color-consistency, Viewers/Legend/filtering, Viewers/Legend/line-chart, Viewers/Legend/scatterplot, Viewers/Legend/structure-rendering, Viewers/Legend/visibility-and-positioning, …)
