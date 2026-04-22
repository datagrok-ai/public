# Test Track — Global Report

**Date**: 2026-04-22
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

**Total**: 187 tests · Run: 138/187 (74%) · Playwright: 108/187 (58%) · Mean Pass: 🟡 77% · Mean Browser: 3m 42s · Mean Spec Gen: 51s · Mean Spec Run: 42s · Mean Total (sum per test): 5m 10s

| Folder | Tests | Run | Playwright | Status | Mean Pass % | Mean Browser | Mean Spec Gen | Mean Spec Run | Mean Total |
|--------|-------|-----|------------|--------|-------------|--------------|---------------|---------------|------------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Bio | 9 | 9/9 (100%) | 1/9 (11%) | 🟡 PARTIAL | 🟡 94% | 2m 45s | 5.0s |  | 2m 46s |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | 🟡 PARTIAL | 🟡 53% |  |  |  |  |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | 🟡 PARTIAL | 🟡 54% | 2m 10s | 45s | 34s | 3m 29s |
| Chem | 14 | 14/14 (100%) | 14/14 (100%) | 🟡 PARTIAL | 🟡 87% | 1m 29s | 32s | 38s | 2m 39s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | 🟡 PARTIAL | 🟡 54% |  |  |  |  |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | 🟡 PARTIAL | 🟢 96% | 3m 10s | 56s | 45s | 4m 52s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | 🟡 PARTIAL | 🟡 67% | 1m 42s | 58s | 17s | 2m 58s |
| General | 11 | 0/11 (0%) | 0/11 (0%) | ⚪ NO DATA |  |  |  |  |  |
| LocalCashing | 1 | 0/1 (0%) | 0/1 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | 🟡 PARTIAL | 🟡 81% |  |  |  |  |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 68% | 28s | 3.0s | 48s | 1m 18s |
| PowerPack | 10 | 10/10 (100%) | 7/10 (70%) | 🟡 PARTIAL | 🟡 56% | 3m 46s | 38s | 24s | 4m 49s |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | 🟡 PARTIAL | 🟡 38% |  |  |  |  |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | 🟡 PARTIAL | 🟡 90% |  |  |  |  |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 68% | 24s | 2.2s | 15s | 37s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Viewers | 49 | 42/49 (86%) | 42/49 (86%) | 🟡 PARTIAL | 🟡 90% | 5m 50s | 1m 7s | 51s | 7m 46s |

## All Tests

**Total**: 187 tests · 🟢 74 PASS / 🟡 45 PARTIAL / 🔴 11 FAIL / 🟡 1 AMBIGUOUS / 🟡 7 SKIP / ⚪ 49 NO RUN · Mean Pass: 🟡 77% · Mean Browser: 3m 42s · Mean Spec Gen: 51s · Mean Spec Run: 42s · Mean Total (sum per test): 5m 10s

| Folder | Test | Status | Pass % | Description | Browser (model+MCP) | Spec Gen (model) | Spec Run (Playwright) | Total (sum) | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Apps | apps | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Apps | tutorials | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [analyze](Bio/analyze-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset typ… | 8m | 5.0s |  | 8m 5s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ — | ⚪ 0.0s |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 80% (4/5) | 4 of 5 steps passed. The Composition/WebLogo viewer opens correctly and properties are accessible. The letter-click sele… |  |  |  |  | 🟡 80% (baseline) | 🟡 80% (baseline) | ⚪ 80 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [convert](Bio/convert-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 4 Bio convert/transform functions work correctly on FASTA data. Get Region extracts a subsequence region. PolyTool >… | 2m |  |  | 2m | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ — | ⚪ — | ⚪ 0.0s |
| Bio | [manage](Bio/manage-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All 3 steps passed. The Manage Monomer Libraries view opens as a full view showing 5 monomer library JSON files (increas… | 30s |  |  | 30s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ — | ⚪ — | ⚪ 0.0s |
| Bio | [msa](Bio/msa-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. MSA dialog opens with correct fields, Alignment Parameters button adds gap penalty inputs as expecte… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 AMBIGUOUS → AMBIGUOUS | 🟡 67% (4/6) | The MSA dialog opens correctly for HELM data and shows MAFFT-based method options (mafft --auto, linsi, ginsi, etc.) ins… |  |  |  |  | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [search](Bio/search-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 4 steps passed. Bio > Search > Subsequence Search opens a filter panel with a Sequence bio substructure filter. Typi… | 30s |  |  | 30s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ — | ⚪ — | ⚪ 0.0s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. Activity Cliffs works correctly with both default parameters (UMAP/Hamming) and custom parameters (t… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | [sequence-space](Bio/sequence-space-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. Sequence Space works correctly with both default (UMAP/Hamming) and custom (t-SNE/Needlemann-Wunsch)… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 90% (4/5) | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. Item… |  |  |  |  | 🟡 90% (baseline) | 🟡 90% (baseline) | ⚪ 90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (0/1) | 1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a single session… |  |  |  |  | 🟡 50% (baseline) | 🟡 50% (baseline) | ⚪ 50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [japanese-in-myfiles](Browse/japanese-in-myfiles-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All 3 steps passed. The file 芹沢 貴之 こんにちは.csv displays correctly with no garbled or incorrect characters. Japanese Kanji … |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | local-deploy | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | package-manager | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [spaces](Browse/spaces-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 24% (8/33) | 8 steps passed, 0 failed, 24 skipped out of 32 total steps. Core CRUD operations (create root/child/nested, rename, add … |  |  |  |  | 🟡 24% (baseline) | 🟡 24% (baseline) | ⚪ 24 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | 🔴 FAIL → FAIL | 🟡 3% (1/30) | 1 step passed, 2 ambiguous, 27 skipped out of 30 total steps. This scenario focuses on UI-only interactions not covered … |  |  |  |  | 🟡 3% (baseline) | 🟡 3% (baseline) | ⚪ 3 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [radar](Charts/radar-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | The Radar viewer reproduced cleanly on dev for both earthquakes.csv (2426 rows) and demog.csv (5850 rows). All 21 Radar … | 1m 30s | 35s | 33s | 2m 38s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 42% (5/12) | The sunburst viewer reproduces structurally on dev — Sunburst can be added to both SPGI and demog, and `hierarchyColumnN… | 3m | 1m | 34s | 4m 34s | 🟡 42% (baseline) | 🟡 42% (baseline) | ⚪ 42 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 20% (1/5) | Setup (open demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduced cleanly on dev. All four test steps are mar… | 2m | 40s | 36s | 3m 16s | 🟡 20% (baseline) | 🟡 20% (baseline) | ⚪ 20 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (3/6) | Smoke coverage only: Scaffold Tree viewer launches from the Chem menu and the magic wand generates a scaffold tree on SP… | 40s | 25s | 50s | 1m 55s | 🟡 50% (baseline) | 🟡 50% (baseline) | ⚪ 50 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Scaffold Tree viewer launches from the Chem → Scaffold Tree menu, and the magic-wand generator produces scaffold nodes f… | 1m 15s | 30s | 40s | 2m 25s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Similarity Search launches from the Chem menu and exposes a viewer that accepts option changes (fingerprint Morgan ↔ Pat… | 25s | 25s | 23s | 1m 13s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Substructure filtering via `grok.chem.searchSubstructure` works on SPGI.csv (3624 rows): benzene substructure yields a b… | 30s | 25s | 24s | 1m 19s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Activity Cliffs computation on SPGI.csv (3624 rows) finishes within 45s and produces a UMAP scatter plot with molecule t… | 30s | 20s | 1m | 1m 50s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) | Calculate Descriptors cannot be exercised on `dev` right now. The Chem top menu fails to open its popup — both through D… | 8m | 1m | 38s | 9m 38s | 🟡 33% (baseline) | 🟡 33% (baseline) | ⚪ 33 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Chemical Space dimensional reduction runs end-to-end on smiles.csv: the dialog opens, OK with defaults produces a Scatte… | 20s | 20s | 56s | 1m 36s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | ChemProp scenario is partially automated: the spec confirms mol1K.sdf opens and the Train Model view is reachable from t… | 1m | 30s | 18s | 1m 48s | 🟡 40% (baseline) | 🟡 40% (baseline) | ⚪ 40 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Elemental Analysis works on dev. The menu path `[name="div-Chem"]` → `Elemental Analysis...` resolves and the dialog ope… | 1m | 30s | 29s | 1m 59s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | 🟢 PASS → PASS | 🟢 100% (2/2) | The filter panel correctly shows a Structure filter for SPGI.csv's Molecule column; clicking the embedded sketch-link op… | 30s | 25s | 21s | 1m 16s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [info-panels](Chem/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | Info panels work correctly on smiles.csv: column-level (Details, Filter, Colors, Style, Chemistry with Rendering/Highlig… | 3m 10s | 1m | 31s | 4m 41s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [mmp](Chem/mmp-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset at the bottom of the view… | 30s | 20s | 1m 12s | 2m 2s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | R-Groups Analysis works on sar_small.csv: MCS auto-populates the sketcher, OK produces a Trellis plot and appends R1–R4 … | 2m 20s | 45s | 59s | 4m 4s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [sketcher](Chem/sketcher-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Sketcher opens via `grok.chem.sketcher(molCol, initialSmiles)` wrapped in `ui.dialog(...).show()`, accepts a typed SMILE… | 30s | 30s | 16s | 1m 16s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Connections | [adding](Connections/adding-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (6/7) | 6 of 7 steps fully passed, Step 5 was partial (TEST button works but actual connection test fails without real credentia… |  |  |  |  | 🟡 93% (baseline) | 🟡 93% (baseline) | ⚪ 93 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [browser](Connections/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 78% (7/9) | 7 of 9 steps passed, 2 ambiguous. Search filtering works correctly and the Context Pane shows all expected tabs (Details… |  |  |  |  | 🟡 78% (baseline) | 🟡 78% (baseline) | ⚪ 78 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | 🟡 9% (1/11) | 1 step passed, 1 failed, 15 skipped. The required `NorthwindTest` MS SQL connection is not present on public.datagrok.ai… |  |  |  |  | 🟡 9% (baseline) | 🟡 9% (baseline) | ⚪ 9 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [delete](Connections/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 steps passed. Both connections were deleted successfully. The confirmation dialog uses a red "DELETE" button (not … |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [edit](Connections/edit-run.md) | 🟢 PASS → PASS | 🟡 86% (6/7) | 6 of 7 steps passed (1 skipped due to missing real credentials). The connection rename, credential modification, and err… |  |  |  |  | 🟡 86% (baseline) | 🟡 86% (baseline) | ⚪ 86 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 FAIL → FAIL | 🔴 0% (0/7) | All 7 steps skipped. This scenario requires a specific Postgres connection at db.datagrok.ai:54327 with superuser creden… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 FAIL → FAIL | 🟡 11% (1/9) | 1 step passed, 1 failed, 7 skipped. This scenario depends on a working Postgres connection to the Northwind database. Th… |  |  |  |  | 🟡 11% (baseline) | 🟡 11% (baseline) | ⚪ 11 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [import-swagger](Connections/import-swagger-run.md) | 🔴 FAIL → FAIL | 🔴 0% (0/7) | All 7 steps skipped. This scenario requires manual interaction: downloading a YAML file to the local machine and drag-dr… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [schema](Connections/schema-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 75% (3/4) | 3 of 4 steps passed, 1 ambiguous. The "Browse schema" context menu option was not found in the current UI, but the schem… |  |  |  |  | 🟡 75% (baseline) | 🟡 75% (baseline) | ⚪ 75 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 86% (6/7) | 6 of 7 steps passed (1 failed). All UI steps worked correctly. The SPARQL connection was created and deleted successfull… |  |  |  |  | 🟡 86% (baseline) | 🟡 86% (baseline) | ⚪ 86 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | Catalog scenario reproduces fully on dev.datagrok.ai. All 6 steps PASS both in MCP and in the Playwright spec (35.8s wal… | 1m 53s | 39s | 38s | 3m 10s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Cyclic Models (PK-PD) scenario reproduces fully on dev.datagrok.ai. The PK-PD library model loads via double-click, Mult… | 1m 2s | 38s | 36s | 2m 16s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Files & Sharing scenario reproduces fully on dev.datagrok.ai. pk.ivp loads via the `DiffStudio:previewIvp` function with… | 2m 31s | 1m 32s | 1m 17s | 5m 20s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (5/6) | Scenario is PARTIAL on dev.datagrok.ai — steps 1–5 pass, step 6 (actually running the fit) does not produce result rows … | 10m 2s | 55s | 1m | 11m 57s | 🟡 83% (baseline) | 🟡 83% (baseline) | ⚪ 83 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | 🟢 PASS → PASS | 🟡 83% (5/6) | Scenario fully reproduces on dev.datagrok.ai. All 6 steps pass in the interactive MCP session and in the Playwright spec… | 1m 31s | 36s | 26s | 2m 33s | 🟡 83% (baseline) | 🟡 83% (baseline) | ⚪ 83 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 scenario steps PASS against dev.datagrok.ai. Edit toggle is reachable via `.d4-ribbon-item .ui-input-bool-switch .… | 4m 30s | 1m 40s | 1m 3s | 7m 13s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Sensitivity Analysis scenario fully reproduces against dev.datagrok.ai. Bioreactor loads from the DiffStudio hub (librar… | 2m 3s | 44s | 39s | 3m 26s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Stages (Acid Production) scenario reproduces fully on dev.datagrok.ai. The library card opens a view named "GA-productio… | 1m 49s | 47s | 24s | 3m | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/linear-regression](EDA/ML methods/linear-regression-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation, c… | 45s | 2.0s | 6.8s | 54s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/pls-regression](EDA/ML methods/pls-regression-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 71% (5/7) | PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-4 c… | 1m | 2.0s | 6.9s | 1m 9s | 🟡 71% (baseline) | 🟡 71% (baseline) | ⚪ 71 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/softmax](EDA/ML methods/softmax-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) | Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns and w… | 10s | 2.0s | 2.6s | 15s | 🟡 33% (baseline) | 🟡 33% (baseline) | ⚪ 33 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/xgboost1](EDA/ML methods/xgboost1-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as Ui… | 5.0s | 2.0s | 2.7s | 9.7s | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/xgboost2](EDA/ML methods/xgboost2-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interactio… | 5.0s | 2.0s | 2.7s | 9.7s | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [anova](EDA/anova-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All 3 scenario steps passed against dev. Dataset opens via JS API in ~1s; ANOVA dialog mounts with sensible defaults (RA… | 1m 30s | 30s | 26s | 2m 26s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | 2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check, where the wording does not spec… | 2m 30s | 2m | 13s | 4m 43s | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 43% (3/7) | 3 of 7 steps passed, 1 failed, 3 were skipped due to the missing prerequisite dataset. The Pareto Front viewer itself is… | 4m | 2m | 32s | 6m 32s | 🟡 43% (baseline) | 🟡 43% (baseline) | ⚪ 43 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 60% (3/5) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 3 PASS / 1 FAIL / 1 SKIP. The dialog path works (menu, F… | 5m | 3m | 1m 7s | 9m 7s | 🟡 60% (baseline) | 🟡 60% (baseline) | ⚪ 60 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🟡 62% (2/4) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 2 PASS / 1 PARTIAL / 1 FAIL. The dialog path (menu, Usin… | 2m | 2m | 15s | 4m 15s | 🟡 62% (baseline) | 🟡 62% (baseline) | ⚪ 62 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| General | api-samples | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
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
| LocalCashing | local-cashing | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [apply](Models/apply-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 4 steps passed. The Apply Model workflow functions correctly end-to-end. The TestDemog model (trained in Train.md) w… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (4/6) | 4 of 6 steps passed; 2 skipped because only 1 model was available (TestDemog was the only model — the second numeric mod… |  |  |  |  | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🟡 28% (5/18) | 2 of 17 sub-steps passed, 13 skipped, 1 failed, 1 ambiguous. The scenario fails entirely due to the Chemprop Docker cont… |  |  |  |  | 🟡 28% (baseline) | 🟡 28% (baseline) | ⚪ 28 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [delete](Models/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps passed. The model deletion workflow works correctly end-to-end. The right-click context menu, confirmation d… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [predictive-models](Models/predictive-models-run.md) | 🟢 PASS → PASS | 🟢 100% (21/21) | All 20 sub-steps passed. The full lifecycle (Train → Apply → Apply on new dataset → Delete) for EDA-based predictive mod… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [train](Models/train-run.md) | 🟢 PASS → PASS | 🟡 90% (9/10) | All 10 steps passed. The Train Model workflow functions correctly on public.datagrok.ai using the built-in EDA engines. … |  |  |  |  | 🟡 90% (baseline) | 🟡 90% (baseline) | ⚪ 90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | browser | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | create | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | delete | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | edit | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [info-panels](Peptides/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are rende… | 17s | 3.0s | 11s | 31s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ -0.1s | ⚪ -0.1s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variability… | 25s | 3.0s | 1m 18s | 1m 46s | 🟡 40% (baseline) | 🟡 40% (baseline) | ⚪ 40 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (4/6) | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clusters… | 18s | 3.0s | 12s | 33s | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ -0.4s | ⚪ -0.4s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Residu… | 50s | 3.0s | 1m 30s | 2m 23s | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 60% (6/10) | Core Add New Column functionality works correctly: creating calculated columns, chaining formulas, and propagating colum… |  |  |  |  | 🟡 60% (baseline) | 🟡 60% (baseline) | ⚪ 60 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 steps passed. Autocomplete works for both functions (typing or Ctrl+Space) and columns ($ symbol). Function select… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps passed. The formula dependency chain (Weight2 -> Weight3 -> Weight4) works correctly. Changing a value in We… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/5) | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on releas… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 4 steps passed. Column names in both ${} (scalar) and $[] (aggregate) syntax are highlighted with the cm-column-name… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All 4 steps passed. Hovering over a function name in the formula editor displays a tooltip with the function signature i… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/6) | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on releas… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps passed on dev.datagrok.ai. The Add New Column dialog opens correctly, displays cleanly with no visual issues… | 22s | 2.0s | 5.3s | 29s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🟡 → 🔴 PARTIAL → FAIL | 🔴 0% (0/10) | The scenario is blocked on infrastructure. Dev has a `NorthwindTest` Postgres connection at `db.datagrok.ai:54322 / nort… | 7m 10s | 1m 15s | 43s | 9m 8s | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | 🔴 +6m 28s | 🔴 +1m 12s | 🟡 43s (new) | 🔴 +8m 23s |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/1) | This scenario file contains only references to GitHub tickets related to formula lines. No executable test steps were de… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [browser](Projects/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 56% (5/9) | 5 of 9 steps passed, 3 skipped, 1 ambiguous. Browse > Dashboards view works correctly: projects are listed, searchable, … |  |  |  |  | 🟡 56% (baseline) | 🟡 56% (baseline) | ⚪ 56 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [complex](Projects/complex-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/13) | All 13 steps skipped. This is the most complex scenario requiring tables from 7+ different sources, drag-and-drop, entit… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/5) | All 5 steps skipped. This scenario requires running a custom JavaScript script with Data Sync enabled, then modifying fi… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (2/4) | 2 of 4 steps passed, 1 skipped, 1 ambiguous. Project deletion works via the API (`grok.dapi.projects.delete()`). The rig… |  |  |  |  | 🟡 50% (baseline) | 🟡 50% (baseline) | ⚪ 50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [opening](Projects/opening-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (5/5) | All 5 steps passed. Projects from the Uploading step are accessible in Browse > Dashboards. Context Panel correctly show… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/4) | All steps skipped. This scenario depends on Projects copy_clone.md (order 5) which was not fully executed. The Link/Clon… |  |  |  |  | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | 🟡 SKIP → SKIP | 🟡 40% (2/5) | 2 of 5 steps passed, 3 skipped. Project preview and opening work. Copy/clone/link operations were not tested because the… |  |  |  |  | 🟡 40% (baseline) | 🟡 40% (baseline) | ⚪ 40 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 57% (8/14) | 8 of 14 steps passed, 6 skipped. Core project creation from local tables, file shares, query results, and join results a… |  |  |  |  | 🟡 57% (baseline) | 🟡 57% (baseline) | ⚪ 57 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
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
| Scripts | [browser](Scripts/browser-run.md) | 🟢 PASS → PASS | 🟡 89% (8/9) | The Scripts Browser scenario passed well. The context pane shows all expected accordions (Details, Script, Run, Activity… |  |  |  |  | 🟡 89% (baseline) | 🟡 89% (baseline) | ⚪ 89 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 92% (10/12) | The Create scenario completed successfully overall. The script `testRscript` was created, parameters configured, saved, … |  |  |  |  | 🟡 92% (baseline) | 🟡 92% (baseline) | ⚪ 92 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [delete](Scripts/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps passed. The delete flow works correctly with a confirmation dialog and immediate removal from the scripts li… |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [edit](Scripts/edit-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. The Edit scenario works correctly — edits are saved persistently and visible on re-open. |  |  |  |  | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | layout | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Core run functionality works: the script can be triggered from context menu with a table selection and from the console.… |  |  |  |  | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All tested steps passed. SPGI.csv opened with TestSchema1 sticky metadata schema pre-configured. The Sticky meta panel i… | 35s | 3.0s | 17s | 55s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | 🟢 PASS → PASS | 🟡 50% (2/4) | Steps 1-2 passed: SPGI.csv opened with TestSchema1 sticky metadata schema, and cloning the table preserves the schema an… | 25s | 3.0s | 21s | 49s | 🟡 50% (baseline) | 🟡 50% (baseline) | ⚪ 50 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All steps passed. The Sticky Meta Schemas browser at `/meta/schemas` shows 20 schemas including TestSchema1. The "NEW SC… | 10s | 3.0s | 7.0s | 20s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | 🟡 20% (1/5) | Step 1 passed: navigated to Databases > Postgres and found CHEMBL connection. Step 2 failed: the "Database meta" section… | 25s | 0.0s |  | 25s | 🟡 20% (baseline) | 🟡 20% (baseline) | ⚪ 20 | ⚪ 0.0s | ⚪ 0.0s | ⚪ — | ⚪ 0.0s |
| Tooltips | actions-in-the-context-menu | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | default-tooltip-visibility | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | edit-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | line-chart---aggregated-tooltip | ⚪ NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | tooltip-properties | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Tooltips | uniform-default-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 14 steps (setup + 13 scenario sections) passed in both the browser-driven MCP run against https://dev.datagrok.ai an… | 7m | 2m | 1m 10s | 10m 10s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -29m | 🟢 -3m | 🔴 +25s | 🟢 -31m 35s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran basic-operations end-to-end against dev. All 31 scenario steps passed in the MCP browser phase (Section 1: structure… | 4m 27s | 9.0s | 50s | 5m 26s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 4m 27s (new) | 🟡 9.0s (new) | 🟡 50s (new) | 🟡 5m 26s (new) |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | Ran chem-and-bio scenario end-to-end against dev. All 11 scenario steps passed in the MCP browser phase (Chem: open spgi… | 2m 50s | 42s | 47s | 4m 19s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 2m 50s (new) | 🟡 42s (new) | 🟡 47s (new) | 🟡 4m 19s (new) |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 scenario steps PASSed on dev. spgi-100.csv loads correctly this time (previous run had to substitute SPGI.csv). C… | 3m 16s | 14s | 55s | 4m 25s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 3m 16s (new) | 🟡 14s (new) | 🟡 55s (new) | 🟡 4m 25s (new) |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed end-to-end on dev: table linking (SELECTION_TO_FILTER and FILTER_TO_FILTER) propagated correctly betw… | 1m 46s | 17s | 35s | 2m 38s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +1m 11s | 🔴 +14s | 🔴 +5.0s | 🔴 +1m 30s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran combined-boolean-filter end-to-end against dev. All 13 numbered scenario steps passed in the MCP browser phase: SEX_… | 2m 37s | 12s | 24s | 3m 13s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 2m 37s (new) | 🟡 12s (new) | 🟡 24s (new) | 🟡 3m 13s (new) |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (14/14) | All 14 steps passed in both the MCP run and the Playwright replay. Expression filter works correctly: 5-rule AND yields … | 1m 14s | 8.0s | 23s | 1m 45s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 1m 14s (new) | 🟡 8.0s (new) | 🟡 23s (new) | 🟡 1m 45s (new) |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 12 steps passed in the MCP run and in the Playwright replay (spec finished in 21.8s). The hierarchical filter correc… | 1m 15s | 21s | 23s | 1m 59s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 1m 15s (new) | 🟡 21s (new) | 🟡 23s (new) | 🟡 1m 59s (new) |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed in the MCP run and in the Playwright replay (spec finished in 8.7s, total wall-clock 11.56s). The tex… | 1m 12s | 20s | 12s | 1m 44s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 1m 12s (new) | 🟡 20s (new) | 🟡 12s (new) | 🟡 1m 44s (new) |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | 🟢 PASS → PASS | 🟢 100% (34/34) | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to reset… | 4m 24s | 40s | 1m 3s | 6m 7s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +2m 59s | 🔴 +35s | 🟢 -9.0s | 🔴 +3m 25s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 69% (5/8) | Color consistency through layout round-trip works — the `.categorical-colors` tag survives save/reload and `R_ONE` stays… | 2m 30s | 35s | 26s | 3m 31s | 🟡 69% (baseline) | 🟡 69% (baseline) | ⚪ 69 | 🔴 +2m | 🔴 +32s | 🔴 +2.0s | 🔴 +2m 34s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 83% (9/12) | Filtering legend updates work end-to-end in the MCP run: numeric filter, categorical filter, layout round-trip, composed… | 3m 10s | 1m 10s | 44s | 5m 4s | 🟡 83% (baseline) | 🟡 83% (baseline) | ⚪ 83 | 🔴 +2m 15s | 🔴 +1m 5s | 🔴 +8.1s | 🔴 +3m 28s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 77% (8/11) | Line chart legend and multi-axis behaviors are mostly correct: 7 legend items for 7 categories, layout round-trip preser… | 2m 10s | 40s | 32s | 3m 22s | 🟡 77% (baseline) | 🟡 77% (baseline) | ⚪ 77 | 🔴 +40s | 🔴 +37s | 🟢 -15s | 🔴 +1m 2s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 58% (7/13) | Categorical legend on scatter plot updates correctly when X axis changes (sub 2) and when the Filter Panel narrows categ… | 4m 15s | 1m 20s | 54s | 6m 29s | 🟡 58% (baseline) | 🟡 58% (baseline) | ⚪ 58 | 🔴 +3m 30s | 🔴 +1m 17s | 🔴 +4.0s | 🔴 +4m 51s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 79% (5/7) | Structure rendering in legends works for Scatter plot, Histogram, Line chart and Pie chart (canvas-based molecule thumbn… | 2m 35s | 40s | 29s | 3m 44s | 🟡 79% (baseline) | 🟡 79% (baseline) | ⚪ 79 | 🔴 +2m | 🔴 +37s | 🔴 +3.0s | 🔴 +2m 40s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 70% (13/20) | Scenario executed end-to-end with a mix of PASS, AMBIGUOUS, and FAIL. Legend display, source-swap, corner positioning, a… | 5m 45s | 1m 30s | 41s | 7m 56s | 🟡 70% (baseline) | 🟡 70% (baseline) | ⚪ 70 | 🔴 +4m 20s | 🔴 +1m 25s | 🟢 -2.7s | 🔴 +5m 42s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | 🟢 PASS → PASS | 🟡 77% (10/13) |  | 7m | 1m | 17s | 8m 17s | 🟡 77% (baseline) | 🟡 77% (baseline) | ⚪ 77 | ⚪ 0.0s | ⚪ 0.0s | ⚪ -0.2s | ⚪ -0.2s |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (82/82) | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color coding… | 3m 3s | 21s | 52s | 4m 16s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +3.0s | 🔴 +11s | 🟢 -8.0s | 🔴 +6.0s |
| Viewers | bar-chart-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 95% (18/19) | 17 of 19 sections passed cleanly; section 8 combined into section 7 in the spec. Section 18 is AMBIGUOUS — `grok.dapi.pr… | 1m 5s | 8.0s | 32s | 1m 45s | 🟢 95% (baseline) | 🟢 95% (baseline) | ⚪ 95 | 🟢 -1m 55s | 🟡 8.0s (new) | 🟢 -15s | 🟢 -2m 2s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selectio… | 25s | 1m | 9.4s | 1m 34s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | 🟢 PASS → PASS | 🟡 67% (8/12) | Core color coding steps passed: Linear (AGE, STARTED), Categorical (SEX, CONTROL) color coding applied and verified via … | 15s | 3.0s | 1m 24s | 1m 42s | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [color-coding-(linked)](Viewers/color-coding-(linked)-run.md) | 🟢 PASS → PASS | 🟡 83% (5/6) | All tested steps passed. Linked color coding works via column tags (`.color-coding-type`, `.color-coding-source-column`)… | 10s | 3.0s | 6.0s | 19s | 🟡 83% (baseline) | 🟡 83% (baseline) | ⚪ 83 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 87% (26/30) | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operation… | 5m | 3.0s | 22s | 5m 26s | 🟡 87% (baseline) | 🟡 87% (baseline) | ⚪ 87 | ⚪ 0.0s | ⚪ 0.0s | ⚪ +0.5s | ⚪ +0.5s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (58/58) | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI interacti… | 18m | 2m |  | 20m | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ — | ⚪ 0.0s |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | 🟡 93% (28/30) | All 14 sections of form-tests-pw.md exercised across 30 steps. 28 PASS, 2 AMBIGUOUS, 0 FAIL in MCP run. Playwright spec … | 18m | 4m | 3m 12s | 25m 12s | 🟡 93% (baseline) | 🟡 93% (baseline) | ⚪ 93 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas elemen… | 18m | 3m | 52s | 21m 52s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 73% (16/22) | Grid tests ran 22 steps (spec softSteps); 17 passed outright and 5 were AMBIGUOUS (Copy/Paste, Column Header Context Men… | 11m | 3m | 1m 18s | 15m 18s | 🟡 73% (baseline) | 🟡 73% (baseline) | ⚪ 73 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | 🟡 94% (15/16) | All 14 heat-map sections exercised across 17 steps. 15 PASS, 1 AMBIGUOUS, 1 SKIP in MCP run. Playwright spec passed full… | 18m | 4m | 49s | 22m 49s | 🟡 94% (baseline) | 🟡 94% (baseline) | ⚪ 94 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (87/94) | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance, l… | 50s | 7.0s | 46s | 1m 43s | 🟡 93% (baseline) | 🟡 93% (baseline) | ⚪ 93 | 🟢 -1m 10s | 🟢 -23s | 🟢 -4.0s | 🟢 -1m 37s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (26/26) | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout sa… | 57s | 8.0s | 1m 43s | 2m 48s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -1m 3s | 🔴 +3.0s | 🟢 -47s | 🟢 -1m 47s |
| Viewers | [map](Viewers/map-run.md) | 🟢 PASS → PASS | 🟡 80% (8/10) | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker size m… | 15s | 3.0s | 9.0s | 27s | 🟡 80% (baseline) | 🟡 80% (baseline) | ⚪ 80 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (16/19) | Matrix Plot tests ran with 15 PASS, 3 AMBIGUOUS, 0 FAIL. The spec executed in 57.7s with all implemented steps passing. … | 20m | 3m | 56s | 23m 56s | 🟡 84% (baseline) | 🟡 84% (baseline) | ⚪ 84 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (8/12) | 9 of 12 steps PASS; 3 SKIP (canvas-based node/edge interactions cannot be automated via DOM). The network diagram viewer… | 8m | 1m 30s | 22s | 9m 52s | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | 🔴 +1m 15s | 🔴 +1m 27s | 🔴 +15s | 🔴 +2m 57s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | All 13 scenario sections (mapped to 12 Playwright softSteps — scale and normalization are combined in the spec) passed d… | 1m 8s | 8.0s | 47s | 2m 3s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -6m 52s | 🔴 +3.0s | 🔴 +9.0s | 🟢 -6m 40s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (81/81) | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appearan… | 40s | 7.0s | 47s | 1m 34s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -3m 20s | 🟢 -3.0s | 🟢 -11s | 🟢 -3m 34s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (17/20) | Pivot Table tests ran with 16 PASS, 2 AMBIGUOUS, 1 SKIP, 0 FAIL. The spec executed in 35.1s with all implemented steps p… | 16m | 3m | 1m 12s | 20m 12s | 🟡 85% (baseline) | 🟡 85% (baseline) | ⚪ 85 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | rendering-structures-on-the-axes | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [row-source](Viewers/row-source-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with all 8… | 4m | 5.0s | 1m 24s | 5m 29s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (20/20) | All 20 sections passed during the MCP run on dev.datagrok.ai. The existing Playwright spec was re-run headed without mod… | 3m 13s | 29s | 52s | 4m 34s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -4m 47s | 🔴 +24s | 🟢 -14s | 🟢 -4m 37s |
| Viewers | scatter-plot-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | statistics-viewer | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | 🟢 PASS → PASS | 🟢 100% (24/24) | 24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to manual… | 4m | 3m | 58s | 7m 58s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (36/36) | All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via `v… | 28m | 4m | 46s | 32m 46s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (70/84) | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range slid… | 3m | 30s | 1m 48s | 5m 18s | 🟡 84% (baseline) | 🟡 84% (baseline) | ⚪ 84 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | viewers-docking | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All 7 MCP scenario steps PASS. The Word Cloud viewer adds via both entry points (Add-Viewer gallery and Toolbox icon), t… | 4m 15s | 1m | 2m 7s | 7m 22s | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +3m 55s | 🔴 +57s | 🔴 +1m 54s | 🔴 +6m 46s |
| Viewers | word-cloud-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | working-with-nan-&-infinity | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |

## Comparison with Previous Reports

Deltas are computed against two baselines pulled from git history of `total-run.md`: `prev1d` (latest commit before today) and `prev7d` (commit closest to today − 7 days, ±3-day window). Signed with `+`/`-`; time deltas use the same format as the values. All status and delta cells carry the Legend icons.

### Totals

**Total (1d)**: Tests Δ **-1** · Run Δ **🔴 -2** · Status **🟡 → 🟡 PARTIAL → PARTIAL** · Browser Δ **🟢 -19s** · Spec Gen Δ **🔴 +2.0s** · Spec Run Δ **⚪ +0.0s** · Total Δ **🟢 -16s**

**Total (7d)**: Spec Gen Δ **🔴 +42s** · Spec Run Δ **🔴 +10s** · Total Δ **🔴 +4m 33s** _(7d-only deltas — count and status deltas live in the 1d row to avoid double-counting.)_

### By Folder

| Folder | Tests Δ | Run Δ | Status | Mean Pass Δ (1d) | Mean Pass Δ (7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|
| Apps | +0 | +0 | ⚪ NO DATA → NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Bio | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 94% (new) | — | ⚪ 0.0s | ⚪ 0.0s | ⚪ — | ⚪ 0.0s |
| Browse | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 53% (new) | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 54% (new) | — | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 87% (new) | — | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Connections | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 54% (new) | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟢 96% (new) | — | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 67% (new) | — | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| General | +0 | +0 | ⚪ NO DATA → NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| LocalCashing | +0 | +0 | ⚪ NO DATA → NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 81% (new) | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Notebooks | +0 | +0 | ⚪ NO DATA → NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 68% (new) | — | ⚪ 0.0s | ⚪ 0.0s | ⚪ -0.1s | ⚪ -0.1s |
| PowerPack | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 56% (new) | — | 🔴 +3m 14s | 🔴 +36s | 🔴 +19s | 🔴 +4m 11s |
| Projects | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 38% (new) | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Queries | +0 | +0 | ⚪ NO DATA → NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 90% (new) | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | +0 | +0 | 🟡 PARTIAL → PARTIAL | 🟡 68% (new) | — | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Tooltips | +1 | +0 | ⚪ NO DATA → NO DATA | — | — | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Viewers | -2 | -2 | 🟡 PARTIAL → PARTIAL | 🟡 90% (new) | — | 🟢 -1m 22s | 🟢 -1.1s | 🟢 -1.9s | 🟢 -1m 23s |

### Per-Test Changes

| Folder | Test | Status | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 80% (baseline) | 🟡 80% (baseline) | ⚪ 80 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 90% (baseline) | 🟡 90% (baseline) | ⚪ 90 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (baseline) | 🟡 50% (baseline) | ⚪ 50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [spaces](Browse/spaces-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 24% (baseline) | 🟡 24% (baseline) | ⚪ 24 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | 🔴 FAIL → FAIL | 🟡 3% (baseline) | 🟡 3% (baseline) | ⚪ 3 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 42% (baseline) | 🟡 42% (baseline) | ⚪ 42 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 20% (baseline) | 🟡 20% (baseline) | ⚪ 20 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (baseline) | 🟡 50% (baseline) | ⚪ 50 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | 🟡 33% (baseline) | 🟡 33% (baseline) | ⚪ 33 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (baseline) | 🟡 40% (baseline) | ⚪ 40 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Connections | [adding](Connections/adding-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (baseline) | 🟡 93% (baseline) | ⚪ 93 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [browser](Connections/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 78% (baseline) | 🟡 78% (baseline) | ⚪ 78 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | 🟡 9% (baseline) | 🟡 9% (baseline) | ⚪ 9 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 FAIL → FAIL | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 FAIL → FAIL | 🟡 11% (baseline) | 🟡 11% (baseline) | ⚪ 11 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [import-swagger](Connections/import-swagger-run.md) | 🔴 FAIL → FAIL | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [schema](Connections/schema-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 75% (baseline) | 🟡 75% (baseline) | ⚪ 75 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 86% (baseline) | 🟡 86% (baseline) | ⚪ 86 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (baseline) | 🟡 83% (baseline) | ⚪ 83 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/pls-regression](EDA/ML methods/pls-regression-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 71% (baseline) | 🟡 71% (baseline) | ⚪ 71 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/softmax](EDA/ML methods/softmax-run.md) | 🔴 FAIL → FAIL | 🟡 33% (baseline) | 🟡 33% (baseline) | ⚪ 33 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/xgboost1](EDA/ML methods/xgboost1-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [ML methods/xgboost2](EDA/ML methods/xgboost2-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 43% (baseline) | 🟡 43% (baseline) | ⚪ 43 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 60% (baseline) | 🟡 60% (baseline) | ⚪ 60 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🟡 62% (baseline) | 🟡 62% (baseline) | ⚪ 62 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🟡 28% (baseline) | 🟡 28% (baseline) | ⚪ 28 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (baseline) | 🟡 40% (baseline) | ⚪ 40 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ -0.4s | ⚪ -0.4s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 60% (baseline) | 🟡 60% (baseline) | ⚪ 60 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | 🟡 SKIP → SKIP | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | 🟡 SKIP → SKIP | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🟡 → 🔴 PARTIAL → FAIL | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | 🔴 +6m 28s | 🔴 +1m 12s | 🟡 43s (new) | 🔴 +8m 23s |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | 🟡 SKIP → SKIP | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [browser](Projects/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 56% (baseline) | 🟡 56% (baseline) | ⚪ 56 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [complex](Projects/complex-run.md) | 🟡 SKIP → SKIP | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 SKIP → SKIP | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (baseline) | 🟡 50% (baseline) | ⚪ 50 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 SKIP → SKIP | 🔴 0% (baseline) | 🔴 0% (baseline) | ⚪ 0 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 57% (baseline) | 🟡 57% (baseline) | ⚪ 57 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 92% (baseline) | 🟡 92% (baseline) | ⚪ 92 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | ⚪ — | ⚪ — | ⚪ — | ⚪ — |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | 🟡 20% (baseline) | 🟡 20% (baseline) | ⚪ 20 | ⚪ 0.0s | ⚪ 0.0s | ⚪ — | ⚪ 0.0s |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -29m | 🟢 -3m | 🔴 +25s | 🟢 -31m 35s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 4m 27s (new) | 🟡 9.0s (new) | 🟡 50s (new) | 🟡 5m 26s (new) |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 2m 50s (new) | 🟡 42s (new) | 🟡 47s (new) | 🟡 4m 19s (new) |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 3m 16s (new) | 🟡 14s (new) | 🟡 55s (new) | 🟡 4m 25s (new) |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +1m 11s | 🔴 +14s | 🔴 +5.0s | 🔴 +1m 30s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 2m 37s (new) | 🟡 12s (new) | 🟡 24s (new) | 🟡 3m 13s (new) |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 1m 14s (new) | 🟡 8.0s (new) | 🟡 23s (new) | 🟡 1m 45s (new) |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 1m 15s (new) | 🟡 21s (new) | 🟡 23s (new) | 🟡 1m 59s (new) |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟡 1m 12s (new) | 🟡 20s (new) | 🟡 12s (new) | 🟡 1m 44s (new) |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +2m 59s | 🔴 +35s | 🟢 -9.0s | 🔴 +3m 25s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 69% (baseline) | 🟡 69% (baseline) | ⚪ 69 | 🔴 +2m | 🔴 +32s | 🔴 +2.0s | 🔴 +2m 34s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 83% (baseline) | 🟡 83% (baseline) | ⚪ 83 | 🔴 +2m 15s | 🔴 +1m 5s | 🔴 +8.1s | 🔴 +3m 28s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 77% (baseline) | 🟡 77% (baseline) | ⚪ 77 | 🔴 +40s | 🔴 +37s | 🟢 -15s | 🔴 +1m 2s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 58% (baseline) | 🟡 58% (baseline) | ⚪ 58 | 🔴 +3m 30s | 🔴 +1m 17s | 🔴 +4.0s | 🔴 +4m 51s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 79% (baseline) | 🟡 79% (baseline) | ⚪ 79 | 🔴 +2m | 🔴 +37s | 🔴 +3.0s | 🔴 +2m 40s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟡 70% (baseline) | 🟡 70% (baseline) | ⚪ 70 | 🔴 +4m 20s | 🔴 +1m 25s | 🟢 -2.7s | 🔴 +5m 42s |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +3.0s | 🔴 +11s | 🟢 -8.0s | 🔴 +6.0s |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 95% (baseline) | 🟢 95% (baseline) | ⚪ 95 | 🟢 -1m 55s | 🟡 8.0s (new) | 🟢 -15s | 🟢 -2m 2s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 87% (baseline) | 🟡 87% (baseline) | ⚪ 87 | ⚪ 0.0s | ⚪ 0.0s | ⚪ +0.5s | ⚪ +0.5s |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 73% (baseline) | 🟡 73% (baseline) | ⚪ 73 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (baseline) | 🟡 93% (baseline) | ⚪ 93 | 🟢 -1m 10s | 🟢 -23s | 🟢 -4.0s | 🟢 -1m 37s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -1m 3s | 🔴 +3.0s | 🟢 -47s | 🟢 -1m 47s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (baseline) | 🟡 84% (baseline) | ⚪ 84 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (baseline) | 🟡 67% (baseline) | ⚪ 67 | 🔴 +1m 15s | 🔴 +1m 27s | 🔴 +15s | 🔴 +2m 57s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -6m 52s | 🔴 +3.0s | 🔴 +9.0s | 🟢 -6m 40s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -3m 20s | 🟢 -3.0s | 🟢 -11s | 🟢 -3m 34s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (baseline) | 🟡 85% (baseline) | ⚪ 85 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🟢 -4m 47s | 🔴 +24s | 🟢 -14s | 🟢 -4m 37s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (baseline) | 🟡 84% (baseline) | ⚪ 84 | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s | ⚪ 0.0s |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | 🟢 PASS → PASS | 🟢 100% (baseline) | 🟢 100% (baseline) | ⚪ 100 | 🔴 +3m 55s | 🔴 +57s | 🔴 +1m 54s | 🔴 +6m 46s |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage is 74% (≥70%) and no folders are fully failing, but 14 folder(s) are PARTIAL:
- Bio, Browse, Charts, Chem, Connections, DiffStudio, EDA, Models, Peptides, PowerPack, Projects, Scripts, StickyMeta, Viewers

### Blocking Issues
- Bio: PARTIAL (Bio/composition-analysis, Bio/pepsea)
- Browse: PARTIAL (Browse/browse-tree-states, Browse/browse, Browse/spaces-(ui-only), Browse/spaces)
- Charts: PARTIAL (Charts/sunburst, Charts/tree)
- Chem: PARTIAL (Chem/Advanced/scaffold-tree, Chem/calculate, Chem/chemprop)
- Connections: PARTIAL (Connections/adding, Connections/browser, Connections/catalogs, Connections/external-provider, Connections/identifiers, Connections/import-swagger, …)
- DiffStudio: PARTIAL (DiffStudio/fitting)
- EDA: PARTIAL (EDA/ML methods/pls-regression, EDA/ML methods/softmax, EDA/ML methods/xgboost1, EDA/ML methods/xgboost2, EDA/multivariate-analysis, EDA/pareto-front-viewer, …)
- Models: PARTIAL (Models/browser, Models/chemprop)
- Peptides: PARTIAL (Peptides/peptide-space, Peptides/peptides, Peptides/sar)
- PowerPack: PARTIAL (PowerPack/AddNewColumn/add-new-column, PowerPack/data-enrichment)
- Projects: PARTIAL (Projects/browser, Projects/deleting, Projects/opening, Projects/uploading)
- Scripts: PARTIAL (Scripts/create, Scripts/run)
- StickyMeta: PARTIAL (StickyMeta/database-meta)
- Viewers: PARTIAL (Viewers/Legend/color-consistency, Viewers/Legend/filtering, Viewers/Legend/line-chart, Viewers/Legend/scatterplot, Viewers/Legend/structure-rendering, Viewers/Legend/visibility-and-positioning, …)
