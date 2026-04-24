# Test Track — Global Report

**Date**: 2026-04-24
**Legend**: 🟢 improvement/pass · 🔴 regression/fail · 🟡 partial/ambiguous/new-missing-delta · ⚪ no change / no data
**Verdict**: Conditionally ready

## Column definitions

- **Pass %** — `(pass + 0.5·partial) / total_steps` from each run's `## Steps` section.
- **Browser** — MCP phase wall-clock = model thinking time + live browser interaction time.
- **Spec Gen** — model-only time to generate the Playwright spec.
- **Spec Run** — Playwright-only spec execution time.
- **Total** (per test) — sum of Browser + Spec Gen + Spec Run for that scenario. **Mean Total** = average of those sums.
- **Pass Δ (1d)** — Pass % change vs. `prev1d` (the most recent committed `total-run.md` strictly before today).
- **Pass Δ (7d)** — Pass % change vs. `prev7d` (the committed `total-run.md` closest to today − 7 days, ±3-day window). Empty when no commit falls in that window or the baseline commit predates Pass %.
- **Trend** — last ≤7 daily Pass % values, oldest → newest, dot-separated. Prefix icon: 🟢 last > first, 🔴 last < first, ⚪ equal or only one point.

## Folder Summary

**Total**: 177 tests · Run: 136/177 (77%) · Playwright: 125/177 (71%) · Mean Pass: 🟡 81% · Mean Browser: 3m 32s · Mean Spec Gen: 54.1s · Mean Spec Run: 45s · Mean Total (sum per test): 4m 48s

| Folder | Tests | Run | Playwright | Status | Mean Pass % | Mean Browser | Mean Spec Gen | Mean Spec Run | Mean Total |
|---|---|---|---|---|---|---|---|---|---|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Bio | 9 | 9/9 (100%) | 9/9 (100%) | 🟡 PARTIAL | 🟢 100% | 2m 43s | 43.2s | 43.4s | 4m |
| Browse | 3 | 2/3 (67%) | 0/3 (0%) | 🟡 PARTIAL | 🟡 70% |  |  |  |  |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | 🟡 PARTIAL | 🟡 54% | 2m 10s | 45s | 34.3s | 3m 29s |
| Chem | 14 | 14/14 (100%) | 14/14 (100%) | 🟡 PARTIAL | 🟡 87% | 1m 37s | 31.8s | 40.4s | 2m 46s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | 🟡 PARTIAL | 🟡 54% |  |  |  |  |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | 🟡 PARTIAL | 🟢 96% | 3m 10s | 56.4s | 45.4s | 4m 52s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | 🟡 PARTIAL | 🟡 67% | 1m 42s | 58s | 17.5s | 2m 58s |
| General | 10 | 0/10 (0%) | 0/10 (0%) | ⚪ NO DATA |  |  |  |  |  |
| LocalCashing | 0 | 0/0 (0%) | 0/0 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | 🟡 PARTIAL | 🟡 55% | 8m 20s | 2m 3s | 1m 45s | 12m 8s |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 68% | 27.5s | 3s | 47.6s | 1m 18s |
| PowerPack | 9 | 9/9 (100%) | 9/9 (100%) | 🟡 PARTIAL | 🟢 97% | 5m 16s | 1m 9s | 29.8s | 6m 54s |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | 🟡 PARTIAL | 🟡 38% |  |  |  |  |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Scripts | 6 | 5/6 (83%) | 5/6 (83%) | 🟡 PARTIAL | 🟡 87% | 1m 17s | 25s | 16s | 1m 58s |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | 🟡 PARTIAL | 🟡 77% | 5m 45s | 1m 26s | 56.2s | 8m 8s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | ⚪ NO DATA |  |  |  |  |  |
| Viewers | 46 | 44/46 (96%) | 44/46 (96%) | 🟡 PARTIAL | 🟡 92% | 4m 8s | 52.3s | 48.6s | 4m 42s |
| WideSmokeTest | 0 | 0/0 (0%) | 0/0 (0%) | ⚪ NO DATA |  |  |  |  |  |

## All Tests

**Total**: 177 tests · 🟢 74 PASS / 🟡 47 PARTIAL / 🔴 11 FAIL / 🟡 4 SKIP / ⚪ 41 NO RUN · Mean Pass: 🟡 81% · Mean Browser: 3m 32s · Mean Spec Gen: 54.1s · Mean Spec Run: 45s · Mean Total (sum per test): 4m 48s

| Folder | Test | Status | Pass % | Description | Browser (model+MCP) | Spec Gen (model) | Spec Run (Playwright) | Total (sum) | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Apps | apps | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Apps | tutorials | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Bio | [analyze](Bio/analyze-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset … | 2m 40s | 20s | 1m 21s | 4m 21s | ⚪ +0% (11/11 → 11/11) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟢 -5m 20s | 🔴 +15s | 🟡 1m 21s (new) | 🟢 -3m 44s |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 100% (5/5) | All 15 checks pass in both browser and Playwright (5 scenario steps × 3 datasets): WebLogo opens from **Bio → Analyze… | 4m 13s | 1m 29s | 41s | 6m 23s | 🟢 +20% (4/5 → 5/5) | 🟢 100% · status PARTIAL → PASS | 🟢 80·80·100 | 🟡 4m 13s (new) | 🟡 1m 29s (new) | 🟡 41s (new) | 🟡 6m 23s (new) |
| Bio | [convert](Bio/convert-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 sub-steps passed in both the MCP run and the Playwright spec. The key changes that closed the two outstanding … | 4m 40s | 1m 10s | 1m 42s | 7m 32s | ⚪ +0% (5/5 → 15/15) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +2m 40s | 🟡 1m 10s (new) | 🟡 1m 42s (new) | 🔴 +5m 32s |
| Bio | [manage](Bio/manage-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All three scenario steps passed in the MCP browser run against dev.datagrok.ai: `HELM.csv` opened (540 rows, HELM → M… | 36s | 25s | 22s | 1m 23s | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +6s | 🟡 25s (new) | 🟡 22s (new) | 🔴 +53s |
| Bio | [msa](Bio/msa-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | The MSA scenario works end-to-end. In the live browser the new "Clusters" column is added via Edit > Add New Column w… | 2m 50s | 40s | 19s | 3m 49s | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 2m 50s (new) | 🟡 40s (new) | 🟡 19s (new) | 🟡 3m 49s (new) |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 → 🟡 AMBIGUOUS → PARTIAL | 🟢 100% (5/5) | Steps 1–5 pass end-to-end against dev: 50-row HELM subset opens with `Macromolecule`/`helm`, Clusters column is added… | 3m 15s | 30s | 25s | 4m 10s | 🟢 +33% (4/6 → 5/5) | 🟢 100% · status AMBIGUOUS → PARTIAL | 🟢 67·67·100 | 🟡 3m 15s (new) | 🟡 30s (new) | 🟡 25s (new) | 🟡 4m 10s (new) |
| Bio | [search](Bio/search-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All four steps passed on the dev server with both the MCP reproduction and the generated Playwright spec. The scenari… | 1m | 25s | 14s | 1m 39s | ⚪ +0% (4/4 → 4/4) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +30s | 🟡 25s (new) | 🟡 14s (new) | 🔴 +1m 9s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets). Activity Cli… | 2m 26s | 1m |  | 3m 26s | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 2m 26s (new) | 🟡 1m (new) | — | 🟡 3m 26s (new) |
| Bio | [sequence-space](Bio/sequence-space-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 scenario steps PASS in both the live grok-browser run and the Playwright spec replay (3 datasets: FASTA, HELM, … | 2m 50s | 30s |  | 3m 20s | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 2m 50s (new) | 🟡 30s (new) | — | 🟡 3m 20s (new) |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 90% (4/5) | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. I… |  |  |  |  | ⚪ +0% (4/5 → 4/5) | 🟡 90% · status unchanged (PARTIAL) | ⚪ 90·90·90 | — | — | — | — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (0/1) | 1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a single sess… |  |  |  |  | ⚪ +0% (0/1 → 0/1) | 🟡 50% · status unchanged (PARTIAL) | ⚪ 50·50·50 | — | — | — | — |
| Browse | package-manager | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Charts | [radar](Charts/radar-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | The Radar viewer reproduced cleanly on dev for both earthquakes.csv (2426 rows) and demog.csv (5850 rows). All 21 Rad… | 1m 30s | 35s | 33s | 2m 38s | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 42% (5/12) | The sunburst viewer reproduces structurally on dev — Sunburst can be added to both SPGI and demog, and `hierarchyColu… | 3m | 1m | 34s | 4m 34s | ⚪ +0% (5/12 → 5/12) | 🟡 42% · status unchanged (PARTIAL) | ⚪ 42·42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 20% (1/5) | Setup (open demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduced cleanly on dev. All four test steps are … | 2m | 40s | 36s | 3m 16s | ⚪ +0% (1/5 → 1/5) | 🟡 20% · status PASS → PARTIAL | ⚪ 20·20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (3/6) | Smoke coverage only: Scaffold Tree viewer launches from the Chem menu and the magic wand generates a scaffold tree on… | 57s | 25s | 49.6s | 2m 12s | ⚪ +0% (3/6 → 3/6) | 🟡 new (prev7d had no entry) | ⚪ 50·50·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Scaffold Tree viewer launches from the Chem → Scaffold Tree menu, and the magic-wand generator produces scaffold node… | 1m 15s | 30s | 40.2s | 2m 25s | ⚪ +0% (3/3 → 3/3) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Similarity Search launches from the Chem menu and exposes a viewer that accepts option changes (fingerprint Morgan ↔ … | 37s | 25s | 24.8s | 1m 27s | ⚪ +0% (3/3 → 3/3) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Substructure filtering via `grok.chem.searchSubstructure` works on SPGI.csv (3624 rows): benzene substructure yields … | 38s | 25s | 29s | 1m 32s | ⚪ +0% (4/4 → 4/4) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Activity Cliffs computation on SPGI.csv (3624 rows) finishes within 45s and produces a UMAP scatter plot with molecul… | 1m 14s | 20s | 1m 4s | 2m 38s | ⚪ +0% (4/4 → 4/4) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | 🟡 33% (1/3) | Calculate Descriptors cannot be exercised on `dev` right now. The Chem top menu fails to open its popup — both throug… | 8m | 1m |  | 9m | ⚪ +0% (1/3 → 1/3) | 🟡 33% · status PASS → FAIL | ⚪ 33·33·33 | ⚪ +0s | ⚪ +0s | 🟡 38s (removed) | 🟢 -38s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Chemical Space dimensional reduction runs end-to-end on smiles.csv: the dialog opens, OK with defaults produces a Sca… | 46s | 20s | 59.8s | 2m 6s | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | ChemProp scenario is partially automated: the spec confirms mol1K.sdf opens and the Train Model view is reachable fro… | 26s | 30s | 19.1s | 1m 15s | ⚪ +0% (2/5 → 2/5) | 🟡 40% · status unchanged (PARTIAL) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Elemental Analysis works on dev. The menu path `[name="div-Chem"]` → `Elemental Analysis...` resolves and the dialog … | 38s | 30s | 30.9s | 1m 39s | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | 🟢 PASS → PASS | 🟢 100% (2/2) | The filter panel correctly shows a Structure filter for SPGI.csv's Molecule column; clicking the embedded sketch-link… | 34s | 25s | 24.2s | 1m 23s | ⚪ +0% (2/2 → 2/2) | 🟢 100% · status PARTIAL → PASS | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [info-panels](Chem/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | Info panels work correctly on smiles.csv: column-level (Details, Filter, Colors, Style, Chemistry with Rendering/High… | 3m 10s | 1m | 31s | 4m 41s | ⚪ +0% (5/5 → 5/5) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [mmp](Chem/mmp-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset at the bottom of the v… | 1m 27s | 20s | 1m 14s | 3m 1s | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status FAIL → PASS | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | R-Groups Analysis works on sar_small.csv: MCS auto-populates the sketcher, OK produces a Trellis plot and appends R1–… | 2m 20s | 45s | 58.7s | 4m 4s | ⚪ +0% (5/5 → 5/5) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [sketcher](Chem/sketcher-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | Sketcher opens via `grok.chem.sketcher(molCol, initialSmiles)` wrapped in `ui.dialog(...).show()`, accepts a typed SM… | 33s | 30s | 19.3s | 1m 22s | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [adding](Connections/adding-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (6/7) | 6 of 7 steps fully passed, Step 5 was partial (TEST button works but actual connection test fails without real creden… |  |  |  |  | ⚪ +0% (6/7 → 6/7) | 🟡 93% · status unchanged (PARTIAL) | ⚪ 93·93·93 | — | — | — | — |
| Connections | [browser](Connections/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 78% (7/9) | 7 of 9 steps passed, 2 ambiguous. Search filtering works correctly and the Context Pane shows all expected tabs (Deta… |  |  |  |  | ⚪ +0% (7/9 → 7/9) | 🟡 78% · status unchanged (PARTIAL) | ⚪ 78·78·78 | — | — | — | — |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | 🟡 9% (1/11) | 1 step passed, 1 failed, 15 skipped. The required `NorthwindTest` MS SQL connection is not present on public.datagrok… |  |  |  |  | ⚪ +0% (1/11 → 1/11) | 🟡 9% · status unchanged (FAIL) | ⚪ 9·9·9 | — | — | — | — |
| Connections | [delete](Connections/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 steps passed. Both connections were deleted successfully. The confirmation dialog uses a red "DELETE" button (n… |  |  |  |  | ⚪ +0% (8/8 → 8/8) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | — | — | — | — |
| Connections | [edit](Connections/edit-run.md) | 🟢 PASS → PASS | 🟡 86% (6/7) | 6 of 7 steps passed (1 skipped due to missing real credentials). The connection rename, credential modification, and … |  |  |  |  | ⚪ +0% (6/7 → 6/7) | 🟡 86% · status unchanged (PASS) | ⚪ 86·86·86 | — | — | — | — |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 FAIL → FAIL | 🔴 0% (0/7) | All 7 steps skipped. This scenario requires a specific Postgres connection at db.datagrok.ai:54327 with superuser cre… |  |  |  |  | ⚪ +0% (0/7 → 0/7), still broken | 🔴 0% · status unchanged (FAIL) | ⚪ 0·0·0 | — | — | — | — |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 FAIL → FAIL | 🟡 11% (1/9) | 1 step passed, 1 failed, 7 skipped. This scenario depends on a working Postgres connection to the Northwind database.… |  |  |  |  | ⚪ +0% (1/9 → 1/9) | 🟡 11% · status unchanged (FAIL) | ⚪ 11·11·11 | — | — | — | — |
| Connections | [import-swagger](Connections/import-swagger-run.md) | 🔴 FAIL → FAIL | 🔴 0% (0/7) | All 7 steps skipped. This scenario requires manual interaction: downloading a YAML file to the local machine and drag… |  |  |  |  | ⚪ +0% (0/7 → 0/7), still broken | 🔴 0% · status unchanged (FAIL) | ⚪ 0·0·0 | — | — | — | — |
| Connections | [schema](Connections/schema-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 75% (3/4) | 3 of 4 steps passed, 1 ambiguous. The "Browse schema" context menu option was not found in the current UI, but the sc… |  |  |  |  | ⚪ +0% (3/4 → 3/4) | 🟡 75% · status unchanged (PARTIAL) | ⚪ 75·75·75 | — | — | — | — |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 86% (6/7) | 6 of 7 steps passed (1 failed). All UI steps worked correctly. The SPARQL connection was created and deleted successf… |  |  |  |  | ⚪ +0% (6/7 → 6/7) | 🟡 86% · status unchanged (PARTIAL) | ⚪ 86·86·86 | — | — | — | — |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | Catalog scenario reproduces fully on dev.datagrok.ai. All 6 steps PASS both in MCP and in the Playwright spec (35.8s … | 1m 53s | 39s | 38s | 3m 10s | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Cyclic Models (PK-PD) scenario reproduces fully on dev.datagrok.ai. The PK-PD library model loads via double-click, M… | 1m 2s | 38s | 36s | 2m 16s | ⚪ +0% (4/4 → 4/4) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Files & Sharing scenario reproduces fully on dev.datagrok.ai. pk.ivp loads via the `DiffStudio:previewIvp` function w… | 2m 31s | 1m 32s | 1m 17s | 5m 20s | ⚪ +0% (4/4 → 4/4) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (5/6) | Scenario is PARTIAL on dev.datagrok.ai — steps 1–5 pass, step 6 (actually running the fit) does not produce result ro… | 10m 2s | 55s | 1m | 11m 57s | ⚪ +0% (5/6 → 5/6) | 🟡 83% · status PASS → PARTIAL | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | 🟢 PASS → PASS | 🟡 83% (5/6) | Scenario fully reproduces on dev.datagrok.ai. All 6 steps pass in the interactive MCP session and in the Playwright s… | 1m 31s | 36s | 26s | 2m 33s | ⚪ +0% (5/6 → 5/6) | 🟡 83% · status unchanged (PASS) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 scenario steps PASS against dev.datagrok.ai. Edit toggle is reachable via `.d4-ribbon-item .ui-input-bool-switc… | 4m 30s | 1m 40s | 1m 3s | 7m 13s | ⚪ +0% (5/5 → 5/5) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Sensitivity Analysis scenario fully reproduces against dev.datagrok.ai. Bioreactor loads from the DiffStudio hub (lib… | 2m 3s | 44s | 39s | 3m 26s | ⚪ +0% (4/4 → 4/4) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | Stages (Acid Production) scenario reproduces fully on dev.datagrok.ai. The library card opens a view named "GA-produc… | 1m 49s | 47s | 24s | 3m | ⚪ +0% (4/4 → 4/4) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | ML methods/linear-regression | 🟢 removed: PASS |  |  |  |  |  |  | — | — | — | — | — | — | — |
| EDA | ML methods/pls-regression | 🟡 removed: PARTIAL |  |  |  |  |  |  | — | — | — | — | — | — | — |
| EDA | ML methods/softmax | 🔴 removed: FAIL |  |  |  |  |  |  | — | — | — | — | — | — | — |
| EDA | ML methods/xgboost1 | 🟡 removed: PARTIAL |  |  |  |  |  |  | — | — | — | — | — | — | — |
| EDA | ML methods/xgboost2 | 🟡 removed: PARTIAL |  |  |  |  |  |  | — | — | — | — | — | — | — |
| EDA | [MLMethods/linear-regression](EDA/MLMethods/linear-regression-run.md) | 🟢 new: PASS | 🟢 100% (3/3) | Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation… | 45s | 2s | 6.8s | 53.8s | 🟢 100% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 100 | 🟡 45s (new) | 🟡 2s (new) | 🟡 6.8s (new) | 🟡 53.8s (new) |
| EDA | [MLMethods/pls-regression](EDA/MLMethods/pls-regression-run.md) | 🟡 new: PARTIAL | 🟡 71% (5/7) | PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-… | 1m | 2s | 6.9s | 1m 9s | 🟡 71% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 71 | 🟡 1m (new) | 🟡 2s (new) | 🟡 6.9s (new) | 🟡 1m 9s (new) |
| EDA | [MLMethods/softmax](EDA/MLMethods/softmax-run.md) | 🔴 new: FAIL | 🟡 33% (1/3) | Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns an… | 10s | 2s | 2.6s | 14.6s | 🟡 33% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 33 | 🟡 10s (new) | 🟡 2s (new) | 🟡 2.6s (new) | 🟡 14.6s (new) |
| EDA | [MLMethods/xgboost1](EDA/MLMethods/xgboost1-run.md) | 🟡 new: PARTIAL | 🟡 67% (2/3) | XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as… | 5s | 2s | 2.7s | 9.7s | 🟡 67% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 67 | 🟡 5s (new) | 🟡 2s (new) | 🟡 2.7s (new) | 🟡 9.7s (new) |
| EDA | [MLMethods/xgboost2](EDA/MLMethods/xgboost2-run.md) | 🟡 new: PARTIAL | 🟡 67% (2/3) | XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interac… | 5s | 2s | 2.7s | 9.7s | 🟡 67% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 67 | 🟡 5s (new) | 🟡 2s (new) | 🟡 2.7s (new) | 🟡 9.7s (new) |
| EDA | [anova](EDA/anova-run.md) | 🟢 PASS → PASS | 🟢 100% (3/3) | All 3 scenario steps passed against dev. Dataset opens via JS API in ~1s; ANOVA dialog mounts with sensible defaults … | 1m 30s | 30s | 26s | 2m 26s | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (2/3) | 2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check, where the wording does not s… | 2m 30s | 2m | 13s | 4m 43s | ⚪ +0% (2/3 → 2/3) | 🟡 67% · status PASS → PARTIAL | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 43% (3/7) | 3 of 7 steps passed, 1 failed, 3 were skipped due to the missing prerequisite dataset. The Pareto Front viewer itself… | 4m | 2m | 32s | 6m 32s | ⚪ +0% (3/7 → 3/7) | 🟡 43% · status unchanged (PARTIAL) | ⚪ 43·43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 60% (3/5) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 3 PASS / 1 FAIL / 1 SKIP. The dialog path works (menu… | 5m | 3m | 1m 7s | 9m 7s | ⚪ +0% (3/5 → 3/5) | 🟡 60% · status PASS → PARTIAL | ⚪ 60·60·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🟡 62% (2/4) | MCP reproduction (phase 2b) on https://dev.datagrok.ai produced 2 PASS / 1 PARTIAL / 1 FAIL. The dialog path (menu, U… | 2m | 2m | 15s | 4m 15s | 🟢 +12% (2/5 → 2/4) | 🟡 62% · status PASS → FAIL | ⚪ 62·50·62 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| General | files-cache | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | first-login | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | inactivity-response | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | login | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | molecule-in-exported-csv | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | network | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | profile-settings | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | startup-time | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | table-manager | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| General | tabs-reordering | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Models | [apply](Models/apply-run.md) | 🟢 → 🔴 PASS → FAIL | 🟡 50% (2/4) | Steps 1 and 2 passed — demog.csv opened and the "Apply predictive model" dialog opened correctly via the top menu. St… | 2m 12s | 30s | 20s | 3m 2s | 🔴 -50% (4/4 → 2/4) | 🟡 50% · status PASS → FAIL | 🔴 100·100·50 | 🟡 2m 12s (new) | 🟡 30s (new) | 🟡 20s (new) | 🟡 3m 2s (new) |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 33% (2/6) | 2 of 6 steps pass fully in both browser and Playwright (step 1: navigation; step 4: Filter Templates). Steps 2, 3, 5,… | 4m 45s | 2m 10s | 51s | 7m 46s | 🔴 -34% (4/6 → 2/6) | 🟡 33% · status unchanged (PARTIAL) | 🔴 67·67·33 | 🟡 4m 45s (new) | 🟡 2m 10s (new) | 🟡 51s (new) | 🟡 7m 46s (new) |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🟡 35% (6/17) | Only 3 of 16 sub-steps pass fully in both browser and Playwright (1.1 open `smiles.csv`, 1.2 open Train view, 2.1 ope… | 13m | 2m | 3m 45s | 18m 45s | 🟢 +7% (5/18 → 6/17) | 🟡 35% · status unchanged (FAIL) | 🟢 28·28·35 | 🟡 13m (new) | 🟡 2m (new) | 🟡 3m 45s (new) | 🟡 18m 45s (new) |
| Models | [delete](Models/delete-run.md) | 🟢 → 🔴 PASS → FAIL | 🟡 20% (1/5) | The delete UI itself could not be exercised on dev because the prerequisite predictive model from train.md/browser.md… | 18m | 4m | 24s | 22m 24s | 🔴 -80% (5/5 → 1/5) | 🟡 20% · status PASS → FAIL | 🔴 100·100·20 | 🟡 18m (new) | 🟡 4m (new) | 🟡 24s (new) | 🟡 22m 24s (new) |
| Models | [predictive-models](Models/predictive-models-run.md) | 🟢 → 🟡 PASS → PARTIAL | 🟢 100% (17/17) | All 17 scenario steps **PASSED** in the interactive MCP browser run (Scenario 1 Train, Scenario 2 Apply, Scenario 3 A… | 5m 50s | 1m 50s | 4m 40s | 12m 20s | ⚪ +0% (21/21 → 17/17) | 🟢 100% · status PASS → PARTIAL | ⚪ 100·100·100 | 🟡 5m 50s (new) | 🟡 1m 50s (new) | 🟡 4m 40s (new) | 🟡 12m 20s (new) |
| Models | [train](Models/train-run.md) | 🟢 PASS → PASS | 🟡 91% (10/11) | All 10 scenario steps passed in the interactive MCP browser run. Both models — classification (TestDemog, Predict pro… | 6m 10s | 1m 50s | 28s | 8m 28s | 🟢 +1% (9/10 → 10/11) | 🟡 91% · status unchanged (PASS) | 🟢 90·90·91 | 🟡 6m 10s (new) | 🟡 1m 50s (new) | 🟡 28s (new) | 🟡 8m 28s (new) |
| Notebooks | browser | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Notebooks | create | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Notebooks | delete | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Notebooks | edit | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Peptides | [info-panels](Peptides/info-panels-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are re… | 17s | 3s | 10.9s | 30.9s | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 40% (2/5) | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variabil… | 25s | 3s | 1m 18s | 1m 46s | ⚪ +0% (2/5 → 2/5) | 🟡 40% · status unchanged (PARTIAL) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | 🟡 1m 18s (new) | 🔴 +1m 18s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (4/6) | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clust… | 18s | 3s | 11.6s | 32.6s | ⚪ +0% (4/6 → 4/6) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Res… | 50s | 3s | 1m 30s | 2m 23s | ⚪ +0% (6/9 → 6/9) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | 🟡 1m 30s (new) | 🔴 +1m 30s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (10/10) | All 10 scenario steps reproduce successfully in the MCP run; 9/10 pass in the Playwright replay. The one FAILED Playw… | 7m 45s | 2m | 32s | 10m 17s | ⚪ +0% (10/10 → 10/10) | 🟡 new (prev7d had no entry) | 🟢 60·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All six autocomplete behaviours PASS in both the MCP run and the Playwright replay. `.cm-tooltip-autocomplete` appear… | 1m 45s | 20s | 8s | 2m 13s | ⚪ +0% (7/7 → 7/7) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All seven sub-steps pass in both the MCP run and the Playwright replay. Dependency propagation across calculated colu… | 2m 15s | 1m | 28s | 3m 43s | ⚪ +0% (7/7 → 7/7) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/functions-sorting](PowerPack/AddNewColumn/functions-sorting-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All five scenario steps PASS in both the MCP run and the Playwright replay (17s, 1 test, 0 failures). The previous ru… | 8m 20s | 2m | 17s | 10m 37s | ⚪ +0% (7/7 → 7/7) | 🟡 new (prev7d had no entry) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps PASS in both the MCP run and the Playwright replay. `${AGE}`, `$[AGE]` and the autocomplete-inserted `${H… | 4m 43s | 1m 13s | 25s | 6m 21s | ⚪ +0% (5/5 → 5/5) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | 🟢 PASS → PASS | 🟢 100% (4/4) | All four steps pass in both MCP and Playwright. The CodeMirror formula editor shows a `.cm-tooltip-hover` on hover wi… | 1m 10s | 15s | 8s | 1m 33s | ⚪ +0% (4/4 → 4/4) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [AddNewColumn/input-functions](PowerPack/AddNewColumn/input-functions-run.md) | 🟢 new: PASS | 🟢 100% (10/10) | All 10 scenario steps pass end-to-end — both in interactive MCP driving and in the Playwright replay (18.8s). A singl… | 3m 30s | 1m | 18.8s | 4m 49s | 🟢 100% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 100 | 🟡 3m 30s (new) | 🟡 1m (new) | 🟡 18.8s (new) | 🟡 4m 49s (new) |
| PowerPack | AddNewColumn/input_functions | 🟢 removed: PASS |  |  |  |  |  |  | — | — | — | — | — | — | — |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All six scenario steps PASS in both the MCP-driven grok-browser run and the Playwright replay (existing spec — not ov… | 1m 54s | 31s | 11s | 2m 36s | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 76% (16/21) | The PowerPack "Enrich column" feature works end-to-end for the primary create/apply/edit/delete flow on a dataframe t… | 16m | 2m | 2m | 20m | ⚪ +0% (16/21 → 16/21) | 🟡 76% · status unchanged (PARTIAL) | 🟢 0·76·76 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Projects | [browser](Projects/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 56% (5/9) | 5 of 9 steps passed, 3 skipped, 1 ambiguous. Browse > Dashboards view works correctly: projects are listed, searchabl… |  |  |  |  | ⚪ +0% (5/9 → 5/9) | 🟡 56% · status unchanged (PARTIAL) | ⚪ 56·56·56 | — | — | — | — |
| Projects | [complex](Projects/complex-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/13) | All 13 steps skipped. This is the most complex scenario requiring tables from 7+ different sources, drag-and-drop, en… |  |  |  |  | ⚪ +0% (0/13 → 0/13), still broken | 🔴 0% · status unchanged (SKIP) | ⚪ 0·0·0 | — | — | — | — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/5) | All 5 steps skipped. This scenario requires running a custom JavaScript script with Data Sync enabled, then modifying… |  |  |  |  | ⚪ +0% (0/5 → 0/5), still broken | 🔴 0% · status unchanged (SKIP) | ⚪ 0·0·0 | — | — | — | — |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 50% (2/4) | 2 of 4 steps passed, 1 skipped, 1 ambiguous. Project deletion works via the API (`grok.dapi.projects.delete()`). The … |  |  |  |  | ⚪ +0% (2/4 → 2/4) | 🟡 50% · status unchanged (PARTIAL) | ⚪ 50·50·50 | — | — | — | — |
| Projects | [opening](Projects/opening-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (5/5) | All 5 steps passed. Projects from the Uploading step are accessible in Browse > Dashboards. Context Panel correctly s… |  |  |  |  | ⚪ +0% (5/5 → 5/5) | 🟢 100% · status unchanged (PARTIAL) | ⚪ 100·100·100 | — | — | — | — |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 SKIP → SKIP | 🔴 0% (0/4) | All steps skipped. This scenario depends on Projects copy_clone.md (order 5) which was not fully executed. The Link/C… |  |  |  |  | ⚪ +0% (0/4 → 0/4), still broken | 🔴 0% · status unchanged (SKIP) | ⚪ 0·0·0 | — | — | — | — |
| Projects | [projects-copy-clone](Projects/projects-copy-clone-run.md) | 🟡 new: SKIP | 🟡 40% (2/5) | 2 of 5 steps passed, 3 skipped. Project preview and opening work. Copy/clone/link operations were not tested because … |  |  |  |  | 🟡 40% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 40 | — | — | — | — |
| Projects | projects-copy_clone | 🟡 removed: SKIP |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 57% (8/14) | 8 of 14 steps passed, 6 skipped. Core project creation from local tables, file shares, query results, and join result… |  |  |  |  | ⚪ +0% (8/14 → 8/14) | 🟡 57% · status unchanged (PARTIAL) | ⚪ 57·57·57 | — | — | — | — |
| Queries | adding | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | browse-&-save-project | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | browse-and-save-project | ⚪ new: NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | browser | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | columns-inspect | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | deleting | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | edit | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | get-all-get-top-100 | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | ms-sql | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | new-sql-query | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | new-visual-query | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | query-layout | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | query-postprocessing | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | transformations | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Queries | visual-query-advanced | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Scripts | [browser](Scripts/browser-run.md) | 🟢 PASS → PASS | 🟡 78% (7/9) | Browser scenario passed in the MCP session — all accordions (Details, Script, Run, Activity, Sharing, Chats, Dev) ren… | 1m 17s | 25s | 16s | 1m 58s | 🔴 -11% (8/9 → 7/9) | 🟡 78% · status unchanged (PASS) | 🔴 89·89·78 | 🟡 1m 17s (new) | 🟡 25s (new) | 🟡 16s (new) | 🟡 1m 58s (new) |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 92% (10/12) | The Create scenario completed successfully overall. The script `testRscript` was created, parameters configured, save… |  |  |  |  | ⚪ +0% (10/12 → 10/12) | 🟡 92% · status unchanged (PARTIAL) | ⚪ 92·92·92 | — | — | — | — |
| Scripts | [delete](Scripts/delete-run.md) | 🟢 PASS → PASS | 🟢 100% (5/5) | All 5 steps passed. The delete flow works correctly with a confirmation dialog and immediate removal from the scripts… |  |  |  |  | ⚪ +0% (5/5 → 5/5) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | — | — | — | — |
| Scripts | [edit](Scripts/edit-run.md) | 🟢 PASS → PASS | 🟢 100% (6/6) | All 6 steps passed. The Edit scenario works correctly — edits are saved persistently and visible on re-open. |  |  |  |  | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | — | — | — | — |
| Scripts | layout | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (6/9) | Core run functionality works: the script can be triggered from context menu with a table selection and from the conso… |  |  |  |  | ⚪ +0% (6/9 → 6/9) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | — | — | — | — |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 3 scenario sections (single-cell add/edit, sticky column behavior, batch edit) were reproduced successfully via M… | 5m 11s | 1m | 44s | 6m 55s | ⚪ +0% (5/5 → 9/9) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +4m 36s | 🔴 +57s | 🔴 +27s | 🔴 +6m |
| StickyMeta | [copy-clone-delete](StickyMeta/copy-clone-delete-run.md) | 🟢 PASS → PASS | 🟢 100% (10/10) | All 10 steps PASS end-to-end, both in the MCP scenario run and in the Playwright replay (spec wall-clock 1m 29s). Ful… | 7m 52s | 2m | 1m 29s | 11m 21s | 🟢 +50% (2/4 → 10/10) | 🟡 new (prev7d had no entry) | 🟢 50·100 | 🔴 +7m 27s | 🔴 +1m 57s | 🔴 +1m 8s | 🔴 +10m 32s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | 🟢 PASS → PASS | 🟢 100% (8/8) | All 8 scenario steps PASS in both the MCP browser reproduction and the generated Playwright spec. TestEntity1 is crea… | 2m 48s | 1m 24s | 59s | 5m 11s | ⚪ +0% (4/4 → 8/8) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +2m 38s | 🔴 +1m 21s | 🔴 +52s | 🔴 +4m 51s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | 🟡 9% (1/11) | The "Database meta" context-panel section is not rendered on `dev.datagrok.ai` for a Postgres DbInfo connection entit… | 7m 10s | 1m 20s | 33s | 9m 3s | 🔴 -11% (1/5 → 1/11) | 🟡 9% · status unchanged (FAIL) | 🔴 20·20·9 | 🔴 +6m 45s | 🔴 +1m 20s | 🟡 33s (new) | 🔴 +8m 38s |
| Tooltips | actions-in-the-context-menu | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Tooltips | default-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Tooltips | default-tooltip-visibility | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Tooltips | edit-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Tooltips | line-chart---aggregated-tooltip | ⚪ removed: NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Tooltips | line-chart-aggregated-tooltip | ⚪ new: NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Tooltips | tooltip-properties | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Tooltips | uniform-default-tooltip | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 14 steps (setup + 13 scenario sections) passed in both the MCP browser run against `https://dev.datagrok.ai` and … | 40s | 10s | 24s | 1m 14s | ⚪ +0% (15/15 → 15/15) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟢 -6m 20s | 🟢 -1m 50s | 🔴 +4s | 🟢 -8m 6s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran basic-operations end-to-end against dev. All 31 scenario steps passed in the MCP browser phase (Section 1: struct… | 4m 27s | 9s | 50s | 5m 26s | ⚪ +0% (26/26 → 13/13) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | Ran chem-and-bio scenario end-to-end against dev. All 11 scenario steps passed in the MCP browser phase (Chem: open s… | 2m 50s | 42s | 47s | 4m 19s | ⚪ +0% (16/16 → 11/11) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | 🟢 PASS → PASS | 🟢 100% (15/15) | All 15 scenario steps PASSed on dev. spgi-100.csv loads correctly this time (previous run had to substitute SPGI.csv)… | 3m 16s | 14s | 55s | 4m 25s | ⚪ +0% (15/15 → 15/15) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed end-to-end on dev: table linking (SELECTION_TO_FILTER and FILTER_TO_FILTER) propagated correctly b… | 1m 46s | 17s | 35s | 2m 38s | ⚪ +0% (9/9 → 9/9) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | Ran combined-boolean-filter end-to-end against dev. All 13 numbered scenario steps passed in the MCP browser phase: S… | 2m 37s | 12s | 24s | 3m 13s | ⚪ +0% (13/13 → 13/13) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (14/14) | All 14 steps passed in both the MCP run and the Playwright replay. Expression filter works correctly: 5-rule AND yiel… | 1m 14s | 8s | 23s | 1m 45s | ⚪ +0% (14/14 → 14/14) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 12 steps passed in the MCP run and in the Playwright replay (spec finished in 21.8s). The hierarchical filter cor… | 1m 15s | 21s | 23s | 1m 59s | ⚪ +0% (12/12 → 12/12) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 steps passed in the MCP run and in the Playwright replay (spec finished in 8.7s, total wall-clock 11.56s). The … | 1m 12s | 20s | 12s | 1m 44s | ⚪ +0% (9/9 → 9/9) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | 🟢 PASS → PASS | 🟢 100% (34/34) | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to re… | 4m 24s | 40s | 1m 3s | 6m 7s | ⚪ +0% (34/34 → 34/34) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 69% (5/8) | Color consistency through layout round-trip works — the `.categorical-colors` tag survives save/reload and `R_ONE` st… | 2m 30s | 35s | 26s | 3m 31s | ⚪ +0% (5/8 → 5/8) | 🟡 new (prev7d had no entry) | ⚪ 69·69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 83% (9/12) | Filtering legend updates work end-to-end in the MCP run: numeric filter, categorical filter, layout round-trip, compo… | 3m 10s | 1m 10s | 44s | 5m 4s | ⚪ +0% (9/12 → 9/12) | 🟡 new (prev7d had no entry) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 77% (8/11) | Line chart legend and multi-axis behaviors are mostly correct: 7 legend items for 7 categories, layout round-trip pre… | 2m 10s | 40s | 32s | 3m 22s | ⚪ +0% (8/11 → 8/11) | 🟡 new (prev7d had no entry) | ⚪ 77·77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 58% (7/13) | Categorical legend on scatter plot updates correctly when X axis changes (sub 2) and when the Filter Panel narrows ca… | 4m 15s | 1m 20s | 54s | 6m 29s | ⚪ +0% (7/13 → 7/13) | 🟡 new (prev7d had no entry) | ⚪ 58·58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 79% (5/7) | Structure rendering in legends works for Scatter plot, Histogram, Line chart and Pie chart (canvas-based molecule thu… | 2m 35s | 40s | 29s | 3m 44s | ⚪ +0% (5/7 → 5/7) | 🟡 new (prev7d had no entry) | ⚪ 79·79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 70% (13/20) | Scenario executed end-to-end with a mix of PASS, AMBIGUOUS, and FAIL. Legend display, source-swap, corner positioning… | 5m 45s | 1m 30s | 41s | 7m 56s | ⚪ +0% (13/20 → 13/20) | 🟡 new (prev7d had no entry) | ⚪ 70·70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | 🟢 PASS → PASS | 🟡 85% (11/13) |  |  |  |  |  | 🔴 -1% (12/14 → 11/13) | 🟡 85% · status unchanged (PASS) | 🟢 77·86·85 | — | — | — | — |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (82/82) | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color cod… | 3m 3s | 21s | 52s | 4m 16s | ⚪ +0% (82/82 → 82/82) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | bar-chart-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 95% (18/19) | 17 of 19 sections passed cleanly; section 8 combined into section 7 in the spec. Section 18 is AMBIGUOUS — `grok.dapi… | 1m 5s | 8s | 32s | 1m 45s | ⚪ +0% (18/19 → 18/19) | 🟢 95% · status unchanged (PARTIAL) | ⚪ 95·95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS | 🟢 100% (11/11) | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selec… | 25s |  |  | 25s | ⚪ +0% (11/11 → 11/11) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | 🟡 1m (removed) | 🟡 9.4s (removed) | 🟢 -1m 9s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | 🟢 PASS → PASS | 🟢 100% (12/12) | All 11 steps passed. The entire test runs on the demog dataset (no SPGI_v2 needed). UI-only steps (Grid Color Coding … |  |  |  |  | ⚪ +0% (12/12 → 12/12) | 🟢 100% · status unchanged (PASS) | 🟢 67·100·100 | — | — | 🟡 29s (removed) | 🟡 29s (removed) |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 87% (26/30) | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operat… |  | 3s | 22.5s | 25.5s | ⚪ +0% (26/30 → 26/30) | 🟡 87% · status unchanged (PARTIAL) | ⚪ 87·87·87 | 🟡 5m (removed) | ⚪ +0s | ⚪ +0s | 🟢 -5m |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (58/58) | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI intera… |  |  |  |  | ⚪ +0% (58/58 → 58/58) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | — | 🟡 2m (removed) | — | 🟡 2m (removed) |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | 🟡 92% (24/26) | All 14 sections of form-tests-pw.md exercised across 30 steps. 28 PASS, 2 AMBIGUOUS, 0 FAIL in MCP run. Playwright sp… |  |  | 3m 12s | 3m 12s | 🔴 -1% (28/30 → 24/26) | 🟡 92% · status unchanged (PASS) | 🔴 93·93·92 | 🟡 18m (removed) | 🟡 4m (removed) | ⚪ +0s | 🟢 -22m |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | 🟢 100% (30/30) | All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas ele… |  |  | 51.5s | 51.5s | ⚪ +0% (36/36 → 30/30) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | 🟡 18m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -21m |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 73% (16/22) | Grid tests ran 22 steps (spec softSteps); 17 passed outright and 5 were AMBIGUOUS (Copy/Paste, Column Header Context … | 11m | 3m | 1m 18s | 15m 18s | ⚪ +0% (16/22 → 16/22) | 🟡 73% · status unchanged (PARTIAL) | ⚪ 73·73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | 🟡 93% (13/14) | All 14 heat-map sections exercised across 17 steps. 15 PASS, 1 AMBIGUOUS, 1 SKIP in MCP run. Playwright spec passed f… |  |  | 48.9s | 48.9s | 🔴 -1% (15/16 → 13/14) | 🟡 93% · status PARTIAL → PASS | 🔴 94·94·93 | 🟡 18m (removed) | 🟡 4m (removed) | ⚪ +0s | 🟢 -22m |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 93% (87/94) | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance… | 50s | 7s | 46s | 1m 43s | ⚪ +0% (87/94 → 87/94) | 🟡 93% · status unchanged (PARTIAL) | ⚪ 93·93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (26/26) | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout… | 57s | 8s | 1m 43s | 2m 48s | ⚪ +0% (26/26 → 26/26) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [map](Viewers/map-run.md) | 🟢 PASS → PASS | 🟡 80% (8/10) | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker siz… | 15s | 3s | 9s | 27s | ⚪ +0% (8/10 → 8/10) | 🟡 80% · status unchanged (PASS) | ⚪ 80·80·80 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (16/19) | Matrix Plot tests ran with 15 PASS, 3 AMBIGUOUS, 0 FAIL. The spec executed in 57.7s with all implemented steps passin… |  |  | 55.9s | 55.9s | ⚪ +0% (16/19 → 16/19) | 🟡 84% · status PASS → PARTIAL | ⚪ 84·84·84 | 🟡 20m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -23m |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 67% (8/12) | 9 of 12 steps PASS; 3 SKIP (canvas-based node/edge interactions cannot be automated via DOM). The network diagram vie… | 8m | 1m 30s | 22s | 9m 52s | ⚪ +0% (8/12 → 8/12) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (13/13) | All 13 scenario sections (mapped to 12 Playwright softSteps — scale and normalization are combined in the spec) passe… | 1m 8s | 8s | 47s | 2m 3s | ⚪ +0% (13/13 → 13/13) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | 🟢 PASS → PASS | 🟢 100% (81/81) | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appea… | 40s | 7s | 47s | 1m 34s | ⚪ +0% (81/81 → 81/81) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (17/20) | Pivot Table tests ran with 16 PASS, 2 AMBIGUOUS, 1 SKIP, 0 FAIL. The spec executed in 35.1s with all implemented step… |  |  | 1m 12s | 1m 12s | ⚪ +0% (17/20 → 17/20) | 🟡 85% · status PASS → PARTIAL | ⚪ 85·85·85 | 🟡 16m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -19m |
| Viewers | [row-source](Viewers/row-source-run.md) | 🟢 PASS → PASS | 🟢 100% (36/36) | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with al… |  | 5s | 1m 24s | 1m 29s | ⚪ +0% (36/36 → 36/36) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 4m (removed) | ⚪ +0s | ⚪ +0s | 🟢 -4m |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | 🟢 PASS → PASS | 🟢 100% (20/20) | All 20 sections passed during the MCP run on dev.datagrok.ai. The existing Playwright spec was re-run headed without … | 3m 13s | 29s | 52s | 4m 34s | ⚪ +0% (20/20 → 20/20) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | scatter-plot-tests | ⚪ NO RUN → NO RUN |  |  |  |  |  |  | — | — | — | — | — | — | — |
| Viewers | [statistics](Viewers/statistics-run.md) | 🟢 PASS → PASS | 🟢 100% (24/24) | All 23 MCP steps passed. The date columns section (STARTED row behavior) was moved to `statistics-tests-ui.md` as a m… | 20m | 4m |  | 24m | ⚪ +0% (24/24 → 24/24) | 🟡 new (prev7d had no entry) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟡 2m (removed) | 🟢 -2m |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | 🟢 PASS → PASS | 🟢 100% (24/24) | 24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to man… |  | 3m | 58s | 3m 58s | ⚪ +0% (24/24 → 24/24) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | — | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | 🟡 PARTIAL → PARTIAL | 🟢 100% (36/36) | All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via… | 28m | 4m | 46s | 32m 46s | ⚪ +0% (36/36 → 36/36) | 🟢 100% · status unchanged (PARTIAL) | ⚪ 100·100·100 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 84% (70/84) | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range s… | 3m | 30s |  | 3m 30s | ⚪ +0% (70/84 → 70/84) | 🟡 84% · status unchanged (PARTIAL) | ⚪ 84·84·84 | ⚪ +0s | ⚪ +0s | 🟡 1m 48s (removed) | 🟢 -1m 48s |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | 🟢 PASS → PASS | 🟢 100% (7/7) | All 7 scenario steps PASS in both the MCP run against https://dev.datagrok.ai and the generated Playwright replay (27… | 2m 17s | 50s | 27s | 3m 34s | ⚪ +0% (7/7 → 7/7) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | 🟢 -1m 58s | 🟢 -10s | 🟢 -1m 40s | 🟢 -3m 48s |
| Viewers | [word-cloud-tests](Viewers/word-cloud-tests-run.md) | 🟡 PARTIAL → PARTIAL | 🟡 85% (39/46) | All 46 scenario steps were exercised against dev.datagrok.ai. In the MCP run, 35 passed, 3 were AMBIGUOUS (visual eff… | 4m 30s | 2m | 34s | 7m 4s | ⚪ +0% (23/27 → 39/46) | 🟡 new (prev7d had no entry) | ⚪ 85·85 | 🟢 -2m | 🔴 +45s | 🔴 +17s | 🟢 -58s |
| Viewers | [working-with-nan-infinity](Viewers/working-with-nan-infinity-run.md) | 🟢 PASS → PASS | 🟢 100% (9/9) | All 9 spec steps PASSED in 1m 24s. NaN and Infinity values in numeric columns are handled gracefully across Scatter P… |  |  | 1m 24s | 1m 24s | ⚪ +0% (9/9 → 9/9) | 🟡 new (prev7d had no entry) | ⚪ 100·100 | 🟡 6m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -9m |

## Comparison with Previous Reports

Deltas are computed against two baselines pulled from git history of `total-run.md`:
`prev1d` (latest commit before today, 2026-04-23) and `prev7d` (commit closest to today − 7 days, 2026-04-15; old format without Pass %).
Signed with `+`/`-`; time deltas use the same format as the values. All
status and delta cells carry the Legend icons.

### Totals

**Total (1d)**: Tests Δ **+0** · Mean Pass Δ **⚪ +0%** · Browser Δ **🟢 -25.7s** · Spec Gen Δ **🟢 -3.4s** · Spec Run Δ **🔴 +4.4s** · Total Δ **🟢 -35.6s**

_7d comparison: prev7d (2026-04-15) predates the Pass % / timing columns, so only status deltas are meaningful for that baseline — reflected in the per-test table below._

### By Folder

| Folder | Tests Δ (1d) | Status | Mean Pass Δ (1d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|
| Apps | +0 | ⚪ NO DATA | — | — | — | — | — |
| Bio | +0 | 🟡 PARTIAL | 🟢 +6% | 🟢 -1.7s | 🔴 +38.2s | 🟡 43.4s (new) | 🔴 +1m 14s |
| Browse | +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — |
| Charts | +0 | 🟡 PARTIAL | ⚪ +0% | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | +0 | 🟡 PARTIAL | ⚪ +0% | ⚪ +0s | ⚪ +0s | 🔴 +0.2s | 🟢 -2.7s |
| Connections | +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — |
| DiffStudio | +0 | 🟡 PARTIAL | ⚪ +0% | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | +0 | 🟡 PARTIAL | 🟢 +1% | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| General | +0 | ⚪ NO DATA | — | — | — | — | — |
| LocalCashing | +0 | ⚪ NO DATA | — | — | — | — | — |
| Models | +0 | 🟡 PARTIAL | 🔴 -26% | 🟡 8m 20s (new) | 🟡 2m 3s (new) | 🟡 1m 45s (new) | 🟡 12m 8s (new) |
| Notebooks | +0 | ⚪ NO DATA | — | — | — | — | — |
| Peptides | +0 | 🟡 PARTIAL | ⚪ +0% | ⚪ +0s | ⚪ +0s | 🔴 +36.4s | 🔴 +42s |
| PowerPack | +0 | 🟡 PARTIAL | ⚪ +0% | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Projects | +0 | 🟡 PARTIAL | ⚪ +0% | — | — | — | — |
| Queries | +0 | ⚪ NO DATA | — | — | — | — | — |
| Scripts | +0 | 🟡 PARTIAL | 🔴 -2% | 🟡 1m 17s (new) | 🟡 25s (new) | 🟡 16s (new) | 🟡 1m 58s (new) |
| StickyMeta | +0 | 🟡 PARTIAL | 🟢 +10% | 🔴 +5m 22s | 🔴 +1m 24s | 🔴 +41.2s | 🔴 +7m 30s |
| Tooltips | +0 | ⚪ NO DATA | — | — | — | — | — |
| Viewers | +0 | 🟡 PARTIAL | ⚪ +0% | 🟢 -2m 3s | 🟢 -24.7s | 🟢 -3.6s | 🟢 -3m 10s |
| WideSmokeTest | +0 | ⚪ NO DATA | — | — | — | — | — |

### Per-Test Changes

| Folder | Test | Status | Pass Δ (1d) | Pass Δ (7d) | Trend (≤7d) | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|---|---|---|---|---|---|---|---|---|---|
| Bio | [analyze](Bio/analyze-run.md) | 🟢 PASS → PASS | ⚪ +0% (11/11 → 11/11) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟢 -5m 20s | 🔴 +15s | 🟡 1m 21s (new) | 🟢 -3m 44s |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | 🟡 → 🟢 PARTIAL → PASS | 🟢 +20% (4/5 → 5/5) | 🟢 100% · status PARTIAL → PASS | 🟢 80·80·100 | 🟡 4m 13s (new) | 🟡 1m 29s (new) | 🟡 41s (new) | 🟡 6m 23s (new) |
| Bio | [convert](Bio/convert-run.md) | 🟢 PASS → PASS | ⚪ +0% (5/5 → 15/15) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +2m 40s | 🟡 1m 10s (new) | 🟡 1m 42s (new) | 🔴 +5m 32s |
| Bio | [manage](Bio/manage-run.md) | 🟢 PASS → PASS | ⚪ +0% (3/3 → 3/3) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +6s | 🟡 25s (new) | 🟡 22s (new) | 🔴 +53s |
| Bio | [msa](Bio/msa-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 2m 50s (new) | 🟡 40s (new) | 🟡 19s (new) | 🟡 3m 49s (new) |
| Bio | [pepsea](Bio/pepsea-run.md) | 🟡 → 🟡 AMBIGUOUS → PARTIAL | 🟢 +33% (4/6 → 5/5) | 🟢 100% · status AMBIGUOUS → PARTIAL | 🟢 67·67·100 | 🟡 3m 15s (new) | 🟡 30s (new) | 🟡 25s (new) | 🟡 4m 10s (new) |
| Bio | [search](Bio/search-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 4/4) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +30s | 🟡 25s (new) | 🟡 14s (new) | 🔴 +1m 9s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 2m 26s (new) | 🟡 1m (new) | — | 🟡 3m 26s (new) |
| Bio | [sequence-space](Bio/sequence-space-run.md) | 🟢 PASS → PASS | ⚪ +0% (6/6 → 6/6) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 2m 50s (new) | 🟡 30s (new) | — | 🟡 3m 20s (new) |
| Browse | [browse](Browse/browse-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/5 → 4/5) | 🟡 90% · status unchanged (PARTIAL) | ⚪ 90·90·90 | — | — | — | — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (0/1 → 0/1) | 🟡 50% · status unchanged (PARTIAL) | ⚪ 50·50·50 | — | — | — | — |
| Charts | [sunburst](Charts/sunburst-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/12 → 5/12) | 🟡 42% · status unchanged (PARTIAL) | ⚪ 42·42·42 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Charts | [tree](Charts/tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (1/5 → 1/5) | 🟡 20% · status PASS → PARTIAL | ⚪ 20·20·20 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/6 → 3/6) | 🟡 new (prev7d had no entry) | ⚪ 50·50·50 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Chem | [calculate](Chem/calculate-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/3 → 1/3) | 🟡 33% · status PASS → FAIL | ⚪ 33·33·33 | ⚪ +0s | ⚪ +0s | 🟡 38s (removed) | 🟢 -38s |
| Chem | [chemprop](Chem/chemprop-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/5 → 2/5) | 🟡 40% · status unchanged (PARTIAL) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Connections | [adding](Connections/adding-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/7 → 6/7) | 🟡 93% · status unchanged (PARTIAL) | ⚪ 93·93·93 | — | — | — | — |
| Connections | [browser](Connections/browser-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (7/9 → 7/9) | 🟡 78% · status unchanged (PARTIAL) | ⚪ 78·78·78 | — | — | — | — |
| Connections | [catalogs](Connections/catalogs-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/11 → 1/11) | 🟡 9% · status unchanged (FAIL) | ⚪ 9·9·9 | — | — | — | — |
| Connections | [external-provider](Connections/external-provider-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (0/7 → 0/7), still broken | 🔴 0% · status unchanged (FAIL) | ⚪ 0·0·0 | — | — | — | — |
| Connections | [identifiers](Connections/identifiers-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (1/9 → 1/9) | 🟡 11% · status unchanged (FAIL) | ⚪ 11·11·11 | — | — | — | — |
| Connections | [import-swagger](Connections/import-swagger-run.md) | 🔴 FAIL → FAIL | ⚪ +0% (0/7 → 0/7), still broken | 🔴 0% · status unchanged (FAIL) | ⚪ 0·0·0 | — | — | — | — |
| Connections | [schema](Connections/schema-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/4 → 3/4) | 🟡 75% · status unchanged (PARTIAL) | ⚪ 75·75·75 | — | — | — | — |
| Connections | [sparql](Connections/sparql-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/7 → 6/7) | 🟡 86% · status unchanged (PARTIAL) | ⚪ 86·86·86 | — | — | — | — |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/6 → 5/6) | 🟡 83% · status PASS → PARTIAL | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | ML methods/linear-regression | 🟢 removed: PASS | — | — | — | — | — | — | — |
| EDA | ML methods/pls-regression | 🟡 removed: PARTIAL | — | — | — | — | — | — | — |
| EDA | ML methods/softmax | 🔴 removed: FAIL | — | — | — | — | — | — | — |
| EDA | ML methods/xgboost1 | 🟡 removed: PARTIAL | — | — | — | — | — | — | — |
| EDA | ML methods/xgboost2 | 🟡 removed: PARTIAL | — | — | — | — | — | — | — |
| EDA | [MLMethods/linear-regression](EDA/MLMethods/linear-regression-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 100 | 🟡 45s (new) | 🟡 2s (new) | 🟡 6.8s (new) | 🟡 53.8s (new) |
| EDA | [MLMethods/pls-regression](EDA/MLMethods/pls-regression-run.md) | 🟡 new: PARTIAL | 🟡 71% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 71 | 🟡 1m (new) | 🟡 2s (new) | 🟡 6.9s (new) | 🟡 1m 9s (new) |
| EDA | [MLMethods/softmax](EDA/MLMethods/softmax-run.md) | 🔴 new: FAIL | 🟡 33% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 33 | 🟡 10s (new) | 🟡 2s (new) | 🟡 2.6s (new) | 🟡 14.6s (new) |
| EDA | [MLMethods/xgboost1](EDA/MLMethods/xgboost1-run.md) | 🟡 new: PARTIAL | 🟡 67% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 67 | 🟡 5s (new) | 🟡 2s (new) | 🟡 2.7s (new) | 🟡 9.7s (new) |
| EDA | [MLMethods/xgboost2](EDA/MLMethods/xgboost2-run.md) | 🟡 new: PARTIAL | 🟡 67% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 67 | 🟡 5s (new) | 🟡 2s (new) | 🟡 2.7s (new) | 🟡 9.7s (new) |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/3 → 2/3) | 🟡 67% · status PASS → PARTIAL | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/7 → 3/7) | 🟡 43% · status unchanged (PARTIAL) | ⚪ 43·43·43 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pca](EDA/pca-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (3/5 → 3/5) | 🟡 60% · status PASS → PARTIAL | ⚪ 60·60·60 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| EDA | [pls](EDA/pls-run.md) | 🔴 FAIL → FAIL | 🟢 +12% (2/5 → 2/4) | 🟡 62% · status PASS → FAIL | ⚪ 62·50·62 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Models | [apply](Models/apply-run.md) | 🟢 → 🔴 PASS → FAIL | 🔴 -50% (4/4 → 2/4) | 🟡 50% · status PASS → FAIL | 🔴 100·100·50 | 🟡 2m 12s (new) | 🟡 30s (new) | 🟡 20s (new) | 🟡 3m 2s (new) |
| Models | [browser](Models/browser-run.md) | 🟡 PARTIAL → PARTIAL | 🔴 -34% (4/6 → 2/6) | 🟡 33% · status unchanged (PARTIAL) | 🔴 67·67·33 | 🟡 4m 45s (new) | 🟡 2m 10s (new) | 🟡 51s (new) | 🟡 7m 46s (new) |
| Models | [chemprop](Models/chemprop-run.md) | 🔴 FAIL → FAIL | 🟢 +7% (5/18 → 6/17) | 🟡 35% · status unchanged (FAIL) | 🟢 28·28·35 | 🟡 13m (new) | 🟡 2m (new) | 🟡 3m 45s (new) | 🟡 18m 45s (new) |
| Models | [delete](Models/delete-run.md) | 🟢 → 🔴 PASS → FAIL | 🔴 -80% (5/5 → 1/5) | 🟡 20% · status PASS → FAIL | 🔴 100·100·20 | 🟡 18m (new) | 🟡 4m (new) | 🟡 24s (new) | 🟡 22m 24s (new) |
| Models | [predictive-models](Models/predictive-models-run.md) | 🟢 → 🟡 PASS → PARTIAL | ⚪ +0% (21/21 → 17/17) | 🟢 100% · status PASS → PARTIAL | ⚪ 100·100·100 | 🟡 5m 50s (new) | 🟡 1m 50s (new) | 🟡 4m 40s (new) | 🟡 12m 20s (new) |
| Models | [train](Models/train-run.md) | 🟢 PASS → PASS | 🟢 +1% (9/10 → 10/11) | 🟡 91% · status unchanged (PASS) | 🟢 90·90·91 | 🟡 6m 10s (new) | 🟡 1m 50s (new) | 🟡 28s (new) | 🟡 8m 28s (new) |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/5 → 2/5) | 🟡 40% · status unchanged (PARTIAL) | ⚪ 40·40·40 | ⚪ +0s | ⚪ +0s | 🟡 1m 18s (new) | 🔴 +1m 18s |
| Peptides | [peptides](Peptides/peptides-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (4/6 → 4/6) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Peptides | [sar](Peptides/sar-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/9 → 6/9) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | 🟡 1m 30s (new) | 🔴 +1m 30s |
| PowerPack | [AddNewColumn/input-functions](PowerPack/AddNewColumn/input-functions-run.md) | 🟢 new: PASS | 🟢 100% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 100 | 🟡 3m 30s (new) | 🟡 1m (new) | 🟡 18.8s (new) | 🟡 4m 49s (new) |
| PowerPack | AddNewColumn/input_functions | 🟢 removed: PASS | — | — | — | — | — | — | — |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (16/21 → 16/21) | 🟡 76% · status unchanged (PARTIAL) | 🟢 0·76·76 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Projects | [browser](Projects/browser-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/9 → 5/9) | 🟡 56% · status unchanged (PARTIAL) | ⚪ 56·56·56 | — | — | — | — |
| Projects | [complex](Projects/complex-run.md) | 🟡 SKIP → SKIP | ⚪ +0% (0/13 → 0/13), still broken | 🔴 0% · status unchanged (SKIP) | ⚪ 0·0·0 | — | — | — | — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | 🟡 SKIP → SKIP | ⚪ +0% (0/5 → 0/5), still broken | 🔴 0% · status unchanged (SKIP) | ⚪ 0·0·0 | — | — | — | — |
| Projects | [deleting](Projects/deleting-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (2/4 → 2/4) | 🟡 50% · status unchanged (PARTIAL) | ⚪ 50·50·50 | — | — | — | — |
| Projects | [project-url](Projects/project-url-run.md) | 🟡 SKIP → SKIP | ⚪ +0% (0/4 → 0/4), still broken | 🔴 0% · status unchanged (SKIP) | ⚪ 0·0·0 | — | — | — | — |
| Projects | [projects-copy-clone](Projects/projects-copy-clone-run.md) | 🟡 new: SKIP | 🟡 40% (baseline) | 🟡 new (prev7d had no entry) | ⚪ 40 | — | — | — | — |
| Projects | projects-copy_clone | 🟡 removed: SKIP | — | — | — | — | — | — | — |
| Projects | [uploading](Projects/uploading-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/14 → 8/14) | 🟡 57% · status unchanged (PARTIAL) | ⚪ 57·57·57 | — | — | — | — |
| Queries | browse-&-save-project | ⚪ removed: NO RUN | — | — | — | — | — | — | — |
| Scripts | [browser](Scripts/browser-run.md) | 🟢 PASS → PASS | 🔴 -11% (8/9 → 7/9) | 🟡 78% · status unchanged (PASS) | 🔴 89·89·78 | 🟡 1m 17s (new) | 🟡 25s (new) | 🟡 16s (new) | 🟡 1m 58s (new) |
| Scripts | [create](Scripts/create-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (10/12 → 10/12) | 🟡 92% · status unchanged (PARTIAL) | ⚪ 92·92·92 | — | — | — | — |
| Scripts | [run](Scripts/run-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (6/9 → 6/9) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | — | — | — | — |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | 🟢 PASS → PASS | ⚪ +0% (5/5 → 9/9) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +4m 36s | 🔴 +57s | 🔴 +27s | 🔴 +6m |
| StickyMeta | [copy-clone-delete](StickyMeta/copy-clone-delete-run.md) | 🟢 PASS → PASS | 🟢 +50% (2/4 → 10/10) | 🟡 new (prev7d had no entry) | 🟢 50·100 | 🔴 +7m 27s | 🔴 +1m 57s | 🔴 +1m 8s | 🔴 +10m 32s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | 🟢 PASS → PASS | ⚪ +0% (4/4 → 8/8) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🔴 +2m 38s | 🔴 +1m 21s | 🔴 +52s | 🔴 +4m 51s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | 🔴 FAIL → FAIL | 🔴 -11% (1/5 → 1/11) | 🟡 9% · status unchanged (FAIL) | 🔴 20·20·9 | 🔴 +6m 45s | 🔴 +1m 20s | 🟡 33s (new) | 🔴 +8m 38s |
| Tooltips | line-chart---aggregated-tooltip | ⚪ removed: NO RUN | — | — | — | — | — | — | — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | 🟢 PASS → PASS | ⚪ +0% (15/15 → 15/15) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟢 -6m 20s | 🟢 -1m 50s | 🔴 +4s | 🟢 -8m 6s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/8 → 5/8) | 🟡 new (prev7d had no entry) | ⚪ 69·69·69 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (9/12 → 9/12) | 🟡 new (prev7d had no entry) | ⚪ 83·83·83 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/11 → 8/11) | 🟡 new (prev7d had no entry) | ⚪ 77·77·77 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (7/13 → 7/13) | 🟡 new (prev7d had no entry) | ⚪ 58·58·58 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (5/7 → 5/7) | 🟡 new (prev7d had no entry) | ⚪ 79·79·79 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (13/20 → 13/20) | 🟡 new (prev7d had no entry) | ⚪ 70·70·70 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | 🟢 PASS → PASS | 🔴 -1% (12/14 → 11/13) | 🟡 85% · status unchanged (PASS) | 🟢 77·86·85 | — | — | — | — |
| Viewers | [box-plot](Viewers/box-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (18/19 → 18/19) | 🟢 95% · status unchanged (PARTIAL) | ⚪ 95·95·95 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | 🟢 PASS → PASS | ⚪ +0% (11/11 → 11/11) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | ⚪ +0s | 🟡 1m (removed) | 🟡 9.4s (removed) | 🟢 -1m 9s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | 🟢 PASS → PASS | ⚪ +0% (12/12 → 12/12) | 🟢 100% · status unchanged (PASS) | 🟢 67·100·100 | — | — | 🟡 29s (removed) | 🟡 29s (removed) |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (26/30 → 26/30) | 🟡 87% · status unchanged (PARTIAL) | ⚪ 87·87·87 | 🟡 5m (removed) | ⚪ +0s | ⚪ +0s | 🟢 -5m |
| Viewers | [density-plot](Viewers/density-plot-run.md) | 🟢 PASS → PASS | ⚪ +0% (58/58 → 58/58) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | — | 🟡 2m (removed) | — | 🟡 2m (removed) |
| Viewers | [form](Viewers/form-run.md) | 🟢 PASS → PASS | 🔴 -1% (28/30 → 24/26) | 🟡 92% · status unchanged (PASS) | 🔴 93·93·92 | 🟡 18m (removed) | 🟡 4m (removed) | ⚪ +0s | 🟢 -22m |
| Viewers | [forms](Viewers/forms-run.md) | 🟢 PASS → PASS | ⚪ +0% (36/36 → 30/30) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | 🟡 18m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -21m |
| Viewers | [grid](Viewers/grid-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (16/22 → 16/22) | 🟡 73% · status unchanged (PARTIAL) | ⚪ 73·73·73 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | 🟢 PASS → PASS | 🔴 -1% (15/16 → 13/14) | 🟡 93% · status PARTIAL → PASS | 🔴 94·94·93 | 🟡 18m (removed) | 🟡 4m (removed) | ⚪ +0s | 🟢 -22m |
| Viewers | [histogram](Viewers/histogram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (87/94 → 87/94) | 🟡 93% · status unchanged (PARTIAL) | ⚪ 93·93·93 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (16/19 → 16/19) | 🟡 84% · status PASS → PARTIAL | ⚪ 84·84·84 | 🟡 20m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -23m |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (8/12 → 8/12) | 🟡 67% · status unchanged (PARTIAL) | ⚪ 67·67·67 | ⚪ +0s | ⚪ +0s | ⚪ +0s | ⚪ +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (17/20 → 17/20) | 🟡 85% · status PASS → PARTIAL | ⚪ 85·85·85 | 🟡 16m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -19m |
| Viewers | [row-source](Viewers/row-source-run.md) | 🟢 PASS → PASS | ⚪ +0% (36/36 → 36/36) | 🟢 100% · status unchanged (PASS) | ⚪ 100·100·100 | 🟡 4m (removed) | ⚪ +0s | ⚪ +0s | 🟢 -4m |
| Viewers | [statistics](Viewers/statistics-run.md) | 🟢 PASS → PASS | ⚪ +0% (24/24 → 24/24) | 🟡 new (prev7d had no entry) | ⚪ 100·100 | ⚪ +0s | ⚪ +0s | 🟡 2m (removed) | 🟢 -2m |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (70/84 → 70/84) | 🟡 84% · status unchanged (PARTIAL) | ⚪ 84·84·84 | ⚪ +0s | ⚪ +0s | 🟡 1m 48s (removed) | 🟢 -1m 48s |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | 🟢 PASS → PASS | ⚪ +0% (7/7 → 7/7) | 🟡 new (prev7d had no entry) | ⚪ 100·100·100 | 🟢 -1m 58s | 🟢 -10s | 🟢 -1m 40s | 🟢 -3m 48s |
| Viewers | [word-cloud-tests](Viewers/word-cloud-tests-run.md) | 🟡 PARTIAL → PARTIAL | ⚪ +0% (23/27 → 39/46) | 🟡 new (prev7d had no entry) | ⚪ 85·85 | 🟢 -2m | 🔴 +45s | 🔴 +17s | 🟢 -58s |
| Viewers | [working-with-nan-infinity](Viewers/working-with-nan-infinity-run.md) | 🟢 PASS → PASS | ⚪ +0% (9/9 → 9/9) | 🟡 new (prev7d had no entry) | ⚪ 100·100 | 🟡 6m (removed) | 🟡 3m (removed) | ⚪ +0s | 🟢 -9m |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage is 77% (136/177). Partial folders: Bio, Browse, Charts, Chem, Connections, DiffStudio, EDA, Models, Peptides, PowerPack, Projects, Scripts, StickyMeta, Viewers. No-data folders: Apps, General, LocalCashing, Notebooks, Queries, Tooltips, WideSmokeTest. No folder is in FAIL state overall, but several contain per-test failures that need attention.

### Blocking / Partial Issues
- Bio (PARTIAL): see per-test rows above.
- Browse (PARTIAL): see per-test rows above.
- Charts (PARTIAL): see per-test rows above.
- Chem (PARTIAL): see per-test rows above.
- Connections (PARTIAL): see per-test rows above.
- DiffStudio (PARTIAL): see per-test rows above.
- EDA (PARTIAL): see per-test rows above.
- Models (PARTIAL): see per-test rows above.
- Peptides (PARTIAL): see per-test rows above.
- PowerPack (PARTIAL): see per-test rows above.
- Projects (PARTIAL): see per-test rows above.
- Scripts (PARTIAL): see per-test rows above.
- StickyMeta (PARTIAL): see per-test rows above.
- Viewers (PARTIAL): see per-test rows above.
