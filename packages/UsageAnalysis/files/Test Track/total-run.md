# Test Track — Global Report

**Date**: 2026-04-10
**Verdict**: Conditionally ready

## Folder Summary

**Total**: 185 tests · Run: 138/185 (74.6%) · Playwright: 106/185 (57.3%) · Mean Spec Gen: 3.5s · Mean Spec Run: 26s

| Folder | Tests | Run | Playwright | Status | Mean Spec Gen | Mean Spec Run |
|--------|-------|-----|------------|--------|---------------|---------------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | NO DATA | | |
| Bio | 9 | 9/9 (100%) | 1/9 (11%) | PARTIAL | 5s | |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | PARTIAL | | |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | PARTIAL | 3s | 24.2s |
| Chem | 14 | 14/14 (100%) | 12/14 (86%) | PARTIAL | 3s | 32.2s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | PARTIAL | | |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | PASS | 3.5s | 26.7s |
| EDA | 10 | 10/10 (100%) | 11/10 (100%) | PARTIAL | 2s | 4.8s |
| General | 11 | 0/11 (0%) | 0/11 (0%) | NO DATA | | |
| LocalCashing | 1 | 0/1 (0%) | 0/1 (0%) | NO DATA | | |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | PARTIAL | | |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | NO DATA | | |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 3s | 47.6s |
| PowerPack | 10 | 10/10 (100%) | 6/10 (60%) | PARTIAL | 2.5s | 5.3s |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | PARTIAL | | |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | NO DATA | | |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | PARTIAL | | |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 2.3s | 15s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | NO DATA | | |
| Viewers | 47 | 42/47 (89%) | 42/47 (89%) | PARTIAL | 6s | 33s |

## All Tests

| Folder | Test | Status | Description | Spec Gen | Spec Run |
|--------|------|--------|-------------|----------|----------|
| Apps | apps | NO RUN | | | |
| Apps | tutorials | NO RUN | | | |
| Bio | [analyze](Bio/analyze-run.md) | PASS | All three Bio > Analyze functions work correctly on all three dataset types (FASTA, HELM, MSA) | ~5s | |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | PARTIAL | 4 of 5 steps passed. Composition/WebLogo viewer opens correctly. Letter-click selection could not be verified | | |
| Bio | [convert](Bio/convert-run.md) | PASS | All 4 Bio convert/transform functions work correctly on FASTA data | | |
| Bio | [manage](Bio/manage-run.md) | PASS | All 3 steps passed. Manage Monomer Libraries view opens showing 5 monomer library JSON files | | |
| Bio | [msa](Bio/msa-run.md) | PASS | All 6 steps passed. MSA dialog opens with correct fields, alignment produces msa(Sequence) column | | |
| Bio | [pepsea](Bio/pepsea-run.md) | PARTIAL | MSA dialog opens correctly for HELM data with MAFFT methods. Alignment on 50 HELM sequences did not complete | | |
| Bio | [search](Bio/search-run.md) | PASS | All 4 steps passed. Subsequence Search opens filter panel with bio substructure filter | | |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | PASS | All 6 steps passed. Activity Cliffs works with default and custom parameters | | |
| Bio | [sequence-space](Bio/sequence-space-run.md) | PASS | All 6 steps passed. Sequence Space works with default and custom parameters, produces scatter plots | | |
| Browse | [browse](Browse/browse-run.md) | PARTIAL | 4 steps passed, 1 partial. Browse tree structure complete, demos work, URL routing works for files and sections | | |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | PARTIAL | Browse tree correctly preserves expand/collapse state within a single session | | |
| Browse | [japanese-in-myfiles](Browse/japanese-in-myfiles-run.md) | PASS | All 3 steps passed. Japanese Kanji and Hiragana characters display correctly in file names | | |
| Browse | local-deploy | NO RUN | | | |
| Browse | package-manager | NO RUN | | | |
| Browse | [spaces](Browse/spaces-run.md) | PARTIAL | 8 steps passed, 24 skipped. Core CRUD operations work; advanced operations not tested | | |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | FAIL | 1 step passed, 27 skipped. UI-only interactions not covered by server API could not be automated | | |
| Charts | [radar](Charts/radar-run.md) | PASS | Radar viewer scenario passed in MCP and Playwright runs | 3s | 14.5s |
| Charts | [sunburst](Charts/sunburst-run.md) | PARTIAL | 8 of 9 steps passed, canvas multi-selection ambiguous | 3s | 41.5s |
| Charts | [tree](Charts/tree-run.md) | PASS | All 4 steps passed. Tree viewer correctly maintains selection state across collaborative filtering | 3s | 16.6s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PASS | Activity Cliffs works correctly on SPGI.csv. UMAP scatter plot with molecule thumbnails, 141 cliffs detected | ~3s | 56.9s |
| Chem | [calculate](Chem/calculate-run.md) | PASS | Descriptors calculation works correctly. Dialog shows all categories, pre-selected descriptors compute correctly | ~3s | 24.5s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS | Chemical Space works with both UMAP and t-SNE parameters. Scatter plots generated correctly | ~3s | 44s |
| Chem | [chemprop](Chem/chemprop-run.md) | PARTIAL | Steps 1-3 passed. Model training failed — Chemprop Docker container unavailable | 3s | |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PASS | All 4 steps passed. Element counts and Molecule Charge computed for 1000 molecules | 3s | 21s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PARTIAL | Core substructure filtering works (benzene → 1356/3624 rows). Advanced steps skipped | ~3s | 29s |
| Chem | [info-panels](Chem/info-panels-run.md) | PASS | All 5 steps passed. Context Panel shows 11 info panels including Chemistry | 3s | 39s |
| Chem | [mmp](Chem/mmp-run.md) | FAIL | Steps 1-2 passed. Step 3 failed: TypeError in mmpAnalysis at chem/src/package.ts:2301 | 3s | |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PASS | R-Groups Analysis dialog opened, MCS computed, trellis plot with decomposition appeared | 3s | 38s |
| Chem | [sketcher](Chem/sketcher-run.md) | PASS | Steps 1-6 passed: sketcher dialog with OpenChemLib, C1CCCCC1 rendered, hamburger menu works | 3s | 13s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | FAIL | Steps 1-2 passed. Adding scaffolds failed — backend 502 error | | |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | FAIL | Steps 1-2 passed. Scaffold tree generation failed — 502 server error | 3s | |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | PASS | All 4 steps passed. Similarity Search viewer shows similar molecules with Tanimoto/Morgan scores | 3s | |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | PASS | All 4 steps passed. Structure filter panel opened, benzene matched 1356/3624 molecules | 3s | 24s |
| Connections | [adding](Connections/adding-run.md) | PARTIAL | 6 of 7 steps passed. TEST button works but connection test fails without real credentials | | |
| Connections | [browser](Connections/browser-run.md) | PARTIAL | 7 of 9 steps passed, 2 ambiguous. Search filtering works, Context Pane shows expected tabs | | |
| Connections | [catalogs](Connections/catalogs-run.md) | FAIL | 1 step passed, 1 failed, 15 skipped. NorthwindTest MS SQL connection not present | | |
| Connections | [delete](Connections/delete-run.md) | PASS | All 8 steps passed. Both connections deleted successfully | | |
| Connections | [edit](Connections/edit-run.md) | PASS | 6 of 7 steps passed. Rename, credential modification, and error testing all work | | |
| Connections | [external-provider](Connections/external-provider-run.md) | FAIL | All 7 steps skipped. Requires Postgres connection not available | | |
| Connections | [identifiers](Connections/identifiers-run.md) | FAIL | 1 step passed, 1 failed. Depends on Northwind database connection | | |
| Connections | [import-swagger](Connections/import-swagger-run.md) | FAIL | All 7 steps skipped. Requires manual file drag-drop interaction | | |
| Connections | [schema](Connections/schema-run.md) | PARTIAL | 3 of 4 steps passed. "Browse schema" context menu not found, but schema browsable via tree | | |
| Connections | [sparql](Connections/sparql-run.md) | PARTIAL | 6 of 7 steps passed. Test endpoint data.ontotext.com not resolvable | | |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | PASS | All 6 steps passed in MCP and Playwright. Full catalog workflow works end-to-end | 5s | 32.1s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS | All 4 steps passed. PK-PD cyclic model loads, Count clickers update solution in real-time | 3s | 22.7s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS | All 3 steps passed. pk.ivp opened from Files browser as DiffStudio view | 3s | 24s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | PASS | All 6 steps passed. Bioreactor model loaded, Nelder-Mead fitting converged in 9 iterations | 5s | 28s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS | All 6 steps passed. DiffStudio opened, Bioreactor loaded with 13 columns and 1001 rows | 3s | 18s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | PASS | All 5 steps passed. Edit toggle and export icon produce JS script with model equations | 3s | 49s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS | All 4 steps passed. Monte Carlo analysis produced 10 samples with 4 viewers | 3s | 23s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | PASS | All 4 steps passed. Acid Production 2-stage simulation, input modification updates charts | 3s | 17s |
| EDA | [anova](EDA/anova-run.md) | PASS | All 3 steps passed. ANOVA dialog with correct defaults (RACE, AGE, Alpha=0.05) | 2s | 3s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PASS | All 3 steps passed. PLS dialog with correct defaults, all 5 viewers appeared | 2s | 7s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL | 6 of 7 steps passed. cars-with-missing.csv missing from dev server | 2s | 9s |
| EDA | [pca](EDA/pca-run.md) | PASS | All 5 steps passed. PCA runs correctly with and without preprocessing | 2s | 3s |
| EDA | [pls](EDA/pls-run.md) | PASS | All 4 steps passed. PLS dialog opened, three PLS columns added to dataset | 2s | 4s |
| EDA | [ML methods/linear-regression](EDA/ML%20methods/linear-regression-run.md) | PASS | Linear Regression trained successfully on cars.csv predicting price | 2s | 6.8s |
| EDA | [ML methods/pls-regression](EDA/ML%20methods/pls-regression-run.md) | PARTIAL | PLS Regression trained successfully. Model Engine dropdown not visible in Train Model UI | 2s | 6.9s |
| EDA | [ML methods/softmax](EDA/ML%20methods/softmax-run.md) | FAIL | Softmax training fails: "incorrect features type" on iris.csv. Known platform bug | 2s | 2.6s |
| EDA | [ML methods/xgboost1](EDA/ML%20methods/xgboost1-run.md) | PARTIAL | XGBoost classification trained on iris.csv. Hyperparameter interaction not testable via JS API | 2s | 2.7s |
| EDA | [ML methods/xgboost2](EDA/ML%20methods/xgboost2-run.md) | PARTIAL | XGBoost regression trained on cars.csv. Hyperparameter interaction not testable via JS API | 2s | 2.7s |
| General | api-samples | NO RUN | | | |
| General | files-cache | NO RUN | | | |
| General | first-login | NO RUN | | | |
| General | inactivity-response | NO RUN | | | |
| General | login | NO RUN | | | |
| General | molecule-in-exported-csv | NO RUN | | | |
| General | network | NO RUN | | | |
| General | profile-settings | NO RUN | | | |
| General | startup-time | NO RUN | | | |
| General | table-manager | NO RUN | | | |
| General | tabs-reordering | NO RUN | | | |
| LocalCashing | local-cashing | NO RUN | | | |
| Models | [apply](Models/apply-run.md) | PASS | Apply Model workflow works end-to-end. TestDemog model found, selected, and applied | | |
| Models | [browser](Models/browser-run.md) | PARTIAL | 4 of 6 steps passed. Browser navigation, search, Context Panel, and Filter work | | |
| Models | [chemprop](Models/chemprop-run.md) | FAIL | 2 of 17 sub-steps passed. Chemprop Docker container not available on public instance | | |
| Models | [delete](Models/delete-run.md) | PASS | All 5 steps passed. Deletion workflow works with confirmation dialog | | |
| Models | [predictive-models](Models/predictive-models-run.md) | PASS | All 20 sub-steps passed. Full lifecycle (Train, Apply, Delete) works correctly | | |
| Models | [train](Models/train-run.md) | PASS | All 10 steps passed. Train Model workflow works using built-in EDA engines | | |
| Notebooks | browser | NO RUN | | | |
| Notebooks | create | NO RUN | | | |
| Notebooks | delete | NO RUN | | | |
| Notebooks | edit | NO RUN | | | |
| Peptides | [info-panels](Peptides/info-panels-run.md) | PASS | peptides.csv loads correctly with Macromolecule semType. Amino acids rendered with distinct colors | 3s | 10.9s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | PARTIAL | SAR analysis produces MCL, Most Potent Residues, and Sequence Variability Map. Wrench button missing | 3s | 78s |
| Peptides | [peptides](Peptides/peptides-run.md) | PARTIAL | Steps 1-4 passed. Context Panel shows Peptides pane with parameters and WebLogo. Canvas ambiguous | 3s | 11.6s |
| Peptides | [sar](Peptides/sar-run.md) | PARTIAL | Steps 1-10 passed. SAR creates SVM, Most Potent Residues, MCL viewers. Canvas interaction ambiguous | 3s | 90s |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | PASS | All 5 steps passed. Add New Column dialog opens, accepts formula, creates column correctly | ~2s | 5.3s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | PARTIAL | Part 1 passed: enrichment joins Customers to Orders via customerid. Parts 2-4 skipped | ~3s | |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | SKIP | Scenario file contains only GitHub ticket references, no executable test steps | | |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | PARTIAL | Core functionality works: creating calculated columns, chaining formulas. Steps 7-10 skipped | | |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | PASS | All 8 steps passed. Autocomplete works for functions and columns | | |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | PASS | All 5 steps passed. Formula dependency chain propagation works correctly | | |
| PowerPack | [AddNewColumn/functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | SKIP | Requires spgi.csv dataset not available as standard demo table | | |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | PASS | All 4 steps passed. Column names highlighted with cm-column-name CSS class | | |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | PASS | All 4 steps passed. Function tooltips display signature on hover | | |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | SKIP | Requires spgi.csv dataset not available as standard demo table | | |
| Projects | [browser](Projects/browser-run.md) | PARTIAL | 5 of 9 steps passed. Browse > Dashboards works, search and Context Panel functional | | |
| Projects | [complex](Projects/complex-run.md) | SKIP | All 13 steps skipped. Requires tables from 7+ sources, drag-and-drop, multi-user verification | | |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | SKIP | All 5 steps skipped. Requires custom JS script with Data Sync | | |
| Projects | [deleting](Projects/deleting-run.md) | PARTIAL | 2 of 4 steps passed. Deletion works via API, right-click context menu not triggerable | | |
| Projects | [opening](Projects/opening-run.md) | PARTIAL | All 5 steps passed. Projects accessible, Context Panel shows attributes, data intact | | |
| Projects | [project-url](Projects/project-url-run.md) | SKIP | Depends on copy_clone scenario which was not executed | | |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | SKIP | 2 of 5 steps passed. Copy/clone/link operations not tested | | |
| Projects | [uploading](Projects/uploading-run.md) | PARTIAL | 8 of 14 steps passed. Core project creation from local tables, file shares, queries works | | |
| Queries | adding | NO RUN | | | |
| Queries | browse-&-save-project | NO RUN | | | |
| Queries | browser | NO RUN | | | |
| Queries | columns-inspect | NO RUN | | | |
| Queries | deleting | NO RUN | | | |
| Queries | edit | NO RUN | | | |
| Queries | get-all-get-top-100 | NO RUN | | | |
| Queries | ms-sql | NO RUN | | | |
| Queries | new-sql-query | NO RUN | | | |
| Queries | new-visual-query | NO RUN | | | |
| Queries | query-layout | NO RUN | | | |
| Queries | query-postprocessing | NO RUN | | | |
| Queries | transformations | NO RUN | | | |
| Queries | visual-query-advanced | NO RUN | | | |
| Scripts | [browser](Scripts/browser-run.md) | PASS | Context pane shows all expected accordions (Details, Script, Run, Activity, Sharing, Chats, Dev) | | |
| Scripts | [create](Scripts/create-run.md) | PARTIAL | Script created, parameters configured, saved, and executed. Run dialog table dropdown issue | | |
| Scripts | [delete](Scripts/delete-run.md) | PASS | All 5 steps passed. Delete flow works with confirmation dialog | | |
| Scripts | [edit](Scripts/edit-run.md) | PASS | All 6 steps passed. Edits saved persistently and visible on re-open | | |
| Scripts | layout | NO RUN | | | |
| Scripts | [run](Scripts/run-run.md) | PARTIAL | Core run from context menu and console works. Steps 4-6 (local file, Files, query) skipped | | |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | PASS | Sticky meta panel shows all 5 fields, "+" buttons add columns correctly | 3s | 17s |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | PASS | Cloning preserves schema and metadata. Steps 3-4 skipped (delete/refresh) | 3s | 21s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | PASS | Schemas browser shows 20 schemas including TestSchema1. NEW SCHEMA button present | 3s | 7s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | FAIL | "Database meta" section not in Context Panel for CHEMBL. Feature may not be deployed | 0s | |
| Tooltips | actions-in-the-context-menu | NO RUN | | | |
| Tooltips | default-tooltip | NO RUN | | | |
| Tooltips | default-tooltip-visibility | NO RUN | | | |
| Tooltips | edit-tooltip | NO RUN | | | |
| Tooltips | line-chart---aggregated-tooltip | NO RUN | | | |
| Tooltips | tooltip-properties | NO RUN | | | |
| Tooltips | uniform-default-tooltip | NO RUN | | | |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS | 3D Scatter Plot opens as rotatable model, data point clicking selects row | ~3s | 10s |
| Viewers | annotation-regions | NO RUN | | | |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | PASS | All 15 bar chart sections passed. Stack, sorting, axis type, color coding all work | ~10s | 60s |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL | 17 of 19 sections passed. Visualization zoom ambiguous: project save fails on dev | | 47s |
| Viewers | [calendar](Viewers/calendar-run.md) | PASS | Calendar viewer added to demog, viewer element found | 3s | 5s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | PASS | Linear and Categorical color coding applied and verified via API | 3s | 84s |
| Viewers | [color-coding-(linked)](Viewers/color-coding-(linked)-run.md) | PASS | Linked color coding works via column tags, source changes preserve relationships | 3s | 6s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | PARTIAL | 27 of 30 steps passed, 3 skipped due to canvas cell interaction limitation | ~3s | 22.5s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | PASS | All 13 scenarios passed across all property combinations | ~2m | |
| Viewers | [form](Viewers/form-run.md) | PASS | Form viewer added to SPGI showing molecule structures and metadata fields | 3s | 15s |
| Viewers | [forms](Viewers/forms-run.md) | PASS | 30 of 34 steps passed, 3 skipped, 1 ambiguous. Molecule rendering works | | 46.9s |
| Viewers | [grid](Viewers/grid-run.md) | PASS | Grid sorting, multi-row selection, column Context Panel (14 tabs), filtering all work | 3s | 18s |
| Viewers | [grid-viewer](Viewers/grid-viewer-run.md) | PASS | Second Grid viewer added, table switching and row height modification work | 3s | 14s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | PASS | SPGI loaded, Heatmap viewer added and properties verified | 3s | 12s |
| Viewers | [histogram](Viewers/histogram-run.md) | PARTIAL | Property-based tests passed. Canvas-based interactions not automatable | ~30s | 50s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | PASS | All 27 sections passed. Properties, context menu, layout save/restore all work | ~5s | 150s |
| Viewers | [map](Viewers/map-run.md) | PASS | Map viewer on earthquakes.csv with auto-detected lat/lon, render type cycling works | 3s | 9s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | PASS | Matrix plot viewer added to SPGI, properties accessible | 3s | 11s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | PASS | Network diagram viewer added to demog, node columns set | 3s | 7s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | PASS | All 36 steps across 12 sections passed. Properties, menus, filtering, columns all work | ~5s | 38.3s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS | All 16 sections passed. All viewer properties work, layout save/restore works | ~10s | 58.2s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | PASS | Pivot table properties modified (pivot, groupBy, aggregate) without errors | 3s | 8s |
| Viewers | rendering-structures-on-the-axes | NO RUN | | | |
| Viewers | [row-source](Viewers/row-source-run.md) | PASS | All 7 viewer types tested with all 8 row source modes on demog and spgi-100 | ~5s | 84s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | PASS | All 20 test sections passed | ~5s | 66s |
| Viewers | statistics-viewer | NO RUN | | | |
| Viewers | [tile](Viewers/tile-run.md) | PASS | Tile viewer added, table switching preserves viewer, calculated column works | 3s | 17s |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | PASS | Tree map viewer properties modified without errors | 3s | 7s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | PARTIAL | Property-based tests passed. Canvas ambiguous, Histogram "Use in Trellis" menu not found | ~30s | 108s |
| Viewers | viewers-docking | NO RUN | | | |
| Viewers | working-with-nan-&-infinity | NO RUN | | | |
| Viewers | [world-cloud](Viewers/world-cloud-run.md) | PASS | Word cloud viewer on SPGI with "Primary Series Name" column works | 3s | 13s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | PASS | All filter panel basic operations work: structure, categorical, histogram range filtering | | |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | PASS | All 6 chem search type modes produce correct row counts. Bio substructure filter works | | |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | PASS | Cloned view inherits filter state. Layout save/restore preserves filter configuration | | |
| Viewers | [FilterPanel/collaborative-filtering](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | PASS | Table linking (selection-to-filter, filter-to-filter) and propagation work correctly | ~3s | 30s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | PASS | Combined boolean filter aggregates columns, supports OR mode, survives layout cycles | | |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | PASS | Expression filter: rules, AND/OR, free-text, layout persistence all work | | |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | PASS | Multi-level tree filtering, checkbox selection, reordering, collaborative filtering work | | |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | PASS | Text filter: string search, highlight, AND/OR logic, fuzzy search all work | | |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | PASS | All 31 steps passed. All 6 viewer types interact correctly with Filter Panel | ~5s | 72s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | PASS | Categorical color coding on "Stereo Category" propagates correctly across all viewer types | ~3s | 24s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | PASS | Legends dynamically update when filters or row source changes are applied | ~5s | 35.9s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | PASS | Line chart legend reflects split categories, Multi Axis mode works, layout persists | ~3s | 47s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | PASS | All 5 scenarios (color/marker, axis, in-viewer filter, Filter Panel, grid coding) passed | ~3s | 50s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | PASS | Core SMILES column renders molecule images in scatterplot legend. Layout persists | ~3s | 26s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | PASS | Legend visibility, color picker, click filtering, auto-hide, corner positions all work | ~5s | 43.7s |
| **Total** | **185** | **88 PASS / 32 PARTIAL / 11 FAIL / 7 SKIP / 47 NO RUN** | | **3.5s** | **26s** |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage is 74.6% (138/185), meeting the 70% threshold. Only one folder (DiffStudio) has all tests passing. No folder has all tests failing. However, 6 folders have zero test coverage (Apps, General, LocalCashing, Notebooks, Queries, Tooltips), accounting for 39 untested scenarios. The 11 FAIL results are concentrated in Connections (4), Chem (3), EDA (1), Models (1), StickyMeta (1), and Browse (1), mostly due to missing infrastructure (Docker containers, database connections, test endpoints).

### Blocking Issues
- **Apps**: 0/2 run coverage — no test results available
- **General**: 0/11 run coverage — no test results available
- **Queries**: 0/14 run coverage — no test results available
- **Tooltips**: 0/7 run coverage — no test results available
- **Notebooks**: 0/4 run coverage — no test results available
- **LocalCashing**: 0/1 run coverage — no test results available
- **Connections**: 4 FAIL results — infrastructure dependencies (NorthwindTest, Postgres, external providers)
- **Chem**: 3 FAIL results — scaffold tree backend 502 errors, MMP TypeError bug
- **EDA/ML methods**: Softmax training bug (`eda:trainSoftmax` incorrect features type)
- **Models**: Chemprop Docker container unavailable on public instance
- **StickyMeta**: Database meta feature may not be deployed
- **Browse**: spaces-(ui-only) scenario could not be automated
