# Test Track — Global Report

**Date**: 2026-04-08
**Verdict**: Not ready

## Folder Summary

| Folder | Tests | Run | Playwright | Status | Mean Spec Gen | Mean Spec Run |
|--------|-------|-----|------------|--------|---------------|---------------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | NO DATA | | |
| Bio | 9 | 9/9 (100%) | 1/9 (11%) | PARTIAL | 5s | |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | PARTIAL | | |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | PARTIAL | 3s | |
| Chem | 14 | 14/14 (100%) | 8/14 (57%) | PARTIAL | 3s | 39.3s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | PARTIAL | | |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | PASS | | |
| EDA | 10 | 5/10 (50%) | 0/10 (0%) | PARTIAL | | |
| General | 11 | 0/11 (0%) | 0/11 (0%) | NO DATA | | |
| LocalCashing | 1 | 0/1 (0%) | 0/1 (0%) | NO DATA | | |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | PARTIAL | | |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | NO DATA | | |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | | |
| PowerPack | 10 | 10/10 (100%) | 6/10 (60%) | PARTIAL | | |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | PARTIAL | | |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | NO DATA | | |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | PARTIAL | | |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | | |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | NO DATA | | |
| Viewers | 66 | 24/66 (36%) | 24/66 (36%) | PARTIAL | 7.4s | 56.4s |
| WideSmokeTest | 23 | 1/23 (4%) | 0/23 (0%) | PARTIAL | | |
| **Total** | **231** | **116/231 (50%)** | **78/231 (34%)** | | **4.9s** | **50.3s** |

## All Tests

| Folder | Test | Status | Description | Spec Gen | Spec Run |
|--------|------|--------|-------------|----------|----------|
| Apps | apps | NO RUN | | | |
| Apps | tutorials | NO RUN | | | |
| Bio | [analyze](Bio/analyze-run.md) | PASS | All three Bio > Analyze functions work correctly on all three dataset types (FASTA, HELM, MSA) | ~5s | |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | PARTIAL | 4 of 5 steps passed. WebLogo viewer opens correctly; letter-click selection not verifiable (canvas-based) | | |
| Bio | [convert](Bio/convert-run.md) | PASS | All 4 Bio convert/transform functions work correctly on FASTA data | | |
| Bio | [manage](Bio/manage-run.md) | PASS | All 3 steps passed. Manage Monomer Libraries view opens showing 5 monomer library JSON files | | |
| Bio | [msa](Bio/msa-run.md) | PASS | All 6 steps passed. MSA dialog opens with correct fields, alignment runs correctly using Kalign 3.3.1 | | |
| Bio | [pepsea](Bio/pepsea-run.md) | AMBIGUOUS | MSA dialog opens for HELM data but alignment computation did not complete within 3+ min timeout | | |
| Bio | [search](Bio/search-run.md) | PASS | All 4 steps passed. Subsequence Search filter works correctly | | |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | PASS | All 6 steps passed. Activity Cliffs works with both default and custom parameters | | |
| Bio | [sequence-space](Bio/sequence-space-run.md) | PASS | All 6 steps passed. Sequence Space works with both default and custom parameters | | |
| Browse | [browse](Browse/browse-run.md) | PARTIAL | 4 steps passed, 1 partial. Item-level selection NOT preserved across page refresh | | |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | PARTIAL | Tree state preserved within session but NOT persisted across page refreshes | | |
| Browse | [japanese-in-myfiles](Browse/japanese-in-myfiles-run.md) | PASS | All 3 steps passed. Japanese characters display correctly | | |
| Browse | local-deploy | NO RUN | | | |
| Browse | package-manager | NO RUN | | | |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | FAIL | 1 passed, 2 ambiguous, 27 skipped. Drag-and-drop not automatable via Chrome DevTools MCP | | |
| Browse | [spaces](Browse/spaces-run.md) | PARTIAL | 8 passed, 24 skipped. Core CRUD works; drag-and-drop and share tests skipped | | |
| Charts | [radar](Charts/radar-run.md) | PASS | Radar viewer opens correctly, properties work for table switching, column selection, style changes | 3s | |
| Charts | [sunburst](Charts/sunburst-run.md) | PARTIAL | 7 of 9 passed. Core functionality works; canvas interaction and old layout steps skipped | 3s | |
| Charts | [tree](Charts/tree-run.md) | PASS | All 4 steps passed. Collaborative filtering works correctly with exact counts | 3s | |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PASS | UMAP embedding produces scatter plot with molecule thumbnails, 141 cliffs detected | ~3s | 56.9s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | PARTIAL | Scaffold Tree opens; complex multi-step interactions require manual testing | | |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | PARTIAL | Scaffold Tree opens but could not generate tree (smiles-50.csv missing, smiles.csv exceeds limit) | | |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | PASS | Similarity Search works correctly with Tanimoto/Morgan fingerprint scores | | |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | PASS | Substructure filtering by benzene correctly reduces rows from 3624 to 1356 | | |
| Chem | [calculate](Chem/calculate-run.md) | PASS | Descriptors calculation works correctly, dialog shows all categories | ~3s | 24.5s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS | Chemical Space works with both UMAP and t-SNE parameters | ~3s | 44.0s |
| Chem | [chemprop](Chem/chemprop-run.md) | SKIP | SDF file loading and ML model training exceed automation limits | | |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PARTIAL | MCP run passed but Playwright spec failed (column count mismatch) | ~3s | 20.9s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PARTIAL | Core substructure filtering works; advanced interactions skipped | ~3s | 29.0s |
| Chem | [info-panels](Chem/info-panels-run.md) | PARTIAL | MCP passed but Playwright spec failed 3 of 9 steps (timing/loading issues) | ~3s | 34.7s |
| Chem | [mmp](Chem/mmp-run.md) | AMBIGUOUS | MMP dialog opens but computation on 20K molecules did not complete within timeout | | |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PASS | All R-Group Analysis steps passed including MCS computation and error handling | ~3s | 66s |
| Chem | [sketcher](Chem/sketcher-run.md) | FAIL | MCP passed core steps but Playwright spec failed entirely (page load timeout) | ~3s | 24.6s |
| Connections | [adding](Connections/adding-run.md) | PARTIAL | 6 of 7 passed. TEST button works but connection test fails without real credentials | | |
| Connections | [browser](Connections/browser-run.md) | PARTIAL | 7 of 9 passed. Search filtering works; "magic wand" filter icon not found in current UI | | |
| Connections | [catalogs](Connections/catalogs-run.md) | FAIL | NorthwindTest MS SQL connection not present on public.datagrok.ai | | |
| Connections | [delete](Connections/delete-run.md) | PASS | All 8 steps passed. Both connections deleted successfully | | |
| Connections | [edit](Connections/edit-run.md) | PASS | 6 of 7 passed. Rename, credential modification, error testing all work correctly | | |
| Connections | [external-provider](Connections/external-provider-run.md) | FAIL | All 7 steps skipped. Requires db.datagrok.ai:54327 with superuser credentials | | |
| Connections | [identifiers](Connections/identifiers-run.md) | FAIL | Depends on working Postgres connection with incorrect credentials | | |
| Connections | [import-swagger](Connections/import-swagger-run.md) | FAIL | Requires manual drag-drop and private API key — not automatable | | |
| Connections | [schema](Connections/schema-run.md) | PARTIAL | 3 of 4 passed. "Browse schema" menu not found but schema browseable via tree | | |
| Connections | [sparql](Connections/sparql-run.md) | PARTIAL | 6 of 7 passed. Test endpoint data.ontotext.com not resolvable from server | | |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | PASS | Full DiffStudio workflow works end-to-end: Library, Save to Model Hub, run from catalog | | |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS | PK-PD cyclic model loaded correctly; modifying Count dynamically updates all charts | | |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS | pk.ivp file opens from Files browser; URL sharing reproduces exact model state | | |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | PASS | Fit view, Process mode, Nelder-Mead fitting all work correctly with 9 iterations | | |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS | Bioreactor example loaded from Library; Multiaxis and Facet tabs work correctly | | |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | PASS | Edit toggle, JS script view, and Model Hub registration all work correctly | | |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS | Sensitivity Analysis produces all 4 viewers (Correlation, PC plot, Scatter, Grid) | | |
| DiffStudio | [stages](DiffStudio/stages-run.md) | PASS | 2-stage model shows distinct visual stages; modifying duration immediately updates charts | | |
| EDA | [anova](EDA/anova-run.md) | PASS | ANOVA dialog works with correct defaults; box plot and F-test interpretation correct | | |
| EDA | ML methods/linear-regression | NO RUN | | | |
| EDA | ML methods/pls-regression | NO RUN | | | |
| EDA | ML methods/softmax | NO RUN | | | |
| EDA | ML methods/xgboost1 | NO RUN | | | |
| EDA | ML methods/xgboost2 | NO RUN | | | |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PASS | MVA/PLS dialog produces all 5 expected viewers | | |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL | Viewer renders correctly but empty column incorrectly appears in dropdowns | | |
| EDA | [pca](EDA/pca-run.md) | PASS | PCA works with and without preprocessing (Center+Scale) | | |
| EDA | [pls](EDA/pls-run.md) | PARTIAL | PLS executes but does not add score columns as scenario expects | | |
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
| Models | [apply](Models/apply-run.md) | PASS | Apply Model workflow works end-to-end; prediction column added correctly | | |
| Models | [browser](Models/browser-run.md) | PARTIAL | 4 of 6 passed. Multi-select and Compare untestable (only 1 model available) | | |
| Models | [chemprop](Models/chemprop-run.md) | FAIL | Chemprop Docker container not available on public.datagrok.ai | | |
| Models | [delete](Models/delete-run.md) | PASS | All 5 steps passed. Model deletion workflow works correctly | | |
| Models | [predictive-models](Models/predictive-models-run.md) | PASS | All 20 sub-steps passed. Full Train/Apply/Delete lifecycle works | | |
| Models | [train](Models/train-run.md) | PASS | All 10 steps passed. Both classification and regression training work | | |
| Notebooks | browser | NO RUN | | | |
| Notebooks | create | NO RUN | | | |
| Notebooks | delete | NO RUN | | | |
| Notebooks | edit | NO RUN | | | |
| Peptides | [info-panels](Peptides/info-panels-run.md) | PASS | All 6 steps passed. Context Panel shows Details and Peptides panels correctly | | |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | PASS | SAR launched successfully; MCL viewer reflects new parameters after recalculation | | |
| Peptides | [peptides](Peptides/peptides-run.md) | PASS | All 6 steps passed. Weblogo chart interaction and row selection work correctly | | |
| Peptides | [sar](Peptides/sar-run.md) | PARTIAL | 12 of 13 passed. Distribution panel had no configurable controls (ambiguous) | | |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | PASS | All 5 steps passed. Formula entry, autocomplete, preview, and history work | | |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | PARTIAL | Core functionality works; datasync persistence steps skipped | | |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | PASS | All 8 steps passed. Autocomplete works for functions and columns | | |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | PASS | All 5 steps passed. Formula dependency chain propagates correctly | | |
| PowerPack | [AddNewColumn/functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | SKIP | Requires spgi.csv dataset not available on server | | |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | PASS | All 4 steps passed. Column names highlighted with cm-column-name CSS class | | |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | PASS | All 4 steps passed. Function tooltips display signature with types | | |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | SKIP | Requires spgi.csv dataset not available on server | | |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | SKIP | Requires Northwind database connection not available | | |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | SKIP | No executable test steps defined (references only) | | |
| Projects | [browser](Projects/browser-run.md) | PARTIAL | 5 of 9 passed. Right-click context menu not triggerable via automation | | |
| Projects | [complex](Projects/complex-run.md) | SKIP | Requires 7+ data sources, drag-and-drop, multi-user verification | | |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | SKIP | Requires server-side file manipulation | | |
| Projects | [deleting](Projects/deleting-run.md) | PARTIAL | 2 of 4 passed. API deletion works; context menu deletion not automatable | | |
| Projects | [opening](Projects/opening-run.md) | PARTIAL | All 5 passed. Projects accessible and open with data intact | | |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | SKIP | Copy/clone/link require multi-step UI not feasible via automation | | |
| Projects | [project-url](Projects/project-url-run.md) | SKIP | Depends on copy_clone scenario which was not executed | | |
| Projects | [uploading](Projects/uploading-run.md) | PARTIAL | 8 of 14 passed. Core project creation works; Spaces and Pivot/Aggregate skipped | | |
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
| Scripts | [browser](Scripts/browser-run.md) | PASS | Context pane shows all expected accordions; script details correct | | |
| Scripts | [create](Scripts/create-run.md) | PARTIAL | Script created and executed; Run dialog didn't auto-select pre-loaded table | | |
| Scripts | [delete](Scripts/delete-run.md) | PASS | All 5 steps passed. Delete flow works with confirmation dialog | | |
| Scripts | [edit](Scripts/edit-run.md) | PASS | All 6 steps passed. Edits saved persistently | | |
| Scripts | layout | NO RUN | | | |
| Scripts | [run](Scripts/run-run.md) | PARTIAL | Core run works from context menu and console; file-based inputs skipped | | |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | PARTIAL | Sticky meta section appears; TestSchema1 not listed due to missing entity type binding | | |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | FAIL | All steps skipped — depends on add-and-edit which didn't fully succeed | | |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | PASS | Entity type and schema created; "Associated with" binding failed via UI | | |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | FAIL | "Database meta" section not found; NorthwindTest connection not present | | |
| Tooltips | actions-in-the-context-menu | NO RUN | | | |
| Tooltips | default-tooltip | NO RUN | | | |
| Tooltips | default-tooltip-visibility | NO RUN | | | |
| Tooltips | edit-tooltip | NO RUN | | | |
| Tooltips | line-chart---aggregated-tooltip | NO RUN | | | |
| Tooltips | tooltip-properties | NO RUN | | | |
| Tooltips | uniform-default-tooltip | NO RUN | | | |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS | 3D Scatter Plot opens correctly; clicking, zoom, Color/Size properties work | ~3s | 10s |
| Viewers | annotation-regions | NO RUN | | | |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | PASS | All 15 sections passed. Stack, sorting, color coding, layout save/restore all work | ~10s | 60s |
| Viewers | bar-chart-tests | NO RUN | | | |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL | Most functionality works; zoom persistence differs from scenario expectation | ~5s | ~104s |
| Viewers | box-plot-tests | NO RUN | | | |
| Viewers | box-plot-tests-pw | NO RUN | | | |
| Viewers | calendar | NO RUN | | | |
| Viewers | color-coding-(linked) | NO RUN | | | |
| Viewers | color-coding | NO RUN | | | |
| Viewers | correlation-plot | NO RUN | | | |
| Viewers | correlation-plot-tests | NO RUN | | | |
| Viewers | correlation-plot-tests-pw | NO RUN | | | |
| Viewers | [density-plot](Viewers/density-plot-run.md) | PASS | Density plot opens with hexagonal bins; axes and zoom work correctly | ~3s | 8s |
| Viewers | density-plot-tests | NO RUN | | | |
| Viewers | density-plot-tests-pw | NO RUN | | | |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | PASS | All filter panel basic operations work correctly | | |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | PASS | All 6 chem search modes and bio substructure filter produce correct row counts | | |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | PASS | All 14 steps passed. Cloned view inherits filter state correctly | | |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | PASS | All 9 steps passed. Table linking filter propagation works correctly | ~3s | 30s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | PASS | All 15 steps passed. Boolean filter with OR mode and layout persistence work | | |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | PASS | Expression filter works in all modes: rules, AND/OR, free-text, layout persistence | | |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | PASS | All 12 steps passed. Multi-level tree filtering with layout save/restore | | |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | PASS | All 9 steps passed. Text filter with AND/OR logic and fuzzy search | | |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | PASS | All 31 steps passed. All viewer types interact correctly with Filter Panel | ~5s | 72s |
| Viewers | form | NO RUN | | | |
| Viewers | forms-tests | NO RUN | | | |
| Viewers | form-tests | NO RUN | | | |
| Viewers | grid | NO RUN | | | |
| Viewers | grid-tests | NO RUN | | | |
| Viewers | grid-viewer | NO RUN | | | |
| Viewers | heatmap | NO RUN | | | |
| Viewers | heat-map-tests | NO RUN | | | |
| Viewers | [histogram](Viewers/histogram-run.md) | PARTIAL | Property-based tests passed; canvas-based interactions not automatable | ~30s | 50s |
| Viewers | histogram-tests | NO RUN | | | |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | PASS | Color coding propagates across all viewer types; layout preserves colors | ~3s | 24s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | PARTIAL | Legends update with filters; in-viewer filter doesn't affect all viewer types | ~3s | 66s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | PASS | Legend reflects split categories; layout and project save/restore work | ~3s | 47s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | PASS | All 38 steps passed. Combined color+marker legends, filtering, layout persistence | ~3s | 50s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | PASS | SMILES renders molecule structures in legend; persists through layout save/restore | ~3s | 26s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | PARTIAL | Legends visible, positions preserved; Color Picker access ambiguous, project save failed | ~3s | 30s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | PASS | All 27 sections passed. Properties, context menu, layout, table switching all work | ~5s | 150s |
| Viewers | line-chart-tests | NO RUN | | | |
| Viewers | map | NO RUN | | | |
| Viewers | matrix-plot | NO RUN | | | |
| Viewers | matrix-plot-tests | NO RUN | | | |
| Viewers | network-diagram | NO RUN | | | |
| Viewers | pc-plot | NO RUN | | | |
| Viewers | pc-plot-tests | NO RUN | | | |
| Viewers | pc-plot-tests-pw | NO RUN | | | |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS | All 16 sections passed. All properties, layout persistence, canvas click work | ~10s | 58.2s |
| Viewers | pie-chart-tests | NO RUN | | | |
| Viewers | pivot-table | NO RUN | | | |
| Viewers | pivot-table-tests | NO RUN | | | |
| Viewers | rendering-structures-on-the-axes | NO RUN | | | |
| Viewers | row-source | NO RUN | | | |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | PASS | All 20 sections passed. Real mouse events for selection; color auto-fallback to JS API | ~5s | 66s |
| Viewers | scatter-plot-tests | NO RUN | | | |
| Viewers | statistics-viewer | NO RUN | | | |
| Viewers | tile | NO RUN | | | |
| Viewers | tree-map-viewer | NO RUN | | | |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | PARTIAL | Property tests passed; canvas interactions ambiguous; "Use in Trellis" not found | ~30s | 108s |
| Viewers | trellis-plot-tests | NO RUN | | | |
| Viewers | viewers-docking | NO RUN | | | |
| Viewers | working-with-nan-&-infinity | NO RUN | | | |
| Viewers | world-cloud | NO RUN | | | |
| WideSmokeTest | apps-and-browse | NO RUN | | | |
| WideSmokeTest | Chem/activity-cliffs | NO RUN | | | |
| WideSmokeTest | Chem/info-panels | NO RUN | | | |
| WideSmokeTest | Chem/r-group-analysis | NO RUN | | | |
| WideSmokeTest | Connections/adding | NO RUN | | | |
| WideSmokeTest | Connections/browser | NO RUN | | | |
| WideSmokeTest | Connections/delete | NO RUN | | | |
| WideSmokeTest | Connections/edit | NO RUN | | | |
| WideSmokeTest | databases | NO RUN | | | |
| WideSmokeTest | [filter-panel](WideSmokeTest/filter-panel-run.md) | PARTIAL | 7 of 10 passed. Layout restore works correctly; scaffold tree not attempted | | |
| WideSmokeTest | layouts | NO RUN | | | |
| WideSmokeTest | Notebooks/browser | NO RUN | | | |
| WideSmokeTest | Notebooks/create | NO RUN | | | |
| WideSmokeTest | Notebooks/delete | NO RUN | | | |
| WideSmokeTest | Notebooks/edit | NO RUN | | | |
| WideSmokeTest | Scripts/browser | NO RUN | | | |
| WideSmokeTest | Scripts/create | NO RUN | | | |
| WideSmokeTest | Scripts/delete | NO RUN | | | |
| WideSmokeTest | Scripts/edit | NO RUN | | | |
| WideSmokeTest | Scripts/run | NO RUN | | | |
| WideSmokeTest | TestPublic/first-login-case | NO RUN | | | |
| WideSmokeTest | uploading | NO RUN | | | |
| WideSmokeTest | Viewers/general-smoke-test | NO RUN | | | |

## Release Readiness

**Verdict**: Not ready

Run coverage is 50% (116/231) — just above the insufficient-data threshold. Six folders have zero test runs (Apps, General, LocalCashing, Notebooks, Queries, Tooltips), and the largest folder (Viewers, 66 tests) has only 36% coverage. While tested areas like DiffStudio (8/8 PASS) and Peptides (4/4 mostly PASS) look solid, the overall picture has too many gaps and failures to confirm release readiness.

### Blocking Issues
- **General** (11 tests): no runs at all
- **Queries** (14 tests): no runs at all
- **Viewers** (66 tests): only 24/66 (36%) run coverage — grid, form, heatmap, pc-plot, map, correlation-plot and many others untested
- **WideSmokeTest** (23 tests): only 1/23 (4%) run coverage
- **Tooltips** (7 tests): no runs at all
- **Notebooks** (4 tests): no runs at all
- **Apps** (2 tests): no runs at all
- **Connections**: 4 of 10 tests FAIL (infrastructure/credential issues)
- **StickyMeta**: 2 of 4 tests FAIL (entity type binding and database meta blocked)
