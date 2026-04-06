# Test Track — Global Report

**Date**: 2026-04-07
**Verdict**: Not ready

## Folder Summary

| Folder | Tests | Run | Playwright | Status | Mean Spec Gen | Mean Spec Run |
|--------|-------|-----|------------|--------|---------------|---------------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | NO DATA | | |
| Bio | 9 | 9/9 (100%) | 0/9 (0%) | PARTIAL | | |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | PARTIAL | | |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | PARTIAL | | |
| Chem | 15 | 15/15 (100%) | 0/15 (0%) | PARTIAL | | |
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
| Viewers | 40 | 18/40 (45%) | 18/40 (45%) | PARTIAL | 3.6s | 44.1s |
| WideSmokeTest | 23 | 1/23 (4%) | 0/23 (0%) | PARTIAL | | |
| **Total** | **202** | **111/202 (55%)** | **58/202 (29%)** | | **3.6s** | **44.1s** |

## All Tests

| Folder | Test | Status | Description | Spec Gen | Spec Run |
|--------|------|--------|-------------|----------|----------|
| Apps | Apps | NO RUN | | | |
| Apps | Tutorials | NO RUN | | | |
| Bio | [Analyze](Bio/Analyze-run.md) | PARTIAL | All three Analyze functions work correctly on all three sample datasets; Activity Cliffs menu fails silently | | |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | PARTIAL | Composition Analysis opens a WebLogo viewer instantly without a dialog; click-to-select not verifiable via automation | | |
| Bio | [Convert](Bio/Convert-run.md) | PARTIAL | 3 of 5 tested steps passed; PolyTool Convert ran but showed no dialog and produced no new column | | |
| Bio | [Manage](Bio/Manage-run.md) | PASS | All 3 steps passed; Manage Monomer Libraries view shows 4 library JSON files with checkboxes | | |
| Bio | [MSA](Bio/MSA-run.md) | PASS | All testable steps passed; MSA dialog opens correctly, msa(Sequence) column created with aligned seqs | | |
| Bio | [Pepsea](Bio/Pepsea-run.md) | PARTIAL | MSA started on HELM subset but did not complete within 2-minute test window | | |
| Bio | [Search](Bio/Search-run.md) | PASS | All 4 steps passed; Subsequence Search correctly filters to 1 matching row | | |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | PARTIAL | Scenario references wrong menu path; Activity Cliffs works via Bio > Analyze > Activity Cliffs | | |
| Bio | [sequence-space](Bio/sequence-space-run.md) | PARTIAL | Sequence Space works correctly; scenario incorrectly specifies menu path as Bio > Search | | |
| Browse | [Browse](Browse/Browse-run.md) | PARTIAL | 4 steps passed, 1 partial; item-level selection not preserved across page refresh | | |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | PARTIAL | Tree state preserved within session but not across page refreshes | | |
| Browse | [japanese-in-myfiles](Browse/japanese-in-myfiles-run.md) | PASS | All 3 steps passed; Japanese characters display correctly throughout the UI | | |
| Browse | local-deploy | NO RUN | | | |
| Browse | package-manager | NO RUN | | | |
| Browse | [Spaces](Browse/Spaces-run.md) | PARTIAL | 8 steps passed, 24 skipped; core CRUD works, drag-and-drop not testable via automation | | |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | FAIL | 1 step passed, 27 skipped; drag-and-drop operations not feasible via browser automation | | |
| Charts | [Radar](Charts/Radar-run.md) | PARTIAL | Viewers render correctly; table switching throws NoSuchMethodError: toLowerCase | | |
| Charts | [Sunburst](Charts/Sunburst-run.md) | PARTIAL | 4 of 9 steps passed; SPGI_v2.csv unavailable in demo files | | |
| Charts | [Tree](Charts/Tree-run.md) | PASS | All 4 steps passed; collaborative filtering works correctly end-to-end | | |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PARTIAL | Analysis completed on SPGI.csv; 141 cliffs found; interactive steps skipped | | |
| Chem | [scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | PASS | All 5 steps passed; Scaffold Tree generation and filtering work correctly | | |
| Chem | [scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | PARTIAL | Core functionality verified (add, filter, edit scaffolds); complex sub-scenarios skipped | | |
| Chem | [similarity-search](Chem/Advanced/similarity-search-run.md) | PASS | All 4 steps passed; property modifications (fingerprint, limit, cutoff) all work | | |
| Chem | [structure-filter](Chem/Advanced/structure-filter-run.md) | PARTIAL | Filter panel opens with all 6 modes; sketcher-based filtering steps skipped | | |
| Chem | [Calculate](Chem/Calculate-run.md) | FAIL | Descriptors function failed to open its dialog via menu and JS API | | |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS | Chemical space analysis completed successfully with UMAP/Tanimoto fingerprints | | |
| Chem | [Chemprop](Chem/Chemprop-run.md) | PARTIAL | Dataset opened; training not attempted due to long computation time on public server | | |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PASS | All steps passed; element counts computed and radar viewer added | | |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PARTIAL | Filter panel opens with structure filter; 9 steps skipped requiring interactive UI | | |
| Chem | [info-panels](Chem/info-panels-run.md) | PARTIAL | Column-level context panel works; cell-level panel could not be triggered programmatically | | |
| Chem | [MMP](Chem/MMP-run.md) | PASS | All 6 steps passed; MMP analysis ran on 20K+ row dataset, all 4 tabs functional | | |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PARTIAL | MCS calculation and R-Groups Analysis produced correct results; repeat-run steps skipped | | |
| Chem | [Sketcher](Chem/Sketcher-run.md) | PARTIAL | Sketcher opens and accepts SMILES input; hamburger menu and clipboard steps skipped | | |
| Connections | [Adding](Connections/Adding-run.md) | PARTIAL | 6 of 7 steps passed; TEST button works but connection test fails without real credentials | | |
| Connections | [Browser](Connections/Browser-run.md) | PARTIAL | 7 of 9 steps passed; search filtering and Context Pane work; magic wand icon not found | | |
| Connections | [Catalogs](Connections/Catalogs-run.md) | FAIL | NorthwindTest MS SQL connection not present on public.datagrok.ai; all steps blocked | | |
| Connections | [Delete](Connections/Delete-run.md) | PASS | All 8 steps passed; both connections deleted successfully | | |
| Connections | [Edit](Connections/Edit-run.md) | PASS | 6 of 7 steps passed; rename, credential modification, and error verification work | | |
| Connections | [external-provider](Connections/external-provider-run.md) | FAIL | All steps skipped; requires specific Postgres connection with superuser credentials | | |
| Connections | [Identifiers](Connections/Identifiers-run.md) | FAIL | Depends on working Postgres Northwind connection; blocked by auth failure | | |
| Connections | [import-swagger](Connections/import-swagger-run.md) | FAIL | All steps skipped; requires manual YAML file download and drag-and-drop | | |
| Connections | [Schema](Connections/Schema-run.md) | PARTIAL | 3 of 4 steps passed; Browse schema context menu option not found in current UI | | |
| Connections | [Sparql](Connections/Sparql-run.md) | PARTIAL | 6 of 7 steps passed; SPARQL connection created and deleted; query execution failed | | |
| DiffStudio | [Catalog](DiffStudio/Catalog-run.md) | PASS | All 6 steps passed; full workflow: load from Library, save to Model Hub, run from catalog | | |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS | All 4 steps passed; PK-PD cyclic model loaded correctly | | |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS | All 3 steps passed; pk.ivp file opened from Files browser as DiffStudio view | | |
| DiffStudio | [Fitting](DiffStudio/Fitting-run.md) | PASS | All 6 steps passed; Fit view and Process mode lookup table work correctly | | |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS | All 6 steps passed; Bioreactor example loaded from Library | | |
| DiffStudio | [Scripting](DiffStudio/Scripting-run.md) | PASS | All 5 steps passed; Edit toggle and JS scripting button work correctly | | |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS | All 4 steps passed; Sensitivity Analysis view opened correctly from ribbon | | |
| DiffStudio | [Stages](DiffStudio/Stages-run.md) | PASS | All 4 steps passed; Acid Production 2-stage model loaded correctly | | |
| EDA | [ANOVA](EDA/ANOVA-run.md) | PASS | All steps passed; ANOVA dialog opens with good defaults and produces correct results | | |
| EDA | ML methods/linear-regression | NO RUN | | | |
| EDA | ML methods/pls-regression | NO RUN | | | |
| EDA | ML methods/Softmax | NO RUN | | | |
| EDA | ML methods/XGBoost1 | NO RUN | | | |
| EDA | ML methods/XGBoost2 | NO RUN | | | |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PASS | All steps passed; MVA/PLS dialog opens correctly with auto-populated fields | | |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL | Viewer renders correctly; auto-selects model column; some interactive steps partial | | |
| EDA | [PCA](EDA/PCA-run.md) | PASS | All steps passed; PCA runs correctly with and without preprocessing | | |
| EDA | [PLS](EDA/PLS-run.md) | PARTIAL | PLS dialog opens correctly; function executes successfully but some steps partial | | |
| General | api-samples | NO RUN | | | |
| General | files-cache | NO RUN | | | |
| General | first-login | NO RUN | | | |
| General | inactivity-response | NO RUN | | | |
| General | login | NO RUN | | | |
| General | molecule-in-exported-csv | NO RUN | | | |
| General | Network | NO RUN | | | |
| General | profile-settings | NO RUN | | | |
| General | startup-time | NO RUN | | | |
| General | table-manager | NO RUN | | | |
| General | tabs-reordering | NO RUN | | | |
| LocalCashing | local-cashing | NO RUN | | | |
| Models | [Apply](Models/Apply-run.md) | PASS | All 4 steps passed; Apply Model workflow functions correctly end-to-end | | |
| Models | [Browser](Models/Browser-run.md) | PARTIAL | 4 of 6 steps passed; only 1 model available so comparison steps skipped | | |
| Models | [Chemprop](Models/Chemprop-run.md) | FAIL | 2 of 17 sub-steps passed; Chemprop Docker container not available on public server | | |
| Models | [Delete](Models/Delete-run.md) | PASS | All 5 steps passed; model deletion workflow works correctly end-to-end | | |
| Models | [predictive-models](Models/predictive-models-run.md) | PASS | All 20 sub-steps passed; full lifecycle Train/Apply/Delete for EDA-based models | | |
| Models | [Train](Models/Train-run.md) | PASS | All 10 steps passed; Train Model workflow functions correctly on public server | | |
| Notebooks | Browser | NO RUN | | | |
| Notebooks | Create | NO RUN | | | |
| Notebooks | Delete | NO RUN | | | |
| Notebooks | Edit | NO RUN | | | |
| Peptides | [info-panels](Peptides/info-panels-run.md) | PASS | All 6 steps passed; peptides dataset opened with color-coded amino acids | | |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | PASS | All 5 steps passed; SAR launched successfully via Bio > Analyze > SAR | | |
| Peptides | [Peptides](Peptides/Peptides-run.md) | PASS | All 6 steps passed; context panel shows weblogo with Activity/Scaling/Clusters | | |
| Peptides | [SAR](Peptides/SAR-run.md) | PARTIAL | 12 of 13 steps passed; full SAR workflow worked with all 4 viewers | | |
| PowerPack | [add-new-column (ANC)](PowerPack/AddNewColumn/add-new-column-run.md) | PARTIAL | Core Add New Column works; chaining formulas and propagation work correctly | | |
| PowerPack | [autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | PASS | All 8 steps passed; autocomplete works for functions and columns | | |
| PowerPack | [formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | PASS | All 5 steps passed; formula dependency chain works correctly | | |
| PowerPack | [functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | SKIP | All steps skipped; requires spgi.csv dataset not available on release server | | |
| PowerPack | [highlight](PowerPack/AddNewColumn/highlight-run.md) | PASS | All 4 steps passed; column name highlighting works in both syntaxes | | |
| PowerPack | [hints](PowerPack/AddNewColumn/hints-run.md) | PASS | All 4 steps passed; function tooltip displays signature and description | | |
| PowerPack | [input_functions](PowerPack/AddNewColumn/input_functions-run.md) | SKIP | All steps skipped; requires spgi.csv dataset not available on release server | | |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | PASS | All 5 steps passed; dialog works with formula entry and autocomplete | | |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | SKIP | All steps skipped; requires Northwind database connection | | |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | SKIP | Scenario contains only GitHub ticket references; no executable test steps | | |
| Projects | [Browser](Projects/Browser-run.md) | PARTIAL | 5 of 9 steps passed; projects listed, searchable, and openable | | |
| Projects | [Complex](Projects/Complex-run.md) | SKIP | All 13 steps skipped; requires tables from 7+ sources and drag-and-drop | | |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | SKIP | All 5 steps skipped; requires custom JS script with Data Sync | | |
| Projects | [Deleting](Projects/Deleting-run.md) | PARTIAL | 2 of 4 steps passed; deletion works via API but UI confirmation dialog ambiguous | | |
| Projects | [Opening](Projects/Opening-run.md) | PARTIAL | All 5 steps passed; Context Panel correctly shows project details | | |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | SKIP | 2 of 5 steps passed; copy/clone/link operations not tested | | |
| Projects | [project-url](Projects/project-url-run.md) | SKIP | All steps skipped; depends on copy_clone scenario not fully executed | | |
| Projects | [Uploading](Projects/Uploading-run.md) | PARTIAL | 8 of 14 steps passed; core project creation works from multiple sources | | |
| Queries | Adding | NO RUN | | | |
| Queries | Browser | NO RUN | | | |
| Queries | browse-&-save-project | NO RUN | | | |
| Queries | columns-inspect | NO RUN | | | |
| Queries | Deleting | NO RUN | | | |
| Queries | Edit | NO RUN | | | |
| Queries | get-all-get-top-100 | NO RUN | | | |
| Queries | ms-sql | NO RUN | | | |
| Queries | new-sql-query | NO RUN | | | |
| Queries | new-visual-query | NO RUN | | | |
| Queries | query-layout | NO RUN | | | |
| Queries | query-postprocessing | NO RUN | | | |
| Queries | Transformations | NO RUN | | | |
| Queries | visual-query-advanced | NO RUN | | | |
| Scripts | [Browser](Scripts/Browser-run.md) | PASS | Context pane shows all expected accordions (Details, Script, Run, Activity, Chats, Dev) | | |
| Scripts | [Create](Scripts/Create-run.md) | PARTIAL | Script created and saved successfully; some parameter configuration steps partial | | |
| Scripts | [Delete](Scripts/Delete-run.md) | PASS | All 5 steps passed; delete flow works with confirmation dialog | | |
| Scripts | [Edit](Scripts/Edit-run.md) | PASS | All 6 steps passed; edits are saved persistently and visible on re-open | | |
| Scripts | Layout | NO RUN | | | |
| Scripts | [Run](Scripts/Run-run.md) | PARTIAL | Core run functionality works from context menu and console | | |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | PARTIAL | Sticky meta section appears in Context Panel; two schemas accessible | | |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | FAIL | All steps skipped due to dependency on Add and edit scenario | | |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | PASS | Entity type and schema created successfully with all field types | | |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | FAIL | Database meta section not found in Context Panel for CHEMBL connection | | |
| Tooltips | actions-in-the-context-menu | NO RUN | | | |
| Tooltips | default-tooltip | NO RUN | | | |
| Tooltips | default-tooltip-visibility | NO RUN | | | |
| Tooltips | edit-tooltip | NO RUN | | | |
| Tooltips | line-chart---aggregated-tooltip | NO RUN | | | |
| Tooltips | tooltip-properties | NO RUN | | | |
| Tooltips | uniform-default-tooltip | NO RUN | | | |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS | 3D Scatter Plot opens correctly; data point clicking selects corresponding row | 3s | 10s |
| Viewers | annotation-regions | NO RUN | | | |
| Viewers | bar-chart | NO RUN | | | |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL | Most functionality works; data properties and filter work correctly | 5s | 104s |
| Viewers | Calendar | NO RUN | | | |
| Viewers | color-coding-(linked) | NO RUN | | | |
| Viewers | color-coding | NO RUN | | | |
| Viewers | correlation-plot | NO RUN | | | |
| Viewers | [density-plot](Viewers/density-plot-run.md) | PASS | Density plot opens correctly showing hexagonal bins; axes changeable | 3s | 8s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | PASS | All filter panel basic operations work correctly on dev server | | |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | PASS | All 6 chem search type modes produce correct row counts | | |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | PASS | All 14 steps passed; cloned view inherits filter state correctly | | |
| Viewers | [FilterPanel/collaborative-filtering](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | PASS | All 9 steps passed; SELECTION_TO_FILTER works correctly for linked tables | | |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | PASS | All 15 steps passed; combined boolean filter supports OR mode correctly | | |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | PASS | Expression filter works in all tested modes; AND/OR switching correct | | |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | PASS | All 12 steps passed; multi-level tree filtering works correctly | | |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | PASS | All 9 steps passed; text filter searches, highlights, and supports multiple terms | | |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | PARTIAL | 27 of 33 steps passed, 3 failed; core filters work but some advanced cases fail | | |
| Viewers | Form | NO RUN | | | |
| Viewers | Grid | NO RUN | | | |
| Viewers | grid-viewer | NO RUN | | | |
| Viewers | Heatmap | NO RUN | | | |
| Viewers | Histogram | NO RUN | | | |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | PASS | All 7 steps passed; categorical coloring propagates across all viewer types | 3s | 24s |
| Viewers | [Legend/Filtering](Viewers/Legend/Filtering-run.md) | PARTIAL | Most filtering scenarios work; legends update properly when filters applied | | |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | PASS | All 10 steps passed; legend correctly reflects split categories | 3s | 47s |
| Viewers | [Legend/Scatterplot](Viewers/Legend/Scatterplot-run.md) | PASS | All 5 scenarios passed (38/38 steps); combined color+marker legends work | 3s | 50s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | PASS | Core column renders molecule structures in scatterplot marker legend | 3s | 26s |
| Viewers | Legend/visibility-and-positioning | NO RUN | | | |
| Viewers | line-chart | NO RUN | | | |
| Viewers | Map | NO RUN | | | |
| Viewers | matrix-plot | NO RUN | | | |
| Viewers | network-diagram | NO RUN | | | |
| Viewers | NX/calc-columns | NO RUN | | | |
| Viewers | NX/Filtering | NO RUN | | | |
| Viewers | NX/formula-lines-and-legend | NO RUN | | | |
| Viewers | NX/Linking | NO RUN | | | |
| Viewers | NX/viewers-for-linked-tables | NO RUN | | | |
| Viewers | pc-plot | NO RUN | | | |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS | All Pie chart functionality works correctly; aggregation, selection, legend | 5s | 84s |
| Viewers | pivot-table | NO RUN | | | |
| Viewers | rendering-structures-on-the-axes | NO RUN | | | |
| Viewers | row-source | NO RUN | | | |
| Viewers | scatter-plot | NO RUN | | | |
| Viewers | scatter-plot-tests | NO RUN | | | |
| Viewers | statistics-viewer | NO RUN | | | |
| Viewers | Tile | NO RUN | | | |
| Viewers | tree-map-viewer | NO RUN | | | |
| Viewers | trellis-plot | NO RUN | | | |
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
| WideSmokeTest | [filter-panel](WideSmokeTest/filter-panel-run.md) | PARTIAL | 7 out of 10 steps passed; scaffold tree step not attempted | | |
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

Run coverage is 55% (111/202), which is below the 70% threshold. Six entire folders have no run data at all (Apps, General, LocalCashing, Notebooks, Queries, Tooltips -- 46 tests total). Among the 111 executed runs, only 1 folder (DiffStudio) achieved a clean PASS across all its scenarios. Multiple folders contain FAIL results: Connections (4 FAILs due to missing server configs), Chem (1 FAIL -- Calculate/Descriptors), Models (1 FAIL -- Chemprop Docker unavailable), StickyMeta (2 FAILs), and Browse (1 FAIL).

### Blocking Issues
- **Apps, General, LocalCashing, Notebooks, Queries, Tooltips**: Zero run coverage -- 46 scenarios never executed
- **Connections**: 4 of 10 runs FAIL (Catalogs, external-provider, Identifiers, import-swagger) -- all require server configurations or credentials not available on public.datagrok.ai
- **Chem/Calculate**: Descriptors dialog fails to open -- possible server-side computation issue
- **Models/Chemprop**: Docker container not available on public server
- **StickyMeta**: copy-clone-delete blocked by add-and-edit dependency; database-meta section missing for CHEMBL connection
- **Browse/spaces-(ui-only)**: FAIL -- nearly all steps require drag-and-drop not feasible via automation
- **Viewers**: Only 18 of 40 scenarios executed (45% coverage); 22 viewers untested
- **WideSmokeTest**: Only 1 of 23 scenarios executed (4% coverage)
