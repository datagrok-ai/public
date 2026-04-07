# Test Track — Global Report

**Date**: 2026-04-07
**Verdict**: Not ready

## Folder Summary

| Folder | Tests | Run | Playwright | Status |
|--------|-------|-----|------------|--------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | NO DATA |
| Bio | 9 | 9/9 (100%) | 0/9 (0%) | PARTIAL |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | FAIL |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | PARTIAL |
| Chem | 14 | 14/14 (100%) | 0/14 (0%) | FAIL |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | FAIL |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | PASS |
| EDA | 10 | 5/10 (50%) | 0/10 (0%) | PARTIAL |
| General | 11 | 0/11 (0%) | 0/11 (0%) | NO DATA |
| LocalCashing | 1 | 0/1 (0%) | 0/1 (0%) | NO DATA |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | FAIL |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | NO DATA |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL |
| PowerPack | 10 | 10/10 (100%) | 6/10 (60%) | PARTIAL |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | PARTIAL |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | NO DATA |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | PARTIAL |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | FAIL |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | NO DATA |
| Viewers | 53 | 19/53 (36%) | 19/53 (36%) | PARTIAL |
| WideSmokeTest | 23 | 1/23 (4%) | 0/23 (0%) | PARTIAL |
| **Total** | **197** | **111/197 (56%)** | **59/197 (30%)** | |

## All Tests

| Folder | Test | Status | Description |
|--------|------|--------|-------------|
| Apps | apps | NO RUN | |
| Apps | tutorials | NO RUN | |
| Bio | [analyze](Bio/analyze-run.md) | PARTIAL | All three Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three sample datasets |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | PARTIAL | Composition Analysis opens a WebLogo viewer instantly without a dialog. The gear icon correctly opens the property panel |
| Bio | [convert](Bio/convert-run.md) | PARTIAL | 3 of 5 tested steps fully passed. Bio > Calculate > Extract Region creates a region-extracted column |
| Bio | [manage](Bio/manage-run.md) | PASS | All 3 steps passed. The Manage Monomer Libraries view opens as a full view and shows 4 monomer libraries |
| Bio | [msa](Bio/msa-run.md) | PASS | All testable steps passed. The MSA dialog opens correctly with method options and Kalign version info |
| Bio | [pepsea](Bio/pepsea-run.md) | PARTIAL | MSA started on the 50-row HELM subset but did not complete within the 2-minute test window |
| Bio | [search](Bio/search-run.md) | PASS | All 4 steps passed. Subsequence Search adds a bioSubstructureFilter panel on the left side |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | PARTIAL | The scenario references "Bio > Search > Sequence Activity Cliffs" which does not exist in the current platform |
| Bio | [sequence-space](Bio/sequence-space-run.md) | PARTIAL | Sequence Space works correctly via Bio > Analyze > Sequence Space. Both runs completed successfully |
| Browse | [browse](Browse/browse-run.md) | PARTIAL | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | PARTIAL | 1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a session |
| Browse | [japanese-in-myfiles](Browse/japanese-in-myfiles-run.md) | PASS | All 3 steps passed. Japanese characters display correctly with no garbled characters |
| Browse | local-deploy | NO RUN | |
| Browse | package-manager | NO RUN | |
| Browse | [spaces](Browse/spaces-run.md) | PARTIAL | 8 steps passed, 0 failed, 24 skipped out of 32 total steps. Core CRUD operations work |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | FAIL | 1 step passed, 2 ambiguous, 27 skipped out of 30 total steps |
| Charts | [radar](Charts/radar-run.md) | PARTIAL | 2 of 3 steps fully passed. Step 3 is PARTIAL — properties need verification |
| Charts | [sunburst](Charts/sunburst-run.md) | PARTIAL | 4 of 9 steps passed, 4 skipped (SPGI_v2.csv unavailable), 1 ambiguous |
| Charts | [tree](Charts/tree-run.md) | PASS | All 4 steps plus setup passed. Tree viewer collaborative filtering works correctly |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PARTIAL | Activity Cliffs analysis completed successfully on SPGI.csv (3624 rows). Found 141 cliffs |
| Chem | [calculate](Chem/calculate-run.md) | FAIL | Descriptors function failed to open its dialog when invoked both via menu and JS API |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS | Chemical space analysis completed successfully with default parameters |
| Chem | [chemprop](Chem/chemprop-run.md) | PARTIAL | Dataset opened successfully. Training was attempted but results are partial |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PASS | All steps passed. Elemental Analysis successfully computed element counts |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PARTIAL | Filter panel opens with structure filter and all 6 filter modes available. 2 steps passed, 9 skipped |
| Chem | [info-panels](Chem/info-panels-run.md) | PARTIAL | Opened smiles.csv successfully, column-level context panel displayed all expected tabs |
| Chem | [mmp](Chem/mmp-run.md) | PASS | All 6 steps passed. MMP analysis ran successfully on 20,267-row dataset |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PARTIAL | Successfully ran R-Groups Analysis on sar_small dataset. MCS calculation found the common core |
| Chem | [sketcher](Chem/sketcher-run.md) | PARTIAL | Sketcher opens correctly and accepts SMILES input. 3 steps passed, 9 skipped |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | PASS | All 5 steps passed. Scaffold Tree generation and filtering work correctly |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | PARTIAL | Core scaffold tree functionality verified: adding, filtering, and editing scaffolds |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | PASS | All 4 steps passed. Similarity Search viewer displays correctly |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | PARTIAL | Filter panel opens correctly with all structure columns detected. 2 steps passed, 5 skipped |
| Connections | [adding](Connections/adding-run.md) | PARTIAL | 6 of 7 steps fully passed, Step 5 was partial |
| Connections | [browser](Connections/browser-run.md) | PARTIAL | 7 of 9 steps passed, 2 ambiguous |
| Connections | [catalogs](Connections/catalogs-run.md) | FAIL | NorthwindTest MS SQL connection is not present on public.datagrok.ai |
| Connections | [delete](Connections/delete-run.md) | PASS | All 8 steps passed. Both connections were deleted successfully |
| Connections | [edit](Connections/edit-run.md) | PASS | 6 of 7 steps passed (1 skipped). Connection rename and credential modification work |
| Connections | [external-provider](Connections/external-provider-run.md) | FAIL | All 7 steps skipped. Requires specific Postgres connection with superuser credentials |
| Connections | [identifiers](Connections/identifiers-run.md) | FAIL | Depends on a working Postgres Northwind connection. The connection fails to establish |
| Connections | [import-swagger](Connections/import-swagger-run.md) | FAIL | All 7 steps skipped. Requires manual drag-drop interaction |
| Connections | [schema](Connections/schema-run.md) | PARTIAL | 3 of 4 steps passed, 1 ambiguous |
| Connections | [sparql](Connections/sparql-run.md) | PARTIAL | 6 of 7 steps passed (1 failed) |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | PASS | Full DiffStudio workflow works end-to-end |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS | All 4 steps passed. PK-PD cyclic model loaded correctly |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS | All 3 functional steps passed |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | PASS | All 6 steps passed. Fit view opened correctly |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS | All 6 steps passed. Diff Studio opened correctly |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | PASS | All 5 steps passed |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS | All 4 steps passed |
| DiffStudio | [stages](DiffStudio/stages-run.md) | PASS | All 4 steps passed |
| EDA | [anova](EDA/anova-run.md) | PASS | All steps passed. ANOVA dialog opens with good defaults |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PASS | All steps passed. MVA/PLS dialog opens correctly |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL | Pareto Front viewer renders correctly but some selection features need verification |
| EDA | [pca](EDA/pca-run.md) | PASS | All steps passed. PCA runs correctly |
| EDA | [pls](EDA/pls-run.md) | PARTIAL | PLS dialog opens correctly. Function executes successfully but some results are partial |
| EDA | ML methods/linear-regression | NO RUN | |
| EDA | ML methods/pls-regression | NO RUN | |
| EDA | ML methods/softmax | NO RUN | |
| EDA | ML methods/xgboost1 | NO RUN | |
| EDA | ML methods/xgboost2 | NO RUN | |
| General | api-samples | NO RUN | |
| General | files-cache | NO RUN | |
| General | first-login | NO RUN | |
| General | inactivity-response | NO RUN | |
| General | login | NO RUN | |
| General | molecule-in-exported-csv | NO RUN | |
| General | network | NO RUN | |
| General | profile-settings | NO RUN | |
| General | startup-time | NO RUN | |
| General | table-manager | NO RUN | |
| General | tabs-reordering | NO RUN | |
| LocalCashing | local-cashing | NO RUN | |
| Models | [apply](Models/apply-run.md) | PASS | All 4 steps passed. Apply Model workflow functions correctly |
| Models | [browser](Models/browser-run.md) | PARTIAL | 4 of 6 steps passed; 2 skipped |
| Models | [chemprop](Models/chemprop-run.md) | FAIL | Chemprop Docker container is not available |
| Models | [delete](Models/delete-run.md) | PASS | All 5 steps passed |
| Models | [predictive-models](Models/predictive-models-run.md) | PASS | All 20 sub-steps passed |
| Models | [train](Models/train-run.md) | PASS | All 10 steps passed |
| Notebooks | browser | NO RUN | |
| Notebooks | create | NO RUN | |
| Notebooks | delete | NO RUN | |
| Notebooks | edit | NO RUN | |
| Peptides | [info-panels](Peptides/info-panels-run.md) | PASS | All 6 steps passed |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | PASS | All 5 steps passed |
| Peptides | [peptides](Peptides/peptides-run.md) | PASS | All 6 steps passed |
| Peptides | [sar](Peptides/sar-run.md) | PARTIAL | 12 of 13 steps passed (1 ambiguous) |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | PASS | All 5 steps passed |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | SKIP | Requires Northwind database connection |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | SKIP | No executable test steps defined |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | PARTIAL | Core Add New Column functionality works correctly |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | PASS | All 8 steps passed |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | PASS | All 5 steps passed |
| PowerPack | [AddNewColumn/functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | SKIP | Requires spgi.csv dataset not available |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | PASS | All 4 steps passed |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | PASS | All 4 steps passed |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | SKIP | Requires spgi.csv dataset not available |
| Projects | [browser](Projects/browser-run.md) | PARTIAL | 5 of 9 steps passed, 3 skipped, 1 ambiguous |
| Projects | [complex](Projects/complex-run.md) | SKIP | All 13 steps skipped — requires tables from 7+ sources |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | SKIP | All 5 steps skipped — requires custom JS script with Data Sync |
| Projects | [deleting](Projects/deleting-run.md) | PARTIAL | 2 of 4 steps passed, 1 skipped, 1 ambiguous |
| Projects | [opening](Projects/opening-run.md) | PARTIAL | All 5 steps passed but context panel details partial |
| Projects | [project-url](Projects/project-url-run.md) | SKIP | Depends on copy_clone which was not fully executed |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | SKIP | 2 of 5 steps passed, 3 skipped |
| Projects | [uploading](Projects/uploading-run.md) | PARTIAL | 8 of 14 steps passed, 6 skipped |
| Queries | adding | NO RUN | |
| Queries | browse-&-save-project | NO RUN | |
| Queries | browser | NO RUN | |
| Queries | columns-inspect | NO RUN | |
| Queries | deleting | NO RUN | |
| Queries | edit | NO RUN | |
| Queries | get-all-get-top-100 | NO RUN | |
| Queries | ms-sql | NO RUN | |
| Queries | new-sql-query | NO RUN | |
| Queries | new-visual-query | NO RUN | |
| Queries | query-layout | NO RUN | |
| Queries | query-postprocessing | NO RUN | |
| Queries | transformations | NO RUN | |
| Queries | visual-query-advanced | NO RUN | |
| Scripts | [browser](Scripts/browser-run.md) | PASS | Scripts Browser scenario passed well |
| Scripts | [create](Scripts/create-run.md) | PARTIAL | Script created and saved, but some parameters had issues |
| Scripts | [delete](Scripts/delete-run.md) | PASS | All 5 steps passed |
| Scripts | [edit](Scripts/edit-run.md) | PASS | All 6 steps passed |
| Scripts | layout | NO RUN | |
| Scripts | [run](Scripts/run-run.md) | PARTIAL | Core run functionality works but some execution paths are partial |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | PARTIAL | Sticky meta section appears in Context Panel. Two schemas available |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | FAIL | All steps skipped due to dependency on scenario 2 |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | PASS | Entity type and schema created successfully |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | FAIL | "Database meta" section not found in Context Panel |
| Tooltips | actions-in-the-context-menu | NO RUN | |
| Tooltips | default-tooltip | NO RUN | |
| Tooltips | default-tooltip-visibility | NO RUN | |
| Tooltips | edit-tooltip | NO RUN | |
| Tooltips | line-chart---aggregated-tooltip | NO RUN | |
| Tooltips | tooltip-properties | NO RUN | |
| Tooltips | uniform-default-tooltip | NO RUN | |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS | All steps passed. 3D Scatter Plot opens correctly |
| Viewers | annotation-regions | NO RUN | |
| Viewers | bar-chart | NO RUN | |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL | Most Box Plot functionality works correctly |
| Viewers | calendar | NO RUN | |
| Viewers | color-coding | NO RUN | |
| Viewers | color-coding-(linked) | NO RUN | |
| Viewers | correlation-plot | NO RUN | |
| Viewers | [density-plot](Viewers/density-plot-run.md) | PASS | All steps passed. Density plot opens correctly |
| Viewers | form | NO RUN | |
| Viewers | grid | NO RUN | |
| Viewers | grid-viewer | NO RUN | |
| Viewers | heatmap | NO RUN | |
| Viewers | histogram | NO RUN | |
| Viewers | line-chart | NO RUN | |
| Viewers | map | NO RUN | |
| Viewers | matrix-plot | NO RUN | |
| Viewers | network-diagram | NO RUN | |
| Viewers | pc-plot | NO RUN | |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS | All Pie chart functionality works correctly |
| Viewers | pivot-table | NO RUN | |
| Viewers | rendering-structures-on-the-axes | NO RUN | |
| Viewers | row-source | NO RUN | |
| Viewers | scatter-plot | NO RUN | |
| Viewers | scatter-plot-tests | NO RUN | |
| Viewers | statistics-viewer | NO RUN | |
| Viewers | tile | NO RUN | |
| Viewers | tree-map-viewer | NO RUN | |
| Viewers | trellis-plot | NO RUN | |
| Viewers | viewers-docking | NO RUN | |
| Viewers | working-with-nan-&-infinity | NO RUN | |
| Viewers | world-cloud | NO RUN | |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | PASS | All filter panel basic operations work correctly |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | PASS | All steps pass. All 6 chem search type modes produce correct row counts |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | PASS | All 14 steps passed. Cloned view correctly inherits filter state |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | PASS | All 9 steps passed |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | PASS | All 15 steps passed |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | PASS | Expression filter works correctly in all tested modes |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | PASS | All 12 steps passed |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | PASS | All 9 steps passed |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | PARTIAL | 27 PASS, 3 FAIL, 2 SKIP, 1 AMBIGUOUS out of 33 steps |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | PASS | All 7 steps passed. Color coding propagates correctly across all viewer types |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | PARTIAL | Most filtering scenarios work correctly. Legends update properly |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | PASS | All 10 steps passed |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | PASS | All 5 scenarios passed (38/38 steps) |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | PASS | All steps passed. Molecule structures render in legend |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | PARTIAL | Legend visibility and positioning features work correctly overall |
| Viewers | NX/calc-columns | NO RUN | |
| Viewers | NX/filtering | NO RUN | |
| Viewers | NX/formula-lines-and-legend | NO RUN | |
| Viewers | NX/legend-backward-compatibility | NO RUN | |
| Viewers | NX/linking | NO RUN | |
| Viewers | NX/viewers-for-linked-tables | NO RUN | |
| WideSmokeTest | apps-and-browse | NO RUN | |
| WideSmokeTest | Chem/activity-cliffs | NO RUN | |
| WideSmokeTest | Chem/info-panels | NO RUN | |
| WideSmokeTest | Chem/r-group-analysis | NO RUN | |
| WideSmokeTest | Connections/adding | NO RUN | |
| WideSmokeTest | Connections/browser | NO RUN | |
| WideSmokeTest | Connections/delete | NO RUN | |
| WideSmokeTest | Connections/edit | NO RUN | |
| WideSmokeTest | databases | NO RUN | |
| WideSmokeTest | [filter-panel](WideSmokeTest/filter-panel-run.md) | PARTIAL | 7 out of 10 steps passed, 0 confirmed failures |
| WideSmokeTest | layouts | NO RUN | |
| WideSmokeTest | Notebooks/browser | NO RUN | |
| WideSmokeTest | Notebooks/create | NO RUN | |
| WideSmokeTest | Notebooks/delete | NO RUN | |
| WideSmokeTest | Notebooks/edit | NO RUN | |
| WideSmokeTest | Scripts/browser | NO RUN | |
| WideSmokeTest | Scripts/create | NO RUN | |
| WideSmokeTest | Scripts/delete | NO RUN | |
| WideSmokeTest | Scripts/edit | NO RUN | |
| WideSmokeTest | Scripts/run | NO RUN | |
| WideSmokeTest | TestPublic/first-login-case | NO RUN | |
| WideSmokeTest | Viewers/general-smoke-test | NO RUN | |

## Release Readiness

**Verdict**: Not ready

Run coverage is at 56% (111/197), which is below the 70% threshold. Multiple folders have FAIL status (Browse, Chem, Connections, Models, StickyMeta), and 6 folders have no test data at all (Apps, General, LocalCashing, Notebooks, Queries, Tooltips). Playwright spec coverage is only 30%.

### Blocking Issues
- **Browse**: FAIL — spaces-(ui-only) scenario failed with 27/30 steps skipped
- **Chem**: FAIL — calculate (Descriptors dialog broken), plus 8 of 14 tests only PARTIAL
- **Connections**: FAIL — catalogs, external-provider, identifiers, import-swagger all failed (missing test infrastructure)
- **Models**: FAIL — chemprop Docker container not available
- **StickyMeta**: FAIL — copy-clone-delete failed (dependency chain), database-meta section not found
- **Apps, General, LocalCashing, Notebooks, Queries, Tooltips**: NO DATA — 0% run coverage
- **WideSmokeTest**: Only 1/23 tests run (4% coverage)
- **Viewers**: Only 19/53 tests run (36% coverage)
