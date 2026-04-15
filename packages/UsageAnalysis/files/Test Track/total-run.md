# Test Track — Global Report

**Date**: 2026-04-15  
**Verdict**: Conditionally ready

## Folder Summary

**Total**: 186 tests · Run: 138/186 (74%) · Playwright: 105/186 (56%) · Mean Spec Gen: 8.8s · Mean Spec Run: 31.6s

| Folder | Tests | Run | Playwright | Status | Mean Spec Gen | Mean Spec Run |
|--------|-------|-----|------------|--------|---------------|---------------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | NO DATA |  |  |
| Bio | 9 | 9/9 (100%) | 1/9 (11%) | PARTIAL | 5.0s |  |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | PARTIAL |  |  |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | PARTIAL | 3.0s | 24.2s |
| Chem | 14 | 14/14 (100%) | 12/14 (86%) | PARTIAL | 3.0s | 32.2s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | PARTIAL |  |  |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | PASS | 3.5s | 26.7s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | PARTIAL | 2.0s | 4.8s |
| General | 11 | 0/11 (0%) | 0/11 (0%) | NO DATA |  |  |
| LocalCashing | 1 | 0/1 (0%) | 0/1 (0%) | NO DATA |  |  |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | PARTIAL |  |  |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | NO DATA |  |  |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 3.0s | 47.6s |
| PowerPack | 10 | 10/10 (100%) | 6/10 (60%) | PARTIAL | 2.5s | 5.3s |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | PARTIAL |  |  |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | NO DATA |  |  |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | PARTIAL |  |  |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 3.0s | 15.0s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | NO DATA |  |  |
| Viewers | 48 | 42/48 (88%) | 42/48 (88%) | PARTIAL | 16.5s | 41.6s |

## All Tests

| Folder | Test | Status | Description | Spec Gen | Spec Run |
|--------|------|--------|-------------|----------|----------|
| Apps | apps | NO RUN |  |  |  |
| Apps | tutorials | NO RUN |  |  |  |
| Bio | [analyze](Bio/analyze-run.md) | PASS | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset ... | 5.0s |  |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | PARTIAL | 4 of 5 steps passed. The Composition/WebLogo viewer opens correctly and properties are accessible. The letter-click s... |  |  |
| Bio | [convert](Bio/convert-run.md) | PASS | All 4 Bio convert/transform functions work correctly on FASTA data. Get Region extracts a subsequence region. PolyToo... |  |  |
| Bio | [manage](Bio/manage-run.md) | PASS | All 3 steps passed. The Manage Monomer Libraries view opens as a full view showing 5 monomer library JSON files (incr... |  |  |
| Bio | [msa](Bio/msa-run.md) | PASS | All 6 steps passed. MSA dialog opens with correct fields, Alignment Parameters button adds gap penalty inputs as expe... |  |  |
| Bio | [pepsea](Bio/pepsea-run.md) | AMBIGUOUS | The MSA dialog opens correctly for HELM data and shows MAFFT-based method options (mafft --auto, linsi, ginsi, etc.) ... |  |  |
| Bio | [search](Bio/search-run.md) | PASS | All 4 steps passed. Bio > Search > Subsequence Search opens a filter panel with a Sequence bio substructure filter. T... |  |  |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | PASS | All 6 steps passed. Activity Cliffs works correctly with both default parameters (UMAP/Hamming) and custom parameters... |  |  |
| Bio | [sequence-space](Bio/sequence-space-run.md) | PASS | All 6 steps passed. Sequence Space works correctly with both default (UMAP/Hamming) and custom (t-SNE/Needlemann-Wuns... |  |  |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | PARTIAL | 1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a single sess... |  |  |
| Browse | [browse](Browse/browse-run.md) | PARTIAL | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. I... |  |  |
| Browse | [japanese-in-myfiles](Browse/japanese-in-myfiles-run.md) | PASS | All 3 steps passed. The file 芹沢 貴之 こんにちは.csv displays correctly with no garbled or incorrect characters. Japanese Kan... |  |  |
| Browse | local-deploy | NO RUN |  |  |  |
| Browse | package-manager | NO RUN |  |  |  |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | FAIL | 1 step passed, 2 ambiguous, 27 skipped out of 30 total steps. This scenario focuses on UI-only interactions not cover... |  |  |
| Browse | [spaces](Browse/spaces-run.md) | PARTIAL | 8 steps passed, 0 failed, 24 skipped out of 32 total steps. Core CRUD operations (create root/child/nested, rename, a... |  |  |
| Charts | [radar](Charts/radar-run.md) | PASS | All steps of the Radar viewer scenario passed both in the MCP browser run and in the Playwright spec execution. The s... | 3.0s | 14.5s |
| Charts | [sunburst](Charts/sunburst-run.md) | PARTIAL | 9 of 9 steps attempted: 8 passed, 1 ambiguous (canvas multi-selection), 1 skipped (old layout). Playwright spec passe... | 3.0s | 41.5s |
| Charts | [tree](Charts/tree-run.md) | PASS | All 4 steps of the Tree viewer collaborative filtering scenario passed. The Tree viewer correctly maintains selection... | 3.0s | 16.6s |
| Chem | [scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | FAIL | Steps 1-2 and the toolbar check passed. The scaffold tree viewer opens correctly with the "empty" state message. Howe... | 3.0s |  |
| Chem | [scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | FAIL | Steps 1-2 passed: SPGI.csv loaded and the Scaffold Tree viewer was added with the "empty" state. However, adding scaf... |  |  |
| Chem | [similarity-search](Chem/Advanced/similarity-search-run.md) | PASS | All 4 steps passed. The Similarity Search viewer opened correctly showing similar molecules with Tanimoto/Morgan simi... | 3.0s |  |
| Chem | [structure-filter](Chem/Advanced/structure-filter-run.md) | PASS | All 4 core steps passed. The Structure filter panel opened with 43 filters for SPGI.csv. Benzene substructure search ... | 3.0s | 24.0s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PASS | Activity Cliffs computation works correctly on SPGI.csv. UMAP embedding produces a scatter plot with molecule thumbna... | 3.0s | 56.9s |
| Chem | [calculate](Chem/calculate-run.md) | PASS | Descriptors calculation works correctly on dev server. Dialog shows all descriptor categories, pre-selected descripto... | 3.0s | 24.5s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS | Chemical Space works correctly with both default (UMAP) and modified (t-SNE) parameters. Scatter plots are generated ... | 3.0s | 44.0s |
| Chem | [chemprop](Chem/chemprop-run.md) | PARTIAL | Steps 1-3 passed: mol1K.sdf loaded correctly, Predictive Model view opened, and columns were configured. However, mod... | 3.0s |  |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PASS | All 4 steps passed. Elemental Analysis computed element counts (H, C, N, O, F, S, Cl, Br, I) and Molecule Charge for ... | 3.0s | 21.0s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PARTIAL | Core substructure filtering works: drawing benzene in the structure filter correctly filters to 1356/3624 rows with h... | 3.0s | 29.0s |
| Chem | [info-panels](Chem/info-panels-run.md) | PASS | All 5 steps passed. The canonical_smiles column Context Panel shows 11 info panels including Chemistry, all expanding... | 3.0s | 39.0s |
| Chem | [mmp](Chem/mmp-run.md) | FAIL | Steps 1-2 passed (dataset opened, MMP dialog opened). Step 3 failed: the MMP analysis function crashed with a `TypeEr... | 3.0s |  |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PASS | All tested steps passed. The R-Groups Analysis dialog opened correctly, MCS computation completed, and the trellis pl... | 3.0s | 38.0s |
| Chem | [sketcher](Chem/sketcher-run.md) | PASS | Core steps 1-6 passed: smiles.csv opened, sketcher dialog displayed with OpenChemLib, C1CCCCC1 entered and rendered, ... | 3.0s | 13.0s |
| Connections | [adding](Connections/adding-run.md) | PARTIAL | 6 of 7 steps fully passed, Step 5 was partial (TEST button works but actual connection test fails without real creden... |  |  |
| Connections | [browser](Connections/browser-run.md) | PARTIAL | 7 of 9 steps passed, 2 ambiguous. Search filtering works correctly and the Context Pane shows all expected tabs (Deta... |  |  |
| Connections | [catalogs](Connections/catalogs-run.md) | FAIL | 1 step passed, 1 failed, 15 skipped. The required `NorthwindTest` MS SQL connection is not present on public.datagrok... |  |  |
| Connections | [delete](Connections/delete-run.md) | PASS | All 8 steps passed. Both connections were deleted successfully. The confirmation dialog uses a red "DELETE" button (n... |  |  |
| Connections | [edit](Connections/edit-run.md) | PASS | 6 of 7 steps passed (1 skipped due to missing real credentials). The connection rename, credential modification, and ... |  |  |
| Connections | [external-provider](Connections/external-provider-run.md) | FAIL | All 7 steps skipped. This scenario requires a specific Postgres connection at db.datagrok.ai:54327 with superuser cre... |  |  |
| Connections | [identifiers](Connections/identifiers-run.md) | FAIL | 1 step passed, 1 failed, 7 skipped. This scenario depends on a working Postgres connection to the Northwind database.... |  |  |
| Connections | [import-swagger](Connections/import-swagger-run.md) | FAIL | All 7 steps skipped. This scenario requires manual interaction: downloading a YAML file to the local machine and drag... |  |  |
| Connections | [schema](Connections/schema-run.md) | PARTIAL | 3 of 4 steps passed, 1 ambiguous. The "Browse schema" context menu option was not found in the current UI, but the sc... |  |  |
| Connections | [sparql](Connections/sparql-run.md) | PARTIAL | 6 of 7 steps passed (1 failed). All UI steps worked correctly. The SPARQL connection was created and deleted successf... |  |  |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | PASS | All 6 steps passed in both MCP and Playwright runs. The full DiffStudio catalog workflow works end-to-end on dev.data... | 5.0s | 32.1s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS | All 4 steps passed in both MCP and Playwright runs. The PK-PD cyclic model loads correctly, Multiaxis/Facet plot tabs... | 3.0s | 22.7s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS | All 3 functional steps passed. The pk.ivp file opened from the Files browser as a DiffStudio view. Modifying Step (0.... | 3.0s | 24.0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | PASS | All 6 steps passed. The Bioreactor model loaded from Library, Fit view opened correctly, Process mode dropdown update... | 5.0s | 28.0s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS | All 6 steps passed. DiffStudio opened correctly, Bioreactor loaded from Library with 13 columns and 1001 rows. Both M... | 3.0s | 18.0s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | PASS | All 5 steps passed. The DiffStudio Edit toggle and `</>` export icon work correctly, producing a JS script with all m... | 3.0s | 49.0s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS | All 4 steps passed. The Sensitivity Analysis view opened correctly from the ribbon. Process mode switching updated FF... | 3.0s | 23.0s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | PASS | All 4 steps passed. The Acid Production (GA-production) model loaded from Library with 2-stage simulation (1st stage ... | 3.0s | 17.0s |
| EDA | [linear-regression](EDA/ML methods/linear-regression-run.md) | PASS | Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation... | 2.0s | 6.8s |
| EDA | [pls-regression](EDA/ML methods/pls-regression-run.md) | PARTIAL | PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-... | 2.0s | 6.9s |
| EDA | [softmax](EDA/ML methods/softmax-run.md) | FAIL | Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns an... | 2.0s | 2.6s |
| EDA | [xgboost1](EDA/ML methods/xgboost1-run.md) | PARTIAL | XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as... | 2.0s | 2.7s |
| EDA | [xgboost2](EDA/ML methods/xgboost2-run.md) | PARTIAL | XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interac... | 2.0s | 2.7s |
| EDA | [anova](EDA/anova-run.md) | PASS | All 3 steps passed successfully. The ANOVA dialog opened correctly with default parameters (Category=RACE, Feature=AG... | 2.0s | 3.0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PASS | All 3 steps passed. The Multivariate Analysis (PLS) dialog opened with correct defaults. After clicking RUN, all 5 ex... | 2.0s | 7.0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL | 6 of 7 steps passed (1 partial due to missing test file). The Pareto Front viewer works correctly: non-numeric column... | 2.0s | 9.0s |
| EDA | [pca](EDA/pca-run.md) | PASS | All 5 steps passed. PCA runs correctly twice — first without preprocessing (raw scores), then with Center+Scale (norm... | 2.0s | 3.0s |
| EDA | [pls](EDA/pls-run.md) | PASS | All 4 steps passed. PLS dialog opened with correct defaults (Predict=price, Using=15 numeric columns, Components=3). ... | 2.0s | 4.0s |
| General | api-samples | NO RUN |  |  |  |
| General | files-cache | NO RUN |  |  |  |
| General | first-login | NO RUN |  |  |  |
| General | inactivity-response | NO RUN |  |  |  |
| General | login | NO RUN |  |  |  |
| General | molecule-in-exported-csv | NO RUN |  |  |  |
| General | network | NO RUN |  |  |  |
| General | profile-settings | NO RUN |  |  |  |
| General | startup-time | NO RUN |  |  |  |
| General | table-manager | NO RUN |  |  |  |
| General | tabs-reordering | NO RUN |  |  |  |
| LocalCashing | local-cashing | NO RUN |  |  |  |
| Models | [apply](Models/apply-run.md) | PASS | All 4 steps passed. The Apply Model workflow functions correctly end-to-end. The TestDemog model (trained in Train.md... |  |  |
| Models | [browser](Models/browser-run.md) | PARTIAL | 4 of 6 steps passed; 2 skipped because only 1 model was available (TestDemog was the only model — the second numeric ... |  |  |
| Models | [chemprop](Models/chemprop-run.md) | FAIL | 2 of 17 sub-steps passed, 13 skipped, 1 failed, 1 ambiguous. The scenario fails entirely due to the Chemprop Docker c... |  |  |
| Models | [delete](Models/delete-run.md) | PASS | All 5 steps passed. The model deletion workflow works correctly end-to-end. The right-click context menu, confirmatio... |  |  |
| Models | [predictive-models](Models/predictive-models-run.md) | PASS | All 20 sub-steps passed. The full lifecycle (Train → Apply → Apply on new dataset → Delete) for EDA-based predictive ... |  |  |
| Models | [train](Models/train-run.md) | PASS | All 10 steps passed. The Train Model workflow functions correctly on public.datagrok.ai using the built-in EDA engine... |  |  |
| Notebooks | browser | NO RUN |  |  |  |
| Notebooks | create | NO RUN |  |  |  |
| Notebooks | delete | NO RUN |  |  |  |
| Notebooks | edit | NO RUN |  |  |  |
| Peptides | [info-panels](Peptides/info-panels-run.md) | PASS | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are re... | 3.0s | 10.9s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | PARTIAL | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variabil... | 3.0s | 1m 18s |
| Peptides | [peptides](Peptides/peptides-run.md) | PARTIAL | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clust... | 3.0s | 11.6s |
| Peptides | [sar](Peptides/sar-run.md) | PARTIAL | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Res... | 3.0s | 1m 30s |
| PowerPack | [add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | PARTIAL | Core Add New Column functionality works correctly: creating calculated columns, chaining formulas, and propagating co... |  |  |
| PowerPack | [autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | PASS | All 8 steps passed. Autocomplete works for both functions (typing or Ctrl+Space) and columns ($ symbol). Function sel... |  |  |
| PowerPack | [formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | PASS | All 5 steps passed. The formula dependency chain (Weight2 -> Weight3 -> Weight4) works correctly. Changing a value in... |  |  |
| PowerPack | [functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | SKIP | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on rel... |  |  |
| PowerPack | [highlight](PowerPack/AddNewColumn/highlight-run.md) | PASS | All 4 steps passed. Column names in both ${} (scalar) and $[] (aggregate) syntax are highlighted with the cm-column-n... |  |  |
| PowerPack | [hints](PowerPack/AddNewColumn/hints-run.md) | PASS | All 4 steps passed. Hovering over a function name in the formula editor displays a tooltip with the function signatur... |  |  |
| PowerPack | [input_functions](PowerPack/AddNewColumn/input_functions-run.md) | SKIP | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on rel... |  |  |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | PASS | All 5 steps passed on dev.datagrok.ai. The Add New Column dialog opens correctly, displays cleanly with no visual iss... | 2.0s | 5.3s |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | PARTIAL | Part 1 (create, save, apply enrichment) passed. The enrichment feature correctly joins the Customers table to Orders ... | 3.0s |  |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | SKIP | This scenario file contains only references to GitHub tickets related to formula lines. No executable test steps were... |  |  |
| Projects | [browser](Projects/browser-run.md) | PARTIAL | 5 of 9 steps passed, 3 skipped, 1 ambiguous. Browse > Dashboards view works correctly: projects are listed, searchabl... |  |  |
| Projects | [complex](Projects/complex-run.md) | SKIP | All 13 steps skipped. This is the most complex scenario requiring tables from 7+ different sources, drag-and-drop, en... |  |  |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | SKIP | All 5 steps skipped. This scenario requires running a custom JavaScript script with Data Sync enabled, then modifying... |  |  |
| Projects | [deleting](Projects/deleting-run.md) | PARTIAL | 2 of 4 steps passed, 1 skipped, 1 ambiguous. Project deletion works via the API (`grok.dapi.projects.delete()`). The ... |  |  |
| Projects | [opening](Projects/opening-run.md) | PARTIAL | All 5 steps passed. Projects from the Uploading step are accessible in Browse > Dashboards. Context Panel correctly s... |  |  |
| Projects | [project-url](Projects/project-url-run.md) | SKIP | All steps skipped. This scenario depends on Projects copy_clone.md (order 5) which was not fully executed. The Link/C... |  |  |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | SKIP | 2 of 5 steps passed, 3 skipped. Project preview and opening work. Copy/clone/link operations were not tested because ... |  |  |
| Projects | [uploading](Projects/uploading-run.md) | PARTIAL | 8 of 14 steps passed, 6 skipped. Core project creation from local tables, file shares, query results, and join result... |  |  |
| Queries | adding | NO RUN |  |  |  |
| Queries | browse-&-save-project | NO RUN |  |  |  |
| Queries | browser | NO RUN |  |  |  |
| Queries | columns-inspect | NO RUN |  |  |  |
| Queries | deleting | NO RUN |  |  |  |
| Queries | edit | NO RUN |  |  |  |
| Queries | get-all-get-top-100 | NO RUN |  |  |  |
| Queries | ms-sql | NO RUN |  |  |  |
| Queries | new-sql-query | NO RUN |  |  |  |
| Queries | new-visual-query | NO RUN |  |  |  |
| Queries | query-layout | NO RUN |  |  |  |
| Queries | query-postprocessing | NO RUN |  |  |  |
| Queries | transformations | NO RUN |  |  |  |
| Queries | visual-query-advanced | NO RUN |  |  |  |
| Scripts | [browser](Scripts/browser-run.md) | PASS | The Scripts Browser scenario passed well. The context pane shows all expected accordions (Details, Script, Run, Activ... |  |  |
| Scripts | [create](Scripts/create-run.md) | PARTIAL | The Create scenario completed successfully overall. The script `testRscript` was created, parameters configured, save... |  |  |
| Scripts | [delete](Scripts/delete-run.md) | PASS | All 5 steps passed. The delete flow works correctly with a confirmation dialog and immediate removal from the scripts... |  |  |
| Scripts | [edit](Scripts/edit-run.md) | PASS | All 6 steps passed. The Edit scenario works correctly — edits are saved persistently and visible on re-open. |  |  |
| Scripts | layout | NO RUN |  |  |  |
| Scripts | [run](Scripts/run-run.md) | PARTIAL | Core run functionality works: the script can be triggered from context menu with a table selection and from the conso... |  |  |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | PASS | All tested steps passed. SPGI.csv opened with TestSchema1 sticky metadata schema pre-configured. The Sticky meta pane... | 3.0s | 17.0s |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | PASS | Steps 1-2 passed: SPGI.csv opened with TestSchema1 sticky metadata schema, and cloning the table preserves the schema... | 3.0s | 21.0s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | PASS | All steps passed. The Sticky Meta Schemas browser at `/meta/schemas` shows 20 schemas including TestSchema1. The "NEW... | 3.0s | 7.0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | FAIL | Step 1 passed: navigated to Databases > Postgres and found CHEMBL connection. Step 2 failed: the "Database meta" sect... |  |  |
| Tooltips | actions-in-the-context-menu | NO RUN |  |  |  |
| Tooltips | default-tooltip-visibility | NO RUN |  |  |  |
| Tooltips | default-tooltip | NO RUN |  |  |  |
| Tooltips | edit-tooltip | NO RUN |  |  |  |
| Tooltips | line-chart---aggregated-tooltip | NO RUN |  |  |  |
| Tooltips | tooltip-properties | NO RUN |  |  |  |
| Tooltips | uniform-default-tooltip | NO RUN |  |  |  |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS | All steps passed. The 3D Scatter Plot opens correctly as a rotatable 3D model, data point clicking selects the corres... | 3.0s | 10.0s |
| Viewers | [basic-operations](Viewers/FilterPanel/basic-operations-run.md) | PASS | All filter panel basic operations work correctly on dev server. Structure/substructure filtering (32 rows for benzene... |  |  |
| Viewers | [chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | PASS | All steps pass. All 6 chem search type modes produce correct row counts matching scenario expectations exactly. Bio s... |  |  |
| Viewers | [cloned-views](Viewers/FilterPanel/cloned-views-run.md) | PASS | All 14 scenario steps passed on dev server. The cloned view correctly inherits filter state from the original view, i... |  |  |
| Viewers | [collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | PASS | All 9 steps passed. Table linking (selection-to-filter and filter-to-filter) works correctly. Filtering propagates be... | 3.0s | 30.0s |
| Viewers | [combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | PASS | All 15 steps passed. The combined boolean filter correctly aggregates boolean columns, supports OR mode filtering, in... |  |  |
| Viewers | [expression-filter](Viewers/FilterPanel/expression-filter-run.md) | PASS | Expression filter works correctly in all tested modes. Adding rules, switching AND/OR, removing rules, free-text mode... |  |  |
| Viewers | [hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | PASS | All 12 steps passed. The hierarchical filter correctly supports multi-level tree filtering with expandable nodes, che... |  |  |
| Viewers | [text-filter](Viewers/FilterPanel/text-filter-run.md) | PASS | All 9 steps passed. The text filter correctly searches string columns, highlights matches, supports multiple search t... |  |  |
| Viewers | [viewers](Viewers/FilterPanel/viewers-run.md) | PASS | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to re... | 5.0s | 1m 12s |
| Viewers | [color-consistency](Viewers/Legend/color-consistency-run.md) | PASS | All 7 steps passed. Categorical color coding on "Stereo Category" propagates correctly across all viewer types (Histo... | 3.0s | 24.0s |
| Viewers | [filtering](Viewers/Legend/filtering-run.md) | PASS | All legend filtering behaviors work correctly. Legends dynamically update when dataframe filters, in-viewer filters, ... | 5.0s | 35.9s |
| Viewers | [line-chart](Viewers/Legend/line-chart-run.md) | PASS | All 10 steps passed. The Line chart legend correctly reflects split categories with diverse colors, Multi Axis mode s... | 3.0s | 47.0s |
| Viewers | [scatterplot](Viewers/Legend/scatterplot-run.md) | PASS | All 5 scenarios passed (38/38 steps). The scatterplot legend correctly displays combined color+marker legends, update... | 3.0s | 50.0s |
| Viewers | [structure-rendering](Viewers/Legend/structure-rendering-run.md) | PASS | All steps passed. The Core column (containing SMILES) renders molecule structure images in the scatterplot marker leg... | 3.0s | 26.0s |
| Viewers | [visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | PASS | All legend visibility and positioning features work correctly. Legends display on viewers with color/split columns, c... | 5.0s | 43.7s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | PASS |  | 1m 0s | 16.8s |
| Viewers | bar-chart-tests | NO RUN |  |  |  |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | PASS | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color cod... | 10.0s | 1m 0s |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL | 17 of 19 sections passed. Section 18 (Visualization zoom) is AMBIGUOUS: project save fails on dev server, and layout ... |  | 47.0s |
| Viewers | [calendar](Viewers/calendar-run.md) | PASS | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selec... | 1m 0s | 9.4s |
| Viewers | [color-coding-(linked)](Viewers/color-coding-(linked)-run.md) | PASS | All tested steps passed. Linked color coding works via column tags (`.color-coding-type`, `.color-coding-source-colum... | 3.0s | 6.0s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | PASS | Core color coding steps passed: Linear (AGE, STARTED), Categorical (SEX, CONTROL) color coding applied and verified v... | 3.0s | 1m 24s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | PARTIAL | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operat... | 3.0s | 22.5s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | PASS | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI intera... | 2m 0s |  |
| Viewers | [form](Viewers/form-run.md) | PASS | All three scenarios passed in both the MCP run and the generated Playwright spec. The spec needed one revision after ... | 5.0s | 34.7s |
| Viewers | [grid-viewer](Viewers/grid-viewer-run.md) | PASS | All 16 steps completed. UI path used where DOM selectors were reliable (Toolbox icon, SAVE dialog, viewer icons). JS ... | 1m 0s | 57.7s |
| Viewers | [grid](Viewers/grid-run.md) | PARTIAL | The Grid scenario was exercised primarily via the JS API due to the canvas-based nature of grid interactions (headers... | 1m 0s | 28.4s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | PARTIAL | Steps 1, 2, 5 passed cleanly. Step 3 (table switching) requires JS-API fallback because the property panel's `<select... | 5.0s | 18.4s |
| Viewers | [histogram](Viewers/histogram-run.md) | PARTIAL | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance... | 30.0s | 50.0s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | PASS | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout... | 5.0s | 2m 30s |
| Viewers | [map](Viewers/map-run.md) | PASS | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker siz... | 3.0s | 9.0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | PASS | Steps 1-2, 5, 7 passed: Matrix plot viewer added to SPGI, properties accessible (xColumnNames, yColumnNames, innerVie... | 3.0s | 11.0s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | PARTIAL | Viewer core functionality (auto-pick, column rebinding, property-driven styling, filter integration, clean close) all... | 3.0s | 7.4s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | PASS | All 36 steps across 12 sections passed. PC Plot properties, context menus, filter interaction, column management, den... | 5.0s | 38.3s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appea... | 10.0s | 58.2s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | PASS | All 10 scenario steps passed. Pivot table viewer behaves correctly on `demog` and `SPGI` — defaults match the documen... | 5.0s | 22.0s |
| Viewers | rendering-structures-on-the-axes | NO RUN |  |  |  |
| Viewers | [row-source](Viewers/row-source-run.md) | PASS | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with al... | 5.0s | 1m 24s |
| Viewers | scatter-plot-tests | NO RUN |  |  |  |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | PASS | All 20 test sections passed on dev.datagrok.ai. Playwright spec executed in 1.1 min (headed mode). Context menu, rect... | 5.0s | 1m 6s |
| Viewers | statistics-viewer | NO RUN |  |  |  |
| Viewers | [tile](Viewers/tile-run.md) | PASS | All 6 scenario sub-steps passed on dev.datagrok.ai. Tile viewer correctly rebinds across dataframes, persists layouts... | 30.0s | 1m 6s |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | PARTIAL | All 9 scenario steps were reproduced. 7 pure PASS, 1 AMBIGUOUS (step 7 — heavy category skew masked shift/ctrl semant... | 2.0s | 7.4s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | PARTIAL | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range s... | 30.0s | 1m 48s |
| Viewers | viewers-docking | NO RUN |  |  |  |
| Viewers | working-with-nan-&-infinity | NO RUN |  |  |  |
| Viewers | [world-cloud](Viewers/world-cloud-run.md) | PASS | All tested steps passed. Word cloud viewer added to SPGI with column set to "Primary Series Name". Properties include... | 3.0s | 13.0s |
| **Total** | **186** | **84 PASS / 35 PARTIAL / 11 FAIL / 48 NO RUN / 1 AMBIGUOUS / 7 SKIP** |  | **8.8s** | **31.6s** |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage is 74%. No FAIL folders, but the following are PARTIAL: Bio, Browse, Charts, Chem, Connections, EDA, Models, Peptides, PowerPack, Projects, Scripts, StickyMeta, Viewers.

### Blocking Issues
- Apps: run coverage 0/2
- General: run coverage 0/11
- LocalCashing: run coverage 0/1
- Notebooks: run coverage 0/4
- Queries: run coverage 0/14
- Tooltips: run coverage 0/7
- Bio: PARTIAL — some scenarios failed or were incomplete
- Browse: PARTIAL — some scenarios failed or were incomplete
- Charts: PARTIAL — some scenarios failed or were incomplete
- Chem: PARTIAL — some scenarios failed or were incomplete
- Connections: PARTIAL — some scenarios failed or were incomplete
- EDA: PARTIAL — some scenarios failed or were incomplete
- Models: PARTIAL — some scenarios failed or were incomplete
- Peptides: PARTIAL — some scenarios failed or were incomplete
- PowerPack: PARTIAL — some scenarios failed or were incomplete
- Projects: PARTIAL — some scenarios failed or were incomplete
- Scripts: PARTIAL — some scenarios failed or were incomplete
- StickyMeta: PARTIAL — some scenarios failed or were incomplete
- Viewers: PARTIAL — some scenarios failed or were incomplete
