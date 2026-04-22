# Test Track — Global Report

**Date**: 2026-04-22
**Verdict**: Conditionally ready

## Folder Summary

**Total**: 187 tests · Run: 138/187 (74%) · Playwright: 107/187 (57%) · Mean Browser: 3m 37.4s · Mean Spec Gen: 50.5s · Mean Spec Run: 41.6s · Mean Total: 5m 4.8s

| Folder | Tests | Run | Playwright | Status | Mean Browser | Mean Spec Gen | Mean Spec Run | Mean Total |
|--------|-------|-----|------------|--------|--------------|---------------|---------------|------------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | NO DATA |  |  |  |  |
| Bio | 9 | 9/9 (100%) | 1/9 (11%) | PARTIAL | 2m 45s | 5s |  | 2m 46.2s |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | PARTIAL |  |  |  |  |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | PARTIAL | 2m 10s | 45s | 34.3s | 3m 29.3s |
| Chem | 14 | 14/14 (100%) | 14/14 (100%) | PARTIAL | 1m 28.6s | 31.8s | 38.4s | 2m 38.8s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | PARTIAL |  |  |  |  |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | PARTIAL | 3m 10.1s | 56.4s | 45.4s | 4m 51.9s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | PARTIAL | 1m 42.5s | 58s | 17.5s | 2m 58s |
| General | 11 | 0/11 (0%) | 0/11 (0%) | NO DATA |  |  |  |  |
| LocalCashing | 1 | 0/1 (0%) | 0/1 (0%) | NO DATA |  |  |  |  |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | PARTIAL |  |  |  |  |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | NO DATA |  |  |  |  |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 27.5s | 3s | 47.8s | 1m 18.2s |
| PowerPack | 10 | 10/10 (100%) | 6/10 (60%) | PARTIAL | 32s | 2.5s | 5.3s | 37.1s |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | PARTIAL |  |  |  |  |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | NO DATA |  |  |  |  |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | PARTIAL |  |  |  |  |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 23.8s | 2.2s | 15s | 37.2s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | NO DATA |  |  |  |  |
| Viewers | 49 | 42/49 (86%) | 42/49 (86%) | PARTIAL | 5m 49.6s | 1m 6.7s | 50.7s | 7m 45.8s |

## All Tests

**Total**: 187 tests · 74 PASS / 46 PARTIAL / 10 FAIL / 1 AMBIGUOUS / 7 SKIP / 49 NO RUN · Mean Browser: 3m 37.4s · Mean Spec Gen: 50.5s · Mean Spec Run: 41.6s · Mean Total: 5m 4.8s

| Folder | Test | Status | Description | Browser | Spec Gen | Spec Run | Total | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|--------|------|--------|-------------|---------|----------|----------|-------|-----------|------------|------------|---------|
| Apps | apps | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Apps | tutorials | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Bio | [analyze](Bio/analyze-run.md) | PASS → PASS | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset... | 8m 0s | 5s |  | 8m 5s | +0s | +0s | — | +0s |
| Bio | [composition-analysis](Bio/composition-analysis-run.md) | PARTIAL → PARTIAL | 4 of 5 steps passed. The Composition/WebLogo viewer opens correctly and properties are accessible. The letter-click s... |  |  |  |  | — | — | — | — |
| Bio | [convert](Bio/convert-run.md) | PASS → PASS | All 4 Bio convert/transform functions work correctly on FASTA data. Get Region extracts a subsequence region. PolyToo... | 2m 0s |  |  | 2m 0s | +0s | — | — | +0s |
| Bio | [manage](Bio/manage-run.md) | PASS → PASS | All 3 steps passed. The Manage Monomer Libraries view opens as a full view showing 5 monomer library JSON files (incr... | 30s |  |  | 30s | +0s | — | — | +0s |
| Bio | [msa](Bio/msa-run.md) | PASS → PASS | All 6 steps passed. MSA dialog opens with correct fields, Alignment Parameters button adds gap penalty inputs as expe... |  |  |  |  | — | — | — | — |
| Bio | [pepsea](Bio/pepsea-run.md) | AMBIGUOUS → AMBIGUOUS | The MSA dialog opens correctly for HELM data and shows MAFFT-based method options (mafft --auto, linsi, ginsi, etc.)... |  |  |  |  | — | — | — | — |
| Bio | [search](Bio/search-run.md) | PASS → PASS | All 4 steps passed. Bio > Search > Subsequence Search opens a filter panel with a Sequence bio substructure filter. T... | 30s |  |  | 30s | +0s | — | — | +0s |
| Bio | [sequence-activity-cliffs](Bio/sequence-activity-cliffs-run.md) | PASS → PASS | All 6 steps passed. Activity Cliffs works correctly with both default parameters (UMAP/Hamming) and custom parameters... |  |  |  |  | — | — | — | — |
| Bio | [sequence-space](Bio/sequence-space-run.md) | PASS → PASS | All 6 steps passed. Sequence Space works correctly with both default (UMAP/Hamming) and custom (t-SNE/Needlemann-Wuns... |  |  |  |  | — | — | — | — |
| Browse | [browse](Browse/browse-run.md) | PARTIAL → PARTIAL | 4 steps passed, 1 partial. Browse tree structure is complete, demos work, URL routing works for files and sections. I... |  |  |  |  | — | — | — | — |
| Browse | [browse-tree-states](Browse/browse-tree-states-run.md) | PARTIAL → PARTIAL | 1 step tested with partial result. The Browse tree correctly preserves its expand/collapse state within a single sess... |  |  |  |  | — | — | — | — |
| Browse | [japanese-in-myfiles](Browse/japanese-in-myfiles-run.md) | PASS → PASS | All 3 steps passed. The file 芹沢 貴之 こんにちは.csv displays correctly with no garbled or incorrect characters. Japanese Kan... |  |  |  |  | — | — | — | — |
| Browse | local-deploy | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Browse | package-manager | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Browse | [spaces](Browse/spaces-run.md) | PARTIAL → PARTIAL | 8 steps passed, 0 failed, 24 skipped out of 32 total steps. Core CRUD operations (create root/child/nested, rename, a... |  |  |  |  | — | — | — | — |
| Browse | [spaces-(ui-only)](Browse/spaces-(ui-only)-run.md) | FAIL → FAIL | 1 step passed, 2 ambiguous, 27 skipped out of 30 total steps. This scenario focuses on UI-only interactions not cover... |  |  |  |  | — | — | — | — |
| Charts | [radar](Charts/radar-run.md) | PASS → PASS | Radar viewer reproduces on dev for earthquakes.csv and demog.csv. Table switching, column toggles, and axis ranges beh... | 1m 30s | 35s | 33s | 2m 38s | +0s | +0s | +0s | +0s |
| Charts | [sunburst](Charts/sunburst-run.md) | PARTIAL → PARTIAL | 9 of 9 steps attempted: 8 passed, 1 ambiguous (canvas multi-selection), 1 skipped (old layout). Playwright spec passe... | 3m 0s | 1m 0s | 34s | 4m 34s | +0s | +0s | +0s | +0s |
| Charts | [tree](Charts/tree-run.md) | PARTIAL → PARTIAL | Setup (demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduces cleanly. The four scenario steps are AMBIGUOUS... | 2m 0s | 40s | 36s | 3m 16s | +0s | +0s | +0s | +0s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PASS → PASS | Activity Cliffs computation works correctly on SPGI.csv. UMAP embedding produces a scatter plot with molecule thumbna... | 30s | 20s | 1m 0s | 1m 50s | +0s | +0s | +0s | +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | PARTIAL → PARTIAL | Smoke coverage only: Scaffold Tree viewer launches from the Chem menu; the magic wand generates a scaffold tree on SPG... | 40s | 25s | 50.3s | 1m 55.3s | +0s | +0s | +0s | +0s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | PASS → PASS | Scaffold Tree viewer launches from Chem → Scaffold Tree; the magic-wand generator produces scaffold nodes within 25s. ... | 1m 15s | 30s | 40.2s | 2m 25.2s | +0s | +0s | +0s | +0s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | PASS → PASS | All 4 steps passed. The Similarity Search viewer opened correctly showing similar molecules with Tanimoto/Morgan simi... | 25s | 25s | 22.7s | 1m 12.7s | +0s | +0s | +0s | +0s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | PASS → PASS | All 4 core steps passed. The Structure filter panel opened with 43 filters for SPGI.csv. Benzene substructure search... | 30s | 25s | 24s | 1m 19s | +0s | +0s | +0s | +0s |
| Chem | [calculate](Chem/calculate-run.md) | FAIL → FAIL | Calculate Descriptors cannot be exercised on dev. The Chem top menu fails to open its popup via both DOM dispatch and ... | 8m 0s | 1m 0s | 38s | 9m 38s | +0s | +0s | +0s | +0s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS → PASS | Chemical Space works correctly with both default (UMAP) and modified (t-SNE) parameters. Scatter plots are generated... | 20s | 20s | 56.3s | 1m 36.3s | +0s | +0s | +0s | +0s |
| Chem | [chemprop](Chem/chemprop-run.md) | PARTIAL → PARTIAL | Steps 1-3 passed: mol1K.sdf loaded correctly, Predictive Model view opened, and columns were configured. However, mod... | 1m 0s | 30s | 18.3s | 1m 48.3s | +0s | +0s | +0s | +0s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PASS → PASS | All 4 steps passed. Elemental Analysis computed element counts (H, C, N, O, F, S, Cl, Br, I) and Molecule Charge for... | 1m 0s | 30s | 28.9s | 1m 58.9s | +0s | +0s | +0s | +0s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PASS → PASS | Filter panel shows a Structure filter for SPGI.csv Molecule column; opening the sketcher, typing c1ccccc1 + Enter + O... | 30s | 25s | 20.8s | 1m 15.8s | +0s | +0s | +0s | +0s |
| Chem | [info-panels](Chem/info-panels-run.md) | PASS → PASS | All 5 steps passed. The canonical_smiles column Context Panel shows 11 info panels including Chemistry, all expanding... | 3m 10s | 1m 0s | 31s | 4m 41s | +0s | +0s | +0s | +0s |
| Chem | [mmp](Chem/mmp-run.md) | PASS → PASS | MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset after ~60s compute. De... | 30s | 20s | 1m 12s | 2m 2s | +0s | +0s | +0s | +0s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PASS → PASS | All tested steps passed. The R-Groups Analysis dialog opened correctly, MCS computation completed, and the trellis pl... | 2m 20s | 45s | 58.7s | 4m 3.7s | +0s | +0s | +0s | +0s |
| Chem | [sketcher](Chem/sketcher-run.md) | PASS → PASS | Core steps 1-6 passed: smiles.csv opened, sketcher dialog displayed with OpenChemLib, C1CCCCC1 entered and rendered,... | 30s | 30s | 16.5s | 1m 16.5s | +0s | +0s | +0s | +0s |
| Connections | [adding](Connections/adding-run.md) | PARTIAL → PARTIAL | 6 of 7 steps fully passed, Step 5 was partial (TEST button works but actual connection test fails without real creden... |  |  |  |  | — | — | — | — |
| Connections | [browser](Connections/browser-run.md) | PARTIAL → PARTIAL | 7 of 9 steps passed, 2 ambiguous. Search filtering works correctly and the Context Pane shows all expected tabs (Deta... |  |  |  |  | — | — | — | — |
| Connections | [catalogs](Connections/catalogs-run.md) | FAIL → FAIL | 1 step passed, 1 failed, 15 skipped. The required `NorthwindTest` MS SQL connection is not present on public.datagrok... |  |  |  |  | — | — | — | — |
| Connections | [delete](Connections/delete-run.md) | PASS → PASS | All 8 steps passed. Both connections were deleted successfully. The confirmation dialog uses a red "DELETE" button (n... |  |  |  |  | — | — | — | — |
| Connections | [edit](Connections/edit-run.md) | PASS → PASS | 6 of 7 steps passed (1 skipped due to missing real credentials). The connection rename, credential modification, and... |  |  |  |  | — | — | — | — |
| Connections | [external-provider](Connections/external-provider-run.md) | FAIL → FAIL | All 7 steps skipped. This scenario requires a specific Postgres connection at db.datagrok.ai:54327 with superuser cre... |  |  |  |  | — | — | — | — |
| Connections | [identifiers](Connections/identifiers-run.md) | FAIL → FAIL | 1 step passed, 1 failed, 7 skipped. This scenario depends on a working Postgres connection to the Northwind database.... |  |  |  |  | — | — | — | — |
| Connections | [import-swagger](Connections/import-swagger-run.md) | FAIL → FAIL | All 7 steps skipped. This scenario requires manual interaction: downloading a YAML file to the local machine and drag... |  |  |  |  | — | — | — | — |
| Connections | [schema](Connections/schema-run.md) | PARTIAL → PARTIAL | 3 of 4 steps passed, 1 ambiguous. The "Browse schema" context menu option was not found in the current UI, but the sc... |  |  |  |  | — | — | — | — |
| Connections | [sparql](Connections/sparql-run.md) | PARTIAL → PARTIAL | 6 of 7 steps passed (1 failed). All UI steps worked correctly. The SPARQL connection was created and deleted successf... |  |  |  |  | — | — | — | — |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | PASS → PASS | All 6 steps passed in both MCP and Playwright runs. The full DiffStudio catalog workflow works end-to-end on dev.data... | 1m 53s | 39s | 38s | 3m 10s | +0s | +0s | +0s | +0s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS → PASS | All 4 steps passed in both MCP and Playwright runs. The PK-PD cyclic model loads correctly, Multiaxis/Facet plot tabs... | 1m 2s | 38s | 36s | 2m 16s | +0s | +0s | +0s | +0s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS → PASS | All 3 functional steps passed. The pk.ivp file opened from the Files browser as a DiffStudio view. Modifying Step (0.... | 2m 31s | 1m 32s | 1m 17s | 5m 20s | +0s | +0s | +0s | +0s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | PARTIAL → PARTIAL | Steps 1-5 pass; step 6 (actually running the fit) does not produce result rows within 2 minutes with recommended swit... | 10m 2s | 55s | 1m 0s | 11m 57s | +0s | +0s | +0s | +0s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS → PASS | All 6 steps passed. DiffStudio opened correctly, Bioreactor loaded from Library with 13 columns and 1001 rows. Both M... | 1m 31s | 36s | 26s | 2m 33s | +0s | +0s | +0s | +0s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | PASS → PASS | All 5 steps passed. The DiffStudio Edit toggle and `</>` export icon work correctly, producing a JS script with all m... | 4m 30s | 1m 40s | 1m 3s | 7m 13s | +0s | +0s | +0s | +0s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS → PASS | All 4 steps passed. The Sensitivity Analysis view opened correctly from the ribbon. Process mode switching updated FF... | 2m 3s | 44s | 39s | 3m 26s | +0s | +0s | +0s | +0s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | PASS → PASS | All 4 steps passed. The Acid Production (GA-production) model loaded from Library with 2-stage simulation (1st stage... | 1m 49s | 47s | 24s | 3m 0s | +0s | +0s | +0s | +0s |
| EDA | [anova](EDA/anova-run.md) | PASS → PASS | All 3 steps passed successfully. The ANOVA dialog opened correctly with default parameters (Category=RACE, Feature=AG... | 1m 30s | 30s | 26s | 2m 26s | +0s | +0s | +0s | +0s |
| EDA | [ML methods/linear-regression](EDA/ML methods/linear-regression-run.md) | PASS → PASS | Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation... | 45s | 2s | 6.8s | 53.8s | +0s | +0s | +0s | +0s |
| EDA | [ML methods/pls-regression](EDA/ML methods/pls-regression-run.md) | PARTIAL → PARTIAL | PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-... | 1m 0s | 2s | 6.9s | 1m 8.9s | +0s | +0s | +0s | +0s |
| EDA | [ML methods/softmax](EDA/ML methods/softmax-run.md) | FAIL → FAIL | Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns an... | 10s | 2s | 2.6s | 14.6s | +0s | +0s | +0s | +0s |
| EDA | [ML methods/xgboost1](EDA/ML methods/xgboost1-run.md) | PARTIAL → PARTIAL | XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as... | 5s | 2s | 2.7s | 9.7s | +0s | +0s | +0s | +0s |
| EDA | [ML methods/xgboost2](EDA/ML methods/xgboost2-run.md) | PARTIAL → PARTIAL | XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interac... | 5s | 2s | 2.7s | 9.7s | +0s | +0s | +0s | +0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PARTIAL → PARTIAL | 2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check). Playwright spec is green — ... | 2m 30s | 2m 0s | 13s | 4m 43s | +0s | +0s | +0s | +0s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL → PARTIAL | 6 of 7 steps passed (1 partial due to missing test file). The Pareto Front viewer works correctly: non-numeric column... | 4m 0s | 2m 0s | 32s | 6m 32s | +0s | +0s | +0s | +0s |
| EDA | [pca](EDA/pca-run.md) | PARTIAL → PARTIAL | MCP produced 3 PASS / 1 FAIL / 1 SKIP. The dialog path works but clicking OK closed the dialog without adding PC1/PC2... | 5m 0s | 3m 0s | 1m 7s | 9m 7s | +0s | +0s | +0s | +0s |
| EDA | [pls](EDA/pls-run.md) | FAIL → FAIL | MCP produced 2 PASS / 1 PARTIAL / 1 FAIL. Dialog path works but PLS validation disables RUN when the Predict column i... | 2m 0s | 2m 0s | 15s | 4m 15s | +0s | +0s | +0s | +0s |
| General | api-samples | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | files-cache | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | first-login | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | inactivity-response | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | login | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | molecule-in-exported-csv | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | network | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | profile-settings | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | startup-time | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | table-manager | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| General | tabs-reordering | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| LocalCashing | local-cashing | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Models | [apply](Models/apply-run.md) | PASS → PASS | All 4 steps passed. The Apply Model workflow functions correctly end-to-end. The TestDemog model (trained in Train.md... |  |  |  |  | — | — | — | — |
| Models | [browser](Models/browser-run.md) | PARTIAL → PARTIAL | 4 of 6 steps passed; 2 skipped because only 1 model was available (TestDemog was the only model — the second numeric... |  |  |  |  | — | — | — | — |
| Models | [chemprop](Models/chemprop-run.md) | FAIL → FAIL | 2 of 17 sub-steps passed, 13 skipped, 1 failed, 1 ambiguous. The scenario fails entirely due to the Chemprop Docker c... |  |  |  |  | — | — | — | — |
| Models | [delete](Models/delete-run.md) | PASS → PASS | All 5 steps passed. The model deletion workflow works correctly end-to-end. The right-click context menu, confirmatio... |  |  |  |  | — | — | — | — |
| Models | [predictive-models](Models/predictive-models-run.md) | PASS → PASS | All 20 sub-steps passed. The full lifecycle (Train → Apply → Apply on new dataset → Delete) for EDA-based predictive... |  |  |  |  | — | — | — | — |
| Models | [train](Models/train-run.md) | PASS → PASS | All 10 steps passed. The Train Model workflow functions correctly on public.datagrok.ai using the built-in EDA engine... |  |  |  |  | — | — | — | — |
| Notebooks | browser | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Notebooks | create | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Notebooks | delete | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Notebooks | edit | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Peptides | [info-panels](Peptides/info-panels-run.md) | PASS → PASS | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are re... | 17s | 3s | 11s | 31s | +0s | +0s | +0s | +0s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | PARTIAL → PARTIAL | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variabil... | 25s | 3s | 1m 18s | 1m 46s | +0s | +0s | +0s | +0s |
| Peptides | [peptides](Peptides/peptides-run.md) | PARTIAL → PARTIAL | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clust... | 18s | 3s | 12s | 33s | +0s | +0s | +0s | +0s |
| Peptides | [sar](Peptides/sar-run.md) | PARTIAL → PARTIAL | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Res... | 50s | 3s | 1m 30s | 2m 23s | +0s | +0s | +0s | +0s |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | PASS → PASS | All 5 steps passed on dev.datagrok.ai. The Add New Column dialog opens correctly, displays cleanly with no visual iss... | 22s | 2s | 5.3s | 29.3s | +0s | +0s | +0s | +0s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | PARTIAL → PARTIAL | Core Add New Column functionality works correctly: creating calculated columns, chaining formulas, and propagating co... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | PASS → PASS | All 8 steps passed. Autocomplete works for both functions (typing or Ctrl+Space) and columns ($ symbol). Function sel... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | PASS → PASS | All 5 steps passed. The formula dependency chain (Weight2 -> Weight3 -> Weight4) works correctly. Changing a value in... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | SKIP → SKIP | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on rel... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | PASS → PASS | All 4 steps passed. Column names in both ${} (scalar) and $[] (aggregate) syntax are highlighted with the cm-column-n... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | PASS → PASS | All 4 steps passed. Hovering over a function name in the formula editor displays a tooltip with the function signatur... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | SKIP → SKIP | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on rel... |  |  |  |  | — | — | — | — |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | PARTIAL → PARTIAL | Part 1 (create, save, apply enrichment) passed. The enrichment feature correctly joins the Customers table to Orders... | 42s | 3s |  | 45s | +0s | +0s | — | +0s |
| PowerPack | [formula-lines](PowerPack/formula-lines-run.md) | SKIP → SKIP | This scenario file contains only references to GitHub tickets related to formula lines. No executable test steps were... |  |  |  |  | — | — | — | — |
| Projects | [browser](Projects/browser-run.md) | PARTIAL → PARTIAL | 5 of 9 steps passed, 3 skipped, 1 ambiguous. Browse > Dashboards view works correctly: projects are listed, searchabl... |  |  |  |  | — | — | — | — |
| Projects | [complex](Projects/complex-run.md) | SKIP → SKIP | All 13 steps skipped. This is the most complex scenario requiring tables from 7+ different sources, drag-and-drop, en... |  |  |  |  | — | — | — | — |
| Projects | [custom-creation-scripts](Projects/custom-creation-scripts-run.md) | SKIP → SKIP | All 5 steps skipped. This scenario requires running a custom JavaScript script with Data Sync enabled, then modifying... |  |  |  |  | — | — | — | — |
| Projects | [deleting](Projects/deleting-run.md) | PARTIAL → PARTIAL | 2 of 4 steps passed, 1 skipped, 1 ambiguous. Project deletion works via the API (`grok.dapi.projects.delete()`). The... |  |  |  |  | — | — | — | — |
| Projects | [opening](Projects/opening-run.md) | PARTIAL → PARTIAL | All 5 steps passed. Projects from the Uploading step are accessible in Browse > Dashboards. Context Panel correctly s... |  |  |  |  | — | — | — | — |
| Projects | [project-url](Projects/project-url-run.md) | SKIP → SKIP | All steps skipped. This scenario depends on Projects copy_clone.md (order 5) which was not fully executed. The Link/C... |  |  |  |  | — | — | — | — |
| Projects | [projects-copy_clone](Projects/projects-copy_clone-run.md) | SKIP → SKIP | 2 of 5 steps passed, 3 skipped. Project preview and opening work. Copy/clone/link operations were not tested because... |  |  |  |  | — | — | — | — |
| Projects | [uploading](Projects/uploading-run.md) | PARTIAL → PARTIAL | 8 of 14 steps passed, 6 skipped. Core project creation from local tables, file shares, query results, and join result... |  |  |  |  | — | — | — | — |
| Queries | adding | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | browse-&-save-project | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | browser | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | columns-inspect | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | deleting | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | edit | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | get-all-get-top-100 | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | ms-sql | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | new-sql-query | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | new-visual-query | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | query-layout | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | query-postprocessing | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | transformations | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Queries | visual-query-advanced | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Scripts | [browser](Scripts/browser-run.md) | PASS → PASS | The Scripts Browser scenario passed well. The context pane shows all expected accordions (Details, Script, Run, Activ... |  |  |  |  | — | — | — | — |
| Scripts | [create](Scripts/create-run.md) | PARTIAL → PARTIAL | The Create scenario completed successfully overall. The script `testRscript` was created, parameters configured, save... |  |  |  |  | — | — | — | — |
| Scripts | [delete](Scripts/delete-run.md) | PASS → PASS | All 5 steps passed. The delete flow works correctly with a confirmation dialog and immediate removal from the scripts... |  |  |  |  | — | — | — | — |
| Scripts | [edit](Scripts/edit-run.md) | PASS → PASS | All 6 steps passed. The Edit scenario works correctly — edits are saved persistently and visible on re-open. |  |  |  |  | — | — | — | — |
| Scripts | [run](Scripts/run-run.md) | PARTIAL → PARTIAL | Core run functionality works: the script can be triggered from context menu with a table selection and from the conso... |  |  |  |  | — | — | — | — |
| Scripts | layout | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | PASS → PASS | All tested steps passed. SPGI.csv opened with TestSchema1 sticky metadata schema pre-configured. The Sticky meta pane... | 35s | 3s | 17s | 55s | +0s | +0s | +0s | +0s |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | PASS → PASS | Steps 1-2 passed: SPGI.csv opened with TestSchema1 sticky metadata schema, and cloning the table preserves the schema... | 25s | 3s | 21s | 49s | +0s | +0s | +0s | +0s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | PASS → PASS | All steps passed. The Sticky Meta Schemas browser at `/meta/schemas` shows 20 schemas including TestSchema1. The "NEW... | 10s | 3s | 7s | 20s | +0s | +0s | +0s | +0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | FAIL → FAIL | Step 1 passed: navigated to Databases > Postgres and found CHEMBL connection. Step 2 failed: the "Database meta" sect... | 25s | 0s |  | 25s | +0s | +0s | — | +0s |
| Tooltips | actions-in-the-context-menu | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | default-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | default-tooltip-visibility | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | edit-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | line-chart---aggregated-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | tooltip-properties | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | uniform-default-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS → PASS | All 14 steps (setup + 13 scenario sections) passed in both the MCP run and the Playwright replay. Spec completes in ~... | 7m 0s | 2m 0s | 50s | 9m 50s | -29m 0s | -3m 0s | +5s | -31m 55s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | PASS → PASS |  | 7m 0s | 1m 0s | 17s | 8m 17s | +0s | +0s | +0s | +0s |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | PASS → PASS | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color cod... | 3m 3s | 21s | 52s | 4m 16s | +3s | +11s | -8s | +6s |
| Viewers | bar-chart-tests | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL → PARTIAL | 17 of 19 sections passed cleanly; section 8 combined into section 7 in the spec. Section 18 is AMBIGUOUS — `grok.dapi.... | 1m 5s | 8s | 32s | 1m 45s | -1m 55s | 8s (new) | -15s | -2m 2s |
| Viewers | [calendar](Viewers/calendar-run.md) | PASS → PASS | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selec... | 25s | 1m 0s | 9.4s | 1m 34.4s | +0s | +0s | +0s | +0s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | PASS → PASS | Core color coding steps passed: Linear (AGE, STARTED), Categorical (SEX, CONTROL) color coding applied and verified v... | 15s | 3s | 1m 24s | 1m 42s | +0s | +0s | +0s | +0s |
| Viewers | [color-coding-(linked)](Viewers/color-coding-(linked)-run.md) | PASS → PASS | All tested steps passed. Linked color coding works via column tags (`.color-coding-type`, `.color-coding-source-colum... | 10s | 3s | 6s | 19s | +0s | +0s | +0s | +0s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | PARTIAL → PARTIAL | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operat... | 5m 0s | 3s | 22s | 5m 25s | +0s | +0s | +0s | +0s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | PASS → PASS | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI intera... | 18m 0s | 2m 0s |  | 20m 0s | +0s | +0s | — | +0s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | PASS → PASS | Ran basic-operations end-to-end against dev. All 31 scenario steps passed in the MCP browser phase (Section 1: struct... | 4m 27s | 9s | 50s | 5m 26s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | PASS → PASS | Ran chem-and-bio scenario end-to-end against dev. All 11 scenario steps passed in the MCP browser phase (Chem: open s... | 2m 50s | 42s | 47s | 4m 19s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | PASS → PASS | All 15 scenario steps PASSed on dev. spgi-100.csv loads correctly this time (previous run had to substitute SPGI.csv)... | 3m 16s | 14s | 55s | 4m 25s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | PASS → PASS | All 9 steps passed end-to-end on dev: table linking (SELECTION_TO_FILTER and FILTER_TO_FILTER) propagated correctly b... | 1m 46s | 17s | 35s | 2m 38s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | PASS → PASS | Ran combined-boolean-filter end-to-end against dev. All 13 numbered scenario steps passed in the MCP browser phase: S... | 2m 37s | 12s | 24s | 3m 13s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | PASS → PASS | All 14 steps passed in both the MCP run and the Playwright replay. Expression filter works correctly: 5-rule AND yiel... | 1m 14s | 8s | 23s | 1m 45s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | PASS → PASS | All 12 steps passed in the MCP run and in the Playwright replay (spec finished in 21.8s). The hierarchical filter cor... | 1m 15s | 21s | 23s | 1m 59s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | PASS → PASS | All 9 steps passed in the MCP run and in the Playwright replay (spec finished in 8.7s, total wall-clock 11.56s). The ... | 1m 12s | 20s | 12s | 1m 44s | +0s | +0s | +0s | +0s |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | PASS → PASS | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to re... | 4m 24s | 40s | 1m 3s | 6m 7s | +0s | +0s | +0s | +0s |
| Viewers | [form](Viewers/form-run.md) | PASS → PASS | All 14 sections of form-tests-pw.md exercised across 30 steps. 28 PASS, 2 AMBIGUOUS, 0 FAIL in MCP run. Playwright spe... | 18m 0s | 4m 0s | 3m 12s | 25m 12s | +0s | +0s | +0s | +0s |
| Viewers | [forms](Viewers/forms-run.md) | PASS → PASS | All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas ele... | 18m 0s | 3m 0s | 51.5s | 21m 51.5s | +0s | +0s | +0s | +0s |
| Viewers | [grid](Viewers/grid-run.md) | PARTIAL → PARTIAL | Grid tests ran 22 steps (spec softSteps); 17 passed outright and 5 AMBIGUOUS (Copy/Paste, Column Header Context Menu, ... | 11m 0s | 3m 0s | 1m 18s | 15m 18s | +0s | +0s | +0s | +0s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | PASS → PASS | All 14 heat-map sections exercised across 17 steps. 15 PASS, 1 AMBIGUOUS, 1 SKIP in MCP run. Playwright spec passed fu... | 18m 0s | 4m 0s | 48.9s | 22m 48.9s | +0s | +0s | +0s | +0s |
| Viewers | [histogram](Viewers/histogram-run.md) | PARTIAL → PARTIAL | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance... | 50s | 7s | 46s | 1m 43s | -1m 10s | -23s | -4s | -1m 37s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | PARTIAL → PARTIAL | Color consistency through layout round-trip works — the `.categorical-colors` tag survives save/reload and `R_ONE` st... | 2m 30s | 35s | 26s | 3m 31s | +0s | +0s | +0s | +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | PARTIAL → PARTIAL | Filtering legend updates work end-to-end in the MCP run: numeric filter, categorical filter, layout round-trip, compo... | 3m 10s | 1m 10s | 44s | 5m 4s | +0s | +0s | +0s | +0s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | PARTIAL → PARTIAL | Line chart legend and multi-axis behaviors are mostly correct: 7 legend items for 7 categories, layout round-trip pre... | 2m 10s | 40s | 32s | 3m 22s | +0s | +0s | +0s | +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | PARTIAL → PARTIAL | Categorical legend on scatter plot updates correctly when X axis changes (sub 2) and when the Filter Panel narrows ca... | 4m 15s | 1m 20s | 54s | 6m 29s | +0s | +0s | +0s | +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | PARTIAL → PARTIAL | Structure rendering in legends works for Scatter plot, Histogram, Line chart and Pie chart (canvas-based molecule thu... | 2m 35s | 40s | 29s | 3m 44s | +0s | +0s | +0s | +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | PARTIAL → PARTIAL | Scenario executed end-to-end with a mix of PASS, AMBIGUOUS, and FAIL. Legend display, source-swap, corner positioning... | 5m 45s | 1m 30s | 41s | 7m 56s | +0s | +0s | +0s | +0s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | PASS → PASS | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout... | 57s | 8s | 1m 43s | 2m 48s | -1m 3s | +3s | -47s | -1m 47s |
| Viewers | [map](Viewers/map-run.md) | PASS → PASS | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker siz... | 15s | 3s | 9s | 27s | +0s | +0s | +0s | +0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | PARTIAL → PARTIAL | Matrix Plot tests ran with 15 PASS, 3 AMBIGUOUS, 0 FAIL. The spec executed in 57.7s with all implemented steps passing... | 20m 0s | 3m 0s | 55.9s | 23m 55.9s | +0s | +0s | +0s | +0s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | PARTIAL → PARTIAL | 9 of 12 steps PASS; 3 SKIP (canvas-based node/edge interactions cannot be automated via DOM). The network diagram vie... | 8m 0s | 1m 30s | 22s | 9m 52s | +1m 15s | +1m 27s | +14.6s | +2m 56.6s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | PASS → PASS | All 13 scenario sections (mapped to 12 Playwright softSteps — scale and normalization are combined in the spec) passe... | 1m 8s | 8s | 47s | 2m 3s | -6m 52s | +3s | +9s | -6m 40s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS → PASS | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appea... | 40s | 7s | 47s | 1m 34s | -3m 20s | -3s | -11s | -3m 34s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | PARTIAL → PARTIAL | Pivot Table tests ran with 16 PASS, 2 AMBIGUOUS, 1 SKIP, 0 FAIL. The spec executed in 35.1s with all implemented steps... | 16m 0s | 3m 0s | 1m 12s | 20m 12s | +0s | +0s | +0s | +0s |
| Viewers | rendering-structures-on-the-axes | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [row-source](Viewers/row-source-run.md) | PASS → PASS | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with al... | 4m 0s | 5s | 1m 24s | 5m 29s | +0s | +0s | +0s | +0s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | PASS → PASS | All 20 sections passed during the MCP run on dev.datagrok.ai. The existing Playwright spec was re-run headed without m... | 3m 13s | 29s | 52s | 4m 34s | -4m 47s | +24s | -14s | -4m 37s |
| Viewers | scatter-plot-tests | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | statistics-viewer | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | PASS → PASS | 24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to manu... | 4m 0s | 3m 0s | 58s | 7m 58s | +0s | +0s | +0s | +0s |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | PARTIAL → PARTIAL | All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via ... | 28m 0s | 4m 0s | 46s | 32m 46s | +0s | +0s | +0s | +0s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | PARTIAL → PARTIAL | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range s... | 3m 0s | 30s | 1m 48s | 5m 18s | +0s | +0s | +0s | +0s |
| Viewers | viewers-docking | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | PASS → PASS | All 7 MCP scenario steps PASS. The Word Cloud viewer adds via both entry points (Add-Viewer gallery and Toolbox icon),... | 4m 15s | 1m 0s | 2m 7s | 7m 22s | +3m 55s | +57s | +1m 54s | +6m 46s |
| Viewers | word-cloud-tests | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | working-with-nan-&-infinity | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |

## Comparison with Previous Report

Delta from the previous `total-run.md` (dated 2026-04-22); signed with `+`/`-`, time deltas use the same format as the values.

### Totals

**Total**: Tests Δ **+0** · Run Δ **+0 (+0%)** · Playwright Δ **+0 (+0%)** · Status **PARTIAL → PARTIAL** · Browser Δ **-28.3s** · Spec Gen Δ **-0.7s** · Spec Run Δ **+0.5s** · Total Δ **-28s**

### By Folder

| Folder | Tests Δ | Run Δ | Playwright Δ | Status | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|--------|---------|-------|--------------|--------|-----------|------------|------------|---------|
| Apps | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Bio | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | — | +0s |
| Browse | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| Charts | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | +0s | +0s |
| Chem | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | +0s | +0s |
| Connections | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| DiffStudio | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | +0s | +0s |
| EDA | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | +0s | +0s |
| General | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| LocalCashing | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Models | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| Notebooks | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Peptides | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | +0s | +0s |
| PowerPack | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | +0s | +0s |
| Projects | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| Queries | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Scripts | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| StickyMeta | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0s | +0s | +0s |
| Tooltips | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Viewers | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | -1m 1.3s | -2s | +1.1s | -1m 0.5s |

### Per-Test Changes

| Folder | Test | Status | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|--------|------|--------|-----------|------------|------------|---------|
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS → PASS | -29m 0s | -3m 0s | +5s | -31m 55s |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | PASS → PASS | +3s | +11s | -8s | +6s |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL → PARTIAL | -1m 55s | 8s (new) | -15s | -2m 2s |
| Viewers | [histogram](Viewers/histogram-run.md) | PARTIAL → PARTIAL | -1m 10s | -23s | -4s | -1m 37s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | PASS → PASS | -1m 3s | +3s | -47s | -1m 47s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | PARTIAL → PARTIAL | +1m 15s | +1m 27s | +14.6s | +2m 56.6s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | PASS → PASS | -6m 52s | +3s | +9s | -6m 40s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS → PASS | -3m 20s | -3s | -11s | -3m 34s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | PASS → PASS | -4m 47s | +24s | -14s | -4m 37s |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | PASS → PASS | +3m 55s | +57s | +1m 54s | +6m 46s |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage holds at 74% (138/187). No folder is outright FAIL and no test changed status this cycle. The only movement is in `Viewers/`: 10 scenarios were re-executed with refreshed timings. Seven got substantially faster (`3d-scatter-plot` dropped 31m 55s after the model-thinking phase was isolated from spec execution; `pc-plot`, `scatter-plot`, `pie-chart`, `line-chart`, `histogram`, `box-plot` all shed ≥1m each), two got slower (`word-cloud` +6m 46s, `network-diagram` +2m 56.6s), and `bar-chart` was essentially flat. `box-plot` also picked up a Spec Gen timing it was missing before. Net effect on the Viewers folder: Mean Browser -1m 1s, Mean Total -1m 0.5s. The wider picture of under-executed folders (Apps, General, LocalCashing, Notebooks, Queries, Tooltips) is unchanged.

### Blocking Issues
- Apps: run coverage 0/2
- General: run coverage 0/11
- LocalCashing: run coverage 0/1
- Notebooks: run coverage 0/4
- Queries: run coverage 0/14
- Tooltips: run coverage 0/7
- Chem: calculate remains FAIL (top menu broken on dev)
- EDA: pls/softmax remain FAIL
- DiffStudio: fitting PARTIAL (fit job produces no rows in 2 min on dev)
- Viewers/Legend: six scenarios remain PARTIAL from the prior cycle — see previous report for specifics
- Bio, Browse, Charts, Connections, EDA, Models, Peptides, PowerPack, Projects, Scripts, StickyMeta, Viewers: PARTIAL
