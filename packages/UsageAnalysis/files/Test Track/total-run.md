# Test Track — Global Report

**Date**: 2026-04-21
**Verdict**: Conditionally ready

## Folder Summary

**Total**: 187 tests · Run: 138/187 (74%) · Playwright: 107/187 (57%) · Mean Browser: 4m 1s · Mean Spec Gen: 52s · Mean Spec Run: 42s · Mean Total: 5m 26s

| Folder | Tests | Run | Playwright | Status | Mean Browser | Mean Spec Gen | Mean Spec Run | Mean Total |
|--------|-------|-----|------------|--------|--------------|---------------|---------------|------------|
| Apps | 2 | 0/2 (0%) | 0/2 (0%) | NO DATA |  |  |  |  |
| Bio | 9 | 9/9 (100%) | 1/9 (11%) | PARTIAL | 2m 45s | 5.0s |  | 2m 46s |
| Browse | 7 | 5/7 (71%) | 0/7 (0%) | PARTIAL |  |  |  |  |
| Charts | 3 | 3/3 (100%) | 3/3 (100%) | PARTIAL | 2m 10s | 45s | 34s | 3m 29s |
| Chem | 14 | 14/14 (100%) | 14/14 (100%) | PARTIAL | 1m 29s | 32s | 38s | 2m 39s |
| Connections | 10 | 10/10 (100%) | 5/10 (50%) | PARTIAL |  |  |  |  |
| DiffStudio | 8 | 8/8 (100%) | 8/8 (100%) | PARTIAL | 3m 11s | 56s | 45s | 4m 52s |
| EDA | 10 | 10/10 (100%) | 10/10 (100%) | PARTIAL | 1m 43s | 58s | 18s | 2m 58s |
| General | 11 | 0/11 (0%) | 0/11 (0%) | NO DATA |  |  |  |  |
| LocalCashing | 1 | 0/1 (0%) | 0/1 (0%) | NO DATA |  |  |  |  |
| Models | 6 | 6/6 (100%) | 6/6 (100%) | PARTIAL |  |  |  |  |
| Notebooks | 4 | 0/4 (0%) | 0/4 (0%) | NO DATA |  |  |  |  |
| Peptides | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 28s | 3.0s | 48s | 1m 18s |
| PowerPack | 10 | 10/10 (100%) | 6/10 (60%) | PARTIAL | 32s | 2.5s | 5.3s | 37s |
| Projects | 8 | 8/8 (100%) | 4/8 (50%) | PARTIAL |  |  |  |  |
| Queries | 14 | 0/14 (0%) | 0/14 (0%) | NO DATA |  |  |  |  |
| Scripts | 6 | 5/6 (83%) | 0/6 (0%) | PARTIAL |  |  |  |  |
| StickyMeta | 4 | 4/4 (100%) | 4/4 (100%) | PARTIAL | 24s | 2.2s | 15s | 37s |
| Tooltips | 7 | 0/7 (0%) | 0/7 (0%) | NO DATA |  |  |  |  |
| Viewers | 49 | 42/49 (86%) | 42/49 (86%) | PARTIAL | 7m 12s | 1m 14s | 53s | 9m 9s |

## All Tests

**Total**: 187 tests · 80 PASS / 40 PARTIAL / 10 FAIL / 1 AMBIGUOUS / 7 SKIP / 49 NO RUN · Mean Browser: 4m 1s · Mean Spec Gen: 52s · Mean Spec Run: 42s · Mean Total: 5m 26s

| Folder | Test | Status | Description | Browser | Spec Gen | Spec Run | Total | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|--------|------|--------|-------------|---------|----------|----------|-------|-----------|------------|------------|---------|
| Apps | apps | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Apps | tutorials | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Bio | [analyze](Bio/analyze-run.md) | PASS → PASS | All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset... | 8m 0s | 5.0s |  | 8m 5s | +0s | +0.0s | — | +0s |
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
| Charts | [radar](Charts/radar-run.md) | PASS → PASS | Radar viewer reproduces on dev for earthquakes.csv and demog.csv. Table switching, column toggles, and axis ranges beh... | 1m 30s | 35s | 33s | 2m 38s | -2m 44s | -15s | +22.2s | -2m 37s |
| Charts | [sunburst](Charts/sunburst-run.md) | PARTIAL → PARTIAL | 9 of 9 steps attempted: 8 passed, 1 ambiguous (canvas multi-selection), 1 skipped (old layout). Playwright spec passe... | 3m 0s | 1m 0s | 34s | 4m 34s | +2m 10s | +57s | -8s | +3m 0s |
| Charts | [tree](Charts/tree-run.md) | PASS → PARTIAL | Setup (demog.csv + Tree viewer + CONTROL/SEX/RACE hierarchy) reproduces cleanly. The four scenario steps are AMBIGUOUS... | 2m 0s | 40s | 36s | 3m 16s | +1m 40s | +37s | +19s | +2m 36s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PASS → PASS | Activity Cliffs computation works correctly on SPGI.csv. UMAP embedding produces a scatter plot with molecule thumbna... | 30s | 20s | 1m 0s | 1m 50s | -20s | +17s | +3s | +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | FAIL → PARTIAL | Smoke coverage only: Scaffold Tree viewer launches from the Chem menu; the magic wand generates a scaffold tree on SPG... | 40s | 25s | 50.3s | 1m 55s | +10s | +22s (new) | +50.3s (new) | +1m 25s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | FAIL → PASS | Scaffold Tree viewer launches from Chem → Scaffold Tree; the magic-wand generator produces scaffold nodes within 25s. ... | 1m 15s | 30s | 40.2s | 2m 25s | -45s | +27s | +40.2s (new) | +22s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | PASS → PASS | All 4 steps passed. The Similarity Search viewer opened correctly showing similar molecules with Tanimoto/Morgan simi... | 25s | 25s | 22.7s | 1m 13s | -10s | +22s | +22.7s (new) | +35s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | PASS → PASS | All 4 core steps passed. The Structure filter panel opened with 43 filters for SPGI.csv. Benzene substructure search... | 30s | 25s | 24.0s | 1m 19s | -10s | +22s | +0s | +12s |
| Chem | [calculate](Chem/calculate-run.md) | PASS → FAIL | Calculate Descriptors cannot be exercised on dev. The Chem top menu fails to open its popup via both DOM dispatch and ... | 8m 0s | 1m 0s | 38s | 9m 38s | +7m 35s | +57s | +14s | +8m 46s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS → PASS | Chemical Space works correctly with both default (UMAP) and modified (t-SNE) parameters. Scatter plots are generated... | 20s | 20s | 56.3s | 1m 36s | -25s | +17s | +12.3s | +4s |
| Chem | [chemprop](Chem/chemprop-run.md) | PARTIAL → PARTIAL | Steps 1-3 passed: mol1K.sdf loaded correctly, Predictive Model view opened, and columns were configured. However, mod... | 1m 0s | 30s | 18.3s | 1m 48s | -1m 0s | +27s | +18.3s (new) | -15s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PASS → PASS | All 4 steps passed. Elemental Analysis computed element counts (H, C, N, O, F, S, Cl, Br, I) and Molecule Charge for... | 1m 0s | 30s | 28.9s | 1m 58s | +30s | +27s | +7.9s | +1m 4s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PARTIAL → PASS | Filter panel shows a Structure filter for SPGI.csv Molecule column; opening the sketcher, typing c1ccccc1 + Enter + O... | 30s | 25s | 20.8s | 1m 15s | +10s | +22s | -8.2s | +23s |
| Chem | [info-panels](Chem/info-panels-run.md) | PASS → PASS | All 5 steps passed. The canonical_smiles column Context Panel shows 11 info panels including Chemistry, all expanding... | 3m 10s | 1m 0s | 31s | 4m 41s | +2m 30s | +57s | -8s | +3m 19s |
| Chem | [mmp](Chem/mmp-run.md) | FAIL → PASS | MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset after ~60s compute. De... | 30s | 20s | 1m 12s | 2m 2s | -1m 20s | +17s | +1m 12s (new) | +9s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PASS → PASS | All tested steps passed. The R-Groups Analysis dialog opened correctly, MCS computation completed, and the trellis pl... | 2m 20s | 45s | 58.7s | 4m 4s | +1m 35s | +42s | +20.7s | +2m 38s |
| Chem | [sketcher](Chem/sketcher-run.md) | PASS → PASS | Core steps 1-6 passed: smiles.csv opened, sketcher dialog displayed with OpenChemLib, C1CCCCC1 entered and rendered,... | 30s | 30s | 16.5s | 1m 16s | +5s | +27s | +3.5s | +35s |
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
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | PASS → PASS | All 6 steps passed in both MCP and Playwright runs. The full DiffStudio catalog workflow works end-to-end on dev.data... | 1m 53s | 39s | 38s | 3m 10s | +1m 8s | +34s | +6s | +1m 48s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS → PASS | All 4 steps passed in both MCP and Playwright runs. The PK-PD cyclic model loads correctly, Multiaxis/Facet plot tabs... | 1m 2s | 38s | 36s | 2m 16s | +32s | +35s | +13s | +1m 20s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS → PASS | All 3 functional steps passed. The pk.ivp file opened from the Files browser as a DiffStudio view. Modifying Step (0.... | 2m 31s | 1m 32s | 1m 17s | 5m 20s | +1m 48s | +1m 29s | +53s | +4m 10s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | PASS → PARTIAL | Steps 1-5 pass; step 6 (actually running the fit) does not produce result rows within 2 minutes with recommended swit... | 10m 2s | 55s | 1m 0s | 11m 57s | +8m 32s | +50s | +32s | +9m 54s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS → PASS | All 6 steps passed. DiffStudio opened correctly, Bioreactor loaded from Library with 13 columns and 1001 rows. Both M... | 1m 31s | 36s | 26s | 2m 33s | +1m 6s | +33s | +8s | +1m 47s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | PASS → PASS | All 5 steps passed. The DiffStudio Edit toggle and `</>` export icon work correctly, producing a JS script with all m... | 4m 30s | 1m 40s | 1m 3s | 7m 13s | +3m 40s | +1m 37s | +14s | +5m 31s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS → PASS | All 4 steps passed. The Sensitivity Analysis view opened correctly from the ribbon. Process mode switching updated FF... | 2m 3s | 44s | 39s | 3m 26s | +1m 33s | +41s | +16s | +2m 30s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | PASS → PASS | All 4 steps passed. The Acid Production (GA-production) model loaded from Library with 2-stage simulation (1st stage... | 1m 49s | 47s | 24s | 3m 0s | +1m 24s | +44s | +7s | +2m 15s |
| EDA | [anova](EDA/anova-run.md) | PASS → PASS | All 3 steps passed successfully. The ANOVA dialog opened correctly with default parameters (Category=RACE, Feature=AG... | 1m 30s | 30s | 26s | 2m 26s | +1m 18s | +28s | +23s | +2m 9s |
| EDA | [ML methods/linear-regression](EDA/ML methods/linear-regression-run.md) | PASS → PASS | Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation... | 45s | 2.0s | 6.8s | 54s | +0s | +0.0s | +0s | +0s |
| EDA | [ML methods/pls-regression](EDA/ML methods/pls-regression-run.md) | PARTIAL → PARTIAL | PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-... | 1m 0s | 2.0s | 6.9s | 1m 9s | +0s | +0.0s | +0s | +0s |
| EDA | [ML methods/softmax](EDA/ML methods/softmax-run.md) | FAIL → FAIL | Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns an... | 10s | 2.0s | 2.6s | 15s | +0s | +0.0s | +0s | +0s |
| EDA | [ML methods/xgboost1](EDA/ML methods/xgboost1-run.md) | PARTIAL → PARTIAL | XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as... | 5.0s | 2.0s | 2.7s | 9.7s | +0s | +0.0s | +0s | +0s |
| EDA | [ML methods/xgboost2](EDA/ML methods/xgboost2-run.md) | PARTIAL → PARTIAL | XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interac... | 5.0s | 2.0s | 2.7s | 9.7s | +0s | +0.0s | +0s | +0s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PASS → PARTIAL | 2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check). Playwright spec is green — ... | 2m 30s | 2m 0s | 13s | 4m 43s | +2m 20s | +1m 58s | +6s | +4m 24s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL → PARTIAL | 6 of 7 steps passed (1 partial due to missing test file). The Pareto Front viewer works correctly: non-numeric column... | 4m 0s | 2m 0s | 32s | 6m 32s | +3m 40s | +1m 58s | +23s | +6m 1s |
| EDA | [pca](EDA/pca-run.md) | PASS → PARTIAL | MCP produced 3 PASS / 1 FAIL / 1 SKIP. The dialog path works but clicking OK closed the dialog without adding PC1/PC2... | 5m 0s | 3m 0s | 1m 7s | 9m 7s | +4m 47s | +2m 58s | +1m 4s | +8m 49s |
| EDA | [pls](EDA/pls-run.md) | PASS → FAIL | MCP produced 2 PASS / 1 PARTIAL / 1 FAIL. Dialog path works but PLS validation disables RUN when the Predict column i... | 2m 0s | 2m 0s | 15s | 4m 15s | +1m 50s | +1m 58s | +11s | +3m 59s |
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
| Peptides | [info-panels](Peptides/info-panels-run.md) | PASS → PASS | All 6 steps passed. The peptides.csv dataset loads correctly with Macromolecule semType detection. Amino acids are re... | 17s | 3.0s | 11s | 31s | +0s | +0.0s | +0s | +0s |
| Peptides | [peptide-space](Peptides/peptide-space-run.md) | PARTIAL → PARTIAL | SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variabil... | 25s | 3.0s | 1m 18s | 1m 46s | +0s | +0.0s | +0s | +0s |
| Peptides | [peptides](Peptides/peptides-run.md) | PARTIAL → PARTIAL | Steps 1-4 passed: peptides.csv loads correctly, the Context Panel shows the Peptides pane with Activity/Scaling/Clust... | 18s | 3.0s | 12s | 33s | +0s | +0.0s | +0s | +0s |
| Peptides | [sar](Peptides/sar-run.md) | PARTIAL → PARTIAL | Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Res... | 50s | 3.0s | 1m 30s | 2m 23s | +0s | +0.0s | +0s | +0s |
| PowerPack | [add-new-column](PowerPack/add-new-column-run.md) | PASS → PASS | All 5 steps passed on dev.datagrok.ai. The Add New Column dialog opens correctly, displays cleanly with no visual iss... | 22s | 2.0s | 5.3s | 29s | +0s | +0.0s | +0s | +0s |
| PowerPack | [AddNewColumn/add-new-column](PowerPack/AddNewColumn/add-new-column-run.md) | PARTIAL → PARTIAL | Core Add New Column functionality works correctly: creating calculated columns, chaining formulas, and propagating co... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/autocomplete](PowerPack/AddNewColumn/autocomplete-run.md) | PASS → PASS | All 8 steps passed. Autocomplete works for both functions (typing or Ctrl+Space) and columns ($ symbol). Function sel... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/formula-refreshing](PowerPack/AddNewColumn/formula-refreshing-run.md) | PASS → PASS | All 5 steps passed. The formula dependency chain (Weight2 -> Weight3 -> Weight4) works correctly. Changing a value in... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/functions_sorting](PowerPack/AddNewColumn/functions_sorting-run.md) | SKIP → SKIP | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on rel... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/highlight](PowerPack/AddNewColumn/highlight-run.md) | PASS → PASS | All 4 steps passed. Column names in both ${} (scalar) and $[] (aggregate) syntax are highlighted with the cm-column-n... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/hints](PowerPack/AddNewColumn/hints-run.md) | PASS → PASS | All 4 steps passed. Hovering over a function name in the formula editor displays a tooltip with the function signatur... |  |  |  |  | — | — | — | — |
| PowerPack | [AddNewColumn/input_functions](PowerPack/AddNewColumn/input_functions-run.md) | SKIP → SKIP | All steps skipped. This scenario requires the spgi.csv dataset which is not available as a standard demo table on rel... |  |  |  |  | — | — | — | — |
| PowerPack | [data-enrichment](PowerPack/data-enrichment-run.md) | PARTIAL → PARTIAL | Part 1 (create, save, apply enrichment) passed. The enrichment feature correctly joins the Customers table to Orders... | 42s | 3.0s |  | 45s | +0s | +0.0s | — | +0s |
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
| Scripts | layout | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Scripts | [run](Scripts/run-run.md) | PARTIAL → PARTIAL | Core run functionality works: the script can be triggered from context menu with a table selection and from the conso... |  |  |  |  | — | — | — | — |
| StickyMeta | [add-and-edit](StickyMeta/add-and-edit-run.md) | PASS → PASS | All tested steps passed. SPGI.csv opened with TestSchema1 sticky metadata schema pre-configured. The Sticky meta pane... | 35s | 3.0s | 17s | 55s | +0s | +0.0s | +0s | +0s |
| StickyMeta | [copy,-clone,-delete](StickyMeta/copy,-clone,-delete-run.md) | PASS → PASS | Steps 1-2 passed: SPGI.csv opened with TestSchema1 sticky metadata schema, and cloning the table preserves the schema... | 25s | 3.0s | 21s | 49s | +0s | +0.0s | +0s | +0s |
| StickyMeta | [create-schema-and-type](StickyMeta/create-schema-and-type-run.md) | PASS → PASS | All steps passed. The Sticky Meta Schemas browser at `/meta/schemas` shows 20 schemas including TestSchema1. The "NEW... | 10s | 3.0s | 7.0s | 20s | +0s | +0.0s | +0s | +0s |
| StickyMeta | [database-meta](StickyMeta/database-meta-run.md) | FAIL → FAIL | Step 1 passed: navigated to Databases > Postgres and found CHEMBL connection. Step 2 failed: the "Database meta" sect... | 25s | 0.0s |  | 25s | +0s | +0.0s | — | +0s |
| Tooltips | actions-in-the-context-menu | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | default-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | default-tooltip-visibility | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | edit-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | line-chart---aggregated-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | tooltip-properties | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Tooltips | uniform-default-tooltip | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS → PASS | All 53 steps passed against dev.datagrok.ai. All 3D Scatter Plot viewer properties set via `viewer.setOptions()`. Canv... | 36m 0s | 5m 0s | 45s | 41m 45s | +34m 5s | +4m 30s | +39s | +39m 14s |
| Viewers | [annotation-regions](Viewers/annotation-regions-run.md) | PASS → PASS |  | 7m 0s | 1m 0s | 17s | 8m 17s | +0s | +0s | +0s | +0s |
| Viewers | [bar-chart](Viewers/bar-chart-run.md) | PASS → PASS | All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color cod... | 3m 0s | 10s | 1m 0s | 4m 10s | +0s | +0s | +0s | +0s |
| Viewers | bar-chart-tests | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [box-plot](Viewers/box-plot-run.md) | PARTIAL → PARTIAL | 17 of 19 sections passed. Section 18 (Visualization zoom) is AMBIGUOUS: project save fails on dev server, and layout... | 3m 0s |  | 47s | 3m 47s | +0s | — | +0s | +0s |
| Viewers | [calendar](Viewers/calendar-run.md) | PASS → PASS | All 11 actions in the Calendar scenario passed on `dev.datagrok.ai`. The viewer correctly renders, tooltips and selec... | 25s | 1m 0s | 9.4s | 1m 34s | +0s | +0s | +0s | +0s |
| Viewers | [color-coding](Viewers/color-coding-run.md) | PASS → PASS | Core color coding steps passed: Linear (AGE, STARTED), Categorical (SEX, CONTROL) color coding applied and verified v... | 15s | 3.0s | 1m 24s | 1m 42s | +0s | +0.0s | +0s | +0s |
| Viewers | [color-coding-(linked)](Viewers/color-coding-(linked)-run.md) | PASS → PASS | All tested steps passed. Linked color coding works via column tags (`.color-coding-type`, `.color-coding-source-colum... | 10s | 3.0s | 6.0s | 19s | +0s | +0.0s | +0s | +0s |
| Viewers | [correlation-plot](Viewers/correlation-plot-run.md) | PARTIAL → PARTIAL | 27 of 30 steps passed, 3 skipped/ambiguous due to canvas-based cell interaction limitation. All property-based operat... | 5m 0s | 3.0s | 22s | 5m 26s | +0s | +0.0s | +0s | +0s |
| Viewers | [density-plot](Viewers/density-plot-run.md) | PASS → PASS | All 13 scenarios passed. The density plot viewer behaves correctly across all tested property combinations. UI intera... | 18m 0s | 2m 0s |  | 20m 0s | +0s | +0s | — | +0s |
| Viewers | [FilterPanel/basic-operations](Viewers/FilterPanel/basic-operations-run.md) | PASS → PASS | All filter panel basic operations work correctly on dev server. Structure/substructure filtering (32 rows for benzene... |  |  |  |  | — | — | — | — |
| Viewers | [FilterPanel/chem-and-bio](Viewers/FilterPanel/chem-and-bio-run.md) | PASS → PASS | All steps pass. All 6 chem search type modes produce correct row counts matching scenario expectations exactly. Bio s... |  |  |  |  | — | — | — | — |
| Viewers | [FilterPanel/cloned-views](Viewers/FilterPanel/cloned-views-run.md) | PASS → PASS | All 14 scenario steps passed on dev server. The cloned view correctly inherits filter state from the original view, i... |  |  |  |  | — | — | — | — |
| Viewers | [FilterPanel/collaborative-filtering-for-linked-tables](Viewers/FilterPanel/collaborative-filtering-for-linked-tables-run.md) | PASS → PASS | All 9 steps passed. Table linking (selection-to-filter and filter-to-filter) works correctly. Filtering propagates be... | 35s | 3.0s | 30s | 1m 8s | +0s | +0.0s | +0s | +0s |
| Viewers | [FilterPanel/combined-boolean-filter](Viewers/FilterPanel/combined-boolean-filter-run.md) | PASS → PASS | All 15 steps passed. The combined boolean filter correctly aggregates boolean columns, supports OR mode filtering, in... |  |  |  |  | — | — | — | — |
| Viewers | [FilterPanel/expression-filter](Viewers/FilterPanel/expression-filter-run.md) | PASS → PASS | Expression filter works correctly in all tested modes. Adding rules, switching AND/OR, removing rules, free-text mode... |  |  |  |  | — | — | — | — |
| Viewers | [FilterPanel/hierarchical-filter](Viewers/FilterPanel/hierarchical-filter-run.md) | PASS → PASS | All 12 steps passed. The hierarchical filter correctly supports multi-level tree filtering with expandable nodes, che... |  |  |  |  | — | — | — | — |
| Viewers | [FilterPanel/text-filter](Viewers/FilterPanel/text-filter-run.md) | PASS → PASS | All 9 steps passed. The text filter correctly searches string columns, highlights matches, supports multiple search t... |  |  |  |  | — | — | — | — |
| Viewers | [FilterPanel/viewers](Viewers/FilterPanel/viewers-run.md) | PASS → PASS | All 31 steps passed. Trellis Plot requires two clicks to apply filter (first selects cell, second applies), Esc to re... | 1m 25s | 5.0s | 1m 12s | 2m 42s | +0s | +0.0s | +0s | +0s |
| Viewers | [form](Viewers/form-run.md) | PASS → PASS | All 14 sections of form-tests-pw.md exercised across 30 steps. 28 PASS, 2 AMBIGUOUS, 0 FAIL in MCP run. Playwright spe... | 18m 0s | 4m 0s | 3m 12s | 25m 12s | +16m 30s | +3m 55s | +2m 37s | +23m 2s |
| Viewers | [forms](Viewers/forms-run.md) | new: PASS | All 15 scenario sections exercised; 36 steps total. 32 PASS, 0 FAIL in MCP run (4 used JS API fallback for canvas ele... | 18m 0s | 3m 0s | 51.5s | 21m 51.5s | 18m 0s (new) | 3m 0s (new) | 51.5s (new) | 21m 51.5s (new) |
| Viewers | [grid](Viewers/grid-run.md) | PARTIAL → PARTIAL | Grid tests ran 22 steps (spec softSteps); 17 passed outright and 5 AMBIGUOUS (Copy/Paste, Column Header Context Menu, ... | 11m 0s | 3m 0s | 1m 18s | 15m 18s | +5m 0s | +2m 0s | +50s | +7m 50s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | PARTIAL → PASS | All 14 heat-map sections exercised across 17 steps. 15 PASS, 1 AMBIGUOUS, 1 SKIP in MCP run. Playwright spec passed fu... | 18m 0s | 4m 0s | 48.9s | 22m 48.9s | +16m 10s | +3m 55s | +30.9s | +20m 35.9s |
| Viewers | [histogram](Viewers/histogram-run.md) | PARTIAL → PARTIAL | Most histogram property-based tests passed successfully. All property setters (bins, split, color, spline, appearance... | 2m 0s | 30s | 50s | 3m 20s | +0s | +0s | +0s | +0s |
| Viewers | [Legend/color-consistency](Viewers/Legend/color-consistency-run.md) | PASS → PASS | All 7 steps passed. Categorical color coding on "Stereo Category" propagates correctly across all viewer types (Histo... | 30s | 3.0s | 24s | 57s | +0s | +0.0s | +0s | +0s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | PASS → PASS | All legend filtering behaviors work correctly. Legends dynamically update when dataframe filters, in-viewer filters,... | 55s | 5.0s | 35.9s | 1m 36s | +0s | +0.0s | -0.1s | -0.1s |
| Viewers | [Legend/line-chart](Viewers/Legend/line-chart-run.md) | PASS → PASS | All 10 steps passed. The Line chart legend correctly reflects split categories with diverse colors, Multi Axis mode s... | 1m 30s | 3.0s | 47s | 2m 20s | +0s | +0.0s | +0s | +0s |
| Viewers | [Legend/scatterplot](Viewers/Legend/scatterplot-run.md) | PASS → PASS | All 5 scenarios passed (38/38 steps). The scatterplot legend correctly displays combined color+marker legends, update... | 45s | 3.0s | 50s | 1m 38s | +0s | +0.0s | +0s | +0s |
| Viewers | [Legend/structure-rendering](Viewers/Legend/structure-rendering-run.md) | PASS → PASS | All steps passed. The Core column (containing SMILES) renders molecule structure images in the scatterplot marker leg... | 35s | 3.0s | 26s | 1m 4s | +0s | +0.0s | +0s | +0s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | PASS → PASS | All legend visibility and positioning features work correctly. Legends display on viewers with color/split columns, c... | 1m 25s | 5.0s | 43.7s | 2m 13.7s | +0s | +0.0s | -0.3s | -0.3s |
| Viewers | [line-chart](Viewers/line-chart-run.md) | PASS → PASS | All 27 scenario sections passed on dev.datagrok.ai. The line chart viewer properties, context menu operations, layout... | 2m 0s | 5.0s | 2m 30s | 4m 35s | +0s | +0.0s | +0s | +0s |
| Viewers | [map](Viewers/map-run.md) | PASS → PASS | Core steps passed: Map viewer added to earthquakes.csv with auto-detected lat/lon, color/size columns set, marker siz... | 15s | 3.0s | 9.0s | 27s | +0s | +0.0s | +0s | +0s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | PASS → PARTIAL | Matrix Plot tests ran with 15 PASS, 3 AMBIGUOUS, 0 FAIL. The spec executed in 57.7s with all implemented steps passing... | 20m 0s | 3m 0s | 55.9s | 23m 55.9s | +19m 40s | +2m 57s | +44.9s | +23m 21.9s |
| Viewers | [network-diagram](Viewers/network-diagram-run.md) | PARTIAL → PARTIAL | Viewer core functionality (auto-pick, column rebinding, property-driven styling, filter integration, clean close) all... | 6m 45s | 3.0s | 7.4s | 6m 55s | +0s | +0.0s | +0s | +0s |
| Viewers | [pc-plot](Viewers/pc-plot-run.md) | PASS → PASS | All 36 steps across 12 sections passed. PC Plot properties, context menus, filter interaction, column management, den... | 8m 0s | 5.0s | 38s | 8m 43s | +0s | +0.0s | +0s | +0s |
| Viewers | [pie-chart](Viewers/pie-chart-run.md) | PASS → PASS | All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appea... | 4m 0s | 10s | 58s | 5m 8s | +0s | +0s | +0s | +0s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | PASS → PARTIAL | Pivot Table tests ran with 16 PASS, 2 AMBIGUOUS, 1 SKIP, 0 FAIL. The spec executed in 35.1s with all implemented steps... | 16m 0s | 3m 0s | 1m 12s | 20m 12s | +15m 15s | +2m 55s | +50s | +18m 50s |
| Viewers | rendering-structures-on-the-axes | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [row-source](Viewers/row-source-run.md) | PASS → PASS | All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with al... | 4m 0s | 5.0s | 1m 24s | 5m 29s | +0s | +0.0s | +0s | +0s |
| Viewers | [scatter-plot](Viewers/scatter-plot-run.md) | PASS → PASS | All 20 test sections passed on dev.datagrok.ai. Playwright spec executed in 1.1 min (headed mode). Context menu, rect... | 8m 0s | 5.0s | 1m 6s | 9m 11s | +0s | +0.0s | +0s | +0s |
| Viewers | scatter-plot-tests | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | statistics-viewer | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | new: PASS | 24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to manu... | 4m 0s | 3m 0s | 58s | 7m 58s | 4m 0s (new) | 3m 0s (new) | 58s (new) | 7m 58s (new) |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | PARTIAL → PARTIAL | All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via ... | 28m 0s | 4m 0s | 46s | 32m 46s | +27m 15s | +3m 58s | +38.6s | +31m 52s |
| Viewers | [trellis-plot](Viewers/trellis-plot-run.md) | PARTIAL → PARTIAL | Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions (bin clicks, range s... | 3m 0s | 30s | 1m 48s | 5m 18s | +0s | +0s | +0s | +0s |
| Viewers | viewers-docking | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | [word-cloud](Viewers/word-cloud-run.md) | PASS → PASS | All tested steps passed. Word cloud viewer added to SPGI with column set to "Primary Series Name". Properties include... | 20s | 3.0s | 13s | 36s | +0s | +0.0s | +0s | +0s |
| Viewers | word-cloud-tests | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | working-with-nan-&-infinity | NO RUN → NO RUN |  |  |  |  |  | — | — | — | — |
| Viewers | ~~grid-viewer~~ | removed: PASS |  |  |  |  |  | 4m 0s (removed) | 1m 0s (removed) | 58s (removed) | 5m 58s (removed) |
| Viewers | ~~tile~~ | removed: PASS |  |  |  |  |  | 5m 0s (removed) | 30s (removed) | 1m 6s (removed) | 6m 36s (removed) |

## Comparison with Previous Report

Delta from the previous `total-run.md` (dated 2026-04-20); signed with `+`/`-`, time deltas use the same format as the values.

### Totals

**Total**: Tests Δ **+0** · Run Δ **+0 (+0%)** · Playwright Δ **+2 (+1%)** · Status **PARTIAL → PARTIAL** · Browser Δ **+2m 16s** · Spec Gen Δ **+42s** · Spec Run Δ **+10s** · Total Δ **+3m 5s**

### By Folder

| Folder | Tests Δ | Run Δ | Playwright Δ | Status | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|--------|---------|-------|--------------|--------|-----------|------------|------------|---------|
| Apps | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Bio | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0.0s | — | +0s |
| Browse | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| Charts | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +22s | +26s | +11s | +59s |
| Chem | +0 | +0 (+0%) | +2 (+14%) | PARTIAL → PARTIAL | +37s | +29s | +6s | +1m 23s |
| Connections | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| DiffStudio | +0 | +0 (+0%) | +0 (+0%) | PASS → PARTIAL | +2m 29s | +52.5s | +18s | +3m 40s |
| EDA | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +1m 24s | +56s | +13.2s | +2m 32s |
| General | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| LocalCashing | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Models | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| Notebooks | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Peptides | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0.0s | +0s | +0s |
| PowerPack | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0.0s | +0s | +0s |
| Projects | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| Queries | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Scripts | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | — | — | — | — |
| StickyMeta | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +0s | +0.0s | +0s | +0s |
| Tooltips | +0 | +0 (+0%) | +0 (+0%) | NO DATA → NO DATA | — | — | — | — |
| Viewers | +0 | +0 (+0%) | +0 (+0%) | PARTIAL → PARTIAL | +4m 12s | +57s | +12s | +5m 12s |

### Per-Test Changes

| Folder | Test | Status | Browser Δ | Spec Gen Δ | Spec Run Δ | Total Δ |
|--------|------|--------|-----------|------------|------------|---------|
| Charts | [radar](Charts/radar-run.md) | PASS → PASS | -2m 44s | -15s | +22.2s | -2m 37s |
| Charts | [sunburst](Charts/sunburst-run.md) | PARTIAL → PARTIAL | +2m 10s | +57s | -8s | +3m 0s |
| Charts | [tree](Charts/tree-run.md) | PASS → PARTIAL | +1m 40s | +37s | +19s | +2m 36s |
| Chem | [activity-cliffs](Chem/activity-cliffs-run.md) | PASS → PASS | -20s | +17s | +3s | +0s |
| Chem | [Advanced/scaffold-tree](Chem/Advanced/scaffold-tree-run.md) | FAIL → PARTIAL | +10s | +22s (new) | +50.3s (new) | +1m 25s |
| Chem | [Advanced/scaffold-tree-functions](Chem/Advanced/scaffold-tree-functions-run.md) | FAIL → PASS | -45s | +27s | +40.2s (new) | +22s |
| Chem | [Advanced/similarity-search](Chem/Advanced/similarity-search-run.md) | PASS → PASS | -10s | +22s | +22.7s (new) | +35s |
| Chem | [Advanced/structure-filter](Chem/Advanced/structure-filter-run.md) | PASS → PASS | -10s | +22s | +0s | +12s |
| Chem | [calculate](Chem/calculate-run.md) | PASS → FAIL | +7m 35s | +57s | +14s | +8m 46s |
| Chem | [chemical-space](Chem/chemical-space-run.md) | PASS → PASS | -25s | +17s | +12.3s | +4s |
| Chem | [chemprop](Chem/chemprop-run.md) | PARTIAL → PARTIAL | -1m 0s | +27s | +18.3s (new) | -15s |
| Chem | [elemental-analysis](Chem/elemental-analysis-run.md) | PASS → PASS | +30s | +27s | +7.9s | +1m 4s |
| Chem | [filter-panel](Chem/filter-panel-run.md) | PARTIAL → PASS | +10s | +22s | -8.2s | +23s |
| Chem | [info-panels](Chem/info-panels-run.md) | PASS → PASS | +2m 30s | +57s | -8s | +3m 19s |
| Chem | [mmp](Chem/mmp-run.md) | FAIL → PASS | -1m 20s | +17s | +1m 12s (new) | +9s |
| Chem | [r-group-analysis](Chem/r-group-analysis-run.md) | PASS → PASS | +1m 35s | +42s | +20.7s | +2m 38s |
| Chem | [sketcher](Chem/sketcher-run.md) | PASS → PASS | +5s | +27s | +3.5s | +35s |
| DiffStudio | [catalog](DiffStudio/catalog-run.md) | PASS → PASS | +1m 8s | +34s | +6s | +1m 48s |
| DiffStudio | [cyclic-models](DiffStudio/cyclic-models-run.md) | PASS → PASS | +32s | +35s | +13s | +1m 20s |
| DiffStudio | [files-and-sharing](DiffStudio/files-and-sharing-run.md) | PASS → PASS | +1m 48s | +1m 29s | +53s | +4m 10s |
| DiffStudio | [fitting](DiffStudio/fitting-run.md) | PASS → PARTIAL | +8m 32s | +50s | +32s | +9m 54s |
| DiffStudio | [open-model](DiffStudio/open-model-run.md) | PASS → PASS | +1m 6s | +33s | +8s | +1m 47s |
| DiffStudio | [scripting](DiffStudio/scripting-run.md) | PASS → PASS | +3m 40s | +1m 37s | +14s | +5m 31s |
| DiffStudio | [sensitivity-analysis](DiffStudio/sensitivity-analysis-run.md) | PASS → PASS | +1m 33s | +41s | +16s | +2m 30s |
| DiffStudio | [stages](DiffStudio/stages-run.md) | PASS → PASS | +1m 24s | +44s | +7s | +2m 15s |
| EDA | [anova](EDA/anova-run.md) | PASS → PASS | +1m 18s | +28s | +23s | +2m 9s |
| EDA | [multivariate-analysis](EDA/multivariate-analysis-run.md) | PASS → PARTIAL | +2m 20s | +1m 58s | +6s | +4m 24s |
| EDA | [pareto-front-viewer](EDA/pareto-front-viewer-run.md) | PARTIAL → PARTIAL | +3m 40s | +1m 58s | +23s | +6m 1s |
| EDA | [pca](EDA/pca-run.md) | PASS → PARTIAL | +4m 47s | +2m 58s | +1m 4s | +8m 49s |
| EDA | [pls](EDA/pls-run.md) | PASS → FAIL | +1m 50s | +1m 58s | +11s | +3m 59s |
| Viewers | [3d-scatter-plot](Viewers/3d-scatter-plot-run.md) | PASS → PASS | +34m 5s | +4m 30s | +39s | +39m 14s |
| Viewers | [form](Viewers/form-run.md) | PASS → PASS | +16m 30s | +3m 55s | +2m 37s | +23m 2s |
| Viewers | [forms](Viewers/forms-run.md) | new: PASS | 18m 0s (new) | 3m 0s (new) | 51.5s (new) | 21m 51.5s (new) |
| Viewers | [grid](Viewers/grid-run.md) | PARTIAL → PARTIAL | +5m 0s | +2m 0s | +50s | +7m 50s |
| Viewers | [heatmap](Viewers/heatmap-run.md) | PARTIAL → PASS | +16m 10s | +3m 55s | +30.9s | +20m 35.9s |
| Viewers | [Legend/filtering](Viewers/Legend/filtering-run.md) | PASS → PASS | +0s | +0.0s | -0.1s | -0.1s |
| Viewers | [Legend/visibility-and-positioning](Viewers/Legend/visibility-and-positioning-run.md) | PASS → PASS | +0s | +0.0s | -0.3s | -0.3s |
| Viewers | [matrix-plot](Viewers/matrix-plot-run.md) | PASS → PARTIAL | +19m 40s | +2m 57s | +44.9s | +23m 21.9s |
| Viewers | [pivot-table](Viewers/pivot-table-run.md) | PASS → PARTIAL | +15m 15s | +2m 55s | +50s | +18m 50s |
| Viewers | [tile-viewer](Viewers/tile-viewer-run.md) | new: PASS | 4m 0s (new) | 3m 0s (new) | 58s (new) | 7m 58s (new) |
| Viewers | [tree-map-viewer](Viewers/tree-map-viewer-run.md) | PARTIAL → PARTIAL | +27m 15s | +3m 58s | +38.6s | +31m 52s |
| Viewers | ~~grid-viewer~~ | removed: PASS | 4m 0s (removed) | 1m 0s (removed) | 58s (removed) | 5m 58s (removed) |
| Viewers | ~~tile~~ | removed: PASS | 5m 0s (removed) | 30s (removed) | 1m 6s (removed) | 6m 36s (removed) |

## Release Readiness

**Verdict**: Conditionally ready

Run coverage holds at 74% (138/187). No folder is outright FAIL, but six new regressions landed this cycle: DiffStudio dropped from PASS to PARTIAL (fitting cannot produce fit rows in 2 min), and five previously-passing scenarios regressed to PARTIAL/FAIL — `Charts/tree`, `Chem/calculate` (FAIL — Chem top menu won't open on dev), `EDA/multivariate-analysis`, `EDA/pca`, `EDA/pls` (FAIL — RUN disabled silently), and `Viewers/matrix-plot`/`pivot-table` (canvas interaction limits). Offsetting gains: `Chem/mmp` and `Chem/Advanced/scaffold-tree-functions` moved FAIL → PASS, `Chem/filter-panel` and `Viewers/heatmap` moved PARTIAL → PASS, Chem Playwright coverage is now 100% (+2 specs), and two new Viewers scenarios (`forms`, `tile-viewer`) replace removed `grid-viewer` and `tile`. Mean wall-clock per scenario jumped because the new `Execute via grok-browser` metric includes model thinking; this is a measurement change, not a throughput regression.

### Blocking Issues
- Apps: run coverage 0/2
- General: run coverage 0/11
- LocalCashing: run coverage 0/1
- Notebooks: run coverage 0/4
- Queries: run coverage 0/14
- Tooltips: run coverage 0/7
- Chem: calculate regressed to FAIL (top menu broken on dev)
- EDA: pls regressed to FAIL (RUN silently disabled when Predict ∈ Using)
- DiffStudio: fitting regressed to PARTIAL (fit job produces no rows in 2 min on dev)
- Bio, Browse, Charts, Connections, EDA, Models, Peptides, PowerPack, Projects, Scripts, StickyMeta, Viewers: PARTIAL — see per-test rows for specifics
