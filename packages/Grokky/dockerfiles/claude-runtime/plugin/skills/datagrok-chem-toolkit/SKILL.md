---
name: datagrok-chem-toolkit
description: Cheminformatics toolkit catalog for small molecules. Covers RDKit/OpenChemLib calculations, descriptors, fingerprints, similarity/substructure search, R-group, MMP, scaffold tree, chemical space, MPO, reactions, and ADMET predictions. Open this skill when the user asks about computing molecular properties, comparing/searching molecules, or predicting ADMET/tox.
---

# datagrok-chem-toolkit

Function catalog for the **structure-and-modeling** half of Datagrok cheminformatics:
Chem, Admetica.

Every function below is callable from a `datagrok-exec` block via
`grok.functions.call('<Package>:<funcName>', {param1, param2, ...})`. For ones
exposed under a top-menu path, the path is shown so you can describe it in plain
language ("Top menu → Chem → ...").

### Param-type conventions used below

Types match the Datagrok function declarations: `dataframe`, `column`, `string`,
`int`, `float`, `bool`, `list<T>`, `object`. A column with a required semType is
written `column<semType>` (e.g. `column<Molecule>`, `column<Macromolecule>`,
`column<numerical>`). Choices for enum-style params are listed inline as
`'a' | 'b' | 'c'`. `?` marks optional, `= x` shows the default.

### Always emit an exec block when the intent is clear

If the user asks for a concrete chemistry action ("generate analogs", "dock these
ligands", "compute IC50", "register this compound", "what is the SMILES of
aspirin") — emit a `datagrok-exec` block calling the relevant function below,
even if some args have to be guessed from context (table columns, sensible
defaults). Don't ask the user to re-state inputs you can already infer; if a
required input is genuinely missing, still emit the block with the best guess
and add a one-line note. Don't answer chemistry questions from memory when a
catalog function exists — call the function.

---

## Cross-package routing — "User wants X → call Y"

| User intent | Function | Package |
|---|---|---|
| Draw / depict a molecule | `Chem:drawMolecule` | Chem |
| Convert SMILES ↔ Molfile ↔ SMARTS ↔ InChI | `Chem:convertMolNotation` | Chem |
| Get InChI / InChI key | `Chem:getInchis`, `Chem:getInchiKeys` | Chem |
| Molecular formula | `Chem:getMolecularFormula` | Chem |
| Canonical SMILES | `Chem:canonicalize` | Chem |
| Validate SMILES (RDKit) | `Chem:validateMolecule`, `Chem:isSmiles`, `Chem:isSmarts` | Chem |
| Basic properties (MW, LogP, HBA, HBD, PSA, ...) | `Chem:addChemPropertiesColumns` or `Chem:getProperties` | Chem |
| Single-molecule cLogP only | `Chem:getCLogP` | Chem |
| Chemical descriptors (RDKit set) | `Chem:getDescriptors` / `Chem:chemDescriptors` | Chem |
| Toxicity risks (OCL: mutagenicity, tumorigenicity, ...) | `Chem:addChemRisksColumns` or `Chem:getToxicityRisks` | Chem |
| Structural alerts (PAINS, BMS, SureChEMBL, ...) | `Chem:structuralAlertsTopMenu` or `Chem:getStructuralAlerts` | Chem |
| Pharmacophore features (donor, acceptor, ...) | `Chem:pharmacophoreFeaturesTopMenu` | Chem |
| Morgan / RDKit / MACCS / AtomPair fingerprints | `Chem:getMorganFingerprints`, `Chem:getFingerprints` | Chem |
| Find similar molecules in a column | `Chem:findSimilar`, `Chem:callChemSimilaritySearch` | Chem |
| Tanimoto similarity matrix | `Chem:similarityMatrixTopMenu` | Chem |
| Substructure search | `Chem:searchSubstructure`, `Chem:SubstructureSearchTopMenu` | Chem |
| Chemical space 2D (UMAP / t-SNE) | `Chem:chemSpaceTopMenu` | Chem |
| Diversity picking | `Chem:callChemDiversitySearch` | Chem |
| BitBIRCH O(N) clustering | `Chem:bitbirchClusteringTopMenu` | Chem |
| Cluster MCS (most-common-substructure per cluster) | `Chem:performClusterMCS` | Chem |
| R-Group decomposition | `Chem:rGroupsAnalysisMenu`, `Chem:rGroupDecomposition` | Chem |
| Matched Molecular Pairs | `Chem:mmpAnalysis` | Chem |
| Activity cliffs | `Chem:activityCliffs` | Chem |
| Scaffold tree (hierarchical) | `Chem:getScaffoldTree`, `Chem:addScaffoldTree` | Chem |
| Elemental analysis (atom counts) | `Chem:elementalAnalysis` | Chem |
| Synthon (REAL-like) space search | `Chem:synthonSearchFunc` | Chem |
| Drug name(s) → SMILES (ibuprofen, aspirin, ...) | `Chem:namesToSmiles` (ChEMBL-backed; **needs Chembl pkg**) | Chem |
| Free text → SMILES (LLM fallback) | `Chem:freeTextToSmiles` | Chem |
| Map IDs incl. InChI key → ChEMBL, ChEMBL → DrugBank, etc. | `Chem:getMapIdentifiers` / `Chem:mapIdentifiersTransform` | Chem |
| MPO score (weighted desirability) | `Chem:mpoCalculate`, `Chem:_mpo` | Chem |
| Apply reaction SMARTS (transformation) | `Chem:transformationReactionsTopMenu` | Chem |
| Two-column reaction enumeration | `Chem:twoComponentReactionTopMenu` | Chem |
| Remove water and salts | `Chem:removeWaterAndSaltsTopMenu` | Chem |
| Deprotect (remove protecting groups) | `Chem:deprotect` | Chem |
| Recalculate 2D coordinates | `Chem:recalculateCoords` | Chem |
| Train Chemprop QSAR model | `Chem:trainChemprop` | Chem |
| Apply Chemprop QSAR model | `Chem:applyChemprop` | Chem |
| ADMET predictions (full panel) | `Admetica:getAdmeProperties` | Admetica |
| ADMET for single molecule | `Admetica:getAdmePropertiesSingle` | Admetica |
| ADMET in Hit Triage pipeline | `Admetica:admeticaHT` | Admetica |

---

## Calling conventions for column-based functions

Most Chem functions operate on a `DG.Column` of molecules and a parent
`DG.DataFrame`. Inside a `datagrok-exec` block:

```js
const mol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);   // first Molecule column
const res = await grok.functions.call('Chem:getProperties', { molecules: mol });
// res is a DG.DataFrame with one column per requested property
```

Function tables below carry a **Tags** column. The values it can hold:

| Tag             | Meaning                                                                                                                                                                |
|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `vectorFunc`    | Returns a `DG.DataFrame`. Usable from the Add-New-Column dialog and inside a `grokky.addCalculatedColumn` formula string. Functions WITHOUT this tag cannot be embedded in a formula — call them via `grok.functions.call`. |
| `join(table)`   | Output declares `action: 'join(table)'` — result columns are auto-appended to the input table. Do not also append yourself.                                            |
| `transform`     | Tagged as a data-sync transform — same behavior as the un-tagged sibling, but tracked in the project history.                                                          |
| `topMenu`       | Has a Top-menu path; otherwise programmatic only.                                                                                                                      |
| `panel`         | Context-panel widget — returns `DG.Widget`, shown in the right side panel.                                                                                             |
| `action`        | Right-click / context-menu action.                                                                                                                                     |

If the molecule notation is unknown, `sourceNotation: DG.chem.Notation.Unknown`
works for SMILES, SMARTS, molblock, and V3K molblock inputs.

If the molecule notation is unknown, `sourceNotation: DG.chem.Notation.Unknown`
works for SMILES, SMARTS, molblock, and V3K molblock inputs.

---

## Chem package — RDKit/OpenChemLib core (`Chem:...`)

The flagship cheminformatics package. RDKit is the engine; OpenChemLib provides
fallback rendering and a few descriptors (logP, logS, drug-likeness, toxicity).

### Four pitfalls to avoid

1. **Do not invent `Chem:CalculateMW`, `Chem:CalculateLogP`, `Chem:CalculateLogS`,
   `Chem:CalculateHBA`, etc.** No such standalone functions exist. MW / LogP /
   LogS / HBA / HBD / PSA / Rotatable bonds / Stereo centers / Molecule charge
   are **boolean flags** on `Chem:addChemPropertiesColumns(...)` and **strings in
   the `selected` list** of `Chem:getProperties(molecules, selected?)`. There IS
   a Python script registered as `Chem:CalculateLogS`, but its own description
   says `DO NOT USE IN FUNCTION PLANNING OR CHAINING` — treat as non-public.
2. **`Chem:convertMolNotation` has no `inchi`.** The `sourceNotation` /
   `targetNotation` enum is `'smiles' | 'cxsmiles' | 'smarts' | 'cxsmarts' |
   'molblock' | 'v3Kmolblock'`. For InChI use `Chem:getInchis` /
   `Chem:getInchiKeys`. For InChI key → ChEMBL ID, use
   `Chem:getMapIdentifiers` / `Chem:mapIdentifiersTransform` (sources
   include `inchi_key`, `chembl`, `pubchem`, `drugbank`, `zinc`, ...).
3. **`Chem:namesToSmiles` is ChEMBL-backed.** It delegates to
   `Chembl:namesToSmiles` and will fail if the Chembl package / DB isn't
   deployed. Still pick it for any "drug name → SMILES" / list-of-names task
   (aspirin, ibuprofen, caffeine, ...). Don't fall back to `Chem:freeTextToSmiles`
   (LLM path, for free-form text fragments). Don't answer from memory.
4. **Docking / generation must route to the right package.**
   "Dock these ligands" → `Docking:getAutodockResults` / `Docking:runAutodock`
   or `Boltz1:docking` — **never** to a `Chem:` function. "Generate analogs" /
   "design new molecules" → `Reinvent4:reinvent` / `Reinvent4:reinventTopMenu`
   — **not** R-group, MMP, or reaction enumeration (those operate on existing
   series). Always emit a `datagrok-exec` block; if a target/seed is missing,
   pass a sensible default (first molecule column, first configured target).

### Rendering & notation

Notation enum (used by `convertMolNotation`, `convertMoleculeNotation`,
`convertNotation`): `'smiles' | 'cxsmiles' | 'smarts' | 'cxsmarts' | 'molblock' | 'v3Kmolblock'`. **No `inchi`** — use `getInchis` / `getInchiKeys` instead.

| Function | Tags | What it does |
|---|---|---|
| `drawMolecule(molStr: string, w?: int, h?: int, popupMenu?: bool)` | — | Returns an `HTMLElement` of a rendered molecule. |
| `canvasMol(x: int, y: int, w: int, h: int, canvas: object, molString: string, scaffoldMolString?: string, options?: object, renderingOptions?: object)` | — | Draws into an existing canvas (with optional scaffold highlight). |
| `convertMolNotation(molecule: string, sourceNotation: <notation>, targetNotation: <notation>)` | — | Returns a single string in the target notation. |
| `convertMoleculeNotation(molecule: column<Molecule>, targetNotation: <notation>, kekulize?: bool = false)` | `vectorFunc` | Notation conversion. Returns a new `DG.Column`. |
| `convertNotation(data: dataframe, molecules: column<Molecule>, targetNotation: <notation>, overwrite?: bool, join?: bool, kekulize?: bool)` | `topMenu` | Top-menu version: appends or overwrites the column. Top menu: Chem → Transform → Convert Notation. |
| `recalculateCoords(table: dataframe, molecules: column<Molecule>, method: 'OCL' \| 'CoordGen' = 'OCL', join: bool = true)` | `topMenu` | Recompute 2D layout. Top menu: Chem → Transform → Recalculate Coordinates. |
| `canonicalize(molecule: string)` | — | → canonical SMILES (string). |
| `getInchis(molecules: column<Molecule>)` / `getInchiKeys(molecules: column<Molecule>)` | `vectorFunc` | Vector InChI / InChI key. |
| `addInchisTopMenu(table: dataframe, molecules: column<Molecule>)` / `addInchisKeysTopMenu(...)` | `topMenu` | Append as a new column. Top menu: Chem → Calculate → To InchI / To InchI Keys. |
| `getMolecularFormula(molecule: string)` | — | Hill formula string (via OCL). |
| `isSmiles(s: string)` / `isSmarts(s: string)` / `validateMolecule(s: string)` | — | RDKit-based validators. `validateMolecule` returns the RDKit error string or `null` if OK. |
| `detectSmiles(col: column, min: int)` | — | Heuristic — tags the column as `SMILES` / `Molecule` if enough cells parse. |
| `getMolFileHandler(molString)` | — | Returns a low-level molfile parser handler. |
| `getRdKitModule()` | — | Returns the raw RDKit-JS module (for advanced calls like `mol.get_descriptors()`). |
| `chemCellRenderer` / `rdKitCellRenderer` / `oclCellRenderer` / `rdKitReactionRenderer` / `rdKitMixtureRenderer` | — | Cell renderers — call only via the grid system, not directly. |
| `editMoleculeCell(cell)` | — | Pops a sketcher dialog bound to a grid cell. |
| `openChemLibSketcher()` | — | Returns the OCL sketcher widget. |

### Properties, descriptors, fingerprints

**Property names** for `getProperties(molecules, selected)` and equivalent flags
on `addChemPropertiesColumns`: `'MW' | 'HBA' | 'HBD' | 'LogP' | 'LogS' | 'PSA' |
'Rotatable bonds' | 'Stereo centers' | 'Molecule charge'`. Pass them in
`selected` as a list of strings, exactly as shown.

**Fingerprint enum**: `'Morgan' | 'RDKit' | 'Pattern' | 'AtomPair' | 'MACCS' | 'TopologicalTorsion'`.

**Toxicity risk names** for `getToxicityRisks(molecules, risks)`: `'mutagenicity' | 'tumorigenicity' | 'irritatingEffects' | 'reproductiveEffects'`.

**Structural alert rule sets** for `getStructuralAlerts(molecules, alerts)` and the per-flag args of `structuralAlertsTopMenu`: `'PAINS' | 'BMS' | 'SureChEMBL' | 'MLSMR' | 'Dundee' | 'Inpharmatica' | 'LINT' | 'Glaxo'`.

| Function | Tags | What it does |
|---|---|---|
| `addChemPropertiesColumns(table: dataframe, molecules: column<Molecule>, MW?: bool = true, HBA?: bool = true, HBD?: bool = true, logP?: bool = true, logS?: bool = true, PSA?: bool = true, rotatableBonds?: bool = true, stereoCenters?: bool = true, moleculeCharge?: bool = false)` | `topMenu` | OCL-based properties → appended as columns. Top menu: Chem → Calculate → Chemical Properties. |
| `getProperties(molecules: column<Molecule>, selected?: list<string>)` | `vectorFunc`, `join(table)` | Vector form. `selected` is a subset of the property names above; omit/empty for all. Returns a DataFrame. |
| `getCLogP(smiles: string)` | — | Single-mol Crippen logP (RDKit). |
| `getChemPropertyFunction(name: string)` | — | Returns a `(smiles) => any` closure for one property — useful for `column.applyFormula`. |
| `getDescriptors(molecules: column<Molecule>, selected?: list<string>)` | `vectorFunc`, `join(table)` | RDKit descriptor set (via Chem docker). Returns DataFrame. Use `chemDescriptorsTree()` for the legal `selected` names. |
| `chemDescriptors(table: dataframe, molecules: column<Molecule>, descriptors: list<string>)` | `topMenu` | Same, appended to `table`. |
| `chemDescriptorsTree()` | — | Returns the tree of available descriptor groups (for UI). |
| `descriptorsDocker()` | `topMenu` | Opens the descriptors picker dialog. Top menu: Chem → Calculate → Descriptors. |
| `calculateDescriptorsTransform(table: dataframe, molecules: column<Molecule>, selected: list<string>)` | `transform` | Same but tagged `transform` for data-sync projects. |
| `getMorganFingerprints(molColumn: column<Molecule>)` | `vectorFunc` | Returns a Column of `DG.BitSet` fingerprints. |
| `getMorganFingerprint(molString: string)` | — | Single-mol → `DG.BitSet`. |
| `getFingerprints(col: column<Molecule>, _metric?: string, fingerprintType?: <fp>)` | — | Returns `{entries, options}` for the dim-reduction preprocessor. |
| `addChemRisksColumns(table: dataframe, molecules: column<Molecule>, mutagenicity?: bool = true, tumorigenicity?: bool = false, irritatingEffects?: bool = false, reproductiveEffects?: bool = false)` | `topMenu` | OCL toxicity risk flags. Top menu: Chem → Calculate → Toxicity Risks. |
| `getToxicityRisks(molecules: column<Molecule>, risks?: list<string>)` | `vectorFunc`, `join(table)` | Vector form. `risks` is a subset of the risk names above. |
| `structuralAlertsTopMenu(table: dataframe, molecules: column<Molecule>, pains: bool = true, bms: bool = false, sureChembl: bool = false, mlsmr: bool = false, dundee: bool = false, inpharmatica: bool = false, lint: bool = false, glaxo: bool = false)` | `topMenu` | Apply rule-based structural alerts. Top menu: Chem → Analyze → Structural Alerts. |
| `getStructuralAlerts(molecules: column<Molecule>, alerts?: list<string>)` | `vectorFunc`, `join(table)` | Vector form. |
| `pharmacophoreFeaturesTopMenu(table: dataframe, molecules: column<Molecule>, donor: bool, acceptor: bool, hydrophobic: bool, aromatic: bool, positive: bool, negative: bool, halogenBond: bool)` | `topMenu` | RDKit pharmacophore family flags. Top menu: Chem → Analyze → Pharmacophore Features. |
| `biochemPropsWidget()` | `topMenu` | Opens the auto-discovery dialog for biochem calculators. Top menu: Chem → Calculate → Biochemical Properties. |

### Similarity, diversity, substructure search

**BitArray similarity metrics** (used by `callChemSimilaritySearch`, `callChemDiversitySearch`, `chemSpaceTopMenu`, `activityCliffs`):
`'Tanimoto' | 'Dice' | 'Asymmetric' | 'Braun-Blanquet' | 'Cosine' | 'Kulczynski' | 'Mc-Connaughey' | 'Rogot-Goldberg' | 'Russel' | 'Sokal' | 'Hamming' | 'Euclidean BitArray'`. UI choosers usually narrow this to `'Tanimoto' | 'Asymmetric' | 'Cosine' | 'Sokal'`.

| Function | Tags | What it does |
|---|---|---|
| `findSimilar(molStringsColumn: column<Molecule>, molString: string, limit?: int = MAX, cutoff?: float = 0.0)` | — | Returns a DataFrame of the most-similar molecules (Morgan/Tanimoto). |
| `callChemSimilaritySearch(df: dataframe, col: column<Molecule>, molecule: string, metricName: <metric>, fingerprint: <fp>, limit: int, minScore: float)` | — | Underlying version with explicit metric and FP type. |
| `getSimilarities(molStringsColumn: column<Molecule>, molString: string)` | — | Returns just the similarity scores as a DataFrame. |
| `getDiversities(molStringsColumn: column<Molecule>, limit?: int)` | — | DataFrame of the most-diverse molecules. |
| `callChemDiversitySearch(col: column<Molecule>, metricName: <metric>, fingerprint: <fp>, limit: int)` | — | Returns indices of diverse picks. |
| `similarityMatrixTopMenu(table: dataframe, molecules: column<Molecule>, symbols: column, fingerprintType?: <fp>)` | `topMenu` | Full pairwise Tanimoto matrix → opens as a new table. Top menu: Chem → Calculate → Similarity Matrix. |
| `searchSubstructure(molStringsColumn: column<Molecule>, molString: string, molBlockFailover: string)` | — | RDKit substructure search → returns a Column wrapping a BitSet of matches. |
| `SubstructureSearchTopMenu(molecules: column<Molecule>)` | `topMenu` | Opens the filter sketcher. Top menu: Chem → Search → Substructure Search. |
| `similaritySearchViewer()` / `diversitySearchViewer()` | — | Returns a `ChemSimilarityViewer` / `ChemDiversityViewer`. Add via `view.addViewer('Chem Similarity Search')`. |
| `similaritySearchTopMenu()` / `diversitySearchTopMenu()` | `topMenu` | Top menu shortcuts. Top menu: Chem → Search. |
| `sortBySimilarity(value: object)` | `action` | Right-click action that sorts the grid by similarity to the picked molecule. |
| `useAsSubstructureFilter(value: object)` | `action` | Right-click action that adds the picked molecule as a substructure filter. |
| `synthonSearchFunc(spaceName: string, molecule: string, maxHits: int = 100, searchType: 'substructure' \| 'similarity' \| 'exact', similarityCutoff?: float = 0.5, includeSynthons?: bool = false)` | — | Search in a synthon (REAL-style) chemical space. |
| `getSynthonSpacesFunc()` | — | Lists installed synthon spaces (filenames in `Chem/files/synthon-data/`). |
| `filterMoleculeDuplicates(molecules: list<string>, molecule: string)` | — | Removes duplicates of `molecule` from a list. |

### Chemical space, R-groups, MMP, activity cliffs, scaffold tree

**Dim-reduction method enum**: `'UMAP' | 't-SNE'`.
**R-group matching strategy**: `'Greedy' | 'GreedyChunks' | 'Exhaustive' | 'NoSymmetrization' | 'GA'`.
**MMP diff type**: `'delta' | 'ratio'`. **MMP activity scaling**: `'none' | 'lg' | '-lg'`.

| Function | Tags | What it does |
|---|---|---|
| `chemSpaceTopMenu(table: dataframe, molecules: column<Molecule>, methodName: 'UMAP' \| 't-SNE', similarityMetric: 'Tanimoto' \| 'Asymmetric' \| 'Cosine' \| 'Sokal' = 'Tanimoto', plotEmbeddings: bool = true, options?: object, preprocessingFunction?: func, clusterEmbeddings?: bool, clusterMCS?: bool)` | `topMenu` | Projects molecules to 2D (UMAP / t-SNE) + scatter plot. Top menu: Chem → Analyze → Chemical Space. |
| `chemSpaceTransform(...)` | `transform` | Transform-tagged version (for data-sync). |
| `getChemSpaceEmbeddings(col: column<Molecule>, methodName: 'UMAP' \| 't-SNE', similarityMetric: <metric>, xAxis: string, yAxis: string, options?: object)` | — | Raw embedding without UI. Returns `ISequenceSpaceResult`. |
| `getChemSimilaritiesMatrix(dim: int, col: column<Molecule>, df: dataframe, colName: string, simArr: list)` | — | Internal — pairwise similarity matrix as columns. |
| `rGroupsAnalysisMenu()` | `topMenu` | Opens the R-group analysis dialog. Top menu: Chem → Analyze → R-Groups Analysis. |
| `rGroupDecomposition(df: dataframe, molColName: string, core: string, rGroupName: string, rGroupMatchingStrategy: 'Greedy' \| 'GreedyChunks' \| 'Exhaustive' \| 'NoSymmetrization' \| 'GA', onlyMatchAtRGroups?: bool)` | — | Programmatic R-group decomposition. Returns `RGroupDecompRes`. |
| `mmpAnalysis(table: dataframe, molecules: column<Molecule>, activities: list<column<numerical>>, diffTypes: list<'delta'\|'ratio'>, scalings: list<'none'\|'lg'\|'-lg'>, fragmentCutoff?: float = 0.4)` | `topMenu` | Matched Molecular Pairs analysis. Top menu: Chem → Analyze → Matched Molecular Pairs. |
| `mmpViewer()` | — | Returns `MatchedMolecularPairsViewer` for manual `view.addViewer`. |
| `activityCliffs(table: dataframe, molecules: column<Molecule>, activities: column<numerical>, similarity: float, methodName: 'UMAP' \| 't-SNE', similarityMetric: <metric>, preprocessingFunction?: func, options?: object, isDemo?: bool, isTest?: bool)` | `topMenu` | Finds pairs of similar molecules with large activity differences. Top menu: Chem → Analyze → Activity Cliffs. |
| `addScaffoldTree()` | `topMenu` | Adds a `ScaffoldTreeViewer` to the current view. Top menu: Chem → Analyze → Scaffold Tree. |
| `scaffoldTreeViewer()` | — | Returns `ScaffoldTreeViewer` for manual use. |
| `scaffoldTreeFilter()` | — | Returns `ScaffoldTreeFilter`. |
| `getScaffoldTree(data: dataframe, ringCutoff?: int, dischargeAndDeradicalize?: bool)` | — | Returns the scaffold tree as a JSON string. |
| `substructureFilter()` | — | Returns an RDKit-based `SubstructureFilter`. |
| `bitbirchClusteringTopMenu(table: dataframe, molecules: column<Molecule>, threshold?: float, fingerprintType?: <fp>)` | `topMenu` | O(N) BitBIRCH clustering → cluster ID column. Top menu: Chem → Calculate → BitBIRCH Clustering. |
| `clusterMCSTopMenu(table: dataframe, molCol: column<Molecule>, clusterCol: column)` | `topMenu` | Most-common substructure per cluster → appends MCS column. Top menu: Chem → Calculate → Cluster MCS. |
| `performClusterMCS(molCol: column<Molecule>, clusterCol: column)` | `vectorFunc` | Vector form. |
| `elementalAnalysis(table: dataframe, molecules: column<Molecule>, radarViewer: bool, radarGrid: bool)` | `topMenu` | Adds atom-count columns; optional radar viewer. Top menu: Chem → Analyze → Elemental Analysis. |
| `runElementalAnalysis(table: dataframe, molecules: column<Molecule>)` | `transform` | Returns the list of added column names. |

### MPO, identifiers, reactions

**MPO aggregation enum** (`WEIGHTED_AGGREGATIONS_LIST`): `'Average' | 'Sum' | 'Product' | 'Geomean' | 'Min' | 'Max'`.

**Identifier source enum** for `mapIdentifiersTransform(fromSource, toSource)` — non-exhaustive: `'smiles' | 'inchi' | 'inchi_key' | 'chembl' | 'pubchem' | 'drugbank' | 'pdb' | 'zinc' | 'chebi' | 'kegg_ligand' | 'hmdb' | 'drugcentral' | 'surechembl' | 'bindingdb' | 'fdasrs' | 'lipidmaps' | 'mcule' | 'gtopdb'` (full list ≈ 35 UniChem sources). The catalog's standard reverse-lookup path (e.g. **InChI key → ChEMBL ID**) is `mapIdentifiersTransform(table, molecules, 'inchi_key', 'chembl')`.

| Function | Tags | What it does |
|---|---|---|
| `_mpo()` | `topMenu` | Opens the MPO Score dialog. Top menu: Chem → Calculate → MPO Score. |
| `mpoCalculate(df: dataframe, columns: list<column>, profileName: string, aggregation: 'Average' \| 'Sum' \| 'Product' \| 'Geomean' \| 'Min' \| 'Max', createDesirabilityColumns?: bool = false)` | — | Computes MPO score + optional per-property desirability columns. |
| `mpoTransformFunction(df: dataframe, profileName: string, aggregation: <agg>, currentProperties: string, silent?: bool)` | `transform` | Wrapper used by data-sync. |
| `mpoProfilesApp(path?: string)` | — | Opens the MPO Profiles app. |
| `getMapIdentifiers()` | `topMenu` | Opens the ID-mapping dialog. Top menu: Chem → Calculate → Map Identifiers. |
| `mapIdentifiersTransform(table: dataframe, molecules: column<Molecule>, fromSource: <id-source>, toSource: <id-source>)` | `transform` | Appends a column with mapped IDs. **This is also the right entry for "InChI key → ChEMBL ID" and similar reverse lookups** (chemblIdToSmilesTs goes only ID→SMILES, never the reverse). |
| `freeTextToSmiles(molfile: string)` | — | LLM-backed: parses free text into SMILES. Returns `string \| null`. Use only when names aren't clean enough for `namesToSmiles`. |
| `namesToSmiles(data: dataframe, names: column<string>)` | `topMenu` | Looks up compound names → canonical SMILES via ChEMBL. **Requires the Chembl package + ChEMBL DB.** Use for any "drug name → SMILES" task (aspirin, ibuprofen, caffeine, ...) — wrap the list of names in a `DG.Column.fromStrings(...)` and call. Top menu: Chem → Transform → Names To Smiles. |
| `deprotect(table: dataframe, molecules: column<Molecule>, fragment: string)` | `topMenu` | Removes the drawn fragment as a protecting group. Top menu: Chem → Transform → Reactions → Deprotect. |
| `removeWaterAndSaltsTopMenu(table: dataframe, molecules: column<Molecule>)` | `topMenu` | Strip water/salts. Top menu: Chem → Transform → Reactions → Remove Water and Salts. |
| `transformationReactionsTopMenu()` | `topMenu` | Run reaction SMARTS over a column of reactants. Top menu: Chem → Transform → Reactions → Transformation. |
| `twoComponentReactionTopMenu()` | `topMenu` | Cross-product reaction between two columns. Top menu: Chem → Transform → Reactions → Two-Component Reaction. |
| `transformationReactionsApp() / twoComponentReactionsApp()` | — | App entry points (browsePath `Chem | Reactions`). |
| `reactionEnumeratorApp()` | — | Forward-reaction library enumeration over building blocks + SMARTS. |
| `beautifyMols(mols: list<string>)` | `vectorFunc` | Returns clean V3K molblocks. |
| `convertToV3KViaOCL(mols: list<string>)` | — | OCL-based V2K → V3K conversion. |
| `convertMolNotationAction(col: column<Molecule>)` | `action` | Right-click action wired to `convertNotation`. |
| `convertMixtureToSmiles(col: column)` | `action` | Right-click action for `ChemicalMixture` columns. |
| `copyAsAction(value: object)` / `copyAsSmiles` / `copyAsMolfileV2000` / `copyAsMolfileV3000` / `copyAsSmarts` / `copyAsImage` | `action` | Right-click "Copy as ..." actions. |

### File I/O

| Function | Tags | What it does |
|---|---|---|
| `importSdf(bytes)` | — | SDF / MOL → DataFrames. Registered as fileHandler for `.sdf`, `.mol`. |
| `importSmi(bytes)` | — | SMI file handler. |
| `importMol2(bytes)` | — | Tripos MOL2 file handler. |
| `importMol(content)` | — | Single MOL file (string). |
| `saveAsSdf()` | — | Export current table as SDF. File exporter: "As SDF…". |

### Chemprop QSAR (deep learning)

| Function | Tags | What it does |
|---|---|---|
| `trainChemprop(df: dataframe, predictColumn: column<numerical>, dataset_type: 'regression' \| 'classification' = 'regression', /* +22 more snake_case args */)` | — | Trains a Chemprop MPN model via the `chem-chemprop` Docker container. Returns the model blob. |
| `applyChemprop(df: dataframe, model: object)` | — | Applies a saved Chemprop model → DataFrame with `outcome` column. |
| `isApplicableNN(df: dataframe, predictColumn: column)` / `isInteractiveNN(df: dataframe, predictColumn: column)` | — | ML framework hooks (not called directly). |

`trainChemprop` has 25 **snake_case** params total (`dataset_type`, `metric`,
`multiclass_num_classes`, `num_folds`, `data_seed`, `split_sizes`, `split_type`,
`activation`, `atom_messages`, `message_bias`, `ensemble_size`,
`message_hidden_dim`, `depth`, `dropout`, `ffn_hidden_dim`, `ffn_num_layers`,
`epochs` (default `50`), `batch_size` (default `64`), `warmup_epochs`,
`init_lr`, `max_lr`, `final_lr`, `no_descriptor_scaling`, ...). Pass at least
`predictColumn`, `dataset_type`, and `epochs`. For the full schema with every
enum (`metric`, `split_type`, `activation`), use the MCP
`get_function('Chem:trainChemprop')` call before invoking — don't guess names.

### Panels (context-panel widgets — show in the right-side panel when a molecule cell is selected)

Useful to call from a `datagrok-exec` block when you want to render a panel
inline. All return `DG.Widget`.

All entries below have the `panel` tag.

| Panel | What |
|---|---|
| `properties(smiles)` | "Chemistry \| Properties" |
| `descriptorsWidget(smiles)` | "Chemistry \| Descriptors" |
| `structure2d(molecule)` / `structure3D(molecule)` | "Structure \| 2D Structure" / "Structure \| 3D Structure" |
| `identifiers(smiles)` | "Structure \| Identifiers" |
| `drugLikeness(smiles)` | "Biology \| Drug Likeness" |
| `toxicity(smiles)` | "Biology \| Toxicity" |
| `structuralAlerts(smiles)` | "Biology \| Structural Alerts" |
| `pharmacophoreFeatures(smiles)` | "Biology \| Pharmacophore Features" |
| `mpoWidget(semValue)` | "Chemistry \| MPO" |
| `mixtureWidget(mixture)` / `mixtureTreeWidget(mixture)` | "Chemistry \| Mixture" / "MixtureTree" |
| `molColumnPropertyPanel(molColumn)` / `molColumnHighlights(molColumn)` | Rendering / highlighting options for a Molecule column. |
| `synthonSubstructureSearchWidget(molecule)` / `synthonSimilaritySearchWidget(molecule)` | Synthon space search panels. |

### Demos (each opens a tutorial layout — use only on user request)

`demoSimilarityDiversitySearch`, `demoMMPA`, `demoRgroupAnalysis`,
`demoMoleculeActivityCliffs`, `demoChemicalSpace`, `demoScaffold`.

### Skipped on purpose (do not advertise to the user)

- `init`, `initChemAutostart` — package init lifecycle, runs automatically.
- Internal helpers: `getContainer`, `getChempropError`, `trainModelChemprop`,
  `applyModelChemprop` (used internally by `trainChemprop` / `applyChemprop`).
- Editors: `SearchSubstructureEditor`, `ChemSpaceEditor`, `ActivityCliffsEditor`,
  `MMPEditor`, `DeprotectEditor`, `mpoProfileEditor` — invoked by the framework.
- `mpoProfilesAppTreeBrowser`, `checkJsonMpoProfile` — UI plumbing.

---

## Admetica package — ADMET predictions (`Admetica:...`)

Docker-backed ADMET model panel covering Absorption, Distribution, Metabolism,
Excretion, and Toxicity endpoints.

**ADMET category enum** for `getModels(property)`: `'Absorption' | 'Distribution' | 'Metabolism' | 'Excretion' | 'Toxicity'` (or omit for all).

| Function | Tags | What it does |
|---|---|---|
| `getModels(property?: 'Absorption' \| 'Distribution' \| 'Metabolism' \| 'Excretion' \| 'Toxicity')` | — | Lists available endpoints. |
| `getAdmeProperties(table: dataframe, molecules: column<Molecule>, props?: list<string>)` | `vectorFunc`, `join(table)` | **Main entry point.** Predicts ADMET for a column → returns a DataFrame appended to `table`. `props` is a list of model names from `getModels()`; omit for all. |
| `getAdmePropertiesSingle(molecule: string)` | — | One molecule → DataFrame of ADMET predictions. |
| `admeticaMenu(table: dataframe, molecules: column<Molecule>, template: string, models: list<string>, addPiechart: bool, addForm: bool)` | `topMenu` | Same as `getAdmeProperties` but with optional pie-chart and form viewers. Top menu: Chem → Admetica → Calculate. |
| `admeticaHT(table: dataframe, molecules: column<Molecule>, absorption: list<string>, distribution: list<string>, metabolism: list<string>, excretion: list<string>)` | — | Hit Triage variant — pass per-category model lists (each from `getModels('<Category>')`). |
| `admeticaWidget(semValue: object)` | `panel` | Per-molecule context panel: "Biology \| Admetica". |
| `runAdmeticaApplication() / admeticaDemo()` | — | Opens the Admetica app / demo (browsePath `Chem`). |

Skipped: `info` (debug), `AdmeticaEditor` (UI editor).

---

## Out of scope

Generative chemistry (Reinvent4), molecular docking (Docking / Boltz1) are intentionally NOT covered by this catalog. Open the matching skill directly if the user asks for those — they're available as registered functions but not described here.

<!-- REMOVED: Reinvent4 package — generative chemistry (`Reinvent4:...`)

Docker-backed REINVENT4 wrapper. Optimizes/derivatizes a seed ligand for a
named target.

| Function | What it does |
|---|---|
| `getFolders()` | Lists available optimization targets (folders in `System:AppData/Reinvent4/...`). |
| `reinvent(ligand: string, optimize: string)` | **Main entry point.** Generates molecules from a seed SMILES (`ligand`). `optimize` is a folder name from `getFolders()`. Returns a DataFrame with `Input_SMILES`, `SMILES`, `Score`. Also writes lineage stickymeta when the `Lineage` schema is installed. |
| `reinventTopMenu(ligand: string, optimize: string)` | Same, then opens result as a new table view. Top menu: Chem → Generate molecules. |
| `runReinvent(ligand: string, optimize: string)` | Lower-level call into the Docker container. Returns raw JSON string. Cached. |

**For "generate N analogs of this ligand" intent:** call `reinvent(ligand, optimize)` where `ligand` is the seed SMILES (read from a molecule column or the user's message). If `optimize` isn't specified, pick the first folder from `getFolders()`. Filter / slice the result DataFrame to N rows after the call.

Skipped: `info`, `ReinventEditor` (UI editor).

---

## Docking package — small-molecule docking (`Docking:...`)

AutoDock Vina via the `autodock` Docker container; renders poses with NGL.

| Function | What it does |
|---|---|
| `getConfigFiles()` | Lists docking targets (folders in `System:AppData/Docking/targets` that contain a `.gpf`). |
| `getAutodockResults(table: dataframe, ligands: column<Molecule>, target: string, poses: int = 10)` | **Main entry point.** Docks each ligand into `target` (folder name from `getConfigFiles()`) with `poses` conformations. Returns a DataFrame (also joined to `table` via `action: join(table)`). Vector-func. |
| `runAutodock(table: dataframe, ligands: column<Molecule>, target: string, poses: int = 10)` | UI wrapper: also adds color coding + sort by binding energy in the current grid. Top menu: Chem → Docking → AutoDock. |
| `dockLigandCached(jsonForm: string, containerId: string)` | Low-level cached call into the container — prefer `getAutodockResults`. |
| `getAutodockService()` | Returns the `IAutoDockService` singleton (for advanced use). |
| `autoDockApp()` / `dockingView(path?: string)` | App entry points (browsePath `Bio`). |
| `autodockWidget(molecule: string)` / `autodockPanel(smiles: string)` / `getAutodockSingle(molecule: string, showProperties?: bool, table?: dataframe)` | Per-molecule context panels and a programmatic single-ligand docking widget. |
| `isApplicableAutodock(molecule: string)` | Returns `true` if a molecule string already carries docking remarks (used by panel registration). |
| `demoDocking()` | Opens the docking demo (Bioinformatics → Docking). |

**For "dock these ligands into the protein target" intent:** call `Docking:getAutodockResults(t, ligands, target, poses)` (or `Boltz1:docking` if the user asked for deep learning). If `target` isn't specified, pick the first folder from `getConfigFiles()`; if `ligands` isn't specified, use the first `column<Molecule>` in `t`. **Never route docking through a `Chem:` function** — Chem has no docking entry point.

Skipped: `info`.

---

## Boltz1 package — Boltz-1 deep-learning structure prediction (`Boltz1:...`)

Boltz-1 (open-source AlphaFold-style) for folding and docking. **For
small-molecule work, the relevant entry point is `docking` — `folding` operates
on macromolecule sequences.**

Small-molecule-relevant:

| Function | What it does |
|---|---|
| `getBoltzConfigFolders()` | Lists available Boltz config sets. |
| `docking(table: dataframe, ligands: column<Molecule>, config: string)` | **Small-molecule docking via Boltz-1.** `config` is a folder name from `getBoltzConfigFolders()`. Returns a DataFrame joined to `table`. Top menu: Chem → Docking → Boltz. |
| `boltzWidget(molecule: string)` | Context panel for a 3D molecule cell that contains Boltz confidence scores. |
| `isApplicableBoltz(molecule: string)` | Returns `true` if the molecule already has Boltz scores. |
| `boltz1App()` | App entry point (browsePath `Bio`). |
| `runBoltz(config: string, msa: string)` | Low-level cached call. |

**For "predict a binding pose via Boltz" / "dock with Boltz" intent:** call `Boltz1:docking(t, ligands, config)` — NOT `Boltz1:folding` (proteins). If `config` isn't specified, pick the first folder from `getBoltzConfigFolders()`. The molecule column param is named **`ligands`** (plural), matching the function signature.

END REMOVED -->

---

## Companion data-side catalog

For ChEMBL, MolTrack, HitTriage, and Curves (compound registration, hit triage
workflows, dose-response analysis, bioactivity DBs), see the sibling skill
**`datagrok-chem-data`**.

## Generic compute / simulation

The `Compute` package exists but contains no small-molecule-specific functions
— it provides the `RichFunctionView` and `Pipeline` scaffolding plus
non-cheminformatics demo scripts (Lorenz attractor, Newton cooling, Lotka-
Volterra, antibody Fab-arm exchange). Skip it unless the user asks for the
RichFunctionView UI shell or compute pipeline framework.
