# Chem's chemical engine: RDKit and Crux

The package setting **Chemical Engine** (`ChemicalEngine`: `RDKit`, the default, or `Crux`) picks the engine behind
Chem's own chemistry. Crux ([crux-core](https://github.com/datagrok-ai/crux-core), Rust compiled to one wasm file,
`vendor/chem-wasm/`) answers RDKit.js's calls with the same names, arguments and return shapes, and adds calls of its
own (descriptors, properties, QED, the scripts' replacements). In Crux mode a feature that has moved runs on Crux, with
its old RDKit path behind it for what Crux cannot do; in RDKit mode it runs as before. Other packages get RDKit objects
from Chem in both modes.

The history of each move, its parity numbers and timings are in `CHANGELOG.md`.

## What lives where

| File | What it holds |
|---|---|
| `engine/chem-engine.ts` | The setting (`requestedEngine`, `isCruxRequested`, the tests' `setChemicalEngine`); the engine in effect (`engineInEffect`: Crux once its module is ready), `engineGeneration` and `onEngineSwitch` for caches and repaints; `byEngineAsync` (a call by the setting's engine); the main-thread helpers `renderingModule`, `byEngineModule`, `withMol`, `cruxReady` |
| `engine/crux-module.ts` | The loader of the Crux file: compiled once (the workers get the compiled module), one main-thread instance (`cruxModuleSync`); a failed load leaves RDKit (one retry after a minute), a wasm trap leaves RDKit for the session (`noteCruxError`) |
| `engine/fallback.ts` | The fallback: `cruxOrRdKit` / `cruxOrRdKitAsync` (the RDKit side when the Crux side throws), `rdkitForRows` (a column job's RDKit side for the rows Crux left), `warnCruxFallback` (logged once per call name, to the console and the platform's log) |
| `engine/worker-jobs.ts` | Column jobs by engine: fingerprints, structural alerts and pharmacophores, notation conversions, InChI; `cruxInput` (a column as the Crux workers get it) |
| `engine/molecule-jobs.ts` | Column jobs by engine: CoordGen coordinates, flattening, beautifying, elemental analysis, the R-group dialog's cores, the SAR matrix's atom counts |
| `engine/analysis-jobs.ts` | Analyses by engine: MMPA fragments, R-groups, MCS, cluster MCS and R-groups, fragment linking; `molzip` (Crux in both modes) |
| `crux/crux-service.ts` | The Crux worker pool (`CruxService`, `getCruxService` in `crux-searches.ts`): `jobInChunks` (a column in chunks across the workers), `idleWorkerJob` (one long call with its own deadline), clusters one per job, typed methods; the search index's segments; time limits, worker restarts |
| `crux/crux.worker.ts` | The worker: the search segments (`build`, `search`, `drop`) and `JOBS`, every job table by name |
| `crux/jobs/*.ts`, `crux/crux-analysis-jobs.ts` | The job tables, one per area (`calculations`, `batch`, `molecule-jobs`, `search-types`, `python-free`, `synthon-search`; the analyses), each job a `WorkerJob` |
| `crux/crux-api.ts` | Crux's calls beyond RDKit.js's; the job shape (`WorkerJobContext`, `WorkerJob`, `throwIfTrap`) |
| `crux/crux-searches.ts`, `crux/crux-search-types.ts` | The substructure search on Crux: Contains / Not contains on the index; Exact, Stereo agnostic, (Not) included in, Similar |
| `crux/crux-inchi.ts` | InChI by the split: InChI's input from Crux, the InChI library in crux-sketch's module |
| `rdkit-service/` | RDKit: the pool (`rdkit-service.ts`, `rdkit-service-worker-client.ts`, `rdkit-service-worker-api.ts`, `../rdkit.worker.ts`) and the worker code (`rdkit-service-worker*.ts`, `query-features.ts`), which the Crux workers also run, on the Crux module |
| `utils/chem-common-rdkit.ts` | RDKit's main-thread module and pool (`getRdKitModule`, `getRdKitService`); drawing (`drawMoleculeToCanvas`, by engine); `canonicalWith`, `readsAsMolecule`, `sameMolecule` (by engine) |
| `utils/mol-creation_rdkit.ts` | Reading molecules and queries with either engine's module (`getMolSafe`, `getQueryMolSafe`, `readsAsQuery`, `clusterMcsSmiles`) |
| `utils/convert-notation-utils.ts` | `_convertMolNotation` (by engine) and `_convertMolNotationByRdKit` (RDKit in both modes) |
| `rendering/rdkit-cell-renderer.ts`, `rendering/mixture-cell-renderer.ts` | The Molecule and mixture cell renderers, by the engine in effect (the "RDKit" names stay: other packages use them); `_fetchMolByEngine` holds the per-molecule fallback |
| `rendering/rdkit-reaction-renderer.ts`, `rendering/reaction-depiction.ts`, `utils/reactions/reaction-drawing.ts` | Reaction drawing: RDKit in both modes until Crux draws reactions |
| `rendering/atom-picker-controller.ts`, `rendering/mol3d-hover-handler.ts` | The 2D↔3D atom picker and hover, by engine (`mapAtomIndices2Dto3DByEngine`) |

RDKit mode loads the Crux file only for what is Crux's in both modes (`molzip`; FCFP and chirality fingerprints, which
Chem's RDKit.js lacks). Crux mode loads it at the package's init; a switch at run time draws with RDKit until it is
ready, then the renderers empty their caches and the grids repaint.

## Moving a feature behind the setting

1. Keep the RDKit path as it is: it stays RDKit mode's path and the fallback. Name RDKit-only code so it can be found
   (`…ByRdKit`, an `rdkit` closure, a function in `rdkit-service/`).
2. Write the Crux path beside it with the same inputs and outputs (the same JSON, columns, types and nulls). Crux
   objects never leave Chem. A helper both paths use takes the module (`RDModule`, RDKit.js's shape) as a parameter
   and never calls `getRdKitModule()` or the RDKit pool itself.
3. Branch once, at the feature's entry, in one of three shapes:
   - **One molecule on the main thread** (a panel, a drawing, a check): `byEngineModule(call, rdKitModule,
     (module, crux) => ...)`, the job written once against the module (`withMol(module, crux, molecule, read)` parses
     the molecule); null (or a throw) from Crux's runs it on RDKit.
   - **A column** (work that would freeze the page and is the same per molecule): a job in a table of `crux/jobs/`
     (`(ctx, chunk, ...args) => result`), sent by `CruxService.jobInChunks`, typed arrays back where they pay; on the
     main thread a `…ByEngine` function: `byEngineAsync(call, crux, rdkit)` (Crux mode: `crux()`, `rdkit()` when it
     fails; RDKit mode: `rdkit()`), with `cruxInput` for the column and `rdkitForRows` for the rows Crux leaves.
   - **One long call** (an embedding, a conformer search): a job sent by `CruxService.idleWorkerJob(op, args,
     timeoutMs)` with its own deadline in its options, behind `byEngineAsync` as a column's.
4. Caches holding engine objects compare `engineGeneration()` and empty themselves when it changes.
5. Tests in the `chemical engine` category (or the feature's own), through `src/tests/engine-utils.ts`:
   `inBothEngines` (the same run in each engine, compared), `inEngine`, `expectNoFallback`, `withCruxUnreadable` (a
   molecule the Crux workers cannot read); its header says what each becomes when RDKit goes.

A whole suite in Crux mode: `setChemicalEngine(ChemicalEngine.Crux)` in code, or the package setting on the server
(`misc/crux/engine-setting.sh set Crux`, `restore` after).

## What follows the setting

Behind every Crux path below, besides what is listed: RDKit for the whole call when Crux fails, and for a molecule
Crux cannot read or a CXSMILES whose query features Crux's reader skips (`cruxReadsCxExtension`).

| Area (where it branches) | Crux mode | RDKit behind it |
|---|---|---|
| Drawing: the cell renderer, `drawMolecule`, `canvasMol`, the mixture renderer (`rendering/rdkit-cell-renderer.ts`, `utils/chem-common-rdkit.ts`) | the Crux module; caches per engine | per molecule, per scaffold Crux cannot read |
| Descriptors (`descriptors/descriptors-engine.ts`) | the tree on the main thread's module, the values on the workers | the chem-chem container |
| Scaffold tree (`getScaffoldTree`, `package.ts`) | a Crux worker, `scaffold_tree.py`'s JSON | the script |
| Substructure search, every type, and the filter's query names (`crux/crux-searches.ts`, `crux/crux-search-types.ts`) | the workers over the column's index segments; Similar by the engine's fingerprints | the RDKit search per query Crux cannot read, per failed search; per row Crux cannot read |
| Fingerprints, similarity, diversity, chemical space (`engine/worker-jobs.ts`) | the workers (`crux/jobs/batch.ts`); caches per engine (`~<column>.Morgan.Crux`) | the RDKit pool; none for FCFP and chirality (Crux's in both modes) |
| Structural alerts, pharmacophore features (`engine/worker-jobs.ts`) | the workers, a chunk per Crux call | the RDKit pool; per SMARTS Crux reads nothing from |
| Notation conversions, `canonicalize`, `grok.chem.convert` (`utils/convert-notation-utils.ts`, `engine/worker-jobs.ts`) | a molecule on the main thread's module, a column on the workers | RDKit's module or pool |
| `isSmiles`, `isSmarts`, `validateMolecule`, `removeDuplicates`, the cell editor's check (`package.ts`, `utils/chem-common-rdkit.ts`) | the main thread's module | RDKit's module (`validateMolecule`'s messages are RDKit's) |
| InChI and InChIKeys (`engine/worker-jobs.ts`, `crux/crux-inchi.ts`) | InChI's input from Crux, the library in crux-sketch's module | the RDKit pool |
| CoordGen, Flatten, `beautifyMols`, elemental analysis (`engine/molecule-jobs.ts`) | the workers (`crux/jobs/molecule-jobs.ts`) | the RDKit pool, or the main thread's loop |
| MMPA, R-groups, MCS, cluster MCS, fragment linking, the MMP viewer's generations (`engine/analysis-jobs.ts`) | the RDKit worker's own code on the workers (`crux/crux-analysis-jobs.ts`) | the RDKit pool per call, molecule, row or cluster |
| An activity cliff's MCS (`analysis/activity-cliffs.ts`) | RDKit's Rascal MCES on a worker | FMCS on the RDKit pool |
| The R-group dialog, the SAR matrix, the Scaffold Tree viewer, the MMP viewer's pictures (`analysis/`, `widgets/scaffold-tree.ts`) | columns on the workers, single structures on the main thread's module | RDKit's module per structure |
| Reactions (`utils/reactions/reaction-engine.ts`) | `get_rxn`, `run_reactants` on the main thread's module | RDKit's module per molecule or run |
| Properties, Drug Likeness, Toxicity (`widgets/properties.ts`, `widgets/drug-likeness.ts`, `widgets/toxicity.ts`) | crux's `chemProperties`, `chemQed`, `chemAlerts` | OpenChemLib (Toxicity Alerts... is Crux's only) |
| 3D coordinates (`Chem:embed3d`, `widgets/structure3d.ts`) | crux's `chemEmbed3d` on an idle worker, its own deadline | the script (none for a timeout) |
| One molecule on the main thread: `getCLogP`, the sketcher's free text, the atom picker, panels (`widgets/`, `rendering/`) | the main thread's module (`byEngineModule`) | RDKit's module |
| The Python scripts replaced (`python-free.ts`: Curate, Mutate, Generate Conformers, Gasteiger charges, Butina, InchiToMol, Calculate logP, Murcko Scaffolds, Filter by Catalogs, Map Identifiers; `widgets/synthon-search.ts`: Synthon Search) | crux's replacements of the scripts (`crux/jobs/python-free.ts`, `crux/jobs/synthon-search.ts`) | the script, renamed with the suffix `Python`, per call |
| `molzip` (`engine/analysis-jobs.ts`) | the workers, in both modes | none: it throws when Crux does not load |

Known differences, both deliberate: Crux is RDKit 2026.03, Chem's RDKit.js is 2024.09, so where they differ (some
canonical SMILES' stereo marks or atom order, atom-pair and torsion fingerprints of a few molecules, an MCS SMARTS'
atom order) Crux gives 2026.03's answer; and the property table's MW is RDKit's average weight where OpenChemLib gave
the monoisotopic mass.

## What other packages see

RDKit mode changes nothing they receive. In Crux mode, the pictures and values of the moved functions come from Crux:

| What they call | Crux mode |
|---|---|
| `getRdKitModule()`, `getRdKitService()` (Bio and the bio library, Helm, SequenceTranslator, Chembl, SureChembl, KetcherSketcher, the BDD bindings) | RDKit, never a Crux object |
| `drawMolecule`, `canvasMol`, the Molecule cell renderer | Crux's picture |
| `chemDescriptors`, `chemDescriptorsTree`, `getScaffoldTree`, `getProperties`, `getChemPropertyFunction` | Crux's: the same names, types and nulls |
| `searchSubstructure`, `filterBySubstructure`, the substructure filter | Crux's search: the same hits (Similar: Crux's fingerprints) |
| `getFingerprints`, `getMorganFingerprints`, `getSimilarities`, `findSimilar`, `getDiversities`, `getStructuralAlerts` | Crux's |
| `getMCS`, `clusterMCS`, `rGroupDecomposition`, MMP analysis | Crux's |
| `convertMolNotation`, `canonicalize`, `grok.chem.convert`, `isSmiles`, `validateMolecule`, `getInchis`, `getInchiKeys`, `removeDuplicates` | Crux's strings |
| `beautifyMols`, `getCLogP`, `embed3d` | Crux's |
| The replaced scripts' functions (`Curate`, `Mutate`, `GenerateConformers`, `ButinaMoleculesClustering`, `InchiToMol`, `CalculateLogP`, `MurckoScaffolds`, `FilterByCatalogs`, `Chem:SynthonSearch`) | Crux's: the scripts' names, inputs and outputs |
| `molzip` | Crux's in both modes |
| every other function | RDKit |

## What stays RDKit, by design or for now

| What | Why |
|---|---|
| `getRdKitModule()` / `getRdKitService()` for other packages | they take RDKit objects; each needs its own move (or a Crux object API) |
| Reaction drawing (`rendering/rdkit-reaction-renderer.ts`, the reaction dialogs' previews: `drawReactionByRdKit` in `utils/reactions/reaction-drawing.ts`, `_convertMolNotationByRdKit`; `validateReactionSmarts` in `utils/reactions/reaction-storage.ts`, which gates the drawing) | crux has no `drawReaction` yet |
| `validateMolecule`'s messages | RDKit's log; Crux keeps none |
| The chem-chem container's other identifiers, Toxicity Risks... (OpenChemLib) | not RDKit |

## Removing RDKit

When Crux becomes the default and RDKit.js goes, the change is deletion, in this order:

1. **Close the gaps first**, each a fallback that would otherwise have nothing behind it: CXSMILES query features in
   Crux's reader (`cruxReadsCxExtension` then always true), reaction drawing, a Crux object API (or moves) for the
   packages that call `getRdKitModule()`, `validateMolecule`'s messages. The platform log's `Crux fallback: <call>`
   warnings show what else RDKit still does in practice.
2. **Make Crux the only engine**: the setting goes (`package.json` `ChemicalEngine`, `chem-engine.ts`' reading of it);
   `isCruxRequested()` and `engineInEffect() === Crux` are true.
3. **Delete every RDKit side**, mechanically:
   - `byEngineAsync(call, crux, rdkit)`, `cruxOrRdKit(call, crux, rdkit)`, `cruxOrRdKitAsync(...)` → `crux()`;
   - the `rdkit` closure of each `…ByEngine` function (and any `if (!isCruxRequested()) return rdkit();`) → gone;
   - `await rdkitForRows(...)` → gone (the rows stay Crux's answer: none);
   - `byEngineModule(call, rdKitModule, job)` → `job(cruxModuleSync(), true)`; `renderingModule` → the Crux module;
   - the `if (crux)` branches inside drawing and conversions (`drawMoleculeWith`, `convertWith`, `withMol`: the
     `crux` flag) → the Crux branch;
   - the fallbacks that are not one call: the cell renderer's per molecule (`_fetchMolByEngine`, `_rendererGetOrCreate`, `_molToSvg`),
     descriptors' container rows, the 3D embedding's (`cruxEmbed3d`), `rdkitCoreMolblock`
     (`analysis/r-group-analysis.ts`), `v3000SmilesByRdKit` (`widgets/scaffold-tree.ts`);
   - each replaced Python script (`*Python`) and its call.
4. **Delete the RDKit-only files and functions** (the table below).
5. Remove RDKit.js itself (`src/RDKit_minimal.js`, its wasm, `rdkit_lib_version.js`, the package init's RDKit load).

| File or area | What goes | What stays |
|---|---|---|
| `engine/chem-engine.ts` | `ChemicalEngine`, `requestedEngine`, `setChemicalEngine`, `isCruxRequested`, `engineInEffect`, the generation and switch listeners, `byEngineAsync`, `renderingModule`, `byEngineModule`, `withMol`'s `crux` flag | `chemSetting`, `cruxReady` (as "the module is loaded") |
| `engine/fallback.ts` | the whole file (`cruxOrRdKit`, `cruxOrRdKitAsync`, `rdkitForRows`, `warnCruxFallback`, the counts and the reporter) | — |
| `engine/crux-module.ts` | the failure-to-RDKit logic (`cruxLoadFailed`, the retry, `noteCruxError`'s switch to RDKit; a load failure becomes an error); `getRdKitWebRoot` (Chem's web root from the package) | the loader |
| `engine/worker-jobs.ts`, `engine/molecule-jobs.ts`, `engine/analysis-jobs.ts` | each function's `rdkit` closure, `byEngineAsync` and `rdkitForRows` calls; `fingerprintsByEngine`'s `engine` parameter and `fingerprintsKey`'s RDKit suffix; `cruxInput`'s CXSMILES skip, `forCrux` and `setCruxUnreadableForTests` (the fallback's tests); `coreMolblocksByCrux`'s `undefined` rows | the Crux sides, renamed without `ByEngine` |
| `crux/crux-searches.ts` | `SubstructureSearchEngine`, `get`/`setSubstructureSearchEngine`, `rdkitMatcher`, `runCruxSearch`'s `rdkitSearch` | the search |
| `crux/crux-search-types.ts` | `rdkitAnswer`, `SegmentState`'s RDKit keys, the `rdkitSearch` parameters | the search types |
| `crux/crux-service.ts`, `crux/crux.worker.ts`, `crux/jobs/`, `crux/crux-analysis-jobs.ts` | the jobs' `rest` / `unread` rows (what the RDKit side did); `cruxReadsCxExtension` and `withinCruxReader` once Crux reads CXSMILES query features | the pool, the worker, the jobs |
| `rdkit-service/` | the pool: `rdkit-service.ts`, `rdkit-service-worker-client.ts`, `rdkit-service-worker-api.ts`, `../rdkit.worker.ts`; the worker code's pool-only parts (the molecule cache, `WORKER_CALL`) | the worker code (`rdkit-service-worker*.ts`, `query-features.ts`): the Crux workers run it on the Crux module (move it under `crux/`) |
| `utils/chem-common-rdkit.ts` | `initRdKitModuleLocal`, `initRdKitService`, `getRdKitModule`, `getRdKitService`, `_rdKitModule`, `_rdKitService`, `moduleInitialized`, `cancelChemOp`, `checkMoleculeValid`, `checkMolEqualSmiles` | drawing, the render options, `canonicalWith`, `readsAsMolecule`, `sameMolecule`, the web root |
| `utils/mol-creation_rdkit.ts` | — (it takes the module) | all of it (rename) |
| `utils/convert-notation-utils.ts` | `_convertMolNotationByRdKit` (with reaction drawing; the substructure filter's last resort for its query's SMARTS becomes `_convertMolNotation`), the CXSMILES check | `_convertMolNotation` on the Crux module |
| `package.ts` | the RDKit load in `init`; `validateMolecule`'s `validationMessageByRdKit` once Crux writes messages; `getScaffoldTree`'s `pythonTree`; `getRdKitModule` once no package calls it (a deprecated alias first) | the registrations; the RDKit module handed to renderers and helpers becomes Crux's |
| `python-free.ts` | each `cruxOrPython` call's script side, `pythonScript`, `scriptOutput`, `PYTHON_SCRIPTS`; the `*Python` scripts; Map Identifiers' container endpoints | the Crux paths (`ChemistryGasteigerPartialCharges` is a Flow node type: a decision) |
| `rendering/` | the cell renderer's RDKit branches (above), `IMolRenderingInfo.engine` and the engine prefix of its cache keys | drawing, the caches; the reaction renderer and depiction until Crux draws reactions |
| `widgets/` | `identifiers.ts`' RDKit side; the OpenChemLib side of drug likeness, properties and toxicity; `structure3d.ts`' script fallback; the substructure filter's `initRdKitService()` and `isUnchanged`'s RDKit side; `synthon-search.ts`' script call; `scaffold-tree.ts`' `v3000SmilesByRdKit` and its `byEngineModule` fallbacks | the panels and viewers |
| `descriptors/descriptors-engine.ts` | the container's tree and values, `fromContainerWhereUnread` | Crux's descriptors |
| `analysis/` | SAR Matrix: `withEngineMol`'s fallback, `canonicalCache`'s `crux` flag, the `engineGeneration()` key prefix, the RDKit warm-up; activity cliffs: `cliffMcs`' FMCS, `drawMolecules`' SMARTS branch, `pairKey`'s prefix; R-groups: `rdkitCoreMolblock`, the cancel race and the `_rdKitService` restart; Deprotect: `addDeprotectedColumn`'s module parameter; MMP: the `rdkit` fields behind `byEngineModule` | the analyses |
| `utils/reactions/`, `utils/reaction-enumeration/` | `EngineReaction`'s RDKit side, `onReactionEngine`'s fallback, `EnumerateOptions.rdkit`, `readOtherwise`; the dialogs' `rdkit` parameters once reactions draw on Crux | `reaction-drawing.ts` and the reaction renderer until then |
| `utils/fingerprint-options.ts` | `needsCrux` (always Crux) | the options |
| `tests/` | the RDKit legs, the parity and fallback tests (`engine-utils.ts`' header says what each helper's callers become) | the Crux tests |
