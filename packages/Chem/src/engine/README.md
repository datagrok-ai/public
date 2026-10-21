# Chem's chemical engine: RDKit and Crux

The package setting **Chemical Engine** (`ChemicalEngine`: `RDKit`, the default, or `Crux`) picks the engine behind
Chem's own chemistry. Crux ([crux-core](https://github.com/datagrok-ai/crux-core), Rust compiled to one wasm file,
`vendor/chem-wasm/`) answers RDKit.js's calls with the same names, arguments and return shapes. In Crux mode, what
has moved to Crux runs on Crux, and every Crux path has the RDKit path behind it. Other packages get RDKit objects
from Chem in both modes (`getRdKitModule()`, `getRdKitService()`); in Crux mode the pictures and values of the moved
functions come from Crux (the last table below).

## What lives where

| File | What it does |
|---|---|
| `chem-engine.ts` | The one reader of the setting (`requestedEngine`, throttled to one read a second, through `chemSetting`, which reads the published Chem's settings from the package-test bundle too) and the tests' override (`setChemicalEngine`); the engine in effect (`engineInEffect`: Crux once its module is ready, else RDKit), a generation number that changes with it (caches compare it), switch listeners (the package repaints the grids), `renderingModule` (the module drawing goes through), `cruxReady` |
| `crux-module.ts` | The lazy loader of `vendor/chem-wasm/chem_wasm_bg.wasm`: fetched and compiled once (`compileCrux`, the compiled Module handed to the Crux workers), instantiated once on the main thread (`loadCrux`, set up like RDKit: `prefer_coordgen(false)`, `use_legacy_stereo_perception(false)`); a load not done within 60 s of its fetch's start (a stalled fetch) has failed; a failure to load is logged once and leaves RDKit until one retry a minute later (a second failure is final for the session); a trap in a Crux call (`WebAssembly.RuntimeError`, `noteCruxError`) leaves RDKit for the session, as the instance may be broken |
| `fallback.ts` | `cruxOrRdKit` / `cruxOrRdKitAsync` and `warnCruxFallback`: a Crux failure or a "not supported by Crux" refusal runs the call on RDKit, logged once per call name with `console.warn` |
| `worker-jobs.ts` | The worker jobs by engine: fingerprints (`fingerprintsByEngine`, and the engine a cache belongs to) and structural alerts (`alertsByEngine`): the Crux pool in Crux mode, the RDKit pool otherwise and behind it (the whole call when Crux fails; otherwise the molecules Crux cannot read, on the RDKit pool) |
| `../crux/crux-service.ts` | The Crux worker pool: molecules in chunks, each worker handed its next chunk when the last is done, each chunk within its own time (30 s plus a per-molecule allowance, from its start); a trap, a timeout or a lost worker restarts the pool (another failure keeps it and its search indexes), and the call falls back to RDKit |

RDKit mode never fetches the file (the Crux substructure search, the `SubstructureSearchEngine` default, does).
In Crux mode the package's init loads it beside RDKit; a switch at run time draws with RDKit until it is ready, then
the renderers empty their caches and the grids repaint.

## What follows the setting

| Area | Crux mode | RDKit behind it |
|---|---|---|
| Molecule cell renderer (`rendering/rdkit-cell-renderer.ts`), `drawMolecule`, the mixture renderer | parse, layout, alignment, highlights and drawing through the Crux module; caches per engine, emptied on a switch | per molecule: what Crux cannot read (`null`), a scaffold it cannot read, a call it refuses |
| `canvasMol` (`utils/chem-common-rdkit.ts`, `drawMoleculeToCanvas`) | the same | per molecule |
| Descriptors (`descriptors/descriptors-engine.ts`: the menu and its editor, the transform, the vector function, the panel, `chemDescriptors`, `chemDescriptorsTree`) | the tree from the main thread's instance, the values from the Crux workers (`crux/crux-service.ts`, each distinct molecule once, in chunks of 1000 across the pool, back as typed arrays) | the chem-chem container, per call and per molecule Crux cannot read (a row of nulls for a molecule given) |
| Scaffold tree (`getScaffoldTree`) | a Crux worker, the same JSON as `scripts/scaffold_tree.py` | the Python script, per call |
| Contains / Not contains search | not this setting: `SubstructureSearchEngine` (default Crux) | RDKit search per query / per unreadable molecule (`crux/crux-searches.ts`) |
| Fingerprints (`chemGetFingerprints` and its cached `~<column>.<type>` columns, the query's fingerprint, similarity and diversity, chemical space and cliffs' preprocessing, BitBIRCH, the SAR matrix's clustering; `engine/worker-jobs.ts`) | the Crux workers run the RDKit worker's own code on the Crux module (the same types and options: Morgan r2/2048 after `remove_hs_in_place`, ...); a cached column belongs to the engine that filled it (`~<column>.Morgan.Crux`, RDKit's `~<column>.Morgan`), and saving one drops the other engine's; options other than the defaults are in its name (`~<column>.Morgan(radius=3)`, `utils/fingerprint-options.ts`); each cache keeps the column's edits since it was read. FCFP and chirality (`useFeatures`, `useChirality`) are Crux's in both modes: Chem's RDKit.js has neither | the RDKit pool, per call and per molecule Crux cannot read (none for FCFP and chirality) |
| Structural alerts and pharmacophore features (`panels/structural-alerts.ts`, `panels/pharmacophore-features.ts`) | the Crux workers, matching as the RDKit worker does, each molecule parsed once | the RDKit pool, per call and per molecule Crux cannot read |
| An activity cliff's MCS, highlighted in its tooltip and panel (`analysis/activity-cliffs.ts`, `cliffMcs`) | RDKit's Rascal MCES (crux's `get_mces`) on a Crux worker (`CruxService.pairUncommonAtomsAndBonds`): each molecule's atoms and bonds outside it, from its own matches; the most bonds the two molecules share, possibly in several fragments, where FMCS finds one connected substructure; no highlight for a pair under its similarity threshold | RDKit's FMCS on the RDKit pool, per call (a Crux failure) |
| The property table (`widgets/properties.ts`: the Properties panel, Chemical Properties' columns, `getProperties`, `getChemPropertyFunction`) | crux's `chemProperties` (a column on the Crux workers, a panel's molecule on the main thread's instance): RDKit's definitions, each with its method and reference (`CRUX_PROPERTY_METHODS`); toxicity and drug-likeness stay OpenChemLib's | OpenChemLib's workers per call; a panel's molecule Crux cannot read: OpenChemLib (a column keeps one method: empty cells) |
| Reactions (`utils/reactions/reaction-engine.ts`: Transformation and Two-Component reactions, Deprotect, Reaction Enumeration) | `get_rxn`, `run_reactants` and their molecules on the main thread's Crux module; drawing a reaction stays RDKit's | RDKit's module per molecule (per run for an enumeration) |
| MMPA fragmentation (`getMmpFrags`: MMP analysis, the SAR matrix's series and clusters; `engine/analysis-jobs.ts`) | the Crux workers run the RDKit worker's `mmpGetFragments` on the Crux module (`crux/crux-analysis-jobs.ts`), the molecules in chunks across the pool | the RDKit pool, per call and per molecule Crux does not read |
| R-group decomposition (`rGroupsMinilib`: the R-group analysis dialog, `rGroupDecomposition`; the SAR matrix's clusters) | a Crux worker runs the RDKit worker's `rGroupAnalysis`, a SAR matrix cluster per job across the pool; a cancel ends the operation at once (the worker's result is dropped) | the RDKit pool, per call (a molecule Crux cannot read, a failure) and per cluster |
| FMCS (`getMCS`: the R-group dialog's MCS; Cluster MCS, `performClusterMCS`; the SAR matrix's anchors) | a Crux worker runs the RDKit worker's `mostCommonStructure`; a cluster per job across the pool, within FMCS's own Timeout of 45 s ('' past it, as the RDKit pool gives a cluster whose worker it restarts after 45 s); Cluster MCS' SMILES from the main thread's Crux module | the RDKit pool, per call and per cluster Crux cannot read |

Fingerprints by Crux are RDKit 2026.03's, Chem's RDKit is 2024.09: on Chem's test molecules (1,402) Morgan, pattern,
MACCS and RDKit fingerprints are identical, atom-pair and topological-torsion ones differ on 17 and 16 molecules, each
of them RDKit's own change between the two releases (RDKit.js 2026.03.6 gives Crux's bits). Structural alerts match
identically. A fingerprint cache belongs to the engine that filled it: the RDKit substructure search keeps its own
pattern fingerprints (the plain-named columns) on RDKit in both modes.

In Crux mode the RDKit pool is not started for what has moved; everything else (the other MCS calls, R-groups, MMPA,
InChI, conversions, reactions, the substructure search's RDKit path, the context panel's one-molecule widgets) is RDKit
in both modes.

Descriptors in Crux mode are `dockerfiles/chem/chem.py` as this package has it, on the RDKit release Crux reproduces:
the same tree (groups, versions, names, order; the descriptions RDKit's docstrings give) and the same values (2D within
1e-12; the 3D ones from Crux's own embedding: chem.py embeds unseeded, so the container's 3D values change from call to
call too). A container image built before commit 30215919ce (GROK-20505: `chem_utils.py`'s curated descriptions, a
"type" key and per-name tags replaced by RDKit's docstrings) shows the old texts in the dialog; the next container
build gives what Crux gives. The scaffold tree is `scripts/scaffold_tree.py`'s (ScaffoldGraphDG 1.1.8, which drops a
molecule RDKit cannot read back from its own aromatic MOL block, e.g. an aromatic `[nH]`; Crux does the same); its
nodes' MOL blocks are Crux's layout.

Settings on existing servers: the server keeps each group's saved settings as one map (the settings editor saves every
property, defaults included) and shows a property the map lacks with the published version's default
(`packages_service.dart`, `getPropertyValues` / `getUserPackagesProperties`); publishing stores nothing. So the new
`SubstructureSearchEngine` default (Crux) reaches only servers where Chem's settings were never saved; elsewhere the
stored value (RDKit, the earlier default) stays until an admin switches it. `ChemicalEngine` is new: RDKit everywhere
until switched.

## What other packages see

RDKit mode (the default) changes nothing they receive. In Crux mode, only the pictures and the values that come
from Crux:

| What they call | RDKit mode | Crux mode |
|---|---|---|
| `getRdKitModule()` (Bio, Helm, SequenceTranslator, Peptides, PowerPack, Chembl, SureChembl, KetcherSketcher, the BDD bindings), `getRdKitService()` and its methods | RDKit | RDKit (never a Crux object) |
| `drawMolecule` / `grok.chem.drawMolecule`, `canvasMol` / `grok.chem.canvasMol`, the Molecule cell renderer, `rdKitCellRenderer` | RDKit's picture | Crux's picture (the same drawing calls; RDKit per molecule Crux cannot do) |
| `chemDescriptors` / `grok.chem.descriptors`, `chemDescriptorsTree`, the descriptors vector function and menu (HitTriage, statistics, Flow) | the chem-chem container | Crux: the same names, types and nulls; the 3D values from Crux's own embedding |
| `getScaffoldTree`, the Scaffold Tree viewer | the Python script | Crux: the same JSON |
| `searchSubstructure`, `filterBySubstructure`, the substructure filter (Bio, DrugBank, PowerPack, Flow, ...) | `SubstructureSearchEngine` (Crux by default): the same hits as RDKit | the same |
| `getMorganFingerprints`, `getFingerprints`, `getSimilarities`, `findSimilar`, `getDiversities` (Bio, Dendrogram, Flow, Chemspace, ...); `getStructuralAlerts`, `runPharmacophoreFeatures` | RDKit's (`getFingerprints`' FCFP and chirality options: Crux's) | Crux's (the same types and options: identical bits on Chem's datasets but where RDKit 2024.09 and 2026.03 differ) |
| `getProperties`, `getChemPropertyFunction`, Chemical Properties (Admetica's MW, HitTriage) | OpenChemLib's | Crux's: RDKit's definitions (MW the average weight, not OpenChemLib's monoisotopic mass) |
| `getMCS`, Cluster MCS (`clusterMCS`), `rGroupDecomposition`, MMP analysis | RDKit | Crux's: RDKit 2026.03's (an MCS SMARTS may name its atoms in another order than Chem's RDKit 2024.09, the same query) |
| every other function (`convertMolNotation`, `canonicalize`, `isSmiles`, `validateMolecule`, InChI, ...) | RDKit | RDKit |

## Moving a feature behind the setting

1. Find its RDKit path and keep it as it is: it stays the RDKit-mode path and the fallback.
2. Write the Crux path beside it with the same inputs and outputs (the same JSON, columns, types and nulls); Crux
   objects never leave Chem.
3. Branch once at the feature's entry: `isCruxRequested()` for a call (a worker job, a calculation), `renderingModule`
   / `engineInEffect` for main-thread drawing; wrap the Crux path in `cruxOrRdKit(Async)` with the call's name.
4. Caches holding engine objects: compare `engineGeneration()` and empty themselves when it changes.
5. Tests in the `chemical engine` category: both engines, the switch, the fallback (`setCruxLoadFailureForTests`),
   Crux against RDKit on Chem's datasets.

Running a whole suite in Crux mode: `setChemicalEngine(ChemicalEngine.Crux)` in code, or the package setting on the
server, which every page then starts with (`misc/crux/engine-setting.sh set Crux`, `restore` after: package tests, the
BDD suite, the tests of the packages that call Chem).
