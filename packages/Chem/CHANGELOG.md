# Chem changelog

## v.next

* [core#1146](https://github.com/datagrok-ai/core/issues/1146): Substructure Search (top menu) opens the sketcher as soon as the filter's sketch link is there instead of after 500 ms; the scaffold tree's load resolves once its nodes are updated; Crux analyses take the next idle worker when it frees instead of polling every 50 ms; demos: no waits before applying a layout, and the activity cliffs demo awaits its scatter plot check
* [core#1146](https://github.com/datagrok-ai/core/issues/1146): Reactions: Removed a 100 ms wait after detecting semantic types, and a table chosen meanwhile no longer gets the previous table's molecule column; R-groups: Removed a 50 ms wait when replacing the latest analysis; MMP: The pair of a clicked row shows on the context panel at once instead of half a second later
* Activity Cliffs: The Activities input defaults to a column named like an activity or potency, else the first float column, instead of the first numeric one (the ID on sar_small)
* Map Identifiers: To Source defaults to chembl instead of the first source in alphabetical order (actor)
* R-Groups Analysis: Fixed "Invalid regular expression" for a column prefix holding a regular expression character, such as `R(`
* [#4403](https://github.com/datagrok-ai/public/issues/4403): MPO: Added a **Coverage Score** option to MPO Score that adds a column with the number of profile parameters measured for each compound
* Substructure search: In Crux mode every search type runs on Crux (the substructure filter, the filter builder's Molecule operators, `searchSubstructure`), so no search starts an RDKit worker: Exact, Stereo agnostic, Included in and Not included in on the Crux workers, Similar by Crux's fingerprints. The hits are RDKit mode's: all of 1,718 Exact, Stereo agnostic, Included in and Not included in searches on ChEMBL 100k and as many on Chem's test datasets; Similar differs only where RDKit's own MACCS and RDKit fingerprints changed from 2024.09 to 2026.03 (5 of 40 searches on ChEMBL 100k). The engine call on ChEMBL 100k over 16 worker threads: a column's first search takes 1.2-3.4 s on Crux, 5.8-6.3 s on RDKit; a later one takes milliseconds on Crux (Not included in: 6 ms; RDKit 525 ms). RDKit answers for a row or a query Crux cannot read and for a search that fails
* Synthon Search: In Crux mode, Synthon Space Search (the Databases | Synthon Search panels, Search Synthons, `Chem:SynthonSearch`) runs in the browser on Crux instead of the `synthon_search.py` script on the server: an exact port of the RDKit 2025.09.3 SynthonSpaceSearch the script runs, with the script's hits in its order, its SMILES and similarities (348 queries over Chem's library and RDKit's synthon spaces; one product's E/Z double bond is written the other way round), the search itself 1.8-2.7x faster. The library, and the fingerprints of the first similarity search, are kept in the browser for the next searches. The script, renamed `SynthonSearchPython`, runs in RDKit mode and for a query Crux does not take (a `#`, which the script reads as SMARTS; a CXSMILES; query atoms)
* Chemical Engine: In Crux mode, Curate, Mutate, Generate Conformers, Butina Cluster, Calculate logP (Chemical Properties' logP), Murcko Scaffolds, Filter by Catalogs, the Gasteiger Partial Charges panel and the InChI converter run in the browser on Crux instead of their Python scripts on the server, and Map Identifiers converts SMILES, InChIs and InChIKeys on Crux instead of the chem-chem container (its ChEMBL and UniChem lookups stay on the server). Each gives the script's (or the container's) results on RDKit 2026.03: Curate the script's SMILES (all but MOL blocks whose hydrogens sit in an S-group), Mutate the script's mutations (drawn as Python's random generator draws them, seeded anew per call), Generate Conformers the script's table (within the spread RDKit itself shows between platforms), Butina the script's clusters; logP, the scaffolds, the catalog flags and the identifiers the scripts' and the container's on Chem's 75,805 file molecules (logP bit for bit; of an InChI's SMILES, 23 written otherwise under the new stereo perception, all but one the same molecule). The functions keep the scripts' names, inputs, outputs and menus, so Flow and saved calls find them; the scripts, renamed with the suffix `Python` (`CuratePython`, `CalculateLogPPython`, ...), run in RDKit mode, when Crux fails and for a CXSMILES whose query features Crux does not read. The Gasteiger panel is the package's in both modes (the script draws it in RDKit mode); the script keeps its name, graphics output and panel role for Flow and shows no pane of its own. Chem's function, end to end on a local stand (warm), the script (RDKit mode) against Crux: Curate of 1,000 molecules 584 ms against 45 ms (every step on, 200 molecules: 3.4 s against 0.8 s), Butina of 1,000 molecules 523 ms against 17 ms, Mutate's 100 random mutations 252 ms against 4 ms (its 4,648 two-step ones 1.95 s against 0.24 s), Generate Conformers (50 requested) 548 ms against 222 ms, the Gasteiger panel's picture 252 ms against 161 ms, an InChI's molecule 100 ms against 1 ms
* Chemical Engine: In Crux mode, Transform | Recalculate Coordinates by CoordGen, Transform | Flatten Molecules, `beautifyMols` (Bio's To Atomic Level) and Elemental Analysis run on the Crux workers, with RDKit behind them per molecule: RDKit.js 2026.03.6's MOL blocks and SMILES (ChEMBL-10k: CoordGen and flat SMILES 10,000 of 10,000, beautified blocks 9,992; RDKit 2026.03 writes another Kekulé form of an aromatic ring in two thirds of the blocks and no valence field on a charged atom) and RDKit mode's element counts and charges (20,009 ChEMBL molecules). Elemental Analysis no longer freezes the page: Chem's function, end to end, is a 32 s main-thread loop for 20,000 SMILES in RDKit mode; in Crux mode the engine call takes 0.75 ms per molecule in one worker, the molecules spread over the workers
* Chemical Engine: In Crux mode, the R-Groups dialog's Core column is written on the Crux workers (its main-thread loop froze the page on a large table) and its MCS and core MOL blocks by Crux; the SAR matrix counts its atoms on the Crux workers and writes its keys, anchors, MW sort, scaffold grouping, templates and depictions by Crux, the engine of the structures its keys are compared with; reactions canonicalize their inputs and Remove Water and Salts its molecules on the Crux workers (RDKit.js 2026.03.6's SMILES)
* Chemical Engine: In Crux mode, `getCLogP` is RDKit 2026.03's Crippen logP (Chem's RDKit.js 2024.09 cuts it to 5 decimals: -0.58179 for -0.5818), and Crux reads the molecule for `isSmarts`, the sketcher's free-text box, the molecule cell editor, Recent and Favorites, the Identifiers panel (it said "Malformed molecule" for a molecule Crux reads and RDKit does not), the Scaffold Tree viewer's sample, nodes and highlights, the diversity viewer's tooltip, OpenChemLib sketcher's Get SMARTS, Deprotect's fragment, Reaction Enumeration's reaction check and the atom picker's 2D→3D mapping, with RDKit behind it per molecule
* Chemical Engine: In Crux mode, structural alerts, pharmacophore features and fingerprints of a column take one Crux call per chunk of molecules instead of one per molecule (for the alerts, per molecule and SMARTS), with the same results (identical on ChEMBL 100k), and the workers no longer pause every 50 molecules. Chem's worker code on ChEMBL 100k, end to end in one node thread, Crux mode now / Crux mode before / RDKit mode: structural alerts 44.7 / 205.3 / 375.6 s, pharmacophore features 8.0 / 9.2 / 20.2 s, Morgan 9.0 / 13.8 / 27.8 s, AtomPair 8.3 / 11.1 / 49.4 s, topological torsions 8.8 / 12.0 / 24.7 s, RDKit 38.5 / 43.0 / 88.6 s, MACCS 45.5 / 48.9 / 90.6 s, pattern 18.9 / 22.2 s (Crux mode)
* Chemical Engine: In Crux mode, RDKit matches an alert SMARTS Crux cannot read (it matched nothing), and a CXSMILES whose query features Crux's reader skips (`rb:`, `m:`, `Q_e` labels, ...) is drawn by RDKit (as a molecule or a scaffold) and gets RDKit's alerts and fingerprints (but FCFP and chirality, Crux's in both modes) and the container's descriptors, with no Crux properties or QED (OpenChemLib's for a panel), as conversions and the search treat it
* Chemical Engine: The first fallback of each call to RDKit per session is also reported to the platform's log (a package warning `Crux fallback: <call>` with the reason), so the switch to Crux can be judged on real usage; nothing is shown to the user
* Sketcher: Fixed Recent and Favorites never being de-duplicated: js-api calls `Chem:removeDuplicates`, which Chem registered as `filterMoleculeDuplicates`
* Molecule cell editor: Fixed OK doing nothing for a SMILES cell whose molecule did not read (the check whether the molecule was edited threw)
* Elemental Analysis: Fixed molecules with enhanced stereo (a V3000 column's chiral molecules: 64 of approved-drugs-100) counting no atoms: RDKit writes their MOL block in V3000, which the count read as empty; it reads V3000 blocks now, in both chemical engines
* Fingerprints: Fixed `getFingerprints` cutting the first and last characters of a fingerprint type that ended, but did not start, with a quote
* Sort by similarity: Fixed its progress text (it said "Installing ScaffoldGraph")
* MMP viewer: Fixed the "Refreshing pairs..." progress bar opening on every refresh; it opens when a refresh takes more than 200 ms
* SAR Matrix: Fixed cells that stayed empty on a large, zoomed-out matrix, however it was scrolled: with more cells on screen than the molecule renderer keeps pictures of (100), each repaint drew again the cells it had dropped, spent its time limit on them and left the rest for the next, which did the same; the matrix keeps its cells' pictures itself
* MPO: Fixed another profile handler being registered on each further `Chem:init` call (HitTriage calls it)
* Molecule rendering: A drawn molecule's menu (right click, or its `⋮`) has Enlarge: the molecule drawn again in a full-screen dialog, for a small one in the context panel or a viewer
* Tests: Added the `chemical engine: column jobs`, `chemical engine: python scripts` and `crux substructure search: search types` categories: each feature moved to Crux, in both engines, against RDKit mode or the script's pinned results, and its fallback; `synthon search`: in Crux mode `Chem:SynthonSearch`'s hits, their order, synthons and similarities (and Search Synthons') are the script's on the category's result files, with no fallback; the script's own run stays skipped (GROK-19863)
* Tests: The engine tests switch engines through one set of helpers (`src/tests/engine-utils.ts`: a Crux leg, an RDKit leg, both, a molecule Crux cannot read), so the RDKit legs and the fallback tests are found by name; the python-script tests call the bundle's functions, so each engine runs what it says (through the platform both ran the stand's engine). Faster: descriptors against the container on 30 molecules (72 before), the scaffold tree against the script on 250 of smiles.csv (750), the InChI comparison on 300 molecules (1,199), a stalled load's limit 1 s (2 s), RDKit's hits for the formerly-RDKit queries reading each molecule once (once per query); `isSmarts` and descriptor timings run with `--benchmark` only; `search: Crux finds what RDKit finds` is gone (`crux substructure search` checks it)
* Sketcher: The Crux sketcher reads a SMILES given as a molblock (`setMolFile`) as the SMILES it is, as OpenChemLib and Ketcher do (a host doing so opened it empty)
* Chemical Properties: the calculator is shown as **Chemical Properties** (its function, `addChemPropertiesColumns`, keeps its name and its declared name `Chemical Properties (OCL)`, so saved HitTriage campaigns, flows and the dialog's history still find it); in Crux mode the calculators dialog names its method as RDKit's definitions, calculated by Crux, and MW is the average molecular weight
* Fingerprints: `getFingerprints` and `chemGetFingerprints` take optional fingerprint options, today's fingerprints by default: Morgan `radius` (2, ECFP4), `size` (2048 bits; every type but MACCS), `useChirality` (Morgan, AtomPair, TopologicalTorsion) and `useFeatures` (FCFP, Morgan). Chemical Space, Activity Cliffs and MCL show them in their preprocessing settings. FCFP and chirality are computed by Crux in both chemical engines: Chem's RDKit.js reads no `useFeatures` (it gives ECFP) and, under the new stereo perception Chem uses, gives L- and D-alanine the same chiral Morgan bits and ignores chirality for atom pairs and torsions; Crux gives RDKit 2026.03's bits (FCFP: RDKit's feature-invariant Morgan, bit for bit). A cached fingerprint column names options other than the defaults (`~smiles.Morgan(radius=3)`), so one option set never reads another's; the defaults keep the names they had, and saved projects' caches stay valid. A column keeps the defaults' cache and the last options asked for
* Fingerprints: Fixed a cached fingerprint column missing an edited row when another cache of the column (another fingerprint type, the substructure search's pattern fingerprints) was read first: the first read took the column's edits for all of them; each cache now keeps its own
* Chemical Engine: In Crux mode, the property table (the Properties panel, Chemical Properties' columns, `getProperties`, `getChemPropertyFunction`) is computed by Crux, in its web workers for a column, instead of OpenChemLib's workers, by RDKit's definitions and identical to RDKit's values: MW the average molecular weight (RDKit's MolWt; OpenChemLib's MW is the monoisotopic mass), HBA Lipinski's N + O (Lipinski et al. 1997), HBD the donor atoms (RDKit's NumHDonors), LogP Wildman-Crippen (1999), LogS ESOL (Delaney 2004, with Walters' 2018 coefficients on RDKit's descriptors), PSA TPSA with S and P (Ertl et al. 2000), rotatable bonds RDKit's strict definition, stereo centres the tetrahedral ones assigned or not, molecule charge the sum of formal charges, the formula in Hill's order (1900). The panel names each property's method and reference on its name. Values differ from OpenChemLib's where the methods differ (on ChEMBL-10k: HBA 100 % equal, HBD 99.3 %, stereo centres 99.8 %, rotatable bonds 77.9 %, PSA 85.1 %; LogP and LogS are other methods). A column of 10,000 molecules takes 2.0 s in one thread where OpenChemLib took 28.3 s
* Chemical Engine: In Crux mode, the Drug Likeness panel gives RDKit's QED (Bickerton et al., Nature Chem. 2012: a desirability from 0 to 1) with the eight properties it weighs, and its button adds a QED column; the Toxicity panel lists the alerts of RDKit's FilterCatalog Brenk (105) and NIH (180) sets the molecule matches, each drawn with its atoms, and per set a button adds the alert counts as a column (`Brenk alerts`, `NIH alerts`). The values and hits are the pinned RDKit's on ChEMBL-10k (crux's `chemQed`, `chemAlerts`). RDKit mode keeps OpenChemLib's score and four risks, and OpenChemLib gives them for a molecule Crux cannot read
* Chemical Engine: Added Chem | Calculate | Toxicity Alerts...: the Brenk and NIH alert counts of each molecule, a column per set (Crux mode; RDKit.js has no FilterCatalog, so in RDKit mode it asks for the Crux engine). Toxicity Risks... stays OpenChemLib's four risks in both modes
* Chemical Engine: In Crux mode, notation conversions (Convert Notation, `convertMolNotation`, `canonicalize`, and so every sketcher's getters through `grok.chem.convert`) are Crux's, a column on the Crux workers: RDKit.js 2026.03.6's strings, the same molecules as RDKit mode's (on ChEMBL 100k 99.13 % of canonical SMILES are the same string as RDKit 2024.09's; the rest write stereo marks or atom order as RDKit 2026.03 does). RDKit converts a molecule Crux cannot read, a CXSMILES whose query features Crux's reader skips (`Q_e`-type labels, `rb:`, `m:`, ...), and for code whose next step is RDKit's (reaction drawing, the substructure filter's query names, an R-group core's checks). `isSmiles` and `validateMolecule` ask Crux first; validation messages stay RDKit's
* Chemical Engine: In Crux mode, the Structural Alerts, Pharmacophores and 2D Structure panes and the MMP viewer's pair alignment, highlights and generated molecules (with a SAR matrix's single-cut cells) are computed by Crux, with RDKit for a molecule Crux cannot read, a SMARTS it does not read and a call that fails
* Chemical Engine: In Crux mode, Save as SDF writes the MOL blocks (and the other molecule columns' SMILES) by Crux, as Convert Notation does, on its web workers: RDKit.js 2026.03.6's
* Chemical Engine: In Crux mode, InChI and InChIKeys (the Identifiers panel, To InchI... and To InchI Keys..., `getInchis`, `getInchiKeys`) are Crux's: Crux fills InChI's input from the molecule as RDKit's MolToInchi fills it, in its web workers for a column, and the InChI library runs in the Crux sketcher's InChI module, which Chem already ships (no new file). RDKit.js 2026.03.6's InChIs and keys (crux-core: 102,120 of 102,120 identical), RDKit for a molecule Crux cannot read
* Chemical Engine: In Crux mode, the 3D Structure panel embeds a molecule drawn in 2D by Crux on its web workers, not by the `SmilesTo3DCoordinates` Python script: RDKit's ETKDGv3 in its all-in-one form without Hs, random starting coordinates above 50 atoms, at most 200 attempts, a UFF clean-up and a 10 s deadline. Large molecules the script never finished embed (a 173-atom macrocyclic peptide in ~0.35 s, Abeta1-42 in ~2.5 s), and one that takes longer than the deadline shows a message instead of freezing the page. The script runs for a molecule Crux cannot read
* Chemical Engine: In Crux mode, the 3D descriptors (the Descriptors 3D group) embed each molecule with RDKit's ETKDGv3 (Hs kept, as the container does), random starting coordinates above 50 atoms, at most 200 attempts and a 5 s deadline per molecule: more molecules get 3D values (ChEMBL-1k 995 of 1,000, 980 before; large molecules 18 of 20, 8 before) and a column of large molecules no longer stalls; the 3D values move from the container's by a median 0.1-6 % per descriptor, 2D values are unchanged
* Chemical Engine: In Crux mode, Crux's calls are faster: aligning molecules to a scaffold (the grid's align-to-scaffold, the MMP viewer's pairs) 2.35x, MCS (Cluster MCS, R-group cores, SAR matrix anchors) 2.2x, Rascal MCES (activity cliffs) 1.76x, R-groups 1.23x, reactions and MMPA 1.15x, every other call measured 1.06-1.17x, with the same results: CoordGen, FMCS, Rascal and the Crux file's own code are built for speed. Against RDKit mode on the same inputs: MCS 3.3-5.6x, MMPA 2.2x, reactions 2.1x, R-groups 1.6-1.8x faster. The Crux file is 2.16 MB gzip, past RDKit_minimal's 1.99 MB
* Chem: `embed3d`, a molecule with 3D coordinates (a MOL block) by the Chemical Engine in effect: Chem's one entry point for 3D coordinates, which the 3D Structure pane and Bio's macromolecule 3D pane use (Crux's embedding with a deadline in Crux mode, the SmilesTo3DCoordinates script in RDKit mode)
* SAR Matrix: cells align to their row's core with RDKit's depictor instead of CoordGen, about 20 times faster per cell (9 ms to 0.4 ms): switching to a 1,000-cell series takes ~0.05 s instead of ~0.8 s (seconds with the browser's developer tools open, where it froze); the cores line up as before
* Chem: `chemicalEngine`, the Chemical Engine setting in effect (RDKit or Crux), for tests and other packages
* Convert Notation: Fixed the RDKit workers keeping molecules they no longer used after a conversion, and a Kekulé conversion (`kekulize`) turning the molecules the substructure search had cached into their Kekulé forms
* Activity Cliffs: an MCS highlight computed by one chemical engine is no longer shown after a switch to the other; each engine computes its own
* Similarity search: in Crux mode, a query molecule whose fingerprint Crux fails to compute gets RDKit's (it failed as malformed)
* Tests: `chemical engine`: conversions, `isSmiles` / `validateMolecule`, the Structural Alerts and Pharmacophores panes and the MMP viewer's joins in Crux mode against RDKit mode on Chem's molecules (SMILES and MOL blocks), and break-it cases: explicit H, query MOL blocks, dative bonds, R-group cores, enhanced stereo, CXSMILES query features, salts, isotopes, radicals, SMARTS, unreadable and over-long strings, a column's chunk boundaries with molecules Crux cannot read, a wasm trap in the panels and conversions, and RDKit's search after a Kekulé conversion
* Chemical Engine: In Crux mode, reactions (Transformation and Two-Component reactions, Deprotect, Reaction Enumeration) run on the Crux module, with RDKit behind it per molecule (per run for an enumeration): the same products as RDKit mode on Chem's reaction tests, its default reactions and the enumeration demo (510 rows; the amino acids' 463). Reaction drawing stays RDKit's in both modes
* Settings: **Chemical Engine** (RDKit, the default, or Crux): in Crux mode the molecule renderer, `drawMolecule` and `canvasMol` draw through Crux (highlights and scaffold alignment included), with RDKit drawing any molecule Crux cannot read or a call it refuses; descriptors (the dialog's tree and the values: the menu, the panel, `chemDescriptors`, the vector function) are calculated by Crux in web workers instead of the chem-chem container, and the Scaffold Tree is built by Crux in a web worker instead of the Python script, each with the server path as its fallback. The Crux file is fetched only in Crux mode (and for `molzip` and FCFP or chirality fingerprints). Other packages get RDKit objects from Chem in both modes (`getRdKitModule`, `getRdKitService`); in Crux mode the pictures and values of the moved functions come from Crux
* Tests: Added BDD scenarios for the Included in / Not included in pair of the substructure card, Align and Highlight in the card's sketcher, and Ketcher picked in the filter card's sketcher (GROK-12581, GROK-12905)
* Chemical Engine: In Crux mode, a molecule whose descriptors Crux cannot calculate (a row of empty values, e.g. a V3000 molblock Crux's reader rejects) gets the container's values, and the fallback is logged, instead of empty cells
* Crux: Chem's Crux engine is one WebAssembly file, `vendor/chem-wasm` (crux-core's `chem-wasm`: rendering, search, descriptors, the scaffold tree), compiled once and shared by the main thread and the Crux workers; it replaces the search-only `src/crux/crux_wasm`
* Crux: Chem's Crux file is 88 KB smaller to download (1,976,315 → 1,888,161 bytes gzip), every result the same; a MOL block read with kekulization off (the retry for a block that cannot be kekulized) keeps the atoms of a bond already marked aromatic as RDKit does (Crux made them aromatic, and a substructure search then missed RDKit's hit)
* Substructure search: Contains / Not contains runs on the engine the **Chemical Engine** setting picks: Crux in Crux mode, RDKit in RDKit mode (the default), with the same hits; the separate **Substructure Search Engine** setting is gone, and a value a server saved for it is not read. In RDKit mode the first search of a column takes RDKit's time again (ChEMBL 100k: 9.8 s where Crux took 0.7 s)
* Substructure search: A Crux search starts no RDKit worker: the Crux workers read a molblock column's molecules themselves (they were converted to SMILES on the RDKit pool first); the RDKit pool starts only for a query Crux cannot express or a search that falls back
* Substructure search: Crux searches every query it reads: it takes the query molecule Chem builds for RDKit (`getQueryMolSafe`) and matches what RDKit matches, so the queries that ran on RDKit's workers (explicit `[#1]` and `!#1` as PAINS writes them, `R<n>`, `r<n>`, `v`, `h`, `x`, isotopes, radicals, atoms that may be an H such as `!C` or `*!O`, a single-or-double ring bond drawn in a sketcher: 130 of the 1,182 structural alerts and drawn queries of Chem's census) run on Crux with RDKit's hits, rings as RDKit defines them (the SSSR of the sanitized molecule). A molblock column's explicit H atoms are kept, as RDKit keeps them (`[#1]`, `D`, `X` and `h` see them), and a query of no atoms matches nothing, as in RDKit (Crux matched every row). A CXSMILES query whose extension holds query features Crux's reader skips (query atom labels such as `Q_e` or `AH_p`, positional variation `m:`, `rb:`, `u:`, `s:`) runs on RDKit
* Substructure search: The RDKit search gives RDKit's hits for ring counts and ring sizes (`R<n>`, `r<n>` and their ranges), hybridization (`^`) and radicals. A molecule its workers did not hold in their cache (past 10,000 molecules per worker, a row edited since the last search, saved fingerprints after the cache moved on) was read unsanitized for speed: its rings from a depth-first search (naphthalene a 10-ring and a 6-ring, so `[c;r10]` hit it and `[c;R2]` took 6 of its atoms), no hybridization, no radicals. A query that asks for them, as its query molecule says, reads such molecules sanitized; every other query reads them as before. Read that way, 6 of Chem's structural alerts hit ChEMBL-10k otherwise than RDKit (BMS 304: 17 rows instead of 8)
* Rendering: The molecule renderer frees the molecules it evicts from its cache (it looked for them in the wrong place and kept them all)
* Chemical Engine: In Crux mode, fingerprints (similarity and diversity search, chemical space, activity cliffs, BitBIRCH, the SAR matrix's clustering, `getMorganFingerprints`, `getFingerprints`) and structural alerts and pharmacophore features are computed by Crux in web workers, with the same types and options; a cached fingerprint column belongs to the engine that filled it, and saving one drops the other engine's. RDKit computes them in RDKit mode, for the molecules Crux cannot read, and wherever Crux fails
* Chemical Engine: In Crux mode, an activity cliff's tooltip and panel highlight what RDKit's Rascal MCES of the two molecules (Crux, in a web worker) leaves out of each, instead of what FMCS leaves out: the most bonds the two share, possibly in several fragments, where FMCS finds one connected substructure. Every cliff of Chem's tests and demo gets one (Rascal's options: similarity threshold 0.3, up to 10,000 matching bond pairs, a 5 s timeout past which Rascal gives its best MCES so far), in milliseconds to seconds where FMCS took 13 s and 88 s (triterpenoid esters) or did not finish (160-atom peptides). RDKit mode and the other MCS calls are unchanged; FMCS does it when Crux fails
* Chemical Engine: In Crux mode, MMPA's fragmentation (MMP analysis, the SAR matrix's series and clusters), R-group decomposition (R-Group Analysis, `rGroupDecomposition`, the SAR matrix's clusters) and FMCS (R-Group Analysis' MCS, Cluster MCS, the SAR matrix's anchors) run on the Crux workers: RDKit's worker code on the Crux module, with the RDKit pool behind it per call, per molecule Crux cannot read and per cluster. The results are RDKit 2026.03's (RDKit.js 2026.03.6 through the same code gives the same fragments of the MMP demo's 20,267 molecules, the same decompositions in all four matching strategies and the same MCS of 128 clusters); where Chem's RDKit 2024.09 differs, RDKit changed between the releases (a fragment keeps the E/Z and chirality marks at its cut; an MCS SMARTS names its atoms in another order, the same query). A cluster's MCS on Crux stops at FMCS's own 45 s timeout and gives none, as RDKit mode gives none for a cluster whose worker it restarts after 45 s. Fragmenting the MMP demo takes 39-42 s in one thread where RDKit takes 55-61 s
* Chem: `molzip`, a package function for Bio and other packages: RDKit's molzip, by Crux in both chemical engines (Chem's RDKit.js has none), joins the pieces of each molecule at their attachment points (`[*:1]` ... `[*:1]` by default; RDKit's MolzipParams by name) and gives the joined molecules' canonical SMILES. In Crux mode the SAR matrix joins cores and fragments with it instead of assembling SMILES strings
* Tests: Added the `crux analyses` category: MMPA fragments, R-group decomposition, FMCS, Cluster MCS, the SAR matrix's clusters and fragment linking from the Crux workers against RDKit mode (no RDKit worker asked, no fallback; FMCS against RDKit 2026.03.6's SMARTS), the fallbacks for a molecule or cluster Crux cannot read, `molzip` in both engines against RDKit 2026.03.6, and the timings
* Tests: Added the `chemical engine` category: RDKit mode never asks for the Crux file; switching both ways and while the file loads; a failed load falls back to RDKit; the renderer in both engines (empty, invalid, SMILES, V2000, V3000, SMARTS, highlights, alignment, very large molecules), its caches per engine, freed evicted molecules, the Crux heap over thousands of renders; descriptors and the scaffold tree from Crux against the chem-chem container and the Python script; the Crux search default against RDKit. `misc/crux/engine-setting.sh` switches a stand's engine settings for whole runs in Crux mode
* [core#874](https://github.com/datagrok-ai/core/issues/874): Scaffold Tree: Fixed an error when the viewer was opened while the current view was not a table view
* Scaffold Tree: Generate is no longer blocked on large tables (it was from 500 molecules): a column of more than 750 molecules is generated from a sample of 750 spread over the whole column, the same rows each time, and Generate's tooltip and the progress say so ("Generate from a sample of 750 of 12,345 molecules")
* [core#875](https://github.com/datagrok-ai/core/issues/875): Analysis: Fixed Activity Cliffs and Chemical Space adding their scatter plot to a different table with the same name when the original view was closed during the run
* Tests: `Crux sketcher: strokes in a sketcher in the context panel leave the current object as it is` and the BDD scenario "Strokes in a molecule column's Filter pane leave the column the current object" (`sketcher/crux-hosts.feature`): js-api's sketcher made each drawing in the Filter pane the current object, so the panel switched at every stroke
* Substructure filter: Fixed a query typed as SMILES turning into SMARTS on a card restored from a layout or a project: its sketcher handed the molecule back as a molblock, which replaced the typed string
* Chemical Engine: In Crux mode, descriptors calculate each distinct molecule of the column once, and the values come from the workers as typed arrays made into columns without converting value by value (the same columns; on 500,000 rows of 1,000 molecules, 20 descriptors take 1.2 s instead of 48 s)
* Tests: Moved the BDD features that clear Filter as you draw (an account setting) or filter while the sketcher is open to the serial lane: in a full run, a query round trip filtered by the previous stroke while another worker had the option cleared
* [#4291](https://github.com/datagrok-ai/public/issues/4291): MPO: Add support for qnum values
* [#4324](https://github.com/datagrok-ai/public/issues/4324): MPO: Fixed Cancel in the property curve settings keeping the changed Score as
* [#4323](https://github.com/datagrok-ai/public/issues/4323): MPO: Fixed Cancel keeping the score column the preview added, and OK leaving the columns of previously picked profiles
* Tests: Added BDD features for MPO qnum scoring by the Score as option ([#4291](https://github.com/datagrok-ai/public/issues/4291)) and for Cancel in MPO Score ([#4323](https://github.com/datagrok-ai/public/issues/4323), [#4324](https://github.com/datagrok-ai/public/issues/4324)); what they cannot claim yet is in `bdd/features/mpo/MISSING.md`
* Tests: Added the `MPO: qnum scoring` package tests (`src/tests/mpo-tests.ts`): the Score as matrix (worst / limit / best) on a qnum column, missing-value strategies, an inverted and a log10 curve, a bound past the curve end, and the smallest exact on a flat curve — the numbers the BDD features defer to a package test (MISSING.md item 6)
* Tests: `sketcher testing: a sketcher fills its dialog`: every sketcher, in a dialog as it opens, in one shown with a size and with the invalid-molecule warning shown, reaches the bottom of the dialog's body and is not covered there (js-api's in-place sketcher left a margin under it)
* Substructure search: A query molblock with an uncharged pentavalent nitrogen (a nitro drawn `N(=O)=O`, an azide `N=N#N`, a pyridine N-oxide `N=O`, from any sketcher) is read sanitized, as one with an H count query is: RDKit makes the group charge-separated, as the molecules searched hold it. Read unsanitized it matched no row (a nitrobenzene: none of smiles.csv's 75). Every other query molblock is read as before
* Substructure search: An H count query on a ring NH in a Kekulé V2000 molblock (OpenChemLib's indole, `hhh` 2) matches its rows: RDKit's V2000 reader gives an atom with an `hhh` no implicit H, so that read sanitized the nitrogen was a radical and its ring not aromatic, and the query matched nothing (its V3000 twin, `HCOUNT`, matched). Such an atom, uncharged, is given the H its query asks for at least (its valence field), when its bonds are single, double or triple and RDKit can then sanitize it
* Tests: `query back and forth: complex queries`: a nitrobenzene and methyl azide drawn with a pentavalent nitrogen in Crux, Ketcher and OpenChemLib, and a nitro, an azide and a pyridine N-oxide so drawn matching the molecules that hold them charge-separated; OpenChemLib's indole with an H count on its NH. `sketcher testing: Ketcher fills the width of its dialog`
* Similarity Search: Fixed Row Source being ignored: Filtered, Selected and FilteredSelected searched every row. The cards come from the rows of the row source and follow the filter or the selection only in the modes that read them; Select all similar selects within the row source (`similarity/diversity: row source` tests)
* Similarity Search, Diversity Search: Row Source defaults to Filtered
* Matched Molecular Pairs: Fixed a second analysis in a page never finishing (its tabs not built, the GPU kept busy) on Intel's Windows driver, and the GPU never comparing the last two molecules (math 1.2.10's WebGPU workgroups and pair ranges); a GPU search that cannot finish falls back to the CPU (`mmpa: gpuPairsAsCpuTwice`)
* [public#4226](https://github.com/datagrok-ai/public/issues/4226): Menu: Similarity Search, Diversity Search, R-Groups Analysis, Scaffold Tree, MPO Score, Chemical Properties and the Transformation and Two-Component reactions are disabled until the current table has a molecule column, and Activity Cliffs, Chemical Space, Matched Molecular Pairs, SAR Matrix and Deprotect unless the current view is such a table, instead of failing or warning
* Sketcher: A SMARTS no syntax tells from SMILES (Daylight's `[NX3][CX3]=[OX1]`, an any bond `~` as its only mark, `[R]`, `A-A`) is set as SMARTS: entered in the field it was set as SMILES, which no sketcher could read, so it drew nothing and filtered nothing; and a CXSMILES with atom labels (`$Q_e$`), an AND group (`&1`) or coordinates is a molecule, not SMARTS: entered, it reached the sketcher as SMARTS, its stereo group and labels lost, and a substructure search by one searched it as SMARTS (`[C@H]` an H count, `C` aliphatic only: 7 rows of smiles.csv where its SMILES finds 64). RDKit is asked only about a text outside SMILES syntax, so a SMILES costs no parse: 100,000 strings of the demo files take some 50 ms, where a parse of each took 20 s (`isSmarts`)
* Copy SMILES, Copy as SMILES: A CXSMILES is copied as RDKit's SMILES, without its extension: one with an OR group or ring flags (`|o1:1|`, `|r|`) was copied as it stood. Copy as SMARTS gives RDKit's SMARTS of a CXSMILES with a stereo group, labels or coordinates, which it gave as it stood
* Tests: `isSmarts`: SMARTS of primitives and CXSMILES, the SMILES-syntax prefilter over Chem's alerts, its pharmacophore features and the demo files, `benchmark.PackageFunctions.isSmarts.100000`, and a test for each caller (the sketcher, `DG.chem.isSmarts`, column highlights, the panels, Copy SMILES, Copy as SMILES and as SMARTS, substructure search, reaction depiction)
* Substructure search: A query molblock with an H count query on an atom (MDL's `hhh`, V3000's `HCOUNT`, as Crux and Ketcher write it) is read sanitized but not kekulized (RDKit's own reading, its rings aromatic as RDKit perceives them; one it cannot sanitize is read unsanitized as before). Read unsanitized, its Kekulé ring stayed Kekulé and matched no aromatic ring, so such a query filtered no rows. Every other query molblock is read as before, unsanitized: sanitizing changed the rows of some (a valence the drawing fixes read as a radical, a pentavalent nitrogen read charge-separated)
* Notation conversion: A molblock converted from SMILES with stereocentres has its chiral flag set: its stereo is absolute, as the SMILES's is. Without it Ketcher took the centres as relative, and the molecule cell editor wrote an edit back with an AND group (`C[C@H](N)C(N)=O |&1:1|`)
* Substructure filter: A query kept in a layout, or restored after its card was reattached, is set back in its notation (a SMARTS or a SMILES entered as such): set as a molblock, the card's sketcher opened empty
* Tests: Added the `query back and forth` categories (complex queries, entered as text, non-realistic drawings, routes): queries of several features drawn by Crux, Ketcher and OpenChemLib, SMARTS, query molfiles and CXSMILES entered in the field, drawings no chemist means, and the routes (reopened and OK unchanged, every sketcher in turn, edited after a reopening, kept in a layout, cleared, the molecule cell editor), each judged by RDKit's matches on mol1K, SPGI or smiles.csv
* Deprotect: Fixed the default Fmoc fragment taking the protected amine's nitrogen with the group (Fmoc-Gly came out as acetic acid, Fmoc-Pro unchanged): its attachment point is now the nitrogen itself, which keeps its hydrogens
* Substructure filter: Fixed a structure card that stopped filtering after Remove All on another view's panel of the same table: the removed card, deactivated first, kept its peers disabled (Use as filter in that view let every row through)
* Cluster MCS: A one-molecule cluster's MCS is the molecule itself with its hydrogens: read as a query, its NH hydrogens were dropped (a urazole came out as the triazolinedione)
* Tests: Added package tests for activity cliffs at the cutoffs 20 and 95, elemental analysis on SMARTS, Chemical Space on an all-empty column, R-groups analyses that find none, BitBIRCH, the similarity matrix and Cluster MCS (`clustering`), the Transform reactions (`Reactions`), the notation round trip over smiles.csv, similarity search on the other fingerprints and the mixture cell renderer; the BDD features they replace keep only their UI gestures
* Substructure filter: Fixed a query marked with a feature RDKit's sanitizing drops (toluene drawn, then its methyl bond marked aromatic) keeping the rows of the query before it, and that query as the filter's state: the filter names a query by the query its search builds (read unsanitized), not by RDKit's sanitized molecule, which took the two for one and skipped the second search
* Substructure search: A V3000 query molblock without explicit hydrogens is read as a query, unsanitized, as a V2000 one is: an aromatic bond outside an aromatic ring is no longer read as single
* Substructure filter: The generic atoms X, XH, M and MH are kept as drawn instead of being rewritten into atom lists: RDKit reads them with their own meanings; the lists lost astatine and the metals beyond 38, and those of XH and MH matched no explicit hydrogen
* Crux substructure search: A query whose bond is single-or-double with a ring or chain topology is searched by RDKit: RDKit writes it `-,=&@`, which a SMARTS parser reads as single anywhere or double in a ring
* Sketcher: OpenChemLib's SMARTS is the query its molfile makes, read unsanitized (an aromatic bond outside an aromatic ring is kept)
* Sketcher: The Crux sketcher's status reports its SMARTS (`smarts`)
* Tests: Added the `query round trip` category: every query feature Crux and Ketcher offer, drawn in the substructure filter on mol1K, SPGI and smiles.csv, filters the rows of its intended SMARTS and comes back the same reopened and through the other sketcher; and the BDD features `filters/query-roundtrip-*` for the report's case on each dataset

* GROK-21063: Chem Space, Activity Cliffs: Fixed an error when the table is closed while the dialog is opening
* Tests: Added the `save as csv` test category (MOLBLOCK columns exported as SMILES read back as molecules); the BDD features no longer repeat the notation permutations and renderer checks the package tests cover
* GROK-21049: Activity Cliffs: Fixed the "Column already exists" error when run alongside another embedding analysis
* Crux substructure search: At most 16 crux workers: each holds a WebAssembly instance, and a page has room for only about 120 WebAssembly memories, which crux's pool and RDKit's (both `hardwareConcurrency - 2`) filled on a 32-thread machine, so that every WebAssembly start after them failed with "Out of memory" (the Crux sketcher, DBSCAN, RDKit). A crux start that failed lets go of the workers it started (each search retried it, and the page piled up hundreds of workers); a second RDKit init (a test category's) keeps the RDKit worker pool instead of leaving the old one running
* RDKit: At most 16 RDKit workers (`hardwareConcurrency - 2` before), as crux's pool: a page has room for only about 120 WebAssembly memories, and on a machine with many threads the RDKit pool would fill it by itself
* Sketcher: Added Crux, the Crux Sketch molecule sketcher on the Crux chemistry engine (WebAssembly), to the choices of the Sketcher setting (crux-sketch 0b22489, engine bceb40f00e97; vendored in `vendor/crux-sketch/`, loaded at its first use); in the substructure filter it opens in query mode
* Sketcher: Fixed the substructure filter losing a structure typed into its sketcher dialog when OK was pressed before the Crux sketcher had loaded ("Filter as you draw" on): Crux, its dialog closed while it loaded, now takes the value the host stored meanwhile and says so, as a sketcher closed after it was ready does
* Sketcher: Fixed a stroke drawn in the Crux sketcher being lost when the page moved focus while it was held (a table view focuses its grid about a second after it opens, which could fall in the middle of a stroke in a molecule cell's editor, and OK then left the cell as it was): only Tab and Shift+Tab cancel a held stroke now
* Sketcher: The Crux sketcher's root is named `cruxSketch` for the platform's automation, and its status reports its canvas, toolbars and label editor, each tool's button (`tool <name>`), and what it holds (the selected atoms and bonds, the tool, whether it holds a query or nothing, the changes it has said); it fires `onRendered` after each drawing it shows. Crux keeps its look under a strict Content-Security-Policy (no `'unsafe-inline'` for styles), which it needs nothing of beyond `'wasm-unsafe-eval'`
* Sketcher: Fixed an edit made in a molecule cell's editor right after its sketcher was ready being undone (OK then left the cell as it was): the editor gives the sketcher the cell's own SMILES back the moment it is ready (`Sketcher.sketcherReady()`), not by polling for it; the substructure filter decides its sketcher's mode once its card is in the page, instead of polling for it for up to 1.5 s
* Sketcher: Fixed a molecule typed into a molecule cell's editor and entered while its sketcher was still loading being ignored on OK (the cell kept its value): the editor gives the sketcher the cell's own SMILES back once it is ready only while the sketcher still holds the cell's molecule, not over a value set meanwhile
* Substructure filter: Fixed a structure filter restored by a layout sometimes leaving the table unfiltered: a filter closed within a second of getting its state (a cloned view's filter panel closed at once) still ran the search its state had scheduled and took the column over, so the filter restored later took that search for its own and never ran it; detaching a filter now cancels its scheduled search
* Sketcher: The 400px minimum width given back to a sketcher dialog opened from a popup no longer applies to a sketcher with no minimum size
* Demo: Fixed Similarity & Diversity Search, Matched Molecular Pairs and Scaffold Tree returning before their layouts were applied, which then landed on whatever view was open next
* Reaction Enumerator: Added "Repeat until it stops" with Max cycles: single-reactant templates are re-applied to their own products within one step, and only the end products are kept (e.g. a triene gives the fully reduced product in step 1)
* Activity Cliffs: The molecules of the pair in the Cliff Details pane carry `role="button"` and `aria-label` ("molecule of row N"), so assistive tools and tests can name them
* Reaction Enumerator: Fixed columns of the reaction templates, building blocks and reagents grids collapsing when another column is resized; columns are now sized once to fixed widths instead of being fitted to the available space
* SAR Matrix: Added building the matrices from core and R-group columns already in the table (Use existing R-groups); attachment points are read in any common notation ([*:n], [n*], R# molfiles, CXSMILES labels), a bridge written in two columns is joined as one piece, and a warning names R-groups that do not rebuild the table's molecules
* SAR Matrix: Fixed predicted structures inverting a stereocenter or double-bond geometry next to an attachment point
* SAR Matrix: A series whose largest connected block has a single row or column now shows its largest block with at least two of each, instead of no matrix or a one-column one
* SAR Matrix: The Core header marks the varied position on the series core, and the context panel shows the core and the row's own R-groups
* Matched Molecular Pairs: Fixed missing values of an integer activity becoming -2147483648 when the activity is log-scaled

* Substructure Search: Added the Crux engine (Substructure Search Engine package property) for Contains / Not contains searches; queries it cannot express run on RDKit
* Substructure Search: Improved Crux parity with RDKit — nitro / N-oxides, perchlorates, elements beyond Rn and aromatic rings RDKit reads only without Kekulize are read as RDKit reads them (Chem datasets: 67 → 0 molecules Crux could not parse, 277 → 7 read differently)
* Substructure Search: Fixed typed queries such as `[OH]`, `[CH3]`, `[N+]`, `c1cc[n+]cc1` finding almost nothing — a query whose SMILES reading is a radical is now read as SMARTS, by both engines
* Substructure Search: Added a fallback to RDKit when a Crux search fails (e.g. a worker runs out of memory); Crux restarts to free its memory
* Improved RDKit start-up: the RDKit wasm is compiled once on the main thread and shared with the workers, and in Crux mode the substructure filter no longer starts the RDKit workers
* Substructure Search: A second run no longer asks for the molecule column, which offered the hidden canonical SMILES column the first search leaves
* Reaction renderer: Replaced inter-step arrows with numbered panels and improved scaffold alignment using cached molecular coordinates
* Substructure filter: Fixed the search progress staying in the task bar at 100% — a search's end closes it, a detached filter terminates every search it started, and the fingerprint precalculation clears its own entry
* GROK-20955: Names To Smiles: Fixed the run failing on a table that already has a canonical_smiles column; the structures go into a column of their own
* GROK-18286: Scaffold Tree: Fixed removing a scaffold under a colored one taking the colors column away, and with it the coloring of a plot colored by that column
* Scaffold Tree: Fixed a removed colored scaffold keeping its color: the removed node was still redrawn once and registered its color again, so the next update brought its highlighting back
* SAR Matrix: Added a Summary tab, the one the viewer opens on: the scale and coverage of the analysis, which R-group, core and measured swap buys the most (with a single cross-validated additive fit ranking every component where the parts are columns), the series to open first, and the untested compounds and virtual analogs worth making; every row lands on its cell
* SAR Matrix: Added the PROTAC SAR Matrix demo: 2,792 patent degraders split into warhead, linker and E3 ligand columns
* SAR Matrix: Fixed a column picked as both the core and an R-group: the dialog drops it from the R-groups and the function refuses it
* SAR Matrix: Docked the viewer as a tab beside the table's grid
* SAR Matrix: The matrix's Columns input picks which R-group runs across the columns; the Vary filter is gone
* SAR Matrix: The analysis rebuilds when rows are added or removed, and follows a renamed molecule, activity, core or R-group column instead of losing it
* SAR Matrix: A proposal the table already holds (under another cut, core spelling, salt form or series) is shown as that compound, untested or excluded, instead of as a new analog; salts are cut to their parent before fragmenting
* SAR Matrix: A project saved before this version rebuilds its matrices once when opened (the saved set now carries the settings that built it)
* Similarity/Diversity search: Each molecule card now announces the row it shows (`data-row`, `name="card-<row>"`), and the viewers report the `search-results` status of `@datagrok-libraries/ml` — `card <row>` hit areas plus `cards` / `current card` / `selected cards` readings — instead of being addressed by position; `isRenderPending` / `onRendered` cover a scheduled or running render
* Generate Conformers: Fixed the run failing with `AttributeError: Cannot set unknown attribute 'maxAttempts'` — ETKDGv3 takes `maxIterations`
* Chem | Calculate | IUPAC Name: Added a Python script that names molecules with openclatura
* Substructure filter: Added the card's structure, search type, fingerprint, similarity cutoff, searching state and whether it is still drawing (the cards below it move until it has) as readings of the filter panel it sits in
* Similarity and Diversity Search: Added a `chem-search` status beside `search-results` — metric, fingerprint, size, row source, header, the set of rows the cards show, card sizes and properties; Similarity adds the target row, cutoff and scores
* Scaffold Tree: Added status readings of its nodes (count, checked, coloured, each node's scaffold, hits and colour), its message and why generation is blocked, hit areas for each node and its icons, and `aria-disabled` on a blocked icon
* Matched Molecular Pairs: Added status readings of its activities, current tab, substitutions, pairs and generated molecules
* Activity cliffs: Added the `only cliffs` reading beside `cliffs`, and named the cliff-count button `button-cliffs`
* R-Groups Analysis: Fixed OK with an empty sketcher removing the latest analysis before reporting that no core was provided
* Added BDD features for the Chem section of TestTrack
* GROK-20753: Added `moleculeFilterOperators` (`meta.role: filterOperators`) — substructure, superstructure, exact, stereo-agnostic and similarity operators for the u2 filter builder
* GROK-20808: MMP: Fixed a crash when no table was selected in the dialog
* GROK-20829: Guard isBitsetStale against the valueless root tree group
* GROK-20824: MPO: Fixed the profile description picking up the title's text styling when styled text was pasted into it
* GROK-20821: MPO: Fixed the dataset list ignoring newly opened tables
* GROK-20823: MPO: Fixed Edit opening duplicate editor tabs and leaving a stale profile preview in the context panel
* GROK-20832: MPO: Fixed Save being disabled on a freshly cloned profile — the copy could not be kept
* GROK-20831: MPO: Fixed truncated profile descriptions not expanding on click after saving a profile
* GROK-20830: MPO: Fixed data-driven Save not using the name and description typed on the tab
* GROK-20805: Chem: Butina cluster doesn't work in a project with datasync

## 1.17.15 (2026-09-02)

* Added SAR Matrix, rendered in a virtualized grid, with SAR Transfer to the grid
* SAR Matrix: Fixed the SAR Transfer tab rebuilding the whole panel (grid flicker, list scroll jumping to the top) on expanding/collapsing a series or selecting a transfer — the list now toggles and swaps in place like the SAR Matrix navigator
* Fixed the package build failing on `TS2610` — `name` is an accessor on the u2 `Component` base, so the viewer overrides it with its own accessor instead of redeclaring it as a property
* Fixed Gasteiger Partial Charges panel on current RDKit — `GetSimilarityMapFromWeights` now requires an explicit `MolDraw2DCairo` drawer (passed by keyword, so pre-2023.09 RDKit still works)

## 1.17.14 (2026-08-13)

* Reaction Enumerator: Added package settings (category "Enumeration") letting admins set per-group default file paths for reaction templates, building blocks and reagents; missing files show a warning balloon and fall back to the bundled ones
* Substructure filter: Added a 'Filter as you draw' checkbox to the sketcher dialog — when unchecked, filtering is deferred until OK is clicked instead of running on every sketch edit
* Added `Filter by Substructure`, `Similarity To` and `Diverse Subset` — table-aware twins of `searchSubstructure`, `getSimilarities` and `getDiversities`, which take a bare column. They declare `(table, column {semType: Molecule}, …)` and return rows or a column added to the table, so a caller has something to carry on with instead of a detached frame or a boxed BitSet
* Added `Apply Reaction` — the non-interactive twin of the Transformation dialog, applying a one-component reaction SMARTS to a molecule column
* Added `To SDF`, which serializes a table to SDF text; the only SDF write path was a zero-argument file exporter reading the current table
* Added `Chemical Space Columns` — Chem Space with typed, enumerated parameters that never plots and returns the X, Y, cluster and cluster-MCS columns it added
* MPO Score by Profile: Refuses to run when a scored property has no column, naming the ones it couldn't resolve — the underlying `computeMpo` only warns and scores over the rest, which returns a plausible number computed from fewer properties than requested
* Added `MPO Score by Profile`, `getMpoProfileNames` and `getMpoProfileProperties`: `mpoCalculate` scores columns that already carry desirability tags and `mpoTransformFunction` wants those tags as a JSON blob, so neither could be driven from a profile name. The new function also takes a property→column mapping, so a profile can score a table whose columns are named differently
* Parameters with a default now declare `initialValue` — the metadata generator reads that, not the TypeScript default, so `aggregation = 'Average'` and friends were registering with no declared default at all
* Chemical Similarity Search: Added choices and defaults for metric, fingerprint, max hits and min similarity, and tagged the query molecule `semType: Molecule`
* Substructure Search: Tagged the query `semType: Molecule` and defaulted the molblock fallback
* Curate: Flipped normalization, reionization, neutralization and main-fragment on by default — with every flag off the function ran and changed nothing
* Fixed display names that were declared via `//friendlyName:`, which the metadata parser ignores: Remove Water and Salts, Run Reaction, Two-Component Reaction
* Remove Water and Salts: Its `molecules` parameter never carried `semType: Molecule` — the decorator options were nested one level too shallow, so it accepted any string column
* Marked genuinely-required parameters non-nullable across the reaction, clustering, scaffold and R-group functions (Two Component Reaction's second table/column, Cluster MCS's cluster column, FindMCS, FindRGroupsWithCore, Murcko, Butina, Curate)
* Docker (chem/chemprop): Cleared reported CVEs — pinned Flask/Werkzeug/gunicorn to fixed releases, upgraded pip/setuptools/wheel (build + runtime tooling), pinned urllib3/idna in chemprop, and added `apt upgrade` to patch base-image OS packages
* Docker (chem/chemprop): Raised security floors (VEX) — torch/torchvision/torchaudio to 2.6.0+cu118 via PyPI wheels (CVE-2025-32434) and removed padelpy's bundled PaDEL-Descriptor jars (log4j 1.2.15, guava 17.0)
* Moved the Chem Playwright E2E suite into the package (playwright/); helpers from @datagrok-libraries/test/src/playwright
* GROK-20239: Chem: OCL renderer: Structures do not appear in labels in some visualisations
* Scaffold Tree: Colors/labels columns are now named `<molColumn> colors/labels (<id>)` using a persisted per-viewer id instead of the editable title, avoiding column-name collisions when viewers share a title
* GROK-20505: Chem: Some descriptors return string values instead of numeric
* InChI / InChI Keys: Moved the conversion to the parallel RDKit service (web workers) instead of the single-threaded module — `To InchI`, `To InchI Keys`, `getInchis` and `getInchiKeys` are now async
* GROK-20486: Fixed a crash in _demoMMPA when the viewer didn't load in time

## 1.17.13 (2026-06-08)

* Chem: MPO: Optimized score calculation for large datasets
* GROK-18818: Chem: MMP: do not allow to choose empty column in dialog
* GROK-15991: Invalid argument(s): Array lengths differ
* Chem: Fixed `waitFor` clearing the wrong timer handle, leaving the polling interval running after the condition was met
* OCL worker: Added missing `break` statements so a drug-likeness or molfile-to-V3K request no longer also runs and posts later cases' results
* OCL worker: Fixed an off-by-one atom loop that could miscalculate molecule charge
* RDKit worker: Fixed the R-group analysis handler posting a stray `undefined` result due to a shadowed variable
* Identifiers: Fixed the UniChem home-page link mapped to the base-ID URL instead of `src_url`
* Map Identifiers: Fixed `column ids.key does not exist` when converting from InChIKey or via the ChEMBL bridge — `_chemMapViaQuery` now renames a cloned input column to `key` to match the converter SQL parameter name
* Chem search: Replaced a floating promise with proper `if`/`.catch()` so substructure-search sub-task rejections are handled
* RDKit service: Worker `terminate()` now rejects in-flight calls instead of leaving callers hung forever; guarded the cluster-MCS timeout path accordingly
* BitBIRCH clustering: Guaranteed the clustering worker is terminated even if the request fails synchronously
* Highlight widget: Unsubscribe `ItemsGrid` event subscriptions when the grid is rebuilt or the widget is detached

## 1.17.12 (2026-06-05)

* GROK-20193: SDF import: Preserve molecule name and comment fields
* MPO: Fixed disabled OK button with no explanation when creating a new profile
* [#3811](https://github.com/datagrok-ai/public/issues/3811): Chem: Added "Export as SVG" action for molecule cells
* Scaffold tree: Improved project loading performance for large datasets and trees — child-node substructure searches are now narrowed to the parent scaffold's matching rows (structure→superstructure relationship), avoiding redundant full-column searches

## 1.17.10 (2026-05-28)

* Chem: Record all generated reaction results in the enumerator

## 1.17.9 (2026-05-26)

* Scaffold tree: Remove the lingering `… colors` fragments column when the last color-coded node is deleted one-by-one (previously only cleaned up on viewer detach or "Drop all trees")
* Chem Space: Guard editor entry against null current dataFrame / missing Molecule column to prevent `NullError: method not found: 'dart' on null` from `setColumnInputTable`
* Chem: Toxicity Risks: Avoid 'Column already exists' crash on repeat invocation by using `getUnusedName` before adding risk columns
* GROK-20108: Docker API: Throw a clear error when the `chem-chem` container is not available, instead of `Cannot read properties of undefined (reading 'id')`
* Scaffold tree: Fixed severe load/scroll slowdown on large dataframes — molecule rendering and substructure searches now run only for on-screen nodes (IntersectionObserver rooted at the actual scroll container) instead of the entire tree
* Scaffold tree: Deduplicated per-node bitset computation (reuses in-flight and completed searches) and added a one-time shared fingerprint warm-up, eliminating redundant full-column searches on load
* Scaffold tree: Debounced scroll-driven substructure searches so fast scrolling no longer spawns a search for every node scrolled past

## 1.17.8 (2026-05-12)

* Scaffold tree: Migrated drop handler to the new `doDrop(args)` signature in `ui.makeDroppable`
* Chem: Scaffold Tree: Fixed `Array.from(undefined)` crash in `getPeerViewers` when the active view is not a TableView
* Reaction enumerator: Reworked UI — replaced the preview with side-by-side reaction templates and building blocks grids that subset the enumeration via row selection; renamed inputs and the full-config button for clarity
* Reaction enumerator: Right pane is now a tabbed view (Reaction templates, Building blocks, Reagents, Preview); preview tab re-runs lazily on activation
* Reaction enumerator: Added reagents-mode enumeration — when a reagents file/column are provided, each step uses exactly one BB or earlier-round product plus reagents in every remaining slot (produces derivatives of each BB across rounds)
* Reaction enumerator: Added a plain-English mode summary under the config buttons describing how the current selections drive the run
* Reaction enumerator: Added test coverage for reagents-mode enumeration (round 1 / multi-round step-pool invariants, empty-reagents fallback, BB-BB combo isolation)

## 1.17.5 (2026-04-15)

### Features:

* Chem: Stereo agnostic search
* Chem: Pharmacophore feature highlighting

### Bug fixes

* GROK-19962: Fixed crash in checkCurrentView when table is opened in Files preview without a registered TableView
* GROK-19994: MMP: Fixed crash when navigating away from MMP view while tab content is lazily rendered
* GROK-19921: Bio, Chem: Make table view active when adding scatter plot for chemical/sequence space and activity cliffs (in case it was switched during function execution)
* GROK-19862: Radar: Viewer type naming issue
* GROK-14952: Chem: Filter Panel and Hamburger menu molecular filters are not sync
* GROK-18530: Chem: Substructure Filter: incorrect filtering in the cloned view when using the "Similar" option
* Chem: Scaffold Tree: Fix double vertical scrollbar
* Chem: MPO: Add EDA package guard and deduplicate profile loading
* Chem: MPO: Preserve profile name and description on dataset change
* Chem: MPO: Fix profile naming and cancel behavior

## 1.17.4 (2026-03-30)

* [#3708](https://github.com/datagrok-ai/public/issues/3708): Chem: Chemical space: Add tags on columns with embeddings

## 1.17.3 (2026-03-26)

* Improved reaction renderer supporting smiles reactions.

## 1.17.2 (2026-03-23)

RDKit render: Render without fixed scaling if rendered outside the global context (for single values like in custom tooltips).

## 1.17.1 (2026-03-17)

* GROK-19867: Chem: Alignment setting (Scaffold column) is not saved to the data-sync project

## 1.17.0 (2026-03-13)

* Chem: MPO: Fix desirability label positioning and rename templates to profiles
* Chem: MPO: Use name derived from the profile
* Chem: MPO: Use mapped columns in calculation
* [#3497](https://github.com/datagrok-ai/public/pull/3624): Chem: MPO: Implement profile management
* Chem: MMP: Fixed broken pairs grid filtration (by current molecule) after switching from Transfromations tab in demo
* GROK-19543: Chem: MMP: Cliffs: broken color icon location
* Chem: MPO: Use line min/max as defaults when not defined in profile
* GROK-19575: Chem | Descriptors: Broken count for selection
* GROK-19569: Chem | MPO: "Add as column" is stuck in updating
* [#3645](https://github.com/datagrok-ai/public/pull/3645): Chem: MPO: Improvements
* GROK-19576: Chem | Descriptors: Disable modification of selected descriptors when added as a column from the panel
* [#3650](https://github.com/datagrok-ai/public/pull/3650): Chem: MPO: Improvements
* GROK-19623: Chem: MPO: Profiles list: Truncate long descriptions
* Chem: Move Generate Conformers to Calculate
* GROK-19620: Chem: MPO: Saving cloned profile replaces parent profile
* GROK-19621: MPO: 'Edit MPO profile' tab including profile`s name
* GROK-19624: Chem | MPO: Deleting a profile should update browse tree without manual refresh
* Chem: MPO: Ensure the OK button only applies changes (not saving to file)
* Chem: MPO: Highlight matching profiles in the profiles input to give users a quick overview of suitable options
* Chem: MPO: Disable OK button when profile is incompatible or mapping is missing
* Chem: MPO: Avoid showing balloons when previewing results
* GROK-19629: Chem: MPO: Enable actions on Browse view
* GROK-19632: Chem: MPO: error on saving profile
* GROK-19627: Chem: MPO: Multiple times cloned profile case
* GROK-19641: Chem: MPO: Close opened/previewed profiles after their deleting from the browse
* Chem: Add pMPO-generated profile to the Chem folder for demo usage
* Chem: MPO: Display warning in dialog footer
* Chem: MPO: Use ObjectHandler for context menu actions and profile preview in the context panel
* Chem: Explicit smiles setting for sketchers
* Chem: MPO: Align option values in the profiles input
* Chem: Add BitBIRCH clustering for O(N) incremental molecule clustering
* Chem: Bitbirch worker
* Chem: Removing explicit Hs for morgan FPS
* Chem: Smiles length limit
* [#3658](https://github.com/datagrok-ai/public/pull/3658): Chem: Scaffold Tree: Add support for custom labels
* Chem: MPO: Fix description column showing clickable cursor without expand behaviour
* Chem: MPO: Fix issue where clicking a profile in the list does not open it in the context panel
* Chem: MPO: Update match emoji and fix best/worst score interactivity
* Chem: Synthon search
* Chem: Scaffold Tree: Support for custom labels (fixes)
* Chem: MPO: Fix issue with Desirability profile detection
* Chem: MPO: Align ui between pMPO and MPO dialogs
* Chem: MPO: Add fallback color section to settings dialog
* Chem: MPO: Prevent context panel from switching the current object when clicks originate from histogram or scores viewers
* Chem: MPO: Add UI hints
* GROK-18490: Chem: MMP: molecule pairs: improve data rendering in context panel
* Chem: Get mol details
* Chem: MPO: Display profile names in the dialog instead of filenames to avoid confusion (in the app we use profile names)
* Chem: Reaction toolkit: Initial UI
* GROK-19715: Chem: MMP: Fragments: correspond empty trellis plot with molecule pairs list
* Chem: Improved reactions editor and app views
* Chem: MPO: Make aggregation to be part of the profile
* Chem: MPO: Remove hints and add descriptive text in front of the input
* Chem: MPO: Fix issues with routing
* Chem: MPO: Fix recalculation and data sync issues
* Chem: Chemprop: Fix docker container discovery
* [#3681](https://github.com/datagrok-ai/public/pull/3681): Chem: MPO: Improvements to the “Create Profile” view
* Chem: MPO: Improve profile creation flow in the dialog
* Chem: MPO: Fix crash when selecting a profile with no matching columns
* GROK-19837: Chem: MMP Demo: clicking on "/\ hERG_pIC50" in Fragments causes errors
* GROK-19838: Chem: Scaffold Tree: Folder style overlaps
* Chem: MPO: Vertical form layout with descriptions under labels
* Chem: MPO: Prevent default override when saving from the MPO dialog
* Chem: MPO: Add profile import and export
* Chem: MPO: Add breadcrumbs for navigation
* Chem: MPO: Move Save button to ribbon and fix dialog opening in design mode by default
* Chem: MPO: Rename “probabilistic” to “data-driven” and rearrange the viewers
* Chem: MPO: Add checkbox for adding the radar plot to the grid
* Chem: MPO: Add reset button and fix onChanged not firing for gaussian and sigmoid on change

## 1.16.14 (2026-01-05)

* Chem: Rendering: Introduce an option to disable post–R-group analysis highlighting

## 1.16.13 (2025-12-23)

* Butina clustering: Expose Python implementation to top-menu

## 1.16.12 (2025-12-22)

* GROK-18406: Chem | Scaffold Tree: Disallow multiple filters on the same column
* Chem: Detector: Extend valid smiles characters to support additional symbols

## 1.16.11 (2025-12-09)

* Add kekulize option to molecule notation conversion functions

## 1.16.10 (2025-12-08)

* R Group analysis: fixed ability to perform analyziz only for certain positions

## 1.16.9 (2025-11-19)

* Validated deprotection feature with more precise reactions.

## 1.16.8 (2025-11-11)

* Support datasync for biochemical functions

## 1.16.7 (2025-10-22)

* Fix compatibility issues with platform version 1.25.3

## 1.16.6 (2025-10-16)

### Bug fixes

* GROK-19062 Similarity Search viewer: incorrect coloring after sorting a dataset
* GROK-18817: Chem | Scaffold tree: Add explanation tooltip why the tree can't be generated
* GROK-18901: Chem | Scaffold tree: Error when loading the layout
* Chem: Fix Distance functions
* GROK-18916: Chem: Remove hint on viewer close
* GROK-18931: Chem | Properties: Molecule charge is NaN
* GROK-18944: Chemprop: train model: 'Creating interactive model preview' step takes too long
* GROK-18945: Chem | Scaffold Tree: Color column not restored when reapplying saved layout
* GROK-18815, GROK-18816, GROK-18820: Top Menu: UI/UX improvements: Fixed names, captions
* [#3497](https://github.com/datagrok-ai/public/pull/3497): Chem: MPO: Update to work correctly with data sync


### Features

* Chem: Mixtures widget in context panel: show components names where available
* Chem: calculate context panel descriptors via docker
* Ability to run Descriptors, Properties, Toxicity risks, Structural alerts using Add new column dialog

## 1.16.5 (2025-09-15)

* Add layout persistant tag for scaffold alignment

## 1.16.4 (2025-09-05)

* Fixed issues with tableview getting for cloned views in functions.

## 1.16.3 (2025-08-26)

### Features

* Removed on demand from docker containers

## 1.16.2 (2025-08-07)

### Bug fixes

* GROK-18622: Implement logic to check and switch TableView if needed
* [#3445](https://github.com/datagrok-ai/public/issues/3445) Structure filter does not realign when coordinates change but structure does not

## 1.16.1 (2025-08-04)

### Bug fixes

* GROK-18609: MMP: Demo throws errors: Cannot read properties of null (reading 'toList')

## 1.16.0 (2025-07-28)

### Features

* Chem: MMP: Added hints
* Chem: Similarity/Diversity: Show error in case no molecule column was found
* Chem: MMP: Collaborative filtering with parent table
* Chem: MMP: Ability to filter fragments/pairs
* Chem: MMP: View results on PC plot, select activity difference function with scaling, in cell pc plots for activities, show trellis filters in popup, added trellis legend
* Chem: Added description to scatter plot for chemical space and activity cliffs
* Chem: Mixtures renderer, mixtures widgets and semantic type
* [#3426](https://github.com/datagrok-ai/public/pull/3426) CPC: Moved to chem, factored into a separate widget file

### Bug fixes

* Chem: Scaffold Tree: Fix bitset not updating on structure change
* Chem: Scaffold Tree: Check if viewer is detached
* GROK-18069: Chem | Scaffold tree: Error when closing
* GROK-18074 Chem: Context Panel: Content tab: df flickers when filtering by Scaffold Tree and Activity Cliffs simultaneously
* Chem: Deprotect: Column already exists
* GROK-18103: Chem | Scaffold Tree: Skip processing if column is empty
* Chem: Correct most diverse structures logic
* GROK-18066: Chem: Context Panel: Structural Alerts: calculate for whole dataset doesn't work for some data
* GROK-18061: Chem: Prevent from running functions which add viewers not on the active table view
* Chem: Scaffold Tree: Exclude fragment column from Select columns dialog
* GROK-18283: Chem | Scaffold Tree: Error when changing scaffold size
* GROK-18284: Chem | Scaffold Tree: Resetting Description > Title creates duplicates
* GROK-18305: Chem | Scaffold tree: Error when uploading tree for non-molecular dataset
* GROK-18269: Chem: Scaffold Tree: Viewer is empty when opening a project
* Chem: MMP: Do not interpret missing values as zero
* GROK-18290 : Substructure filter: multiple identical filters issue
* [#3346](https://github.com/datagrok-ai/public/issues/3346): Export as SDF: Fixed none values are written as numbers ant not as empty strings
* GROK-18388: Chem | Scaffold tree filter: Error after Close All
* [#3399](https://github.com/datagrok-ai/public/issues/3399): Do not allow to modify Scaffold colors from Scatterplot or other viewers
* Chem: do not allow multiple substructure filters
* GROK-18391: Chem: the project opens with different filtering results before and after clearing the cache
* GROK-18409: Chem: R-group analysis: column prefix in dialog not working
* GROK-18405: Chem | Scaffold tree: Incorrect filtering after layout application
* GROK-18417: Chem: Fix column parent in structural alerts
* GROK-18419: Chem: MMP: activity list is empty
* GROK-18495: Chem: Fix errors on current row unselection
* GROK-18477: Chem: MMP: Update legend when changing activities on tp
* GROK-18491: Chem: MMP: Fragments: "Full screen" trellis plot icon is not working
* GROK-18479: Chem: MMP: Generation: should be no empty values
* GROK-18493: Chem: MMP: Fragments: trellis legend optimization in case of big amount of legend items
* GROK-18494: Chem: MMP: Properties: changing Molecules name cause no effect and triggers "reading 'toList'" errors in console
* GROK-18570: Chem | Structural Alerts: Broken styles
* GROK-18530 Chem: Substructure Filter: incorrect filtering in the cloned view when using the "Similar" option
* GROK-18567: Chem | Chemprop: Cross validation fails for mol1k.csv

## 1.15.1 (2025-04-23)

### Features

* [#3302](https://github.com/datagrok-ai/public/issues/3302): Multiple scaffold trees: Naming/coloring should be saved in projects

### Bug fixes

* GROK-17929: Chem: Similarity sorting checks
* GROK-17958: Chem: some demos overwrite Datagrok tab name

## 1.15.0 (2025-03-29)

### Features

* [#3273](https://github.com/datagrok-ai/public/issues/3273): Scaffold tree: add the ability to apply Scaffold tree colors to Scatterplot or other viewers
* GROK-17731: Chem | Scaffold Tree: Remember previous hide/show state
* GROK-17867: Chem: Demo improvements
* Add caching to Gasteiger charges
* [#3302](https://github.com/datagrok-ai/public/issues/3302): Multiple scaffold trees: add the ability to apply Scaffold tree colors to Scatterplot
* MMP: Open scatter plot context menu when right click on cliffs sp, separate generations scatter plots for each activity, cliffs tab scatter plot improvements, trellis filters improvements, use column descriptions instead of custom tooltips for non molecule columns, added help url
* Chem: Cluster MCL Implementation
* [#3314](https://github.com/datagrok-ai/public/pull/3314): Chem: Scaffold Tree: Demo improvements
* Chem: Parallel beautification

### Bug fixes

* GROK-17761: Chem: Demo: Chem: Activity cliffs: Viewer initialization function failed
* GROK-17791: Chem: Scaffold tree: Coloring the scaffold tree element triggers the grid scroll to the end
* Chem: Similarity/DIversity search: some progress bar related fixes
* GROK-17892: Chem | Scaffold tree: Errors on closing viewer
* GROK-17867: Chem: Demo: Added setTimeout for showing help

## 1.14.0 (2025-02-24)

### Features

* New improved MMP implementation
* Integration with Reinvent
* Copy molecule as image
* Context actions for chemical columns

### Bug fixes

* GROK-17483: Chem: Structure filter is not reset on reseting the filtering in some cases
* GROK-17461: Scaffold tree: Error when uploading the saved tree
* GROK-17511: Chem | ChemProp: Error when training a model
* GROK-17555: Filter Panel: Reset filter causes errors in some cases
* GROK-17471: Chem: Incorrect filter UI when adding filter via 'Use as filter' or 'Add filter'
* GROK-17597: Chem: Wrong filtering when opening the project

## 1.13.1 (2025-01-06)

### Bug fixes

* Fix RDKit cell renderer isolating MonomerHoverData.getSubstruct
* Fix RDKit cell renderer isolating ISubstructProvider.getSubstruct

## 1.13.0 (2024-10-25)

### Bug Fixes

* [#3055](https://github.com/datagrok-ai/public/issues/3055): Page is freezing after applying layout with multiple
  scaffold trees with large number of structures

## 1.12.4 (2024-10-18)

### Bug fixes

* Fix rdkit cell renderer for viewers tooltip, enable test
* Invalidate fingerprints column only in case column data (not metadata) has been changed

## 1.12.3 (2024-10-01)

### Bug fixes

* GROK-16730: Demo: Activity cliffs: Cannot read properties of null

## 1.12.2 (2024-09-25)

### Bug fixes

* Fix rdkit cell renderer error in scatter plot tooltip, add test

## 1.12.1 (2024-09-25)

### Bug fixes

* Expect GPU to be unavailable on some browsers
* Fix MMPa crushing without errors

## 1.11.6 (2024-09-20)

### Bug Fixes

* [#3055](https://github.com/datagrok-ai/public/issues/3055): Page is freezing after applying layout with multiple
  scaffold trees with large number of structures

## 1.11.5 (2024-08-27)

### Bug fixes

* Add tableView to runChemSpace for reduce dimensionality ui options
* [#3004](https://github.com/datagrok-ai/public/issues/3004): Chem: Scaffold tree is based on wrong table by default in
  some cases
* GROK-15938: Chem | Scaffold Tree: Array length differ
* Chem: MMPA: Fix demo

## 1.11.4 (2024-08-23)

### Features

* Chem: MMPA: filter fragments by significance
* Chem: MMPA: project and layout saving
* Chem: MMPA: demo loaded not calculated
* Chem: Removed r-groups settings download from init

### Bug Fixes

* Chem: MMPA: Fixed property panel on cliff line click
* [#2983](https://github.com/datagrok-ai/public/issues/2983): Chem: fix tooltip styles for structure columns
* GROK-16493: MMPA: Fixed Summary viewer doesn't show the values

## 1.11.3 (2024-08-12)

### Features

* Chem: MMPA: Work with multucomponent entities

## 1.11.2 (2024-08-12)

### Features

* Chem: Chemprop: Remove redundant calls in the core and move it fully to package
* Chem: MMPA: Add meaningfull dataset for demo
* Chem: MMPA: Restrictions and CPU fallback
* Chem: MMPA: Compensate filtering by cutoff
* Chem: MMPA: Add reflection for selection
* Chem: MMPA: Add pairs filter, expose filters object so it can be used in any tab, minor improvements
* Chem: Update demos in order to be shown in the Browse

### Bug Fixes

* Chem: MMPA: Fixes for cliffs tab and lines improvements
* [#2964](https://github.com/datagrok-ai/public/issues/2964): Scaffold Tree: 'Selection to filter' reset immediately for
  some data if scaffold tree is present

## 1.11.1 (2024-07-25)

### Bug Fixes

* Removed unnecessary code

## 1.11.0 (2024-07-24)

*Dependency: datagarok-api >= 1.20.0*

### Bug Fixes

* GROK-16328: Chem: R-Group analysis: Broken styles
* GROK-16329: Chem: R-Group analysis: Error when none R groups were found
* GROK-16338: Chem: Chemical space: Reseting marker type to circles
* Chem: Scaffold Tree: Fix representation when the viewers are stacked
* Chem: Chemprop: Fix applying the model to different datasets
* Chem: Fix substructure filter style
* Chem: Do not render molecules with invisible width

## 1.10.0 (2024-06-12)

### Features

* Webgpu acceleration for chem space, activity cliffs, molecule generations
* GROK-15299: Chem: Ability to save results of Chem functions to projects with and without data sync
* GROK-15454: Chem: Similarity/Diversity search viewers: Row source property
* GROK-15508: Chem: MMP: Generated - show what cases are actually generated
* GROK-15501: Chem: MMPA: UI/UX improvements
* GROK-15502: Chem: MMPA: Ability to swith activity on/off
* GROK-15507: Chem: MMP: Control of fragments size
* Chem: Integrate with Chemprop
* Chem: Chemprop upgrade from version "1.6.7" to version "2.0" (reduced memory usage, improved performance)

### Bug Fixes

* GROK-15285: Chem: MMP: Error on re-run
* GROK-15321: Chem: MMP: Error while running Close All after running MMP

## 1.9.2 (2024-04-12)

### Features

* R Group analysis improvements: ability to select 'Only match at R groups', saving R group user settings, small UI
  improvements

## 1.9.1 (2024-04-12)

### Features

* R Group analysis: in browser R group analysis using RDKit JS, highlight of r groups, ability to define enumerated R
  groups in core
* Similarity/Diversity search: ability to recalculate results on filtering

## 1.9.0 (2024-03-21)

### Features

* [#2495](https://github.com/datagrok-ai/public/issues/2495): Chem: MMP - initial vesrion
* [#2540](https://github.com/datagrok-ai/public/issues/2540): Chem | Substructure search: Add an option to toggle the
  realignment of molecules by molecular filter
* [#2525](https://github.com/datagrok-ai/public/issues/2525): Chem | Scaffold Tree: Remove the intermediate dialog
  before the color picker
* [#2504](https://github.com/datagrok-ai/public/issues/2504): Filter panel: The question mark hovering doesn’t show the
  filtering of the scaffold tree viewer

### Bug Fixes

* GROK-15147: Chem | Scaffold Tree: Highlighting disappears in some cases
* GROK-15240: Chem | RGroups: Handling malformed data
* GROK-15213: Chem | RGroups: Undo latest analysis doesn't work when run on multiple dataframes
* GROK-15187: Chem | RGroups: Fixed memory consumption
* GROK-15171: Chem: Substructure filter doesn't apply correctly after removing rows and disable

## 1.8.13 (2024-02-29)

### Bug Fixes

* Warm up substructure search web workers only once when first filter is created on column

## 1.8.12 (2024-02-15)

### Bug Fixes

* Fixed multiple substructure search calls with the same structure (in case multiple layouts are applied to cloned views
  at once)

## 1.8.11 (2024-02-08)

### Bug Fixes

* [#2642](https://github.com/datagrok-ai/public/issues/2642): Filtering done using viewers is unexpectedly reset on
  applying other filters in some cases
* [#2639](https://github.com/datagrok-ai/public/issues/2639): Structure in filter is unexpectedly changed when there are
  two views in some cases
* [#2646](https://github.com/datagrok-ai/public/issues/2646): Structure highlighting is not removed when filters panel
  is closed

## 1.8.10 (2024-01-19)

### Features

* Molecule labels on scatter plot zoom in Chemical Space and Activity Cliffs

### Bug Fixes

* [#2626](https://github.com/datagrok-ai/public/issues/2626): Unexpected warning on adding child node to the scaffold
  tree if parent contains 'H' atom
* [#2628](https://github.com/datagrok-ai/public/issues/2628): Structure filter is not applied in some cases when there
  are two views opened
* [#2629](https://github.com/datagrok-ai/public/issues/2629): Structure filter filters out all molecules in the presence
  of filter-by-zoom scatterplot and multiple views
* GROK-14542: Filter Panel | Structure filter: incorrect behavior on deleting a filter in a cloned view

## 1.8.9 (2024-01-09)

### Features

* Update MinimalLib to 1.2.15 featuring support for RGroups (support for R-group decomposition)
* Chem: Add ability to specify fingerprint types and cluster embeddings for chemical space

### Bug Fixes

* GROK-14383: Chem | Scaffold Tree: The selection doesn't change on folder node click
* GROK-14411: Demo: Exception while closing Scaffold Tree
* GROK-14409: Chem | R-Group Analysis: Error with prefix
* GROK-14375: Chem | Activity cliffs: Implement scatter plot lines renderer
* [#2612](https://github.com/datagrok-ai/public/issues/2612): Chem | Scaffold tree: Structure part containing H atom is
  replaced with '?' in some cases

## 1.8.8 (2023-12-07)

### Bug Fixes

* Chem: Fix demos for scaffold tree, chemical databases and similarity/diversity search.

## 1.8.7 (2023-12-05)

### Features

* [#2550](https://github.com/datagrok-ai/public/issues/2550): Chem: Scaffold Tree: Counter is not updated when applying
  the filters
* [#2525](https://github.com/datagrok-ai/public/issues/2525): Chem: Scaffold Tree: Highlighting improvements:
  * More specific fragments should be drawn on top
  * New icons for different coloring states

### Bug Fixes

* [#2536](https://github.com/datagrok-ai/public/issues/2536): Chem: Substructure Filter: Not terminated on fiter reset
* [#2553](https://github.com/datagrok-ai/public/issues/2553): Chem | Scaffold tree: Realignment issue when the scaffold
  color is set
* GROK-14305: Chem | Scaffold Tree: Error when deleting rows in grid

## 1.8.6 (2023-11-17)

### Bug Fixes

* [#2525](https://github.com/datagrok-ai/public/issues/2525): Chem: Scaffold tree: Highlighting improvements:
  * Add a tooltip for the checkbox
  * Add a tooltip for NOT
* [#2526](https://github.com/datagrok-ai/public/issues/2526): Chem: Molecules realignment improvements:
  * Scaffold tree shouldn't reorient molecules in the grid
* [#2533](https://github.com/datagrok-ai/public/issues/2533): RdKitService stays unresponsive as terminateFlag is not
  reset after substructure search has completed

## 1.8.5 (2023-11-14)

### Features

* [#2459](https://github.com/datagrok-ai/public/issues/2459): Chem: Scaffold Tree: Improvements:
  * Add resetting the scaffold tree viewer (on the reset filter click)
  * Add a tooltip for removing orphans

### Bug Fixes

* [#2511](https://github.com/datagrok-ai/public/issues/2511): Scaffold tree: highlighting is not updated when structure
  is edited
* [#2512](https://github.com/datagrok-ai/public/issues/2512): Filtering options (contains, exact, included in etc.) are
  not synchronised for structure filters on different tabs in some cases

## 1.8.4 (2023-11-07)

### Features

* [#2459](https://github.com/datagrok-ai/public/issues/2459): Chem: Scaffold Tree: Improvements:
  * The “trash can” icon is red and stands out from the general style; make a blue one from the same collection
  * Add Expand All/Collapse All icons for easy navigation
  * Save the state of the tree (what was selected, what was expanded) when saving
* [#2420](https://github.com/datagrok-ai/public/issues/2420): Chem | Scaffold tree coloring:
  * use the trash iscon instead of '-'

### Bug Fixes

* [#2474](https://github.com/datagrok-ai/public/issues/2474): Chem | Scaffold tree: Colour is not inherited by child
  node properly
* [#2476](https://github.com/datagrok-ai/public/issues/2476): Chem | Scaffold tree: Rows are shown as filtered for a
  while after applying structure even if checkbox is not selected
* [#2448](https://github.com/datagrok-ai/public/issues/2448): Some structures are displayed incorrectly when highlighted
  if they are in SMILES format
* [#2473](https://github.com/datagrok-ai/public/issues/2473): Structure filter: 'current value > use as filter' breaks
  filtering in some cases

## 1.8.3 (2023-10-26)

### Features

* [#2421](https://github.com/datagrok-ai/public/issues/2421): Chem | Molecular Search options:
  * Add some hint to explain the new options to users
  * Fix the ability to collapse the controls after switching from Similarity to Contains
  * Make the icon more noticeable
* [#2420](https://github.com/datagrok-ai/public/issues/2420): Chem | Scaffold Tree coloring:
  * The custom coloring should override the default one
  * Coloring for child structures - it's very inconvenient that you can't click the icon again to return to the custom
    color
  * Add the icon for removing the specific node of the scaffold tree

### Bug Fixes

* [#2399](https://github.com/datagrok-ai/public/issues/2399): Scaffold tree viewer: coloring can't be saved to the
  Layout
* [#2450](https://github.com/datagrok-ai/public/issues/2450): Scaffold tree imported from file is not restored on
  applying saved layout
* [#2448](https://github.com/datagrok-ai/public/issues/2448): Some structures are displayed incorrectly when highlighted
  if they are in SMILES format

## 1.8.2 (2023-10-22)

### Bug Fixes

* GROK-14149: Add progress indicator to sequence/chem space and remove jumping points

## 1.8.1 (2023-10-16)

### Bug Fixes

* [#2399](https://github.com/datagrok-ai/public/issues/2399): Scaffold Tree: Coloring can't be saved to the Layout
* [#2400](https://github.com/datagrok-ai/public/issues/2400): Molecular search: Three dots infinite loading

## 1.8.0 (2023-08-09)

### Features

* GROK-13817: Highlighting multiple substructures with different colors inside one molecule structure
* [#2355](https://github.com/datagrok-ai/public/issues/2355): Integrate Scaffold Tree with color-coded fragments
* GROK-13993: Fast Chemical space on large datasets (> 20k)  using sparse matrix
* GROK-13994: Implemented filter by superstructure, exact structure, similarity score
* GROK-13966: Exposed fingerprints options for similarity, diversity search

## 1.7.2 (2023-09-05)

### Bug Fixes

* Chem: Scaffold Tree: Checkbox shouldn't be set, when group is expanded
* Chem: Scaffold Tree: Fix the behaviour of allowGenerate property

## 1.7.1 (2023-08-31)

### Features

* GROK-13571: Chem | Ability to terminate substructure search if substructure has been changed

### Bug Fixes

* GROK-13327: Chem | Substructure Search: two identical panels open on the Filter Panel
* GROK-13791: Chem | Chemical space (using t-SNE) fails on smiles dataset
* [#2135](https://github.com/datagrok-ai/public/issues/2135):
  * The structure rendering is too small.
* [#2322](https://github.com/datagrok-ai/public/issues/2322): Properties panel is unexpectedly reset on changing viewer
  properties if there is a scaffold tree filter in filters panel
* GROK-13848: Chem: Substructure search results flickering

## 1.7.0 (2023-08-09)

### Features

* GROK-13172: Chem | implement substructure search using preliminary filtration by pattern fingerprints

## 1.6.22 (2023-08-07)

### Bug Fixes

* GROK-13713: Chem | Incorrect molecule rendering

## 1.6.21 (2023-08-02)

*Dependency: datagarok-api >= 1.16.0*

### Features

* Calculate drug likeness, toxicity and alerts for whole table from widgets
* Color Coding for toxicity
* Scaffold Tree improvements:
  * [#2154](https://github.com/datagrok-ai/public/issues/2154): Scaffold Tree: harmonization.

### Bug Fixes

* GROK-13586: _chemFindSimilar fails with 'Cannot read properties of null (reading 'rows')'
* [#2135](https://github.com/datagrok-ai/public/issues/2135):
  * The counts and controls are partially hidden
  * When change the drawing of a scaffold, then selection (checkbox) of all other scaffolds gets reset (to unchecked)
* [#2139](https://github.com/datagrok-ai/public/issues/2139): Scaffold tree stops working after adding an invalid
  structure

## 1.6.20 (2023-07-21)

This release focuses on improvements and bug fixes.

*Dependency: datagarok-api >= 1.14.0*

### Features

* Set default values in all dialogs where appropriate.
* Unified the layout of the search panel results.
* Work with DBs APIs: Added links to sources, similarity search viewer style, reproduce for ChemicalSpace, Enamine.
* **Elemental Analysis**: the resulting column name now includes the specific column for which calculations were
  conducted.
* **R-Groups Analysis**: a new option to choose between searching for MCS exact atoms or MCS exact bonds.
* **Similarity Search**: a new **Follow Current Row** setting prevent recalculation when the current row is changed.
* For proper handling of properties and rendering, we now check for smarts and molecular fragments separately.
* Ability to copy data from the **Descriptors** and **Properties** tabs on the **Context Pane**.
* Moved **Descriptors** and **Fingerprints** from the  **Context Pane**  to the **Top Menu** ( **Chem** > **Calculate
  **).
* Added **Substructure Search** to the **Top Menu** ( **Chem** > **Search**)
* Modified the tooltip for dialog and drag-n-drop to prevent it from overlapping with the data.
* Implemented RDKit rendering for Chembl, ChemblAPI, PubChem, and DrugBank databases if OCL is used currently.
* UI polishing: harmonized input field names, repositioned elements in dialogs, added tooltips, and organized the top
  menu items into groups: Calculate, ADME/Tox, Search, Analyze, and Transform.
* Scaffold Tree improvements:
  * [#1730](https://github.com/datagrok-ai/public/issues/1730): Implemented Scaffold Tree integration to the **Filters
    Panel**.
  * [#1998](https://github.com/datagrok-ai/public/issues/1998): Now a scaffold tree shows a confirmation dialog before
    dropping all trees on the **Clear** icon click.
  * Added the **Allow Generate** property for the viewer in order to control autogeneration.

### Bug Fixes

* GROK-13105: Substructure search doesn't work after similarity search.
* GROK-13118: Activity cliffs selects non-numeric column as activity.
* GROK-13123: Substructure Search: error when the Filter Panel is opened with substructure filtering.
* [#1492](https://github.com/datagrok-ai/public/issues/1492): Elemental analysis: malformed data handling.
* GROK-11898: Orientation for smiles in Structural Alerts.
* GROK-12115: Hamburger menu closing while switching a sketcher.
* GROK-12905: Sketcher is not opening from the **Filter Panel**.
* GROK-12929: Scripts don't work if called from the package.
* GROK-12933: Drug likeness: set score precision.
* GROK-12961: Elemental Analysis: `Unsupported operation: NaN.round()` error on some data.
* GROK-12962: Similarity Search: doesn't work on malformed data.
