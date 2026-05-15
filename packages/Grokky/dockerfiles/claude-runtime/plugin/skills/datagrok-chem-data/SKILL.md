---
name: datagrok-chem-data
description: Bioactivity and small-molecule data-side catalog. Covers ChEMBL search/lookup, MolTrack compound and batch registration, HitTriage / HitDesign / PeptiHit workflows, and Curves dose-response analysis (IC50/EC50 fitting). Open this skill when the user asks about querying ChEMBL, registering or looking up compounds in MolTrack, running a hit-triage campaign, or fitting/inspecting dose-response curves.
---

# datagrok-chem-data

Function catalog for the **data and workflow** half of Datagrok cheminformatics:
Chembl, MolTrack, HitTriage, Curves.

All functions are callable from a `datagrok-exec` block via
`grok.functions.call('<Package>:<funcName>', {param1, param2, ...})`. For the
data-side packages, results are usually `DG.DataFrame` objects produced by SQL
queries or by Docker services (MolTrack).

For the structure-and-modeling side (Chem, Admetica, Reinvent4, Docking,
Boltz1), see the sibling skill **`datagrok-chem-toolkit`**.

### Param-type conventions used below

Types match Datagrok function declarations: `dataframe`, `column`, `string`,
`int`, `float`, `bool`, `list<T>`, `object`. A column with a required semType
is written `column<semType>` (e.g. `column<Molecule>`). Enum values are listed
inline as `'a' | 'b' | 'c'`. `?` marks optional, `= x` shows the default.

### Namespace casing

The package is referenced as **`Chembl:`** in `grok.functions.call(...)` and
`grok.data.query(...)`. It's branded as "ChEMBL" in UI text, but the function
namespace is namespace-exact lowercase-after-first-letter: `Chembl:foo`,
never `ChEMBL:foo` or `chembl:foo`.

### Always emit an exec block when the intent is clear

If the user asks for a concrete data-side action ("register this compound",
"fit IC50 from this data", "what is the SMILES of aspirin", "list ChEMBL
compounds for organism X") — emit a `datagrok-exec` block calling the relevant
function below. Don't answer chemistry/data questions from memory when a
catalog function exists. If a required input (SMILES, payload, mapping JSON)
is missing, emit the block with a clearly-marked placeholder or sensible
default and add a one-line note — don't bail out into a conversational reply.

### Tag legend (used in function tables below)

Same tags as `datagrok-chem-toolkit`'s "Calling conventions" section, plus a
few data-side specific ones. Common values that appear in the **Tags** column:

| Tag                  | Meaning                                                                                                                                  |
|----------------------|------------------------------------------------------------------------------------------------------------------------------------------|
| `vectorFunc`         | Returns a `DG.DataFrame`. Usable inside a `grokky.addCalculatedColumn` formula and the Add-New-Column dialog. Without it, call via `grok.functions.call`. |
| `transform`          | Tagged as a data-sync transform — same behavior as the un-tagged sibling but tracked in project history.                                 |
| `topMenu`            | Has a Top-menu path.                                                                                                                     |
| `app` / `appTreeBrowser` | App entry point / its tree-browser sidebar hook.                                                                                     |
| `panel`              | Context-panel widget — `DG.Widget` shown in the right side panel.                                                                        |
| `demo`               | Demo function — opens a tutorial layout.                                                                                                 |
| `fileViewer` / `fileHandler` | Registered for a specific file extension.                                                                                        |
| `cellRenderer`       | Custom cell renderer.                                                                                                                    |
| `converter`          | Type converter (`meta.role: 'converter'` with an `inputRegexp`).                                                                         |
| `curveConverter`     | Curves package — registers a curve-format converter (`3dx`, `compact-dr`, `pzfx`).                                                       |
| `hitTriageDataSource` / `pepTriageDataSource` | HitTriage / PepTriage data source registered for the campaign-ingest picker.                                |

---

## Cross-package routing — "User wants X → call Y"

| User intent | Function | Package |
|---|---|---|
| ChEMBL substructure search | `Chembl:chemblSubstructureSearchPanel` (panel) or query `Chembl:patternSubstructureSearch` | Chembl |
| ChEMBL similarity search | `Chembl:chemblSimilaritySearchPanel` (panel) or query `Chembl:patternSimilaritySearch` | Chembl |
| ChEMBL ID → SMILES | `Chembl:chemblIdToSmilesTs` | Chembl |
| Get ChEMBL molregno for SMILES | `Chembl:chemblMolregno` | Chembl |
| Pull ChEMBL compounds by organism | `Chembl:getChemblCompoundsByOrganism` | Chembl |
| Open MolTrack | `MolTrack:molTrackApp` | MolTrack |
| Register one compound / batch | `MolTrack:registerMolTrackProperties` (props), `MolTrack:registerBulk` (bulk CSV) | MolTrack |
| Look up compound by corporate ID | `MolTrack:getCompoundByCorporateId` | MolTrack |
| Look up batch by corporate ID | `MolTrack:getBatchByCorporateId` | MolTrack |
| Advanced compound search with structured filter | `MolTrack:advancedSearch` | MolTrack |
| Generic entity search (compounds / batches / assays / runs / results) | `MolTrack:search` | MolTrack |
| List compound / batch / assay schema | `MolTrack:fetchCompoundProperties`, `MolTrack:fetchBatchProperties`, `MolTrack:fetchSchema`, `MolTrack:fetchDirectSchema` | MolTrack |
| Register assay / assay results | `MolTrack:registerAssays` | MolTrack |
| Retrieve all entities of a scope | `MolTrack:retrieveEntity` | MolTrack |
| Start a Hit Triage campaign | `HitTriage:hitTriageApp` | HitTriage |
| Start a Hit Design campaign | `HitTriage:hitDesignApp` | HitTriage |
| Start a PeptiHit campaign | `HitTriage:peptiHitApp` | HitTriage |
| Get demo molecule dataset (100 / 5000 / N) | `HitTriage:demoFileIngest`, `HitTriage:demoFileIngest1`, `HitTriage:demoFileIngest2` | HitTriage |
| Convert raw plate readouts into fit-curve cells | `Curves:dataToCurves` or `Curves:dataToCurvesTopMenu` | Curves |
| Compute IC50 / Hill / AUC / Top / Bottom from dose-response data | `Curves:dataToCurves` (fit) → `Curves:addStatisticsColumn` (extract) | Curves |
| Pull a stat from existing fit cells | `Curves:addStatisticsColumn` | Curves |
| Aggregate fit stats across series | `Curves:addAggrStatisticsColumn` | Curves |
| Open a `.pzfx` GraphPad Prism file | `Curves:previewPzfx` | Curves |
| Compute Minimum Significant Ratio | script `Curves:CalculateMSR` (see Curves Scripts) | Curves |

---

## Chembl package — local ChEMBL search and lookup (`Chembl:...`)

Thin TS layer over SQL queries against a local ChEMBL PostgreSQL connection.
The heavy lifting is in `Chembl/queries/` — these TS functions just expose
useful entry points.

### Programmatic functions

| Function | Tags | What it does |
|---|---|---|
| `chemblIdToSmilesTs(id: string = 'CHEMBL1185')` | `converter` | CHEMBL ID (e.g. `'CHEMBL1185'`) → canonical SMILES. Registered as a `converter` for the `(CHEMBL[0-9]+)` regex — string fields auto-resolve. **Only goes ID → SMILES**; for InChI key / other-ID → ChEMBL ID, use `Chem:mapIdentifiersTransform` from the toolkit skill. |
| `chemblMolregno(table: dataframe, molecules: column<Molecule>)` | — | Appends a `CHEMBL molregno` column to the table by looking up each SMILES in the ChEMBL DB. |
| `getChemblCompoundsByOrganism(maxNumberOfMolecules: int = 1000, organism: string = 'Shigella')` | `hitTriageDataSource` | Pulls compounds active against targets of the given organism. Returns DataFrame. |
| `getChemblCompounds(maxNumberOfMolecules: int = 1000)` | `hitTriageDataSource` | Pulls a generic ChEMBL compound sample. |
| `namesToSmiles(names: list<string>)` (SQL query) | — | The underlying SQL query that powers `Chem:namesToSmiles`. Returns a DataFrame with a `canonical_smiles` column. Invoke via `grok.functions.call('Chembl:namesToSmiles', {names: ['aspirin','ibuprofen',...]})` or via `grok.data.query`. |

### Panels (context-panel widgets)

All entries below have the `panel` tag.

| Panel | What |
|---|---|
| `chemblSubstructureSearchPanel(mol)` | "Databases \| ChEMBL \| Substructure Search (Internal)". |
| `chemblSimilaritySearchPanel(mol)` | "Databases \| ChEMBL \| Similarity Search (Internal)". |
| `chemblSearchWidgetLocalDb(mol: string, substructure: bool = false)` | The underlying widget — call directly to embed in custom UI. |

### Demos

`demoDatabasesChembl()` — opens the "Database Queries" demo. Top menu: Cheminformatics → Database Queries.

### Skipped

- `init` — registers the CHEMBL ID handler at startup.

---

## MolTrack package — small-molecule registration system (`MolTrack:...`)

Docker-backed compound and batch registration with dynamic properties, assay
schema, and a structured search DSL. Backed by the `moltrack` Docker container
and a `moltrack` PostgreSQL connection.

**Entity scope enum** (used by `registerBulk`, `search`, `retrieveEntity`): `'compounds' | 'batches' | 'assays' | 'assay_runs' | 'assay_results'`.

### App and views

| Function | Tags | What it does |
|---|---|---|
| `molTrackApp(path?: string)` | `app` | Opens the MolTrack app. `path` (URL) routes to sub-views: `'Compound' \| 'Batch' \| 'Register' \| 'Search' \| 'Saved Searches' \| 'Bulk' \| 'Assay' \| 'Schema'`. URL params `corporate_compound_id` / `corporate_batch_id` jump to a specific entity. browsePath `Chem`. |
| `molTrackAppTreeBrowser(appNode: object, browseView: object)` | `appTreeBrowser` | Wires the tree-browser sidebar (auto). |
| `initDB()` | — | Initializes all MolTrack schemas + assay data — call once after first deploy. |

### Lookups and metadata

| Function | Tags | What it does |
|---|---|---|
| `getCompoundByCorporateId(corporateValue: string)` | — | Compound JSON by corporate compound ID. |
| `getBatchByCorporateId(corporateValue: string)` | — | Batch JSON by corporate batch ID. |
| `fetchCompoundProperties()` / `fetchBatchProperties()` | — | JSON describing all properties for the `compound` / `batch` scope. Cached monthly. |
| `fetchSchema()` | — | JSON of all dynamic fields. |
| `fetchDirectSchema()` | — | JSON of all static (direct-table-column) fields. |

### Registration

| Function | Tags | What it does |
|---|---|---|
| `registerMolTrackProperties(jsonPayload: string)` | — | Registers compound/batch property definitions. `jsonPayload` is a JSON string defining fields `name`, `value_type`, `entity_type`, `nullable`, etc. |
| `registerAssays(assayPayload: string)` | — | Register assay schema (JSON string with `name`, properties, result types). |
| `registerBulk(csvFile: file, scope: 'compounds' \| 'batches' \| 'assays' \| 'assay_runs' \| 'assay_results', mapping: string, errorHandling: 'reject_all' \| 'reject_row')` | — | Bulk-register rows from CSV. `mapping` is a JSON string `{ csvColumn: moltrackProperty }`. `errorHandling`: `'reject_all'` aborts on any invalid row; `'reject_row'` skips. |

**For "register this compound" intent:** even if the payload isn't fully spelled out, emit an exec block that calls `MolTrack:registerMolTrackProperties` with a placeholder `jsonPayload` string (e.g. `{"name":"...","value_type":"...","entity_type":"compound","nullable":false}`) or `MolTrack:registerBulk` with a placeholder CSV/mapping. Do NOT bail into a conversational reply — surface the function call with TODOs.

### Search

| Function | Tags | What it does |
|---|---|---|
| `search(query: string, entityEndpoint: 'compounds' \| 'batches' \| 'assays' \| 'assay_runs' \| 'assay_results')` | — | Generic search — `query` is a JSON string with `level`, `output`, `filter`, `output_format`. |
| `advancedSearch(outputFields: list<string>, filter: object)` | — | **Most useful programmatic entry.** Returns a DataFrame of compounds. `outputFields` are fully-qualified (`compounds.canonical_smiles`, `compounds.details.<dynamic_prop>`). `filter` is a structured boolean tree. **Operators:** string `'=' \| '!=' \| 'IN' \| 'STARTS WITH' \| 'ENDS WITH' \| 'LIKE' \| 'CONTAINS'`; numeric `'<' \| '>' \| '<=' \| '>=' \| 'RANGE'`; datetime `'BEFORE' \| 'AFTER' \| 'ON'`; molecular (only for `compounds.structure`) `'IS SIMILAR'` (needs `threshold`), `'IS SUBSTRUCTURE OF'`, `'HAS SUBSTRUCTURE'`. Compose with `{operator: 'AND' \| 'OR', conditions: [...]}`. |
| `retrieveEntity(scope: 'compounds' \| 'batches' \| 'assays' \| 'assay_runs' \| 'assay_results')` | — | Returns *all* entities of a scope as a flattened DataFrame. |
| `getMoltrackPropPanelById(id: semantic_value<Grok ID>)` | `panel` | "Databases \| MolTrack" context panel for a `Grok ID` cell. |

### Skipped

- `init` — bootstraps MolTrack queries, semantic types, and schema once per session.

### Queries shipped in `queries/`

- `setup.sql` — schema-init queries; called via `init` (`checkDBInitialized`, etc.).
- `search.sql` — backing query used by the search UI.

---

## HitTriage package — small-mol hit triage workflows (`HitTriage:...`)

Multi-step campaign apps for triaging hits, designing new molecules, and the
peptide-flavored variants.

### Apps (each returns a `DG.ViewBase`)

All entries below have the `app` tag.

| Function | What it does |
|---|---|
| `hitTriageApp()` | **Hit Triage** — end-to-end small-mol triage pipeline. browsePath `Chem`. |
| `hitDesignApp()` | **Hit Design** — design + score new analogs. browsePath `Chem`. |
| `peptiHitApp()` | **PeptiHit** — peptide variant. Internally inits Bio. browsePath `Peptides`. |
| `pepTriageApp()` | **PepTriage** — peptide triage. Internally inits Bio. browsePath `Peptides`. |

### Tree-browser hooks (auto-registered)

`hitTriageAppTreeBrowser`, `hitDesignAppTreeBrowser`, `peptiHitAppTreeBrowser`,
`pepTriageAppTreeBrowser` — populate the campaigns tree in the sidebar.

### Demo data sources

| Function | Tags | What it does |
|---|---|---|
| `demoFileIngest()` | `hitTriageDataSource` | 100-molecule demo dataset. |
| `demoFileIngest1()` | `hitTriageDataSource` | 5000-molecule demo dataset. |
| `demoFileIngest2(numberOfMolecules: int)` | `hitTriageDataSource` | Variable-size demo dataset. |
| `demoPeptideSequences(peptideCount: int)` | `pepTriageDataSource` | Demo HELM peptide sequences for PeptiHit. |
| `demoFileSubmit(df: dataframe, molecules: string)` | — | Demo submit handler (logs row count). |

### Misc

| Function | Tags | What it does |
|---|---|---|
| `registerMoleculesToViD()` | — | Pushes all campaign molecules into MolTrack-managed ViD. |
| `hitDesignVidPanel(vid: semantic_value<HIT_DESIGN_VID>)` | `panel` | "Hit Design V-iD" context panel — shows campaigns for a V-iD. |
| `gasteigerCellRenderer()` | `cellRenderer` | PNG cell renderer for Gasteiger charge images. |
| `htPackageSettingEditor(propList: object)` | — | Hit Triage package settings editor (UI). |

### Skipped (UI plumbing)

Tree-browser callbacks above are not invoked directly by user code.

### Queries shipped in `queries/`

- `locks.sql` — concurrency / lock tracking queries.
- `vid.sql` — campaign V-iD lookups: `getCampaignsByVid`, `getMoleculeByVid`
  (used by the `hitDesignVidPanel` panel).

---

## Curves package — dose-response curves (`Curves:...`)

Stores curve data in a single cell (one `DG.Column` of fit chart JSON) and
renders it as a per-cell curve. Includes converters from XML 3DX, compact-DR,
and GraphPad Prism `.pzfx` formats.

**Fit-statistic property enum** (used by `addStatisticsColumn` and `addAggrStatisticsColumn`'s `propName`): keys from the curve `FitStatistics` type — `'rSquared' | 'auc' | 'interceptX' | 'interceptY' | 'slope' | 'top' | 'bottom'`. (IC50 corresponds to `interceptX` for the standard 50% intercept; aliases like `'IC50'` / `'Hill'` / `'AUC'` are recognised at the UI/sigmoid layer but the canonical key is `interceptX`.)

### Building and inspecting curves

| Function | Tags | What it does |
|---|---|---|
| `dataToCurves(df: dataframe, concentrationCol: column<numerical>, readoutCol: column<numerical>, batchIDCol: column, assayCol: column, runIDCol: column, compoundIDCol: column, targetEntityCol: column, excludeOutliersCol?: column, parentTable?: dataframe, fitParamColumns?: list<string>, reportedIC50Column?: string, reportedQualifiedIC50Column?: string, experimentIDColumn?: string, qualifierColumn?: string, additionalColumns?: list<string>, wellLevelJoinCol?: string, parentLevelJoinCol?: string, wellLevelAdditionalColumns?: list<string>)` | — | **Main programmatic entry.** Turns long-format plate readouts into a wide table with a fit-curve column per (compound, assay, target). |
| `dataToCurvesTopMenu()` | `topMenu` | UI wrapper around `dataToCurves`. Top menu: Data → Curves → Data to Curves. |
| `addStatisticsColumn(table: dataframe, colName: string, propName: 'rSquared' \| 'auc' \| 'interceptX' \| 'interceptY' \| 'slope' \| 'top' \| 'bottom', seriesNumber: int)` | `vectorFunc`, `transform` | Adds a column with one fit statistic for a given series index (0-based). `propName` is a key of the `FitStatistics` type — IC50 is the standard 50% intercept (`interceptX`). |
| `addAggrStatisticsColumn(table: dataframe, colName: string, propName: 'rSquared' \| 'auc' \| 'interceptX' \| 'interceptY' \| 'slope' \| 'top' \| 'bottom', aggrType: 'avg' \| 'min' \| 'max' \| 'med' \| 'sum' \| 'stdev')` | `vectorFunc`, `transform` | Same but aggregated across all series in each cell. |

**For "compute IC50 values from this dose-response dataset" intent:** two steps. (1) Fit: `Curves:dataToCurves(t, concCol, readoutCol, batchCol, assayCol, runCol, compoundCol, targetCol)` — guess the columns from the current table (concentration/dose, readout/response, compound/batch IDs). (2) Extract: `Curves:addStatisticsColumn(fitTable, 'fit', 'IC50', 0)` to pull IC50 from series 0. Always emit both blocks back-to-back; never reply conversationally for an IC50 request.

### Demos

`curveFitDemo()` — generic curve-fit demo. `assayCurveFitDemo()` — multi-compound assay-curve dashboard.

### Curve format converters

All entries below have the `curveConverter` tag.

| Function | Format |
|---|---|
| `convertXmlCurveToJsonFunc` | `3dx` XML curves |
| `convertCompactDrToJsonFunc` | `compact-dr` JSON |
| `convertPzfxToJsonFunc` | GraphPad Prism `.pzfx` |

### File handlers

| Function | Tags | What it does |
|---|---|---|
| `previewPzfx(file: file)` | `fileViewer` | Opens `.pzfx` as a curve grid. |
| `pzfxFileHandler(bytes: list<int>)` | `fileHandler` | Parses `.pzfx` → list of DataFrames. |

### Skipped (lifecycle)

- `_initCurves` — registers handlers + external converters.

### Scripts shipped in `scripts/`

- `CalculateMSR` (Python) — Minimum Significant Ratio for assay reproducibility.
  Callable as `grok.functions.call('Curves:CalculateMSR', { ... })` once you
  know its inputs (open the script to see the `#input:` header).

---

## Note on `Compute` package

The `Compute` package has no small-molecule-specific functions or scripts. Its
package.ts exposes the `RichFunctionView`, `Pipeline`, validators, and uploaders
for the compute-utils framework; its `scripts/` are physics / dynamics demos
(Lorenz attractor, Lotka-Volterra, Newton cooling, antibody Fab-arm exchange,
maximal-filter-volume in R, bioreactor utilities). If a small-molecule task
seems to need Compute, you almost certainly want a different package — usually
one of the ones cataloged here or in `datagrok-chem-toolkit`.
