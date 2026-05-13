# CLAUDE.md

Guidance for Claude Code working on the Chembl package.

## Purpose

Datagrok plugin that exposes the [ChEMBL](https://www.ebi.ac.uk/chembl/) bioactivity database (and a
companion [UniChem](https://www.ebi.ac.uk/unichem/) connection) to platform users. The bulk of the
package is **SQL queries** against a Postgres ChEMBL instance with the **RDKit cartridge** installed,
plus a small TypeScript layer that wires those queries into context panels, info panels, the
`CHEMBL_ID` semantic type, the database explorer, the HitTriage data-source picker, and the demo app.

Auth / credentials: built-in. The three connections in `connections/` point at the public Datagrok
demo Postgres at `db.datagrok.ai:54325/54326` with hard-coded credentials — no secrets to set.

## Architecture

- **`src/package.ts`** — entry point. `PackageFunctions` holds the registered functions (decorator
  style, regenerated into `package.g.ts` by `grok api`): the `init` autostart that registers the
  CHEMBL ID handler, two info-panel widgets (`Substructure Search (Internal)`, `Similarity Search
  (Internal)`) backed by the shared `chemblSearchWidgetLocalDb`, two `hitTriageDataSource` functions
  (`Chembl Compounds`, `Chembl targets by organism`), the `chemblMolregno` HitTriage column
  function, the `chemblIdToSmilesTs` converter (regexp `(CHEMBL[0-9]+)` → `Molecule`), and the
  `Database Queries` demo. Also exports two free functions used by the widgets:
  `chemblSubstructureSearch` and `chemblSimilaritySearch`, both of which round-trip through
  `Chem:getRdKitModule` to canonicalise input before calling the cartridge queries
  `patternSubstructureSearch` / `patternSimilaritySearch`.
- **`src/handlers.ts`** — `registerChemblIdHandler()` boots a `DBExplorer` from `explorerConfig` and
  attaches a custom renderer that draws SMILES cells as molecules. Called once from `init`.
- **`src/explorer-config.ts`** — declarative `DBExplorerConfig` for `@datagrok-libraries/db-explorer`.
  Defines entry points (`CHEMBL_ID`, `molregno`), join shortcuts (e.g. `molecule_dictionary` →
  `compound_structures` brings `canonical_smiles`/`standard_inchi`), explicit cross-schema FKs from
  `rdk.mols`/`rdk.fps` to `public.molecule_dictionary`, the column subset shown for
  `molecule_dictionary`, and the per-table "preview" column shown in `headerNames`.
- **`src/demo.ts`** — `_demoDatabasesChembl()` builds a two-tab view (input form + raw SQL) that
  reruns `FracClassificationWithSubstructure` whenever an input changes.
- **`src/tests/`** — `cartridge.ts` (smoke-tests the three search queries) and `converters.ts`
  (golden values for every converter query). All marked `stressTest: true`.
- **`queries/*.sql`** — where most of the product lives. See "Queries" below.
- **`connections/*.json`** — `Chembl` (Postgres provider, used by everything that needs `Choices` /
  `batchMode` / parameterised inputs), `ChemblSql` (PostgresDart provider, used **only** by the
  converters in `converters.sql`), `Unichem` (separate DB, used by one unit-test query).
- **`enrichments/*.json`** — declarative DBExplorer enrichments (extra columns pulled in on the fly
  when a CHEMBL ID / molregno is materialised). One file per logical enrichment, each specifying a
  key table/column and a join plan against `Chembl:Chembl`. Not queries — don't confuse with
  `queries/`. Edit when you want an extra field to appear alongside existing CHEMBL drill-downs.

## Queries

Each `.sql` file is a flat list of `--name`-delimited queries. Adding a new query = appending a
block; no TS wiring is needed unless TS code calls it directly. Group conventions:

| File | What it holds |
|---|---|
| `cartridge.sql` | RDKit-cartridge queries: `patternSimilaritySearch`, `patternSimilaritySearchWithThreshold` (uses `--meta.batchMode: true` to set `rdkit.tanimoto_threshold` before the SELECT), `patternSubstructureSearch`, plus the simple `ChemblNumberOfStructures`, `ChemblMolregNoBySmiles`, `StructuresByOrganism`. |
| `converters.sql` | All ID-conversion queries (ChEMBL ↔ SMILES ↔ InChI ↔ InChIKey ↔ molregno ↔ name). Two are registered as type converters via `--meta.role: converter` + `--meta.inputRegexp`. |
| `queries.sql` | Ad-hoc browse/search queries shown to end users via the friendly-name menu (`Browse | …`, `Search | …`, `Misc | …`). Includes the `MolregnoInfo` / `ChemblInfo` info-panel widgets (registered by `--tags: panel, widget`). |
| `browser.sql` | Internal `_cb…`-prefixed queries that back the ChEMBL browser UI. |
| `suggestions.sql` | Autocomplete sources — referenced from other queries as `{suggestions: Chembl:organisms}`, `{choices: Query("…")}`, etc. |
| `cartridge.sql` extras | — |
| `activityDetailsForTarget.sql`, `bioactivityForBacterialTargets.sql`, `pkData.sql` | Stand-alone parameterised queries; each has a sibling `.layout` file that the platform applies to results. |

A few SQL conventions worth knowing:

- The bind-parameter syntax is `@name`. For RDKit functions that need a Postgres positional cast
  (e.g. `mol_from_smiles($1::cstring)`), `@pattern` and `$1` refer to the same first parameter —
  see `patternSimilaritySearch`.
- Multi-statement queries (e.g. `set_config(...)` followed by `SELECT`) require
  `--meta.batchMode: true` plus a `--batch` separator between statements.
- Friendly names use the prefix as a menu category: `Browse | …`, `Search | …`, `Converters | …`,
  `Misc | …`, `Suggestions | …`. Names beginning with `_` (e.g. `_cbAllChemblStructures`) are
  internal queries hidden from the database tree.
- Layout files (`<queryName>.layout`) are picked up by name match — keep them in sync when renaming.

## Glossary — domain concepts → code

| Concept | Code / SQL | Meaning | Key relationships |
|---|---|---|---|
| **CHEMBL ID** | semType `CHEMBL_ID`; column `molecule_dictionary.chembl_id`; matches regex `CHEMBL\d+` | Public, stable compound identifier (e.g. `CHEMBL1185`). What end-users paste. | Maps 1:1 to `molregno` via `molecule_dictionary`. |
| **molregno** | semType `molregno`; integer PK on `molecule_dictionary.molregno` | Internal compound registration number. The join key for almost every cross-table query. | Foreign key from `compound_structures`, `compound_records`, `activities`, `rdk.mols`, `rdk.fps`, `molecule_synonyms`, … |
| **canonical_smiles / standard_inchi / standard_inchi_key** | columns on `compound_structures` | Three canonical structure encodings keyed by `molregno`. | Source of truth for the converter queries. |
| **max_phase** | column on `molecule_dictionary` (numeric 0–4) | Highest clinical phase a compound has reached. `max_phase = 4` ≡ approved drug (ChEMBL convention, not in the column name). | Filter on this for approved-compound queries. |
| **RDKit cartridge** | schema `rdk` (`rdk.mols`, `rdk.fps`), GUC `rdkit.tanimoto_threshold`, operators `@>` (substructure), `%` and `<%>` (Tanimoto), functions `morganbv_fp`, `mol_from_smiles`, `get_mfp2_neighbors` | Postgres extension that adds chemical types/operators. Substructure & similarity searches go through it. | `rdk.mols` / `rdk.fps` join back to `molecule_dictionary` on `molregno` (declared in `explicitReferences` in `explorer-config.ts`). |
| **target / assay / activity** | tables `target_dictionary`, `assays`, `activities`, `target_components`, `component_sequences` | Bioactivity backbone. A target has assays; assays produce activities for compounds. | Used by `StructuresByOrganism`, `compound activity details for all targets containing @protein`, `bioactivityForBacterialTargets`, etc. |
| **drug_mechanism / research_companies / molecule_synonyms** | tables of the same names | Used together (`molregno` → `res_stem_id` → `country/company`) in commercial-origin search queries. | See `QueryBySubstructure`, `MolregnoInfo`. |
| **FRAC classification** | `frac_classification`, `molecule_frac_classification` | Fungicide Resistance Action Committee taxonomy. Four-level hierarchy (level1..level4_description). | Powers the `FracClassification*` queries. |
| **Chembl vs ChemblSql connection** | `connections/chembl.json` (`dataSource: Postgres`) vs `connections/chemblsql.json` (`dataSource: PostgresDart`) | Same DB, two providers. The Postgres (Java) provider supports `Choices`, `Suggestions`, `batchMode`. The PostgresDart provider is faster for the high-volume dataframe-in/dataframe-out converter queries. | New `--input` queries with `Choices` or `batchMode` should use `Chembl`. New row-by-row converters that take/return a dataframe should use `ChemblSql`. |
| **CHEMBL ID handler** | `DBExplorer` from `@datagrok-libraries/db-explorer`, configured by `explorerConfig` | Lets the user click a CHEMBL ID anywhere and explore the schema starting from that row. | Initialised once in `init`; the molecule renderer is added on top. |

## Conventions specific to this package

- New queries: pick the file by topic (`cartridge.sql` for RDKit, `converters.sql` for ID
  conversion, `queries.sql` for user-facing browse/search, `browser.sql` for internal `_cb…`
  queries). Add `--description`, `--friendlyName` (with `Browse | …` / `Search | …` / `Misc | …`
  prefix), and a default value for every `--input`. Match the connection convention above.
- Multi-statement SQL (e.g. `SET LOCAL`, `set_config`) must use `--meta.batchMode: true` and a
  `--batch` line between statements; otherwise the provider rejects it.
- Substructure / similarity TS code must canonicalise the input via `Chem:getRdKitModule` before
  hitting the cartridge — see `chemblSubstructureSearch` / `chemblSimilaritySearch` in `package.ts`.
  Mirror that pattern instead of passing user input straight to the SQL.
- The `Chem` package is a **devDependency only** — it's available at runtime on every server that
  loads Chembl, but never `import` from it; call `grok.functions.call('Chem:…')` instead.
- Functions in `package.ts` use the **TypeScript decorator** form (`@grok.decorators.func`,
  `@grok.decorators.panel`, `@grok.decorators.param`). `grok api` regenerates `package.g.ts` and
  `package-api.ts` from these — don't edit those files. The repo-level rule about JSDoc `//name:`
  comments still applies for `package-test.ts` and for scripts/queries.
- When adding a new entry point or join to the database explorer, edit
  `src/explorer-config.ts` (not `handlers.ts`). Cross-schema joins (e.g. `rdk.*` →
  `public.molecule_dictionary`) must be declared in `explicitReferences` because the cartridge
  tables don't have real Postgres FKs.
