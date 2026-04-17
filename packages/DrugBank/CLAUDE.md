# DrugBank

Context panels and a drug-name converter backed by the DrugBank Open Structures dataset
(~11k approved/investigational drugs). When a cell's semantic type is `Molecule`, the
panels run substructure or similarity search against the local dataset and render up
to 20 matching molecule tiles; clicking a tile opens the drug on go.drugbank.com.

The dataset (`files/drugbank-open-structures.d42`) is loaded once at init from a binary
Datagrok dataframe and kept in memory — there is no remote API call. Chemistry operations
delegate to the Chem package (`grok.chem.searchSubstructure`, `grok.chem.findSimilar`).

## Architecture

- `src/package.ts` — `PackageFunctions` class with `@grok.decorators`. `initDrugBank`
  reads the `.d42` file into module-scoped `dbdf`, `synonymsCol`, `moleculeCol`. The
  three exposed functions are two panels (`#panel` + `meta.role: widgets`) and the
  `db:<name>` converter (`meta.role: converter`, `meta.inputRegexp: (db\:.+)`).
- `src/widgets.ts` — `searchWidget` builds the panel (header + grid of up to 20 tiles,
  each drawn with `grok.chem.drawMolecule`; similarity panel adds a score line).
  `drugNameMoleculeConvert` does a case-insensitive `includes()` scan over `SYNONYMS`
  and returns the first matching `molecule` value. Also defines `SEARCH_TYPE` and
  `COLUMN_NAMES` (the d42 column names — `DRUGBANK_ID`, `COMMON_NAME`, `molecule`,
  `SYNONYMS`, plus `score` on the similarity result).
- `src/searches.ts` — `searchSubstructure` sets `dbdf.filter` from the Chem bitset and
  returns `dbdf` itself (not a copy); `findSimilar` returns the Chem result dataframe
  with `DRUGBANK_ID` and `COMMON_NAME` joined in via the `index` column.
- `src/tests/drugbank-tests.ts` — smoke tests only (no assertions); constants in `const.ts`.

## Dataset columns (in `drugbank-open-structures.d42`)

| Column        | Type   | Notes                                            |
|---------------|--------|--------------------------------------------------|
| `molecule`    | string | SMILES / molfile — the chem search target        |
| `DRUGBANK_ID` | string | e.g. `DB00945` — used in the drugbank.com URL    |
| `COMMON_NAME` | string | Shown in tooltips                                |
| `SYNONYMS`    | string | Delimited list; the converter does a substring match (not exact) |

## Conventions

- The panel functions keep `dbdf` at module scope and mutate its `filter` in place
  (substructure search). Don't refactor this to return copies without checking that
  downstream grid-view code still works.
- The `db:<name>` converter is invoked via the platform's input converter pipeline
  (matched by `meta.inputRegexp`); it is not called directly from TS.
- `grok publish` must be run against a server that already has the Chem package
  installed — the chem calls will fail at runtime otherwise.
