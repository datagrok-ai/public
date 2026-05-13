# SureChembl

Patent search over a locally-deployed SureCHEMBL Postgres database. The database lives in a package-managed Docker container (RDKit cartridge enabled); the package exposes substructure and similarity search as both data functions and info-panel widgets that appear when a `Molecule`-typed value is selected anywhere in the platform.

## Architecture

```
[ Context panel widget ]                 [ Direct call: grok.functions.call ]
            \                                       /
             \                                     /
        sureChembl{Substructure,Similarity}SearchWidget
                     |
                     v
              patentSearch()  --(debounced inputs)-->  runSearch()
                                                          |
                                                          v
                            sureChembl{Substructure,Similarity}Search
                                                          |
                                          (RDKit normalizes the query)
                                                          |
                                                          v
                                  searchPatentBy{Substructure,Similarity}
                                            (queries/queries.sql)
                                                          |
                                                          v
                                   SureChembl connection -> Surechembl container
```

- The two `*Widget` panels are pure UI; they delegate to the two data functions through `grok.functions.call(`${_package.name}:${fnName}`, ...)` where `fnName` is read off the static method (`PackageFunctions.sureChemblXSearch.name`) inside `patentSearch`. Widgets only pass `SearchType`, never a function-name string — IDE renames of the static methods stay in sync. Going through `grok.functions.call` (not a direct call) is required so the per-function `meta.cache` annotations apply.
- The two data functions are thin wrappers around the shared `patentDataSearch(molecule, limit, searchType, threshold?)` helper. The helper normalizes the query molecule with RDKit (cached via `getRdkitModule()`), then dispatches through `grok.data.query` to the corresponding SQL query.

## The query-normalization contract (load-bearing)

This is the only non-obvious bit of the package and the part most likely to break under naive edits. The branch lives at the single `searchType === SearchType.substructure ? mol.get_smarts() : mol.get_smiles()` line inside `patentDataSearch` — touching it changes both flows at once:

| `searchType`                | RDKit method called | What is sent to SQL    | Why                                                          |
|-----------------------------|---------------------|------------------------|--------------------------------------------------------------|
| `SearchType.substructure`   | `mol.get_smarts()`  | SMARTS string          | Postgres query casts to `qmol` for `m@>` substructure match  |
| `SearchType.similarity`     | `mol.get_smiles()`  | Canonical SMILES       | `get_mfp2_neighbors(@pattern)` expects a molecule, not a query|

Do not unify these two outputs. SMARTS sent to the similarity query, or SMILES sent to the substructure query, will silently return wrong/empty results — not an error.

`patentDataSearch` enforces a `MAX_SMILES_LENGTH` (5000) limit on non-molblock input before invoking RDKit; `grok.chem.isMolBlock(molecule)` is the gate.

## Stale-result guard

`patentSearch` keeps a per-search-type token `${molecule}|${limit}|${threshold}` in `currentSearches`. `runSearch` writes the new token before launching the query and the resolve/reject handlers compare against it before mutating the DOM. If you add a new input that affects results, include it in this token or fast-typing users will see results from a superseded query.

## Wiring to the container and DB

- Connection: [connections/surechembl.json](connections/surechembl.json) — `name: "SureChembl"`, server resolves via the `${Surechembl<DockerContainer>}` placeholder, so spinning up the container automatically points the connection at it. The two SQL queries reference this connection by name (`--connection: SureChembl`); renaming the connection requires updating both queries.
- Container: [dockerfiles/Dockerfile](dockerfiles/Dockerfile) — `datagrok/demo_db_surechembl:2024.1`, on-demand, 30-min idle timeout, 1 CPU / 2 GB. Tests must `await ensureContainerRunning('surechembl', CONTAINER_TIMEOUT)` before any query call; without it the first call races the cold-start.
- Cache: both data functions and both SQL queries declare `meta.cache: all` with daily invalidation `0 0 * * *`. Keep these aligned — invalidating only the TS side leaves stale rows on the server.

## SQL result shape

The `fields` column is a comma-joined string-aggregation of the integer codes in `schembl_document_chemistry.field`:

| code | label       |
|------|-------------|
| 1    | DESCRIPTION |
| 2    | CLAIMS      |
| 3    | ABSTRACT    |
| 4    | TITLE       |
| 5    | IMAGES      |
| 6    | ATTACHMENTS |

`updateSearchPanel` reads the column verbatim — it does not parse the joined string.

The similarity query also emits a `similarity` column; substructure does not. `updateSearchPanel` branches on `table.col('similarity')` to decide whether to render the score, so the two queries' column shapes must stay in sync apart from this one column.

## Dead code in `src/`

- [src/download-patents.ts](src/download-patents.ts) and [src/surechembl-api.ts](src/surechembl-api.ts) are not referenced from `package.ts`. The single former call site in `updateSearchPanel` (around [src/package.ts:122-140](src/package.ts#L122-L140)) is commented out. Treat them as not-yet-shipped scaffolding for downloading patent PDFs from `surechembl.org`, not as live infrastructure. Do not assume changes to either file affect runtime behavior.

## Conventions specific to this package

- Self-reference functions as `${_package.name}:funcName`, never as a hardcoded `Surechembl:` literal. The package's `friendlyName` is `SureChEMBL` but the registered name is `Surechembl` — the placeholder avoids drift if either is ever renamed.
- Tests live under [src/tests/](src/tests/) and are aggregated through `package-test.ts`; they hit the real container, not a mock.
