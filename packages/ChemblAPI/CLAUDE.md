# ChemblAPI

## Purpose

Thin HTTP client to EBI's public **ChEMBL** and **UniChem** REST APIs. Exposes similarity/substructure
search widgets and info panels for molecules, plus a UniChem `inchikey → external IDs` lookup.

**Not to be confused with the sibling `Chembl` package** ([../Chembl/](../Chembl/)), which talks to a
ChEMBL PostgreSQL mirror via SQL queries and a data connection. This package has no `queries/`,
no `connections/`, and no database dependency — all lookups go out to `https://www.ebi.ac.uk`.

## Architecture

- [src/package.ts](src/package.ts) — `PackageFunctions` class. Only decorated, registered Datagrok
  functions live here: the two molecule info panels, `GetCompoundsIds`, and `Chembl Get by Id`.
- [src/utils.ts](src/utils.ts) — internal helpers, not registered with the platform:
  `getData` (XML fetch + parse), the two thin wrappers `chemblSubstructureSearch` /
  `chemblSimilaritySearch`, `getSmiles`, and `chemblSearchWidget` (the grid UI the panels render).
  The `SEARCH_TYPE` and `ELEMENTS` enums and the `BASE_URL`/`WIDTH`/`HEIGHT` constants also live here.
- [src/package-test.ts](src/package-test.ts) — test entry point, wires `runTests`.
- [src/tests/searches.ts](src/tests/searches.ts) — `getData` similarity/substructure tests against a
  fixed SMILES (`O=C1CN=C(c2ccccc2N1)C3CCCCC3`), both expect exactly **20** rows.
- [src/tests/chem-panels.ts](src/tests/chem-panels.ts) — UniChem `GetCompoundsIds` test against a fixed
  InChI key.
- [swaggers/chembl.json](swaggers/chembl.json) — OpenAPI description of the ChEMBL REST endpoints.
  The platform reads the swagger and exposes each operation as a data query callable via
  `grok.data.query('ChemblAPI:<OperationName>', {...})`. Not imported by the TS code, but
  load-bearing at runtime — `getById` depends on it (see exposed-functions table).

## Exposed functions

Declared in [src/package.ts](src/package.ts); typed wrappers regenerated into `src/package.g.ts` by
`grok api`.

| Function | Role | What it does |
|---|---|---|
| `Databases \| ChEMBL \| Substructure Search API` | `panel` | Molecule info panel — renders `chemblSearchWidget` in substructure mode. |
| `Databases \| ChEMBL \| Similarity Search API` | `panel` | Molecule info panel — renders `chemblSearchWidget` in similarity mode. |
| `GetCompoundsIds` | plain | UniChem lookup: InChI key → list of `{src_id, src_compound_id}`. |
| `Chembl Get by Id` | plain | Calls the swagger-generated query `ChemblAPI:MoleculeJson` (= ChEMBL `/molecule.json`) for a given ChEMBL ID. |

`chemblSearchWidget` is **not** registered — it's an internal UI helper in [src/utils.ts](src/utils.ts)
that the two panels call directly. `getData(searchType, smiles, score?)` is the HTTP/XML workhorse
behind the widget; also internal and exported only for tests.

## Glossary

- **ChEMBL ID** — `CHEMBLnnnnn`. `Chembl Get by Id` prepends `CHEMBL` if the caller passes just digits.
- **`SEARCH_TYPE`** — `'substructure'` or `'similarity'`; used as a URL path segment against
  `https://www.ebi.ac.uk/chembl/api/data/{searchType}/{smiles}[/score]`.
- **`ELEMENTS`** — XML tag names ChEMBL returns (`molecule_chembl_id`, `canonical_smiles`,
  `molecule_properties`, `similarity`).
- **Similarity score** — ChEMBL returns `0–100`; `getData` divides by 100 and stores as float.
- **Panel vs widget** — the two panels exist so they can be surfaced separately in the context
  panel. Both delegate to the same `chemblSearchWidget`.

## Conventions and gotchas

1. **External HTTP must go through `grok.dapi.fetchProxy`** (repo rule; CORS would otherwise break).
   Already correct throughout this package — do not swap in raw `fetch`.
2. **ChEMBL responses are XML, parsed with `DOMParser`.** Fields come out of
   `xmlDoc.getElementsByTagName(...)` — there is no JSON path.
3. **Result count is hard-capped at 20** in two places in [src/utils.ts](src/utils.ts): `getData`'s
   `rowCount = Math.min(molecules.length, 20)` and the widget's own `molCount = Math.min(table.rowCount, 20)`.
   Changing the limit requires touching both.
4. **Runtime dependency on the `Chem` package.** `getSmiles` calls `Chem:convertMolNotation`. If
   `Chem` is not deployed on the target server, the widget short-circuits to
   `"Molecule string is malformed"`.
5. **The `MoleculeJson` query comes from the swagger, not from a `queries/` folder.** `Chembl Get by Id`
   calls `grok.data.query('${_package.name}:MoleculeJson', {molecule_chembl_id__exact: id})`. There is
   no `queries/` directory in this package — the platform auto-registers a query per swagger operation,
   so `MoleculeJson` ↔ the `/molecule.json` path in [swaggers/chembl.json](swaggers/chembl.json).
   When grepping for these query names you won't find SQL or JS for them; look in the swagger.
   Note also that no other package in `public/packages/` calls `Chembl Get by Id` today (it's only
   reachable through `grok.functions.call`).

## Build and test

Standard Datagrok package layout — see [../CLAUDE.md](../CLAUDE.md) for `npm run build`, `grok publish`,
`grok test`, linking, and webpack externals. Nothing here deviates from the defaults.
