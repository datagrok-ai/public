# PubChemApi

Thin REST client for the PubChem PUG API (`https://pubchem.ncbi.nlm.nih.gov/rest/pug`).
Exposes three info panels on `Molecule` columns (substructure / similarity / identity
search against PubChem), two semantic-type converters (`pubchem:<CID>` → SMILES,
InChIKey → SMILES), and a `GetIupacName` function.

## Architecture

- `src/package.ts` — function registrations via `@grok.decorators` on static methods of
  `PackageFunctions`. Panels use `meta.role: widgets`; converters use
  `meta.role: converter` + `meta.inputRegexp`. `grok api` regenerates `package.g.ts` /
  `package-api.ts` from these.
- `src/pubchem.ts` — REST wrappers. Searches go through the unified `search()` function
  which hits PubChem's async `fast*` endpoints (`fastsimilarity_2d`, `fastidentity`,
  `fastsubstructure`) and handles three possible response shapes: `Waiting.ListKey`
  (async — then polled via `getListById`, up to 30× at 500 ms intervals), a bare
  PC_Compound array (fastidentity), or `{PC_Compounds: [...]}` (fastsimilarity_2d).
  Identity returns full records with `.props`; similarity/substructure are flattened
  via `flattenPcCompound` to `{CID, CanonicalSMILES, ConnectivitySMILES}`. Thin wrappers
  `similaritySearch` / `identitySearch` / `substructureSearch` forward to `search()`.
  `getBy` is a synchronous identifier lookup (e.g. InChIKey → CID). `getCompoundInfo`
  hits the richer `pug_view/data/compound/{id}` endpoint consumed by the info panel.
- `src/widget.ts` — all UI. `getSearchWidget` renders the grid of molecule cards
  (CID tooltip + open-in-PubChem click + "open as table" icon). `buildInfoPanel`
  renders a flat info panel: GHS hazard icons if the record has a `Primary Hazards`
  section, then a single `ui.tableFromMap` with Name / CID /
  Formula / MW / CAS / IUPAC / SMILES / InChIKey pulled from the `pug_view` tree via
  two helpers — `sectionAt(record, path)` walks the `Section[]` tree by `TOCHeading`
  (e.g. `['Names and Identifiers', 'Computed Descriptors', 'SMILES']`), and
  `readString(record, path)` extracts the first stringifiable `Value` under it
  (`StringWithMarkup[0].String` or `Number[0]`, optionally suffixed with `Unit`).
  `extractHazardIcons` reuses `sectionAt` to reach the icon `Markup[]` array.
  A single footer link opens the full record on pubchem.ncbi.nlm.nih.gov. The 2D
  structure is intentionally omitted — it's already visible in the grid. Similarity
  flow re-ranks PubChem hits locally with `grok.chem.findSimilar` (limit 20, cutoff
  0.75).
- `src/utils.ts` — shared types (`pubChemSearchType`, `pubChemIdType`, `anyObject`,
  `paramsType`), `urlParamsFromObject`, and `getSmiles` (a normalizer that calls
  `Chem:convertMolNotation`).
- `src/constants.ts` — PubChem base URLs (`pubChemBaseURL`, `pubChemRest`, `pubChemPug`)
  and the `COLUMN_NAMES` enum used by `widget.ts` and `flattenPcCompound`.
- `src/tests/const.ts` — canned test molecule strings (SMILES, MOL2000, MOL3000,
  SMARTS, EMPTY) used by the test suite.

## Glossary

- **CID** — PubChem Compound ID (integer).
- **InChIKey** — 27-char hashed InChI, regex `[A-Z]{14}-[A-Z]{10}-N`. Matched by the
  `inchiKeysToSmiles` converter.
- **`pubchem:<CID>`** — user-facing identifier literal. The `pubChemToSmiles` converter
  regex (`^\s*pubchem\s*:\s*[0-9]+\s*$`, case-insensitive) routes these to SMILES.
- **ListKey** — handle returned by async PubChem searches; client polls until `Waiting`
  disappears from the response.
- **`PubChemApi` connection** (referenced in `@grok.decorators.func` `connection:`
  metadata) — auto-registered from `swaggers/pubchem.json` (Swagger 2.0, host
  `pubchem.ncbi.nlm.nih.gov/rest/pug`). There is no local `connections/` folder —
  update the swagger file to change the connection or its endpoints.

## Conventions

- **Always use `grok.dapi.fetchProxy`** for PubChem URLs — never raw `fetch` (CORS will
  fail).
- **Normalize mol inputs via `getSmiles`** (`Chem:convertMolNotation`) before hitting
  PubChem. This makes `@datagrok/chem` a runtime dependency even though `package.json`
  lists it under `devDependencies` (it's resolved by name via `grok.functions.call`).
- **Always `encodeURIComponent` SMILES before embedding in a URL path segment.** SMILES
  contain `#` (triple bond), `[`, `]`, `@`, `=`, `(`, `)`, `/`, `\`, `+`, etc. — all of
  which the proxy will reject if left raw. Partial escaping (e.g. only `#`) breaks on
  real-world inputs like `[H][C@]1(O[C@H]...)`.
- **Expect either `CanonicalSMILES` or `ConnectivitySMILES`** in PubChem property
  responses — always fall back between the two (see `pubChemToSmiles` and `widget.ts`
  `moleculesCol` selection).
- Query strings are built with `urlParamsFromObject`, not `URLSearchParams`.
