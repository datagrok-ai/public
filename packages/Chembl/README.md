# ChEMBL

ChEMBL is a [package](https://datagrok.ai/help/develop/#packages) for the [Datagrok](https://datagrok.ai)
platform that integrates the [ChEMBL](https://www.ebi.ac.uk/chembl/) bioactivity database and the
companion [UniChem](https://www.ebi.ac.uk/unichem/) database. The package ships as a thin
TypeScript layer on top of a large collection of SQL queries against a Postgres ChEMBL instance with
the [RDKit cartridge](https://rdkit.org/docs/Cartridge.html) installed.

The bundled connections point at the public Datagrok demo ChEMBL database, so the package works
out of the box — no credentials required.

## Features

* **Substructure and similarity search** against the local ChEMBL DB via the RDKit cartridge,
  exposed both as a context panel (`Databases | ChEMBL | Substructure Search`,
  `Databases | ChEMBL | Similarity Search`) and as standalone functions.
* **Chemical identifier converters** — round-trip between ChEMBL ID, SMILES, InChI, InChI Key,
  molregno, and compound name. The `chemblIdToSmilesTs` converter is wired to the
  `CHEMBL[0-9]+` regexp, so any recognised ChEMBL ID anywhere in the platform can be expanded
  into a molecule.
* **CHEMBL ID handler and database explorer** — clicking a ChEMBL ID opens a
  [DBExplorer](https://github.com/datagrok-ai/public/tree/master/libraries/db-explorer) view
  that starts from `molecule_dictionary` and lets the user walk through related tables
  (structures, activities, assays, targets, synonyms, mechanisms, …). Declarative
  [enrichments](https://github.com/datagrok-ai/public/tree/master/packages/Chembl/enrichments)
  add extra columns on the fly.
* **Info panels** (`MolregnoInfo`, `ChemblInfo`) that surface compound metadata in the context
  panel based on the `CHEMBL_ID` / `molregno` semantic types.
* **HitTriage data sources** — `Chembl Compounds`, `Chembl targets by organism`, and the
  `Chembl molregno` column function, for seeding HitTriage campaigns from ChEMBL.
* **Ad-hoc browse/search queries** (`Browse | …`, `Search | …`, `Misc | …`) covering FRAC
  classification, bioactivity for bacterial targets, PK data, activity details for a target,
  and more. Each has a matching `.layout` file where applicable.
* **ChEMBL browser UI** backed by the internal `_cb…` queries.
* **Demo app** — `Database Queries` (registered at `Cheminformatics | Database Queries`)
  demonstrates parameterised queries against ChEMBL.

## Directory layout

* [connections](https://github.com/datagrok-ai/public/tree/master/packages/Chembl/connections)
  — `Chembl` (Postgres, used for most queries), `ChemblSql` (PostgresDart, used by the
  dataframe-in/dataframe-out converters), and `Unichem`.
* [queries](https://github.com/datagrok-ai/public/tree/master/packages/Chembl/queries) —
  `cartridge.sql` (RDKit substructure/similarity), `converters.sql` (ID conversion),
  `queries.sql` (user-facing browse/search + info-panel widgets), `browser.sql` (internal
  browser queries), `suggestions.sql` (autocomplete sources), plus standalone parameterised
  queries (`activityDetailsForTarget.sql`, `bioactivityForBacterialTargets.sql`, `pkData.sql`)
  with sibling `.layout` files.
* [enrichments](https://github.com/datagrok-ai/public/tree/master/packages/Chembl/enrichments)
  — declarative DBExplorer enrichments that pull extra columns when a ChEMBL ID or molregno is
  materialised.
* [src](https://github.com/datagrok-ai/public/tree/master/packages/Chembl/src) —
  `package.ts` (function registrations), `handlers.ts` (CHEMBL ID handler),
  `explorer-config.ts` (DBExplorer entry points, joins, and column subsets), `demo.ts`
  (demo app).

See also:

* [Grok API](https://datagrok.ai/help/develop/packages/js-api)
* [Packages](https://datagrok.ai/help/develop/#packages)
* [Data Connection](https://datagrok.ai/help/access/#data-connection)
* [Data Query](https://datagrok.ai/help/access/#data-query)
* [Info Panels](https://datagrok.ai/help/explore/data-augmentation/info-panels)
* [Semantic Types](https://datagrok.ai/help/govern/catalog/semantic-types)
