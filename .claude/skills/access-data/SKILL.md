---
name: access-data
version: 0.2.1
description: |
  Pull rows from an external SQL database into a Datagrok `DataFrame` —
  declare a connection, store credentials, author a parameterized query,
  and call it from TypeScript. For plugin authors whose code needs data
  from outside the browser. Produces `connections/<name>.json` +
  `queries/<name>.sql` (registering a `<Package>:<Query>` function) and
  the TS that resolves the query into a `DataFrame`.
  Use when asked to "pull rows from a database into a dataframe", "wire
  a parameterized SQL query into my plugin", or "store DB credentials
  for a Datagrok connection". For REST endpoints see
  `fetch-external-rest`; for bundled CSVs see `load-bundled-csv`; for
  user/server files see `read-server-file`.
triggers:
  - pull rows from postgres into dataframe
  - wire a parameterized sql query
  - register a sql query as a package function
  - store database credentials securely
  - declare a datagrok data connection
allowed-tools:
  - Read
  - Edit
  - Write
  - Bash
harness-authored: true
---

# access-data

## When to use

Your package needs rows from an external SQL database materialized as a
`DataFrame` — with `grok publish` wiring the connection, credentials
kept out of git, and one JS call to run the query. The spine is four
steps: declare the connection, transfer credentials, author the query,
call it from TypeScript.

## Prerequisites

- Package scaffold (`grok create <Name>`); paths are relative to the
  package root. Produce one with `create-package` if absent.
- Imports: `grok` from `datagrok-api/grok`, `_package` from
  `./package`. API key from `<GROK_HOST>/u` for step 2.

## Steps

1. **Declare the data connection in `connections/<name>.json`.**
   Required: `parameters` (data-source-specific) and `dataSource`
   (exact connector string — `Postgres`, `PostgresDart`, `MariaDB`,
   `MS SQL`, `Snowflake`, `Files`, ...; lowercase forms like
   `postgres` / `mssql` don't match). `name` defaults to the filename.
   Secrets go under `credentials.parameters` — never inside
   `parameters` or `connString`. The article example (line 41) uses
   `PostgresDart`; the live `Chembl` reference uses `Postgres` — both
   are valid, the connector registry accepts either. Scaffold with
   `grok add connection <name>`. Knowledge: `DG-FACT-400`,
   `DG-FACT-401`, `DG-FACT-409`.
   ```json
   {
     "#type": "DataConnection", "name": "Chembl",
     "parameters": {"server": "db.datagrok.ai", "port": 54325,
                    "db": "chembl", "cacheResults": true},
     "credentials": {"parameters": {"login": "", "password": ""}},
     "dataSource": "Postgres"
   }
   ```
   Expected: `grok publish` registers the connection as
   `<Package>:<Name>` (case-insensitive at lookup). The `#type` field
   is the deserializer discriminator used by the Chembl reference; the
   article omits it from its example but the live JSON includes it.
   Source:
   `{{ lattice.harness.help_develop_root }}/how-to/db/access-data.md:29-114`.
   Reference: `packages/Chembl/connections/chembl.json` (real shape
   including `#type`).

2. **Transfer credentials after the first deploy.** POST them once
   published so they route through Datagrok's credentials store. The
   dotted path targets one connection; `.../for/<Package>` (no dot)
   targets the package itself. Knowledge: `DG-FACT-409`.
   ```bash
   curl -X POST "$GROK_HOST/api/credentials/for/$PACKAGE.$CONNECTION" \
     -H "Authorization: $API_KEY" -H "Content-Type: application/json" \
     -d '{"login":"abc","password":"123"}'
   ```
   Expected: `200 OK`. Read back via `await _package.getCredentials()`
   → `c.parameters['<key>']` (NOT `c.openParameters`, the redacted
   display copy). Source:
   `{{ lattice.harness.help_develop_root }}/how-to/db/access-data.md:81-116`.
   Reference: `js-api/src/entities/misc.ts:163-166` (the
   `getCredentials` typing and `parameters` vs `openParameters` split).

3. **Author the parameterized SQL query in `queries/<basename>.sql`.**
   Each query is a `--`-prefixed header block followed by SQL,
   terminated by `--end`. `--end` is OPTIONAL when the file holds a
   single query — the parser only needs it to separate consecutive
   queries. Headers: `--name`, `--friendlyName`, `--description`,
   `--connection: <Package>:<Connection>` (omit if the package has
   only one), `--input: <type> <name> [= <default>]` (reference in SQL
   as `@name`). Keep one query per file when attaching a sibling
   `<basename>.js` post-process or `<basename>.layout` — those bind by
   basename, not by query name. Header annotation details (input
   types, friendlyName, description) follow the standard
   function-parameter syntax referenced from the article (line 165).
   Scaffold with `grok add query <name>`. Knowledge: `DG-FACT-402`,
   `DG-FACT-403`.
   ```sql
   --name: ProteinClassification
   --connection: Chembl
   --input: int level = 1
   select * from protein_classification where class_level = @level
   ```
   Expected: `grok publish` registers a function
   `<Package>:ProteinClassification`. Source:
   `{{ lattice.harness.help_develop_root }}/how-to/db/access-data.md:120-213`
   (header grammar at 132-156; sibling `.js`/`.layout` binding at
   189-213; annotation deferral at 165).
   Reference: `packages/Chembl/queries/queries.sql:1-6`.

4. **Execute the query from TypeScript.** `grok.data.query` is async
   and returns the result coerced to `T` (default `DataFrame`). Pass
   only `(queryName, parameters)` — the 3rd `adHoc` argument is
   JSDoc-marked `@deprecated`; there is no 4th polling-interval
   argument (verified against the JS-API signature, not the article).
   Knowledge: `DG-FACT-404`.
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const df = await grok.data.query(
     `${_package.name}:ProteinClassification`, {level: 1});
   grok.shell.addTableView(df);
   ```
   Expected: a `TableView` opens with the query result. For a scalar
   return type use `grok.data.query<number>(...)`. Source:
   `{{ lattice.harness.help_develop_root }}/how-to/db/access-data.md:169-185`.
   Reference: `js-api/src/data.ts:255-262` (signature, deprecation
   marker, positional-arg count).

## Common failure modes

- **`grok.data.query(...)` warns about the 3rd argument.** `adHoc`
  is `@deprecated`; an apparent 4th polling-interval argument doesn't
  exist in the signature. Drop both — pass only `(queryName,
  parameters)` (`DG-FACT-404`).
- **`dataSource` set to `postgres` / `mssql` — connection silently
  fails.** Server matches the string exactly; use `Postgres` (or
  `PostgresDart`), `MS SQL` (with space), `MySQL` (`DG-FACT-401`).
- **Login/password committed in `parameters` or `connString`.** That
  block is world-readable; the credentials store is bypassed. Move
  secrets under `credentials.parameters` and POST per step 2
  (`DG-FACT-400`, `DG-FACT-409`).
- **Sibling `<basename>.js` or `<basename>.layout` never applied.**
  The `.sql` file holds multiple queries; binding is by basename, not
  query name. Split into one query per file (`DG-FACT-403`).
- **Sharing the connection didn't share the queries.** Sharing flows
  one way: query → connection auto-shares; the reverse does NOT.
  Share each query explicitly (`DG-FACT-410`).

## Verification

- After step 1, the file is valid JSON and `grok publish` accepts it
  without a "missing dataSource" error.
- After step 2, `curl` returns `200 OK` and
  `await _package.getCredentials()` resolves with the posted keys.
- After step 3, the `<Package>:ProteinClassification` function appears
  in the function registry after `grok publish`.
- After step 4, running `<Package>:ProteinClassification()` in the
  Datagrok console opens a `TableView` with the expected rows.

## See also

- Source: `{{ lattice.harness.help_develop_root }}/how-to/db/access-data.md`
  (mirror: `docs/_internal/articles-mirror/how-to/db/access-data.md`).
- Knowledge (`docs/_internal/knowledge/knowledge-graph.md`):
  `DG-FACT-400` – `DG-FACT-404`, `DG-FACT-409`, `DG-FACT-410`.
- Related skills: `create-package` (scaffolds `connections/` and
  `queries/`), `manage-credentials` (general form of step 2),
  `fetch-external-rest` (REST via `grok.dapi.fetchProxy`),
  `load-bundled-csv` (`_package.webRoot` + `grok.data.loadTable`),
  `read-server-file` (`OpenServerFile` / `grok.dapi.files.*`).
