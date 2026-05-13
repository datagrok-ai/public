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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# access-data

## When to use

Your package needs rows from an external SQL database materialized as a
`DataFrame` — with `grok publish` wiring the connection, credentials
kept out of git, and one JS call to run the query. The spine is four
steps: declare the connection, transfer credentials, author the query,
call it from TypeScript.

## Steps

1. **Declare the data connection in `connections/<name>.json`.**
   Scaffold with `grok add connection <name>`. Secrets go under
   `credentials.parameters`, never inside `parameters` or `connString`.
   See `DG-FACT-400` (file shape), `DG-FACT-401` (exact `dataSource`
   strings — `Postgres`/`PostgresDart`/`MariaDB`/`MS SQL`/...; lowercase
   forms don't match), `DG-FACT-409` (credentials routing).
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
   `<Package>:<Name>` (case-insensitive at lookup).

2. **Transfer credentials after the first deploy.** POST them once
   published so they route through Datagrok's credentials store. The
   dotted path targets one connection; `.../for/<Package>` (no dot)
   targets the package itself (see `DG-FACT-409`).
   ```bash
   curl -X POST "$GROK_HOST/api/credentials/for/$PACKAGE.$CONNECTION" \
     -H "Authorization: $API_KEY" -H "Content-Type: application/json" \
     -d '{"login":"abc","password":"123"}'
   ```
   Expected: `200 OK`. Read back via `await _package.getCredentials()`
   → `c.parameters['<key>']` (NOT `c.openParameters`).

3. **Author the parameterized SQL query in `queries/<basename>.sql`.**
   Scaffold with `grok add query <name>`. See `DG-FACT-402` (header
   grammar: `--name`, `--friendlyName`, `--description`, `--connection`,
   `--input`; `@name` references; `--end` separator) and `DG-FACT-403`
   (one query per file when attaching sibling `<basename>.js` post-
   process or `<basename>.layout` — binding is by basename).
   ```sql
   --name: ProteinClassification
   --connection: Chembl
   --input: int level = 1
   select * from protein_classification where class_level = @level
   ```
   Expected: `grok publish` registers `<Package>:ProteinClassification`.

4. **Execute the query from TypeScript.** Pass only `(queryName,
   parameters)` — see `DG-FACT-404` for signature (3rd `adHoc` is
   `@deprecated`; no 4th argument exists).
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const df = await grok.data.query(
     `${_package.name}:ProteinClassification`, {level: 1});
   grok.shell.addTableView(df);
   ```
   Expected: a `TableView` opens with the query result. For a scalar
   return type use `grok.data.query<number>(...)`.

## Common failure modes

- 3rd/4th arg to `grok.data.query(...)` — drop them (`DG-FACT-404`).
- `dataSource` set to lowercase `postgres`/`mssql` — exact match required (`DG-FACT-401`).
- Credentials in `parameters`/`connString` — move under `credentials.parameters` and POST per step 2 (`DG-FACT-400`, `DG-FACT-409`).
- Sibling `<basename>.js`/`.layout` not applied — split multi-query files (`DG-FACT-403`).
- Sharing connection didn't share queries — sharing is one-way query→connection (`DG-FACT-410`).

## See also

- Source: `help/develop/how-to/db/access-data.md`
  (mirror: `docs/_internal/articles-mirror/how-to/db/access-data.md`).
- Knowledge (`docs/_internal/knowledge/knowledge-graph.md`):
  `DG-FACT-400` – `DG-FACT-404`, `DG-FACT-409`, `DG-FACT-410`.
- Related skills: `create-package` (scaffolds `connections/` and
  `queries/`), `manage-credentials` (general form of step 2),
  `fetch-external-rest` (REST via `grok.dapi.fetchProxy`),
  `load-bundled-csv` (`_package.webRoot` + `grok.data.loadTable`),
  `read-server-file` (`OpenServerFile` / `grok.dapi.files.*`).
