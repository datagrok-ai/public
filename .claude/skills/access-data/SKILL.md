---
name: access-data
description: Pull data into a Datagrok package via package connections, queries, REST proxies, and file shares
---

# access-data

## When to use

Your package needs to load a DataFrame from somewhere outside the browser
— a SQL database, a REST endpoint, a file share, or a CSV bundled with
the package itself. Triggers: "connect to Postgres from my package",
"call a parameterised query", "fetch from an API without CORS", "open a
CSV from /Home", "ship a connection.json".

## Prerequisites

- A package scaffold (`grok create <Name>`). Paths below are relative to
  the package root.
- `datagrok-api` available — entry points are `grok.data.query`,
  `grok.functions.call`, `grok.dapi.fetchProxy`, `grok.data.loadTable`,
  `grok.dapi.files.*` (knowledge `DG-FACT-022`/`023`/`026`/`028`/`032`).
- A user API key from `<GROK_HOST>/u` for POSTing credentials
  (knowledge `DG-FACT-039`).

## Steps

1. **Add a database connection JSON.**
   Drop a file under `connections/<name>.json`. `name` is optional (defaults
   to filename); `dataSource` and `parameters` are required; sensitive
   values go in a separate `credentials.parameters` block (knowledge
   `DG-FACT-033`, `DG-FACT-040`).
   ```json
   {
     "name": "Chembl",
     "parameters": {"server": "$GROK_DB_SERVER", "db": "chembl_24"},
     "credentials": {"parameters": {"login": "", "password": ""}},
     "dataSource": "PostgresDart",
     "description": "CHEMBL db",
     "tags": ["demo", "chem"]
   }
   ```
   Expected: a JSON file at `connections/chembl.json` that
   `grok publish` will register. See
   `packages/Chembl/connections/chembl.json` for a real example. For
   allowed `dataSource` values see knowledge `DG-FACT-034`.

2. **Transfer credentials after deploy.**
   Don't commit secrets — POST them after publishing the package so they
   route through the credentials store (knowledge `DG-FACT-039`).
   ```bash
   curl -X POST "$GROK_HOST/api/credentials/for/$PACKAGE_NAME.$CONNECTION_NAME" \
     -H "Authorization: $API_KEY" \
     -H "Content-Type: application/json" \
     -d '{"login":"abc","password":"123"}'
   ```
   Expected: `200 OK`. Get `$API_KEY` from your profile page (e.g.
   `https://public.datagrok.ai/u`).

3. **Author a parameterised SQL query.**
   Place the file under `queries/<name>.sql`. One query per file so a
   sibling `<name>.js` post-process script and `<name>.layout` can
   bind unambiguously (knowledge `DG-FACT-035`, `DG-FACT-037`).
   ```sql
   --name: ProteinClassification
   --connection: Chembl
   --input: int level = 1
   select * from protein_classification where class_level = @level
   ```
   Expected: `grok publish` registers a function
   `<PackageName>:ProteinClassification` discoverable in the platform
   console. Omit `--connection` if the package has only one connection;
   omit `--end` if the file has only one query.

4. **Run the query from JavaScript.**
   Use `grok.data.query` with the namespaced name. The third positional
   argument (`adHoc`) is `@deprecated` — do not pass it (knowledge
   `DG-FACT-022`, drifts `DG-FACT-DRIFT-009`/`010`).
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const df = await grok.data.query(
     `${_package.name}:ProteinClassification`,
     {level: 1},
   );
   grok.shell.addTableView(df);
   ```
   Expected: a TableView opens with query results. `grok.functions.call`
   is the equivalent generic form (knowledge `DG-FACT-023`).

5. **Proxy a REST request through Datagrok.**
   Bypass CORS and gain server-side caching with
   `grok.dapi.fetchProxy`. Same shape as `fetch`; pass `maxAge` (seconds)
   to cache GET/HEAD responses (knowledge `DG-FACT-026`).
   ```typescript
   const resp = await grok.dapi.fetchProxy(
     'https://jsonplaceholder.typicode.com/posts',
     {method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({title: 'hello'})},
   );
   const data = await resp.json();
   ```
   Expected: `resp.ok === true`; `data` is the parsed JSON.

6. **Load a CSV bundled with the package.**
   Static assets resolve under `_package.webRoot` (knowledge
   `DG-FACT-027`). `grok.data.loadTable` handles the URL → DataFrame step
   (knowledge `DG-FACT-028`).
   ```typescript
   import {_package} from './package';

   const df = await grok.data.loadTable(
     `${_package.webRoot}data-samples/test.csv`);
   grok.shell.addTableView(df);
   ```
   Expected: a TableView opens with the contents of
   `<package>/data-samples/test.csv`.

7. **Open a file from a file share.**
   File-share paths follow the `<NAMESPACE>:<CONNECTION>/<path>` grammar
   (knowledge `DG-FACT-025`). `OpenServerFile` returns a list of tables;
   take `[0]` (knowledge `DG-FACT-024`).
   ```typescript
   const tables = await grok.functions.eval(
     `OpenServerFile("${grok.shell.user.login}:Home/data.csv")`);
   grok.shell.addTableView(tables[0]);
   ```
   Expected: the user's `/Home/data.csv` opens. For a package-bundled
   share, swap the prefix for `<PackageName>:<ConnectionName>/...`.

8. **Read/write file-share contents from code.**
   To register a shared directory, ship a `Files` connection with `dir`
   and `indexFiles` (knowledge `DG-FACT-038`); use the same JSON shape
   as step 1 with `"dataSource": "Files"`. `grok.dapi.files.*` accepts
   a path, a `FileInfo`, or a connection GUID; all methods are async
   (knowledge `DG-FACT-031`, `DG-FACT-032`).
   ```typescript
   const home = `${grok.shell.user.login}:Home`;
   await grok.dapi.files.writeAsText(`${home}/note.txt`, 'hello');
   const text = await grok.dapi.files.readAsText(`${home}/note.txt`);
   ```
   Expected: the file lands in the user's `/Home`; the read returns
   `'hello'`.

## Common failure modes

- **`grok.data.query(...)` extra args are silently dropped or
  deprecated.** The 4th polling-interval arg in older docs does not
  exist; the 3rd `adHoc` arg is `@deprecated` (drifts
  `DG-FACT-DRIFT-009`/`010`). Fix: only pass `(name, params)`.
- **Credentials end up in the connection JSON and get checked in.**
  Login/password belong in `credentials.parameters` and must NEVER be
  embedded in `connString` (knowledge `DG-FACT-040`). Fix: split them
  out and POST credentials per step 2.
- **`OpenServerFile(...)` returns an empty result, no tables open.**
  Path grammar is `<NAMESPACE>:<CONNECTION>/<path>` and is
  case-insensitive on the connection (knowledge `DG-FACT-025`,
  `DG-FACT-033`). Fix: confirm the connection name and that the file
  exists under the connection's `dir`.
- **CORS error from a direct `fetch`.** Direct browser fetch hits
  cross-origin policies; the platform offers a proxy. Fix: switch to
  `grok.dapi.fetchProxy` (knowledge `DG-FACT-026`).
- **Wrong `openTable` overload.** `grok.data.openTable(id)` takes a
  GUID; `grok.data.files.openTable(path)` takes a path (knowledge
  `DG-FACT-030`). Fix: pick the form matching your input.
- **Query parser swallows two queries as one.** Without `--end` the
  parser cannot tell where one query stops (knowledge `DG-FACT-035`).
  Fix: add `--end` between queries, or split into separate `.sql` files
  (recommended — see step 3).
- **Sharing a connection didn't share the queries.** Sharing flows
  asymmetrically: query → connection auto-shares, connection → query
  does NOT (knowledge `DG-FACT-041`). Fix: share each query explicitly
  or share at the package/space level.

## Verification

- After step 4, run `<PackageName>:ProteinClassification()` in the
  Datagrok console; a TableView opens with the expected row count.
- After step 5, `resp.status === 200` and `await resp.json()` parses.
- After step 7, `tables.length >= 1` and `tables[0].rowCount > 0`.
- After step 8, `await grok.dapi.files.exists(<path>) === true`.

## See also

- Source articles:
  - `help/develop/how-to/db/access-data.md`
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts `DG-FACT-022`
    through `DG-FACT-041` and drifts `DG-FACT-DRIFT-009..011`.
- Related skills:
  - `build-an-app` (apps consume query results to render views).
  - `user-settings-storage` (cache last-used query inputs per user).
