---
name: db-in-plugin
version: 0.1.0
description: |
  Give a Datagrok package its own persistent SQL tables that live in the
  platform's shared Postgres instance, isolated by schema. For plugin
  authors who need a small CRUD backend for user-generated data
  (annotations, locks, lookup dictionaries) without standing up a
  separate database server. Produces a `databases/<name>/*.sql`
  migration directory and a `queries/*.sql` file that reads/writes
  through the auto-registered `<Pkg>:<name>` connection.
  Use when asked to "persist user annotations server-side", "add a
  small CRUD store to my package", or "give my plugin its own SQL
  tables without managing a server".
triggers:
  - persist annotations server-side
  - small crud backend for a package
  - package-owned sql tables
  - server-side storage for user-generated rows
  - shared sql tables for users
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# db-in-plugin

## When to use

Your package needs a place to put rows that survive page reloads and
are visible to every authorized user — a notes table, a lock table, a
seed dictionary — but you don't want to ship a sidecar database
container. The platform hosts your tables as a Postgres schema in its
own metadata DB and auto-registers the connection (knowledge
`DG-FACT-419`, `DG-FACT-420`). For an own-engine/own-extensions case,
use `db-in-docker` instead.

## Prerequisites

- A Datagrok server you can publish to with `--release` (debug
  publishes do NOT apply migrations — knowledge `DG-FACT-424`).
- Query-authoring conventions from `access-data` — this skill consumes
  the tables via `grok.data.query(...)` (knowledge `DG-FACT-425`).

## Steps

1. **Pick a schema name and create the directory.** The dir name becomes the schema, connection (`<Pkg>:<name>`), and namespace prefix (`DG-FACT-419`). Short lowercase identifier, no hyphens.
   ```bash
   mkdir -p packages/<PackageName>/databases/<name>
   ```

2. **Write the initial migration.** Create `0000_init.sql` (4-digit zero-padded — `DG-FACT-421`). Use `id SERIAL PRIMARY KEY` (the article's `floor(random()*1000)` collides — `DG-FACT-DRIFT-DBPLUGIN-003`).
   ```sql
   CREATE TABLE <name>.notes (
       id SERIAL PRIMARY KEY,
       author TEXT NOT NULL,
       body TEXT NOT NULL,
       created_at TIMESTAMPTZ DEFAULT NOW()
   );

   GRANT ALL PRIVILEGES ON SCHEMA <name> TO CURRENT_USER;
   GRANT ALL PRIVILEGES ON ALL TABLES   IN SCHEMA <name> TO :LOGIN;
   GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA <name> TO :LOGIN;
   ```
   `:LOGIN` is substituted at deploy time; schema-wide GRANT covers SERIAL sequences and future tables (`DG-FACT-423`, `DG-FACT-DRIFT-DBPLUGIN-002`).

3. **Write a query against the new tables.** Queries live under `queries/`, not `databases/`. Use the auto-registered `<PackageName>:<name>` connection (`DG-FACT-425`).
   ```sql
   --name: insertNote
   --connection: <PackageName>:<name>
   --input: string author
   --input: string body
   INSERT INTO <name>.notes (author, body)
   VALUES (@author, @body)
   RETURNING id;
   --end

   --name: getNotes
   --connection: <PackageName>:<name>
   SELECT id, author, body, created_at FROM <name>.notes
   ORDER BY created_at DESC;
   --end
   ```
   `--end` separates queries in a multi-query file. One-query files when you need a sibling `.js`/`.layout` (`DG-FACT-403`).

4. **Call the queries from the package** — same dispatcher as any other Datagrok query (`DG-FACT-425`).
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const id = await grok.data.query<number>(
     `${_package.name}:insertNote`,
     {author: grok.shell.user.login, body: 'hello'},
   );
   const df = await grok.data.query(`${_package.name}:getNotes`);
   grok.shell.addTableView(df);
   ```
   Expected: `insertNote` returns the new row's `id`; `getNotes` opens
   a TableView showing all rows.

5. **Release-publish to apply migrations.** Debug publishes skip migrations (`DG-FACT-424`); applied files tracked by NAME in lex order (`DG-FACT-421`, `DG-FACT-DRIFT-DBPLUGIN-004`).
   ```bash
   cd packages/<PackageName>
   grok publish <server> --release
   ```

6. **Add a forward-only follow-up migration.** Bump the prefix; only additive changes are safe (`DG-FACT-422`). Never edit a published migration — see `DG-FACT-DRIFT-DBPLUGIN-004`.
   ```sql
   -- 0001_add_tags.sql
   ALTER TABLE <name>.notes ADD COLUMN IF NOT EXISTS tags TEXT[];
   ```

## Common failure modes

- **`permission denied for relation <name>.<table>`** — single-table GRANT instead of schema-wide; misses sequences/future tables (`DG-FACT-DRIFT-DBPLUGIN-002`).
- **`duplicate key`** — `floor(random()*1000)` collides after ~37 rows. Use `SERIAL` (`DG-FACT-DRIFT-DBPLUGIN-003`).
- **Edited migration doesn't reach installs** — tracked by filename. Add a new file with next prefix (`DG-FACT-DRIFT-DBPLUGIN-004`).
- **Migration didn't run** — published without `--release` (`DG-FACT-424`).
- **Filenames sort wrong** — mixed prefix widths; zero-pad (`DG-FACT-421`).

## See also

- Source articles:
  - `help/develop/how-to/db/db-in-plugin.md`
  - `help/develop/how-to/db/db-in-docker.md` (alternative)
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-419..425`, drifts `DG-FACT-DRIFT-DBPLUGIN-001..004`. Also
  `DG-FACT-403` (query basename binding), `DG-FACT-424` (debug vs
  release publish).
- Related skills: `access-data` (query authoring + credentials),
  `db-in-docker` (own DB engine), `publish-packages` (release
  lifecycle).
