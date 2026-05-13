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

1. **Pick a schema name and create the directory.**
   The directory name becomes the Postgres schema, the connection
   (`<Pkg>:<name>`), and the namespace prefix in query headers
   (knowledge `DG-FACT-419`). Use a short lowercase identifier — no
   hyphens, no mixed case (production: `hitdesign`, `plts`, `todo`).
   ```bash
   mkdir -p packages/<PackageName>/databases/<name>
   ```
   Expected: an empty directory. No JSON file is required — the
   platform scans `databases/` on release.

2. **Write the initial migration.**
   Create `0000_init.sql`. Use a 4-digit zero-padded prefix; mixing
   widths inside one directory sorts wrong (knowledge `DG-FACT-421`).
   Use `id SERIAL PRIMARY KEY` — the article's `floor(random()*1000)`
   pattern collides after ~37 rows (drift `DG-FACT-DRIFT-DBPLUGIN-003`).
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
   Expected: a single file at
   `packages/<PackageName>/databases/<name>/0000_init.sql`. The
   `:LOGIN` token is substituted at deploy time by the runtime user;
   the schema-wide GRANT (not single-table) covers SERIAL sequences
   and future tables in the same migration directory (knowledge
   `DG-FACT-423`, drift `DG-FACT-DRIFT-DBPLUGIN-002`).

3. **Write a query against the new tables.**
   Queries live under `queries/`, NOT `databases/`. Reference the
   connection as `<PackageName>:<name>` — the platform auto-registers
   it at install time (knowledge `DG-FACT-425`).
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
   Expected: a file under `packages/<PackageName>/queries/notes.sql`.
   `--end` is required between queries in a multi-query file. If a
   query needs a sibling `.js` post-process or `.layout`, put it in
   its own one-query file (knowledge `DG-FACT-403`).

4. **Call the queries from the package.**
   The same dispatcher as any other Datagrok query — there is no
   plugin-DB-specific JS API (knowledge `DG-FACT-425`).
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

5. **Release-publish to apply migrations.**
   Debug publishes do not run migrations — only `--release` does
   (knowledge `DG-FACT-424`). Migrations are applied in lexicographic
   filename order; the platform tracks applied files by NAME, not
   content (knowledge `DG-FACT-421`, drift `DG-FACT-DRIFT-DBPLUGIN-004`).
   ```bash
   cd packages/<PackageName>
   grok publish <server> --release
   ```
   Expected: the build succeeds and the server log shows the
   migration was applied. The new connection appears under
   `Browse > Databases > <PackageName>:<name>`.

6. **Add a forward-only follow-up migration.**
   Bump the prefix; only additive changes are safe (knowledge
   `DG-FACT-422`). Never edit or rename a published migration — the
   platform sees the filename as already-applied and skips it, so the
   new SQL never reaches existing installs (drift
   `DG-FACT-DRIFT-DBPLUGIN-004`).
   ```sql
   -- 0001_add_tags.sql
   ALTER TABLE <name>.notes ADD COLUMN IF NOT EXISTS tags TEXT[];
   ```
   Expected: a new file `databases/<name>/0001_add_tags.sql`. Re-run
   `grok publish <server> --release` to roll it out.

## Common failure modes

- **`permission denied for relation <name>.<table>`** at first query.
  The migration ran as the install user but never granted access to
  the runtime user. Fix: end every migration with the schema-wide
  `GRANT ... TO :LOGIN` lines from step 2 — the article's single-table
  `GRANT ALL ON TABLE ...` form silently misses sequences and future
  tables (drift `DG-FACT-DRIFT-DBPLUGIN-002`).
- **`duplicate key value violates unique constraint` on the primary
  key.** You followed the article's `floor(random()*1000)::int` pattern;
  it collides after ~37 inserts. Fix: use `id SERIAL PRIMARY KEY`
  (drift `DG-FACT-DRIFT-DBPLUGIN-003`).
- **Edited an already-published migration; the change isn't visible
  on installed clients.** The platform tracks applied migrations by
  filename. Fix: revert the edit, add a new file with the next prefix
  (drift `DG-FACT-DRIFT-DBPLUGIN-004`).
- **`grok publish` succeeded but the migration didn't run.** You
  published without `--release` — debug mode deploys per-user and
  skips migrations. Fix: republish with `--release` (knowledge
  `DG-FACT-424`).
- **Filenames sort wrong (`10_x.sql` before `2_x.sql`).** Mixed prefix
  widths. Fix: zero-pad consistently (`0010_x.sql`, `0002_x.sql`) —
  knowledge `DG-FACT-421`.

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
