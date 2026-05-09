---
name: db-in-plugin
description: Add a custom Postgres database to a Datagrok package using the platform-managed instance
---

# db-in-plugin

## When to use

Your package needs persistent server-side storage and an empty Postgres
schema is enough — CRUD app, lookup tables, app state. Triggers: "give my
plugin its own DB", "store records that survive page reload", "build a
compound registry / lock table / dictionary inside Datagrok". For a
non-Postgres engine, a Postgres extension RDS doesn't ship, or seeded data
baked into an image, use `db-in-docker` instead (knowledge `DG-FACT-052`).

## Prerequisites

- A package scaffold (`grok create <Name>`). Paths below are relative to
  the package root.
- Knowledge of SQL DDL — every file under `databases/<Name>/` is plain
  Postgres SQL the platform applies in lexicographical order.
- Familiarity with package query files (`queries/*.sql`) — see the
  `access-data` skill (knowledge `DG-FACT-035`, `DG-FACT-036`).
- `grok publish ... --release` access for the target Datagrok instance —
  dev publishes don't propagate migrations to other users (knowledge
  `DG-FACT-055`).

## Steps

1. **Create the database directory.**
   The directory name `<Name>` becomes BOTH the connection name AND the
   Postgres schema name (knowledge `DG-FACT-051`). Use lowercase to match
   every canonical package (`hitdesign`, `plts`, `moltrack`, `biologics`).
   No `connections/<name>.json` file is needed — the platform
   auto-registers the connection at deploy (knowledge `DG-FACT-050`).
   ```bash
   mkdir -p databases/compounds
   ```
   Expected: `databases/compounds/` directory at the package root.

2. **Author the initial migration `0000_init.sql`.**
   Always namespace into a schema named after the directory — never use
   `public`. The plugin DB lives in the SAME Postgres instance as
   Datagrok's metadata and other plugins; unprefixed objects collide
   (knowledge `DG-FACT-052`, drift `DG-FACT-DRIFT-019`).
   ```sql
   CREATE SCHEMA IF NOT EXISTS compounds;

   CREATE TABLE compounds.list (
       id SERIAL PRIMARY KEY,
       smiles VARCHAR(512) NOT NULL,
       created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
   );

   GRANT ALL PRIVILEGES ON SCHEMA compounds TO CURRENT_USER;
   GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA compounds TO :LOGIN;
   GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA compounds TO :LOGIN;
   ```
   Expected: `databases/compounds/0000_init.sql` exists. The `:LOGIN`
   placeholder is mandatory — Datagrok substitutes the auto-issued
   Postgres role at apply-time; without these GRANTs queries get
   `permission denied` at runtime (knowledge `DG-FACT-056`, drift
   `DG-FACT-DRIFT-017`).

3. **Add forward-only migrations as the schema evolves.**
   Use a 4-digit zero-padded prefix + underscore + descriptive name. Files
   apply in lexicographical filename order — `0000_init.sql` before
   `0001_add_optional_column.sql` (knowledge `DG-FACT-054`). Migrations
   are forward-only: ADD COLUMN works; DROP COLUMN, ALTER TYPE that loses
   data, or DROP TABLE break previously-deployed instances. There is NO
   rollback (knowledge `DG-FACT-053`).
   ```sql
   -- databases/compounds/0001_add_supplier.sql
   ALTER TABLE compounds.list
       ADD COLUMN IF NOT EXISTS supplier VARCHAR(255);
   ```
   Expected: filename sorts AFTER `0000_init.sql`; new column visible
   after the next `--release` publish.

4. **Author queries against the connection.**
   Place SQL query files under `queries/`. Use the fully-qualified
   `--connection: <Package>:<Name>` grammar — capitals in `<Package>`
   match `package.json`, `<Name>` matches the directory name (knowledge
   `DG-FACT-057`). The article tutorial uses the looser `--connection:
   Compounds` short form, but every canonical package
   (`HitTriage:hitdesign`, `Plates:plts`) uses the namespaced form (drift
   `DG-FACT-DRIFT-018`).
   ```sql
   --name: insertCompound
   --connection: CompoundRegistrator:compounds
   --input: string smiles
   INSERT INTO compounds.list (smiles) VALUES (@smiles)
   RETURNING id, smiles;
   --end

   --name: getCompounds
   --connection: CompoundRegistrator:compounds
   SELECT id, smiles, created_at FROM compounds.list ORDER BY id DESC;
   --end
   ```
   Expected: `queries/compounds.sql` exists; query namespace shows up
   under `Functions → Queries → CompoundRegistrator` after publish.

5. **Call the queries from package code.**
   Use `grok.data.query('<Package>:<QueryName>', params)`. Returns a
   `Promise<DataFrame>` — see `access-data` knowledge `DG-FACT-022`.
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as ui from 'datagrok-api/ui';

   const smilesInput = ui.input.molecule('Molecule');
   const addButton = ui.button('Add Compound', async () => {
     await grok.data.query(
       `${_package.name}:insertCompound`,
       {smiles: smilesInput.value!},
     );
     grok.shell.info('Compound added.');
     const df = await grok.data.query(`${_package.name}:getCompounds`);
     grok.shell.addTableView(df);
   });
   ```
   Expected: clicking the button inserts a row and re-renders the table
   view with the new entry.

6. **Publish in `--release` mode to apply the migration.**
   Each new migration requires a release publish to propagate to other
   users on the shared server (knowledge `DG-FACT-055`). Dev publishes
   only update your own session.
   ```bash
   webpack
   grok publish <host> --release
   ```
   Expected: publish exits `0`. The connection appears under **Browse →
   Databases → CompoundRegistrator:compounds** with the new column
   visible in the schema browser.

## Common failure modes

- **`permission denied for table compounds.list` at query time.** The
  migration created the table but never granted to `:LOGIN`. The
  platform-issued login is a different role from the schema owner
  (knowledge `DG-FACT-056`, drift `DG-FACT-DRIFT-017`). Fix: append
  `GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA <name> TO :LOGIN;` (and
  the same for `ALL SEQUENCES`) at the end of the migration that creates
  them. Mirror the pattern from
  `packages/HitTriage/databases/hitdesign/0000_init.sql:18-20`.
- **Migration silently doesn't apply.** Filename prefix doesn't sort
  after the previous one, or `--release` was skipped (knowledge
  `DG-FACT-054`, `DG-FACT-055`). Fix: re-check lexicographical order
  (`ls databases/<name>/` should print files in apply order); republish
  with `--release`.
- **Migration breaks existing deployments.** Used `DROP COLUMN`, `ALTER
  TYPE` that loses data, or `DROP TABLE` — the platform has no rollback
  and rejects backwards-incompatible changes (knowledge `DG-FACT-053`).
  Fix: write a forward-compatible migration (add a new column, dual-write
  during transition, deprecate the old column without dropping it).
- **Query 404s with `connection not found`.** `--connection` value
  doesn't match `<Package>:<directoryName>` (knowledge `DG-FACT-057`,
  drift `DG-FACT-DRIFT-018`). Fix: use the fully-qualified form —
  `--connection: CompoundRegistrator:compounds`, not `--connection:
  Compounds`.
- **Schema collision with another plugin or platform metadata.** Tables
  created in `public` instead of a namespaced schema (drift
  `DG-FACT-DRIFT-019`). Fix: prefix every DDL object with `<name>.` and
  add `CREATE SCHEMA IF NOT EXISTS <name>;` at the top of `0000_init.sql`.

## Verification

- `ls databases/<name>/` lists all migrations in apply order.
- `grok publish <host> --release` exits `0`.
- **Browse → Databases** lists `<Package>:<name>`; clicking it opens the
  schema browser and the table you just created is visible.
- A read query (`grok.data.query('<Package>:getCompounds')`) returns a
  `DataFrame` without permission errors.
- A write query (`grok.data.query('<Package>:insertCompound', {...})`)
  returns the inserted row's `id` and a follow-up read sees the new row.

## See also

- Source articles:
  - `help/develop/how-to/db/db-in-plugin.md`
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts `DG-FACT-050`
    through `DG-FACT-057` and drifts `DG-FACT-DRIFT-016..019`.
- Related skills:
  - `access-data` (sibling — covers query file headers, parameter
    grammar, `grok.data.query` JS-API).
  - `db-in-docker` (alternative — when a non-Postgres engine, a Postgres
    extension RDS doesn't ship, or seeded data baked into an image is
    required).
