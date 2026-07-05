# Inventory

Inventory is the reference app for the **hybrid / dynamic-schemas + scale story** of
entity-mapped domain schemas (`databases/inventory/schema.json`): plugin-owned relational
tables inside the Datagrok database with platform row/column security, Datlas-managed CRUD,
batch ingestion, and typed TypeScript clients generated from the manifest.

The schema:

* `items` — `table`-mode security, business key on `sku`. The relational core
  (`sku`, `name`, `quantity`, `location`) is complemented by per-department **jsonb property
  schemas**: *chemistry* (`cas_number`, `hazard_class`), *procurement* (`unit_cost`,
  `reorder_point`), *quality* (`inspection_status`). Different groups can be granted
  different schemas and thereby see different columns of the same table.
* `stock_movements` — `master` mode delegating security to the item, with **`audit: false`**:
  the high-churn example where the movement log itself is the history.

The app (`Apps | Inventory`) dogfoods `src/generated/db.ts` — the typed clients that
`grok api` generates from the manifest (run `grok api` in this package to regenerate).

## Demo walkthrough

1. **Import.** Open `Apps | Inventory`, click **Import...**, and pick a CSV or Parquet file
   with `sku,name,quantity,location` columns (extra columns matching jsonb properties, e.g.
   `cas_number`, also load). Rows are **batch-upserted by SKU**: re-importing an updated
   file updates existing items in place and inserts only the new SKUs — the import summary
   reports the inserted/updated/skipped counts. Parquet files are converted client-side via
   the Arrow package (install it for Parquet support).
2. **Adjust.** Select an item in the grid and click **Adjust stock...** — enter a signed
   delta and a reason. The adjustment updates `quantity` with **optimistic concurrency**
   (the update carries the row version and retries on a 409 conflict), and logs a
   `stock_movements` row. Quantities cannot go below zero (`min: 0` is engine-enforced).
3. **Aggregate.** Switch to the **Stock by location** tab: a `sum(quantity) group by location`
   aggregation computed server-side over the rows and columns the current user can see.
4. **Per-group columns.** Create the department groups and add users to them:

   ```
   # verify the schemas are registered and create the Chemists/Procurement groups
   grok s functions run 'Inventory:SetupInventoryDemo()'

   # add users to the departments (idempotent)
   grok s groups add-members Chemists alice
   grok s groups add-members Procurement bob
   ```

   Then grant row and column visibility to the groups. Domain security checks **direct**
   permission rows on the table and property-schema entities (the shape the package
   deployer writes for the publisher). Granting them to further groups has no platform
   API yet, so this step is a server-admin SQL statement for now (run against the
   Datagrok database, e.g. via the `System:Datagrok` connection):

   ```sql
   -- View on the items table (row visibility) + a department schema (column visibility)
   insert into permissions (id, entity_id, user_group_id, permission_id)
   select uuid_generate_v4(), s.id, g.id, '34da1550-e870-11e6-9cb3-825892686412' -- View
   from (select id from domain_tables where name = 'items'
         union select id from entity_property_schemas where name = 'chemistry') s,
        groups g
   where g.name = 'Chemists'
     and not exists (select 1 from permissions p where p.entity_id = s.id
       and p.user_group_id = g.id
       and p.permission_id = '34da1550-e870-11e6-9cb3-825892686412');
   ```

   (Repeat with `'procurement'` for Procurement.) Relational columns belong to the table's
   auto-generated *core* schema, granted to Everyone at deploy; the department schemas stay
   default-closed, so `alice` sees `cas_number`/`hazard_class` but not `unit_cost`, and
   `bob` the other way around — the same grid, different columns. Filters and sorts
   referencing a column the user cannot see are rejected server-side.

## Tests

`grok test` runs the CI suite (`src/tests/inventory-tests.ts`): batch CSV upsert count
assertions (the phase-2 milestone "stock sync upserts by business key"), Parquet import via
Arrow (falls back to asserting the clear no-Arrow error when Arrow is not deployed),
optimistic-concurrency conflict + retry, aggregation totals, per-schema column
attribution, and the demo group setup.

See `core/docs/features/db-schemas/ARCHITECTURE.md` (§9.4) in the Datagrok monorepo for the design.
