---
name: data-enrichments
version: 0.1.0
description: |
  Ship a declarative "one-click join from a database" action with a
  Datagrok plugin as a JSON file under `enrichments/`. When a user
  opens the context panel on a column whose values match the file's
  key column, PowerPack renders an Apply button that pulls in the
  listed fields from the joined tables. For package authors who want
  users to attach lookup columns without writing a query.
  Use when asked to "let users attach company info next to customer
  ids", "expose a postgres join in the column context panel", or
  "ship a one-click lookup with my plugin".
triggers:
  - attach lookup columns from a database
  - one-click join in the column context panel
  - pull related fields next to a chembl id
  - bundle a column-side database action with a plugin
  - expose a join from a connection alongside a column
  - add company info next to a customer id
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# data-enrichments

## When to use

Your package ships (or relies on) a Datagrok connection, and you want
users to one-click "join in these extra columns" whenever a column with
a matching key value (e.g. `customerid`, `chembl_id`) appears in any
open table. The user never writes SQL — they click the column header,
the context panel lists your enrichment, they hit **Apply**, and the
fields land in the dataframe.

## Prerequisites

- A package scaffold (`grok create <Name>`). Enrichments live at
  `<PackageRoot>/enrichments/`, one JSON file per enrichment
  (`DG-FACT-122`).
- A Datagrok connection already deployed — shipped by this package
  under `connections/<name>.json` or by another installed package.
  The enrichment references it by nqName `<Package>:<ConnectionName>`
  (`DG-FACT-123`).
- **PowerPack must be installed on the target Datagrok instance.** The
  enrichments runtime, the column context-panel UI, and the
  `runEnrichment` function all live in PowerPack (`DG-FACT-448`,
  `packages/PowerPack/src/db-explorer.ts`). Without it, the JSON is
  shipped but never surfaced.
- Release-mode publish — `grok publish <host> --release`. Dev publishes
  do not propagate plugin assets to other users.

## Steps

1. **Identify the trigger column and the verbatim connection nqName.**
   Open **Browse → Databases** in the target instance, find the
   connection, copy its nqName as displayed (every shipping enrichment
   uses first-letter-only-capitalized — `Dbtests:PostgresTest`,
   `Chembl:Chembl` — not `DBTests:`,
   `packages/DBTests/enrichments/company_info.json:3`). Note
   `keySchema.keyTable.keyColumn` for the column users will click.
   `keyColumn` must be `STRING | BIG_INT | INT | FLOAT`; PowerPack
   rejects `BOOL` and `DATE_TIME` inline with "Cannot use column of
   <type> type for enrichment" (`DG-FACT-127`,
   `packages/PowerPack/src/db-explorer.ts:396-398`).

2. **Create the `enrichments/` folder at the package root.**
   The platform discovers enrichments by directory name; no manifest
   needed (`DG-FACT-122`).
   ```bash
   mkdir -p enrichments
   ```
   Expected: `enrichments/` sits next to `package.json`.

3. **Author the bundled JSON.**
   Six top-level keys plus `joins[]` (`DG-FACT-123`). File name does
   not affect runtime; shipping packages use `<key_table>_<purpose>.json`
   (`packages/Chembl/enrichments/compound_id_properties.json`).
   ```json
   {
     "name": "Add company info",
     "connection": "Dbtests:PostgresTest",
     "keySchema": "public",
     "keyTable": "orders",
     "keyColumn": "customerid",
     "fields": [
       "public.orders.customerid",
       "public.customers.companyname",
       "public.customers.contacttitle"
     ],
     "joins": [
       {
         "leftTableName": "public.orders",
         "rightTableName": "public.customers",
         "joinType": "left",
         "leftTableKeys": ["customerid"],
         "rightTableKeys": ["customerid"]
       }
     ]
   }
   ```
   Use the BUNDLED shape (with `connection`, no `keyDb`). The runtime
   `Enrichment` interface PowerPack uses for user-saved enrichments
   omits `connection` and adds `keyDb` instead; the shapes are not
   interchangeable (`DG-FACT-DRIFT-DE-001`,
   `packages/PowerPack/src/db-explorer.ts:755-763`).

4. **List every result column in `fields[]` as 3-segment selectors.**
   Each entry is `<schema>.<table>.<column>`; both `keyTable` columns
   and joined-table columns may appear (`DG-FACT-126`). Include
   `keyTable.keyColumn` itself so the result row carries the lookup
   key. No alias / rename syntax — result column names are the
   unqualified `<column>`. See
   `packages/Chembl/enrichments/compound_records_docs.json:7-17` (9
   columns from one joined table).

5. **Bring every right-side table in via `joins[]`.**
   Five keys per entry — `leftTableName`, `rightTableName` (fully
   qualified `schema.table`), `joinType`, `leftTableKeys`,
   `rightTableKeys` (`DG-FACT-124`). Key arrays must be the same
   length; positional pairs form the condition. Every shipping
   enrichment uses single-element arrays — default to one key per
   side. `joinType` ∈ {`left`, `inner`, `right`} per the article; all
   20 shipping enrichments use `left` (`DG-FACT-125`) — default to
   `left`. For multi-hop chains, add one entry per hop in dependency
   order.

6. **Publish in `--release` mode.**
   The platform loads enrichment files when the plugin is published
   (`DG-FACT-122`). Dev publishes only update your own session.
   ```bash
   webpack
   grok publish <host> --release
   ```
   Expected: publish exits `0`; build output lists the new
   `enrichments/` files.

## Common failure modes

- **Nothing appears in the column context panel after publish.**
  Either PowerPack is not installed on the target instance
  (`DG-FACT-122`), or the focused column is `BOOL` / `DATE_TIME` and
  PowerPack renders the inline rejection message instead
  (`DG-FACT-127`). Fix: install PowerPack; pick a
  `STRING` / `INT` / `BIG_INT` / `FLOAT` column whose values match
  `keyColumn`.
- **`connection` not found at apply time.** Wrong package casing in the
  nqName (`Dbtests:…` not `DBTests:…`), wrong connection name, or the
  upstream package isn't deployed. Fix: copy the nqName verbatim from
  **Browse → Databases**.
- **A field column never lands in the result.** A `fields[]` entry
  points at a table no `joins[]` entry brings in (`DG-FACT-124`). Fix:
  add a join that connects the missing table to one already reachable
  from `keyTable`, or drop the orphan field.
- **Join silently returns empty rows.** `leftTableKeys` and
  `rightTableKeys` differ in length, or the named columns don't exist
  on the named tables — the platform produces an empty join, not an
  error (`DG-FACT-124`). Fix: align array lengths; verify each key
  column exists via `information_schema.columns`.
- **Bundled file copied from a user-saved enrichment fails to load.**
  User-saved enrichments under `System:AppData/PowerPack/enrichments/`
  use the `keyDb` shape with no `connection` (`DG-FACT-DRIFT-DE-001`).
  Plugin-bundled JSON must use `connection` and omit `keyDb`. Fix:
  rewrite the file in bundled shape.

## See also

- Source: `help/develop/how-to/packages/data-enrichments.md`. In-app
  counterpart: `help/access/databases/databases.md` § Data enrichment
  (this skill covers only the publish-with-plugin path).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — `DG-FACT-122`
  (folder + PowerPack requirement), `DG-FACT-123` (top-level keys),
  `DG-FACT-124` (join entry shape), `DG-FACT-125` (`joinType` enum),
  `DG-FACT-126` (`fields[]` grammar), `DG-FACT-127` (allowed key-column
  types), `DG-FACT-448` (PowerPack runtime owner),
  `DG-FACT-DRIFT-DE-001` (bundled vs runtime shape).
- Reference packages: `packages/DBTests/enrichments/` (3 files, simplest
  single-hop) and `packages/Chembl/enrichments/` (17 files, covers
  every shape this skill ships).
- Related skills: `access-data` (connection nqName grammar),
  `db-in-plugin` (ship the connection alongside the enrichment).
