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
   Open **Browse → Databases**, copy the nqName as displayed (e.g.
   `Dbtests:PostgresTest`, `Chembl:Chembl` — first-letter-only-cap).
   `keyColumn` must be `STRING|BIG_INT|INT|FLOAT` — PowerPack rejects
   `BOOL`/`DATE_TIME` inline (`DG-FACT-127`).

2. **Create the `enrichments/` folder at the package root.**
   ```bash
   mkdir -p enrichments
   ```

3. **Author the bundled JSON.** Six required keys plus `joins[]`
   (`DG-FACT-123`). File name does not affect runtime.
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
   Use BUNDLED shape (with `connection`, no `keyDb`) — the runtime
   shape PowerPack reads for user-saved enrichments is different and
   not interchangeable (`DG-FACT-DRIFT-DE-001`).

4. **List every result column in `fields[]` as 3-segment selectors.**
   `<schema>.<table>.<column>`; include `keyTable.keyColumn` so the
   result carries the lookup key. No alias/rename syntax (`DG-FACT-126`).

5. **Bring every right-side table in via `joins[]`.** Five keys per
   entry (`DG-FACT-124`); key arrays must be same length, positional
   pairs form the condition. Default `joinType: left` (`DG-FACT-125`,
   `inner`/`right` documented but unused in shipping packages). For
   multi-hop chains, one entry per hop in dependency order.

6. **Publish in `--release` mode.** Dev publishes don't propagate
   plugin assets to other users (`DG-FACT-122`).
   ```bash
   webpack
   grok publish <host> --release
   ```

## Common failure modes

- Nothing in column context panel — PowerPack missing, or column is `BOOL`/`DATE_TIME` (`DG-FACT-122`, `DG-FACT-127`).
- `connection` not found — wrong package casing in nqName; copy verbatim from **Browse → Databases**.
- Field never lands — `fields[]` entry refers to a table no `joins[]` brings in (`DG-FACT-124`).
- Empty join silently returned — `leftTableKeys`/`rightTableKeys` length mismatch or wrong column names (`DG-FACT-124`).
- Bundled file copied from user-saved enrichment fails — user-saved uses `keyDb`, bundled uses `connection`; rewrite (`DG-FACT-DRIFT-DE-001`).

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
