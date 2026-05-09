---
name: data-enrichments
description: Ship a column enrichment with a Datagrok plugin by bundling a JSON file under enrichments/
---

# data-enrichments

## When to use

Your package needs to give users a one-click "join extra columns from a
database" action on a specific column — e.g. given a `customerid`, pull
`companyname` and `contacttitle`; given a ChEMBL `chembl_id`, attach
molecular weight and ALogP. Triggers: "ship enrichments with my plugin",
"add Enrich… button content", "expose lookup data in the column context
panel". For interactive enrichment authoring (the in-app "Enrich…" UI),
see `help/access/databases/databases.md` § Data enrichment — this skill
covers only the publish-with-plugin path.

## Prerequisites

- A package scaffold (`grok create <Name>`); enrichments live at
  `<PackageRoot>/enrichments/` (knowledge `DG-FACT-122`).
- A Datagrok DB connection that already exists — shipped by this package
  under `connections/` or by another deployed package. The enrichment
  references it by nqName `<Package>:<Name>` (knowledge `DG-FACT-123`).
- **PowerPack must be installed on the target Datagrok instance** —
  enrichments only surface in the column context panel when PowerPack
  is present (knowledge `DG-FACT-122`). No PowerPack → file is loaded
  but no UI.
- `grok publish ... --release` access; dev-mode publishes don't
  propagate plugin assets to other users.

## Steps

1. **Pick the lookup target — connection, key table, key column.**
   The enrichment fires when a user opens the column context panel on a
   column whose values match `keyColumn`. `connection` is the nqName of
   an existing Datagrok connection (`Chembl:Chembl`,
   `Dbtests:PostgresTest`). `keySchema` is the schema of `keyTable`
   (always `public` in the canonical examples — knowledge
   `DG-FACT-123`). `keyColumn` type must be `STRING`, `BIG_INT`, `INT`,
   or `FLOAT`; `BOOL` and `DATE_TIME` are rejected at runtime
   (knowledge `DG-FACT-127`).
   ```bash
   grok config show <host>          # confirm the alias resolves
   ```
   Expected: working alias for the target server. Verify the connection
   nqName in **Browse → Databases**; the displayed casing is the casing
   the JSON must use.

2. **Create the `enrichments/` folder at the package root.**
   The platform discovers enrichments by directory name; one JSON file
   per enrichment, no manifest required (knowledge `DG-FACT-122`).
   ```bash
   mkdir -p enrichments
   ```
   Expected: `enrichments/` directory next to `package.json`.

3. **Author the bundled-shape JSON.**
   Six required keys plus the `joins` array (knowledge `DG-FACT-123`).
   File name does not affect runtime behavior; shipping packages use
   `<key_table>_<purpose>.json`
   (`packages/DBTests/enrichments/company_info.json`,
   `packages/Chembl/enrichments/compound_id_properties.json`).
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
   Expected: `enrichments/<file>.json` parses. Use the BUNDLED shape
   (with `connection`, no `keyDb`) — the runtime `Enrichment` type
   PowerPack uses for user-saved enrichments has `keyDb` instead and is
   not interchangeable (knowledge `DG-FACT-DRIFT-048`).

4. **Fill `fields[]` with 3-segment selectors.**
   Each entry is `<schema>.<table>.<column>` — both source-table and
   joined-table columns may appear; include `keyTable.keyColumn` itself
   so the result row carries the lookup key (knowledge `DG-FACT-126`).
   No alias / rename syntax — result column name is the unqualified
   `<column>`.

5. **Define `joins[]` so every right-side table is reachable from `keyTable`.**
   Each entry has `leftTableName`, `rightTableName` (fully-qualified
   `schema.table`), `joinType`, `leftTableKeys`, `rightTableKeys`
   (knowledge `DG-FACT-124`). `joinType` ∈ {`left`, `inner`, `right`};
   every shipping enrichment uses `left` — default to `left` and only
   deviate when an inner / right join is intentional (knowledge
   `DG-FACT-125`). The two key arrays must be the same length;
   positional pairs form the join condition. Multi-hop chains use
   multiple `joins[]` entries listed in dependency order.

6. **Publish in `--release` mode.**
   The platform loads enrichment files when the plugin is published
   (knowledge `DG-FACT-122`). Dev publishes only update your own
   session.
   ```bash
   webpack
   grok publish <host> --release
   ```
   Expected: publish exits `0`; build output lists the new
   `enrichments/` files.

## Common failure modes

- **Enrichment never appears in the column context panel.** Either
  PowerPack is not installed on the target instance (knowledge
  `DG-FACT-122`), or the column the user clicked is `BOOL` /
  `DATE_TIME` and PowerPack renders "Cannot use column of <type> type
  for enrichment" inline (knowledge `DG-FACT-127`). Fix: install
  PowerPack; pick a `STRING` / `INT` / `BIG_INT` / `FLOAT` column.
- **`connection` not found at apply time.** Wrong package casing in
  the nqName (every shipping enrichment uses the
  first-letter-only-capitalized form — `Dbtests:…` not `DBTests:…`),
  wrong connection name, or the upstream package isn't installed.
  Fix: copy the nqName verbatim from **Browse → Databases**.
- **Field column missing from the result.** A `fields[]` entry
  references a table that no `joins[]` entry brings in (knowledge
  `DG-FACT-124`, `DG-FACT-126`). Fix: add a `joins[]` entry that
  joins the missing table to one already present, or drop the orphan
  field selector.
- **Join silently returns empty rows.** `leftTableKeys` and
  `rightTableKeys` differ in length, or named columns don't exist on
  the named tables — the platform produces an empty join rather than
  an error (knowledge `DG-FACT-124`; no production composite-key
  example exists). Fix: align array lengths; verify each key column
  via `information_schema.columns`.
- **Mixed bundled and runtime shapes.** Hand-writing the file with
  `keyDb` instead of `connection` after copying from a user-saved file
  under `System:AppData/PowerPack/enrichments/…` (knowledge
  `DG-FACT-DRIFT-048`). Fix: bundled JSON in `enrichments/` uses
  `connection`; never use `keyDb` in shipped files.

## Verification

- `jq . enrichments/<file>.json` exits `0`.
- Every `fields[]` entry is 3-segment:
  `jq -r '.fields[] | select((split(".") | length) != 3)' enrichments/<file>.json`
  prints nothing.
- Every join has equal-length key arrays:
  `jq -e '.joins | all((.leftTableKeys | length) == (.rightTableKeys | length))' enrichments/<file>.json`
  exits `0`.
- After `grok publish <host> --release`: open a table whose schema
  contains `keyColumn`, click that column, and the column context panel
  shows the enrichment under its `name` with an Apply button. Clicking
  Apply attaches the columns listed in `fields[]`.

## See also

- Source articles:
  - `help/develop/how-to/packages/data-enrichments.md`
  - `help/access/databases/databases.md` (§ Data enrichment — the in-app
    counterpart to publishing).
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts `DG-FACT-122`
    (folder layout + PowerPack requirement), `DG-FACT-123` (top-level
    keys), `DG-FACT-124` (join entry shape), `DG-FACT-125` (joinType
    enum), `DG-FACT-126` (fields selector grammar), `DG-FACT-127`
    (allowed column types), `DG-FACT-128` (runtime storage —
    informational, not part of the publish workflow), drift
    `DG-FACT-DRIFT-048` (bundled vs runtime shape divergence).
- Related skills:
  - `access-data` — connection nqName grammar.
  - `db-in-plugin` — sibling for shipping a connection alongside
    enrichments when the target DB is the platform-managed Postgres.
