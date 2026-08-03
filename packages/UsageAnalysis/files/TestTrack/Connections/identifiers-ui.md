---
feature: connections
target_layer: playwright
coverage_type: edge
priority: p2
realizes: []
realized_as: []
related_bugs: []
---

# Identifiers — manual UI checks

This is the **manual companion** to `identifiers.md`. Covers the part that
`identifiers.test.ts` cannot exercise because the grid is canvas-rendered:
the auto-test verifies the column's semantic type via the JS API
(`tv.dataFrame.col(...).semType`), but the manual `customerid` blue
highlighting in the rendered grid cells lives entirely in canvas pixels —
unreachable from the DOM.

The autotest covers (when `DG_PG_PASSWORD` is set):

- Right-click `test_postgres` → **Configure Identifiers...**, fill Schema, fill
  identifier (CUSTOMER_ID / customers / customerid / `[A-Z]{5}`), SAVE
- Reload, expand `test_postgres > Schemas > public > customers`, **Get All**,
  verify `customerid` carries `semType = "CUSTOMER_ID"` server-side
- Right-click connection → remove identifiers config, reload, verify column
  no longer carries the type

## Pre-conditions

- `test_postgres` exists with valid credentials (`DG_PG_PASSWORD` set during
  the autotest run, or the connection edited manually with valid creds)
- Identifier configured per `identifiers.md` step 5

## Steps

1. Open **Browse > Databases > Postgres > test_postgres > Schemas > public > customers** → **Get All**
2. Wait for the grid to render
3. Locate the `customerid` column

## What to look for

- Values in the `customerid` column are rendered with **blue** text (the
  Datagrok identifier-highlight style)
- Click the `customerid` column header → the **Context Panel > Details**
  pane shows `Semantic type: CUSTOMER_ID`
- After **Remove identifiers config** + page reload, the same column renders
  with default text colour (no blue) and no `Semantic type` row in the panel
