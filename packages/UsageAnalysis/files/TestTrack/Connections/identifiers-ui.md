---
feature: connections
realizes: [views.connections]
priority: p2
target_layer: manual-only
coverage_type: edge
manual_only_reason: |
  The customerid blue highlighting in the rendered grid cells is drawn on
  canvas, so it must be checked by eye.
related_bugs: []
---

# Identifiers — manual UI checks

This is the **manual companion** to `identifiers.md`. Covers the part the
autotest cannot exercise: the `customerid` blue highlighting in the rendered
grid cells is drawn on canvas, so it must be checked by eye.


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

## Cleanup

- The final check already removes the identifiers config; if the run stopped
  early, right-click `test_postgres` → **Configure Identifiers...** and remove
  the CUSTOMER_ID entry so `identifiers.md` step 5 can re-create it cleanly on the next run.
- Close the opened `customers` table view.

---
{
  "order": 2,
  "datasets": []
}
