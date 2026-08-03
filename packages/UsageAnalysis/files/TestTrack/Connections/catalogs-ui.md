---
feature: connections
target_layer: playwright
coverage_type: edge
priority: p2
realizes: []
realized_as: []
related_bugs: []
---

# Catalogs — manual UI checks

This is the **manual companion** to `catalogs.md`. The autotest
(`catalogs.test.ts`) is fully gated behind `test.skip(true)` because no
catalog-capable connection (MS SQL / NorthwindTest) exists on the dev test
account, so there is no fixture to drive the UI against. The whole scenario
delegates to this checklist.

## Pre-conditions

- A Datagrok account with at least one catalog-capable connection (e.g. the
  shared MS SQL `NorthwindTest`).

## Steps — Catalog browsing

1. Open **Browse > Databases**
2. Expand the catalog-capable connection (e.g. `MS SQL > NorthwindTest`)
3. Verify the tree shows a **Catalogs** node with the double-database icon
   (instead of a flat **Schemas** node)
4. Expand **Catalogs** — list of available DB catalogs renders
5. Expand one catalog — its schemas render under it
6. Expand a schema — its tables render

## Steps — Catalog preview

7. Click a catalog node so it becomes the current object
8. On the **Context Panel** (right) check that a preview area renders
9. Verify the catalog's name is displayed correctly in the panel header

## Steps — Catalog meta properties

10. Select a catalog
11. On the Context Panel locate the **Comment** and **LLM comment** meta inputs
12. Type a comment, click away to a sibling node, then click the catalog again — the comment must persist
13. Set an LLM comment the same way — verify it renders correctly

## Steps — Catalog schema view

14. Right-click a catalog node → **Browse**
15. A schema view opens listing all schemas inside the catalog
16. Right-click a catalog node → **Open as table**
17. A table view opens with all DB objects from the catalog

## What to look for

- Catalogs tree icons are correct (double-database for catalogs, single for schemas)
- No red error balloons during expansion / preview / comment edits
- Catalog comment/LLM-comment inputs persist after click-away/re-select
- Both **Browse** and **Open as table** open without errors
