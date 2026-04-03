# Schema — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Browse > Databases | PASS | PASSED | Databases view loaded showing all DB providers |
| 2 | Expand the Postgres provider | PASS | PASSED | Clicked expand triangle on Postgres tree node; CHEMBL, Datagrok, Northwind, Starbucks, SureCHEMBLDocker, test_postgres, World connections visible |
| 3 | Right-click Northwind and select Browse schema | AMBIGUOUS | N/A | Right-click context menu on Northwind does NOT contain "Browse schema" option; menu shows: Browse, New Query..., New Visual Query..., Delete..., Edit..., Rename..., Clone..., Clear cache, Browse queries, Test connection, Share..., Copy, Add to favorites. Used tree expansion instead: Northwind → Schemas → public → tables |
| 4 | Check you can interact with structures as with DB Tables | PASS | PASSED | Expanded Northwind → Schemas → public → visible tables: categories, customers, employees, orders, products, etc. Right-clicked `customers` → context menu shows: Get All, Get Top 100, New SQL Query..., New Visual Query... — standard DB table operations |

## Summary

3 of 4 steps passed, 1 ambiguous. The "Browse schema" context menu option was not found in the current UI, but the schema can be browsed via the tree (Northwind → Schemas → public → tables). Tables in the schema have proper DB interaction menus.

## Retrospective

### What worked well
- Northwind tree node expands to show Schemas subtree correctly
- public schema contains all expected Northwind tables (categories, customers, employees, orders, etc.)
- Table context menus show DB operations: Get All, Get Top 100, New SQL Query..., New Visual Query...

### What did not work
- "Browse schema" context menu item does not exist in the current UI — may have been removed or renamed to "Browse"
- Clicking "Browse" from context menu opens the queries browser, not a schema diagram view

### Suggestions for the platform
- Add "Browse schema" back to the connection context menu, or rename the existing option to clarify it opens a schema view
- Provide a dedicated schema diagram viewer (like dbdiagram or ER diagram) accessible from the connection context menu

### Suggestions for the scenario
- Update step 3: "Right-click the Northwind connection and expand Schemas in the tree" — the context menu option "Browse schema" no longer exists
- Add a precondition: verify which Postgres connection has accessible Northwind DB
