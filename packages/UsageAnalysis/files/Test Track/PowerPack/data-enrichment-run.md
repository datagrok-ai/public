# Data Enrichment — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1.1 | Navigate to Databases > Postgres > NorthwindTest | PASS | N/A | Clicked Databases > Postgres > NorthwindTest in Browse tree |
| 1.2 | Create SQL query for Orders table | PASS | N/A | Created "orders VQ" with `select * from orders`, saved |
| 1.3 | Run query, open Orders table | PASS | N/A | 830 rows, 14 columns |
| 1.4 | Click customerid column | PASS | N/A | Used `grok.shell.o = col` to select column |
| 1.5 | Find Enrich section in context panel | PASS | N/A | NorthwindTest > Enrich showed existing enrichments |
| 1.6 | Click + Add enrichment, add Customers table join | PASS | N/A | Navigated icon-database-tables > northwind > public > customers |
| 1.7 | Save enrichment with name | PASS | N/A | Named "test enrichment", clicked SAVE |
| 1.8 | Click enrichment to apply — columns added | PASS | N/A | Columns grew from 14 to 24 (10 customer cols added) |
| 1.9 | Create visual query | SKIP | N/A | Not attempted |
| 1.10 | Edit enrichment | SKIP | N/A | Not attempted |
| 2.x | Multiple enrichments | SKIP | N/A | Not attempted |
| 3.x | Persistence across projects/layouts | SKIP | N/A | Not attempted |
| 4.x | Visibility for different users | SKIP | N/A | Requires multi-user login |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~42s |
| Spec file generation | ~3s |
| Spec script execution | FAILED — Dart submenu dispatchEvent doesn't work in Playwright |

## Summary

Part 1 (create, save, apply enrichment) passed. The enrichment feature correctly joins the Customers table to Orders via customerid, adding 10 new columns. Parts 2-4 (multiple enrichments, persistence, multi-user) were skipped.

## Retrospective

### What worked well
- Enrichment dialog correctly shows joined table preview with all columns
- Database schema browser navigation works via nested menus (with MCP hover)
- Saved enrichment appears immediately in the Enrich section
- Clicking enrichment name applies it and adds columns to the grid

### What did not work
- Dart submenu expansion requires MCP hover, not JS dispatchEvent
- Pressing Escape closed the enrichment dialog instead of just the table menu
- First click on enrichment appeared to do nothing (async column fetch)
- ENRICH button was disabled after adding table join — reason unclear

### Suggestions for the platform
- Add name= attributes to enrichment links in context panel
- Escape should close only the deepest open popup, not the parent dialog
- Show loading indicator when enrichment columns are being fetched

### Suggestions for the scenario
- Separate SQL and visual query creation into distinct steps
- Step 4 (multi-user) needs specific credentials
- Clarify single-click (apply) vs double-click (edit) for enrichments
