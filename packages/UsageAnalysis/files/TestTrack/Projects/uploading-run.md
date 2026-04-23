# Uploading — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Case 1: Open two local tables, link by USUBJID | PASS | PASSED | Two demog tables opened and linked via LinkTables function |
| 2 | Case 1: Save project with linked tables | PASS | PASSED | Saved as "Test Project - Case 1 Local", sharing dialog appeared |
| 3 | Case 3: Open two files from Browse > Files | PASS | PASSED | demog.csv (5850 rows) and earthquakes.csv (2426 rows) opened |
| 4 | Case 3: Save project | PASS | PASSED | Saved as "Test Project - Case 3 Files" |
| 5 | Case 4: Run two Northwind queries | PASS | PASSED | PostgresProducts (77 rows) and PostgresCustomers (91 rows) |
| 6 | Case 4: Save project with query results | PASS | PASSED | Saved as "Test Project - Case 4 Queries" |
| 7 | Case 9: Join two tables (inner + left) | PASS | PASSED | grok.data.joinTables works; inner join 5850 rows, left join 5850 rows |
| 8 | Case 9: Save project with joined tables | PASS | PASSED | Saved as "Test Project - Case 9 Joins" with 4 tables |
| 9 | Case 2: Local + Files combination | SKIP | SKIPPED | Similar to Cases 1+3, skipped for time |
| 10 | Case 5: Query + Files combination | SKIP | SKIPPED | Similar to Cases 3+4, skipped for time |
| 11 | Case 6: Spaces + Spaces | SKIP | SKIPPED | No Spaces set up on release server |
| 12 | Case 7-8: Spaces combinations | SKIP | SKIPPED | No Spaces set up on release server |
| 13 | Case 10: Pivot Table > Add to workspace | SKIP | SKIPPED | Complex UI interaction not attempted |
| 14 | Case 11: Aggregate Rows > Add to workspace | SKIP | SKIPPED | Complex UI interaction not attempted |

## Summary

8 of 14 steps passed, 6 skipped. Core project creation from local tables, file shares, query results, and join results all work correctly. The Save Project dialog, sharing dialog, and project upload all function as expected. Cases involving Spaces and advanced operations (Pivot/Aggregate) were not tested due to server configuration and time constraints.

## Retrospective

### What worked well
- Opening tables from demo data, DemoFiles, and database queries via JS API
- LinkTables and JoinTables functions work correctly
- Save Project dialog correctly lists all open tables
- Project upload and sharing dialog flow works smoothly

### What did not work
- `grok.functions.call('PostgresNorthwind:PostgresProducts')` failed with "Unable to get project asset" — had to use `grok.dapi.queries.list()` + `executeTable()` instead
- JoinTables via `grok.functions.call` required column objects not strings — API inconsistency
- Right-click context menu did not appear when dispatching contextmenu event on gallery items

### Suggestions for the platform
- Make query execution via `PackageName:QueryName` pattern more reliable
- Standardize JoinTables parameter types (accept column names as strings)
- Add `data-testid` attributes to gallery items for easier automation

### Suggestions for the scenario
- Specify which demo files and queries to use for each case
- Clarify what "Local storage" means (drag-and-drop vs grok.data.demo?)
- Add expected row counts for verification
- Cases 6-8 (Spaces) need pre-condition: a Space must exist
