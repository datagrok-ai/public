# Add sticky metadata to a single cell — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI dataset | PASS | PASSED | SPGI.csv opened from System.DemoFiles; 3,624 rows, 88 columns |
| 2 | Select a row, click on cell under column Structure | PASS | PASSED | Row 0 selected, Structure column current cell set via df.currentRowIdx/currentCol API |
| 3 | Open Context Panel → Sticky Meta. Find TestSchema1 | FAIL | FAILED | Sticky meta accordion section appeared after Structure cell selected (semtype=Molecule); TestSchema1 NOT listed — schema was saved without entity type association (entity type picker UI failed in scenario 1); only `apisamples-int-meta-ntg6l` and `int-meta` schemas shown |
| 4 | Fill metadata: rating=5, notes="test note", verified=true, review_date=current date | SKIP | SKIP | Skipped — TestSchema1 not present in panel |
| 5 | Save the schema | SKIP | SKIP | Skipped — no schema to save |

**Add sticky column (sub-scenario)**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | In SPGI table, add sticky column via Context Panel | PASS | PASSED | Clicked "+" in Sticky meta pane; `int-meta` column appeared in grid header (sticky column added) |
| 2 | Verify column appears with header marker (circle) | AMBIGUOUS | AMBIGUOUS | Column "int-meta" appeared with dot marker (•) in header; circle indicator present |
| 3 | Sort/filter by sticky column | SKIP | SKIP | Not attempted — no real values set in the column |
| 4 | Remove sticky column | SKIP | SKIP | Not attempted |

**Batch edit metadata on multiple rows**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Select multiple rows, open Sticky Meta > Edit for all | SKIP | SKIP | TestSchema1 not available; batch edit not tested |

## Summary

The Sticky meta section correctly appears in the Context Panel when a Molecule-typed column cell is selected. Two schemas (`apisamples-int-meta-ntg6l` and `int-meta`) were visible for the molecule. A sticky column was successfully added to the grid view. However, TestSchema1 was not listed because the "Associated with" entity type picker failed in scenario 1 — the schema was saved without entity type, so it didn't appear for molecule cells.

## Retrospective

### What worked well
- Sticky meta accordion section automatically appears for Molecule semtype columns
- Sticky column "int-meta" was added and appears in the grid with a dot marker in the header
- The feature mechanism (schema-per-entity-type) works correctly for pre-configured schemas

### What did not work
- TestSchema1 not listed in Sticky meta panel — entity type association missing due to UI limitation in scenario 1
- The `setAllValues` API returned an error ("undefined") when called programmatically
- No inputs/form fields appeared inside the sticky meta pane for editing cell values (pane shows schema structure but no editable form)

### Suggestions for the platform
- When opening the "New Schema" dialog, "Associated with" should show all entity types in a searchable dropdown, not a breadcrumb navigator
- When a schema has no entity type, it should clearly appear in the panel for ALL cells (or clearly document the behavior)
- The Sticky meta pane should show input fields even when no metadata is set (not just schema structure)

### Suggestions for the scenario
- Add prerequisite: "TestSchema1 must have an entity type with matchBy=Molecule or semtype=Molecule"
- Clarify whether batch edit requires the context menu or is done via the Sticky Meta panel
- Add a step to verify the blue marker visually (take screenshot)
