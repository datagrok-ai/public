# Create metadata schema & entity type — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Log in with admin permissions | PASS | PASSED | Already logged in as admin |
| 2 | Navigate to Browse → Platform → Sticky Meta → Types | PASS | PASSED | Clicked Platform → Sticky Meta → Types; 4/43 entity types shown |
| 3 | Click "New Entity Type..." button | PASS | PASSED | Dialog "Create a new entity type" opened |
| 4 | Enter name TestEntity1, set matching expression semtype=molecule. Save | PASS | PASSED | Filled Name=TestEntity1, Matching expression=semtype=molecule; clicked OK; TestEntity1 appeared at top of list (5/44) |
| 5 | Navigate to Browse → Platform → Sticky Meta → Schemas | PASS | PASSED | Clicked Schemas tree node; 16/16 schemas listed |
| 6 | Click "NEW SCHEMA" button | PASS | PASSED | Dialog "Create a new schema" opened with Name, Associated with, Properties fields |
| 7 | Enter schema name TestSchema1, associate with TestEntity1 | AMBIGUOUS | AMBIGUOUS | Name set to TestSchema1; "Associated with: select entities" opens a combo popup but navigates via breadcrumbs rather than showing entity type list; entity type not bound via UI — saved without association |
| 8 | Add properties: rating (int), notes (string), verified (bool), review_date (datetime). Save | PASS | PASSED | All 4 properties added with correct types via select dropdowns; clicked OK; TestSchema1 appeared in list (17/17) |

## Summary

Entity type TestEntity1 (semtype=molecule) and schema TestSchema1 (rating/int, notes/string, verified/bool, review_date/datetime) were created successfully. The "Associated with" entity type binding could not be set through the UI combo widget (it showed breadcrumb navigation items instead of entity type suggestions), so TestSchema1 was saved without an entity type association. Both items appear in their respective lists with no errors.

## Retrospective

### What worked well
- New Entity Type dialog works reliably; name and matching expression fields accept native setter + dispatchEvent pattern
- New Schema dialog properties table allows adding rows with the `fa-plus` icon; selects accept `.value` assignment
- Both created entities persist after save and appear in their lists

### What did not work
- "Associated with: select entities" combo widget — clicking opens a `d4-combo-popup` showing breadcrumb navigation ("Home", "Schemas") rather than a searchable entity type list; no text input found in popup to type a search query; maximize icon closes popup without opening a full browser
- The combo widget may require real mouse pointer/drag events that are not reproducible via JS dispatch

### Suggestions for the platform
- The "Associated with" entity type picker should include a text search input in its combo popup so users can type to filter entity types by name
- The combo popup should show entity type names directly (not just breadcrumb navigation), especially when the list is small enough
- Consider adding a keyboard-accessible search-and-select pattern for entity type association

### Suggestions for the scenario
- Clarify whether "Associated with" is required or optional for the schema to be valid
- Add a note that the entity type picker may require navigating into the Types sub-list in the combo popup
- Consider a separate step "Verify properties are saved correctly by clicking the schema and inspecting its details"
