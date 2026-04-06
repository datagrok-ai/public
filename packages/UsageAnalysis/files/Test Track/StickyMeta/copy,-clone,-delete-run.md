# Copy / clone / move objects with metadata — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: FAIL

## Steps

**Copy / clone / move sub-scenario**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI table with sticky meta added | SKIP | SKIP | Prerequisite: no metadata was successfully added in scenario 2 (TestSchema1 entity type not associated) |
| 2 | Verify metadata on cloned table | SKIP | SKIP | No metadata to verify |
| 3 | Verify metadata on new view opened | SKIP | SKIP | No metadata to verify |
| 4 | Save as project and reopen | SKIP | SKIP | No metadata to verify |
| 5 | Move project to Space and open from there | SKIP | SKIP | No metadata to verify |
| 6 | Import/export | SKIP | SKIP | No metadata to verify |

**Delete metadata sub-scenario**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | On cell with metadata, open Context Panel → Sticky Meta | SKIP | SKIP | No metadata set on any cell |
| 2 | Delete fields rating and notes, save | SKIP | SKIP | No metadata to delete |
| 3 | Refresh and verify removal | SKIP | SKIP | No metadata to verify |

**Persistency sub-scenario**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Add metadata to objects | SKIP | SKIP | No metadata set in previous scenario |
| 2 | Refresh browser tab, logout/login | SKIP | SKIP | No metadata to verify |
| 3 | Check metadata on same objects | SKIP | SKIP | No metadata to verify |

## Summary

All steps skipped due to dependency on scenario 2 (Add and edit). No sticky metadata was successfully added to any cell because TestSchema1 was not properly associated with an entity type, causing it not to appear in the Sticky meta Context Panel. The entire scenario chain (add → copy/persist → delete) could not be tested.

## Retrospective

### What worked well
- N/A — scenario blocked by prerequisite failure

### What did not work
- Full scenario blocked by entity type association failure in scenario 1
- Sequential test dependency: if scenario 1 fails to set up the schema correctly, scenarios 2, 3 are all unrunnable

### Suggestions for the platform
- Sticky Meta scenarios should be more independent — each should set up its own data fixtures
- Consider a simpler default schema (e.g., one note field, no entity type required) for basic testing

### Suggestions for the scenario
- Add setup/teardown fixtures so scenarios are independently executable
- "Copy/clone" sub-scenario should explicitly state that a pre-existing metadata cell is required and provide steps to create it if missing
- Separate the three sub-scenarios into three individual scenario files with their own prerequisites
