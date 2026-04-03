# Spaces (UI only) — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 2a | Update > Copy / Rename copy | AMBIGUOUS | Context menu "Copy" submenu shows clipboard operations (ID, Grok name, Markup, URL) only. No "Duplicate space" found. If intent is copying a space entity, feature may be absent or accessible differently |
| 2b | Update > Rename link | SKIP | Could not test — no linked entities available in space |
| 2c | Update > Rename entity in space (Browse Tree) | SKIP | No entities in space during test |
| 2d | Update > Rename entity in space (view) | SKIP | No entities in space during test |
| 3a | Drag and Drop > From local storage (single file) to space | SKIP | Not tested — drag-and-drop not reliably automatable |
| 3b | Drag and Drop > From local storage (group) to space | SKIP | Not tested |
| 3c | Drag and Drop > Between spaces (Link) | SKIP | Not tested |
| 3d | Drag and Drop > Between spaces (Copy) | SKIP | Not tested |
| 3e | Drag and Drop > Between spaces (Move) | SKIP | Not tested |
| 3f | Drag and Drop > Group of files between spaces | SKIP | Not tested |
| 3g | Drag and Drop > From Dashboards (link/move, preserve table linking) | SKIP | Not tested |
| 3h | Drag and Drop > From Files (file/group/folder — Link/Copy/Move) | SKIP | Not tested |
| 3i | Drag and Drop > From Packages | SKIP | Not tested |
| 3j | Drag and Drop > Alt + dragging (default = Link) | SKIP | Not tested |
| 4a | Read > Preview entities (file, project, child space, linked file, linked project, linked function) | SKIP | No entities in space |
| 4b | Read > Open file or project | SKIP | No entities in space |
| 4c | Read > Search in view | PASS | Search box filters in real time; "Child" query returned ClaudeChildRenamed correctly (tested in Spaces.md run) |
| 5a | Share > Linked entity | SKIP | No linked entities |
| 5b | Share > Entity | SKIP | No entities |
| 5c | Share > Root space | PASS | Share dialog accessible via right-click → "Share..."; shows user/group/email input, owner "Full access", Advanced editor link |
| 5d | Share > Child space | SKIP | Child space deleted before this test |
| 6a | Delete > Linked space | SKIP | No linked space created during test |
| 7a | Forbidden > Moving functions from packages to spaces | SKIP | Not tested — requires drag interaction |
| 7b | Forbidden > Dragging root space into its child | SKIP | Not tested |
| 7c | Forbidden > Dragging child to root Spaces level | SKIP | Not tested (awaiting clarification per scenario) |

## Summary

This scenario covers the subset of Spaces functionality not validated by server autotests. Most drag-and-drop scenarios (section 3) and entity-level operations (sections 2, 4, 5 entity sub-items, 6, 7) were skipped due to either no content being present in the test space or the inherent limitations of browser automation for drag-and-drop. Search (4c) and Share root space (5c) were confirmed working.

## Retrospective

### What worked well
- Share root space dialog works correctly (same as Spaces.md)
- Search in space view filters correctly

### What did not work
- The vast majority of scenarios in this file require drag-and-drop, which is not reliably testable via DevTools automation
- "Copy / Rename copy" remains ambiguous — clipboard copy is available but duplication is not obvious

### Suggestions for the platform
- Provide a JS API for adding linked entities to spaces, enabling automated testing without drag-and-drop
- Consider adding a "Paste" or "Add link" button in the space view to allow adding entities without drag

### Suggestions for the scenario
- This scenario file is essentially a "manual testing checklist" — add explicit note that drag-and-drop sections require manual or Playwright-native testing
- Add a prerequisite step: populate the space with at least one linked project and one linked file before testing sections 2, 4, 5, 6
- Clarify "Copy / Rename copy" — specify whether this means duplicating a space or copying an entity within a space
