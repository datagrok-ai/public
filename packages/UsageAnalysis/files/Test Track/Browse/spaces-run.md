# Spaces — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1.1 | Create root space | PASS | Created "ClaudeTestSpace" via + button in Spaces view |
| 1.2 | Create child space | PASS | Created "ChildSpace1" under ClaudeTestSpace via bottom + button |
| 1.3 | Multiple levels of nesting | PASS | Created "GrandChild1" under ChildSpace1 (3 levels deep) |
| 2.1 | Rename root | PASS | Renamed ClaudeTestSpace → ClaudeTestRenamed via context menu > Rename |
| 2.2 | Rename child | PASS | Context menu > Rename confirmed working on child spaces |
| 2.3 | Add to favorites | PASS | Context menu "Add to favorites" works, space appears in Favorites |
| 2.4 | Remove from favorites | SKIP | Not tested in this run |
| 2.5 | Copy / Rename copy | SKIP | Context menu "Copy" copies identifiers, no duplicate-space action found |
| 2.6 | Link / Rename link | SKIP | Requires entities inside the space |
| 2.7 | Rename entity in space (Browse Tree and view) | SKIP | Requires entities inside the space |
| 3.1 | DnD from local storage to space | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 3.2 | DnD between spaces (Link/Copy/Move) | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 3.3 | DnD group of files between spaces | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 3.4 | DnD from Dashboards (link/move) | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 3.5 | DnD from Files (file/folder, link/copy/move) | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 3.6 | DnD from Packages | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 3.7 | Alt + dragging (default action is link) | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 4.1 | Preview entities (file, project, child space, etc.) | SKIP | No entities added to space for preview |
| 4.2 | Open file or project | SKIP | No entities in space to open |
| 4.3 | Search in view | PASS | Search box present and functional in Spaces view |
| 5.1 | Share linked entity | SKIP | Skipped to avoid affecting other users on public server |
| 5.2 | Share entity | SKIP | Skipped to avoid affecting other users on public server |
| 5.3 | Share root space | SKIP | Skipped to avoid affecting other users on public server |
| 5.4 | Share child space | SKIP | Skipped to avoid affecting other users on public server |
| 6.1 | Delete linked files/projects/functions from space | SKIP | No linked entities created |
| 6.2 | Delete linked space | SKIP | No linked spaces created |
| 6.3 | Delete space (root and child) | PASS | All test spaces deleted successfully, cleanup complete |
| 6.4 | Delete files/projects in space | SKIP | No files/projects in space |
| 7.1 | Moving functions from packages to spaces | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 7.2 | Dragging root space into its child | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 7.3 | Dragging child to root Spaces level | SKIP | Drag-and-drop not feasible via Chrome DevTools MCP |
| 7.4 | Creating duplicate space names | SKIP | Not tested in this run |
| 7.5 | Using empty name | SKIP | Not tested in this run |

## Summary

8 steps passed, 0 failed, 24 skipped out of 32 total steps. Core CRUD operations (create root/child/nested, rename, add to favorites, search, delete) all work correctly. Most skipped steps require drag-and-drop interactions (not feasible via Chrome DevTools MCP) or entities within spaces. Share tests were skipped to avoid impacting other users on the public server.

## Retrospective

### What worked well
- Space creation via both + buttons (root and child) works reliably
- Rename dialog properly pre-fills current name and responds to keyboard input
- Delete confirmation dialog works with proper warning text
- Add to favorites via context menu confirmed working
- Tree navigation and node selection works reliably via JS evaluation

### What did not work
- `fill()` tool doesn't trigger Datagrok's Dart-based change listeners — OK button stays disabled. Fixed by using `press_key Control+A` then `type_text` which sends real keyboard events
- Drag-and-drop operations cannot be tested via Chrome DevTools MCP — covers ~40% of the scenario
- Dialog inputs don't use standard `<input>` elements — snapshot-based `uid` approach needed

### Suggestions for the platform
- Add `data-testid` attributes to dialog inputs and buttons for automated testing
- Add a JS API for space management (create, rename, delete, move entities) to enable programmatic testing
- The "Copy" context menu item copies identifiers — consider renaming to "Copy identifier" to avoid confusion with "Duplicate"

### Suggestions for the scenario
- Split drag-and-drop steps into a separate scenario marked as "manual only" or provide JS API alternatives
- Add pre-condition: "Ensure you have a test space with sample entities before testing Read/Share/Delete entity steps"
- Clarify what "Copy / Rename copy" means — is it duplicating the space or copying entities?
- Consider adding JS API commands as alternatives for drag-and-drop steps
