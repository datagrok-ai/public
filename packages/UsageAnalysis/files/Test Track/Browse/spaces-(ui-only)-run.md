# Spaces (UI only) — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 2.1 | Copy / Rename copy | AMBIGUOUS | No "Duplicate Space" action found in context menu |
| 2.2 | Rename link | SKIP | Requires linked entities in space |
| 2.3 | Rename entity in space (Browse Tree and view) | SKIP | Requires entities in space |
| 3.1 | DnD from local storage (group and single file) | SKIP | Drag-and-drop not feasible via browser automation |
| 3.2 | DnD between spaces: Link | SKIP | Drag-and-drop not feasible via browser automation |
| 3.3 | DnD between spaces: Copy | SKIP | Drag-and-drop not feasible via browser automation |
| 3.4 | DnD between spaces: Move | SKIP | Drag-and-drop not feasible via browser automation |
| 3.5 | DnD group of files: Link | SKIP | Drag-and-drop not feasible via browser automation |
| 3.6 | DnD group of files: Copy | SKIP | Drag-and-drop not feasible via browser automation |
| 3.7 | DnD group of files: Move | SKIP | Drag-and-drop not feasible via browser automation |
| 3.8 | DnD from Dashboards (link and move) | SKIP | Drag-and-drop not feasible via browser automation |
| 3.9 | DnD from Files: file/group Link | SKIP | Drag-and-drop not feasible |
| 3.10 | DnD from Files: file/group Copy | SKIP | Drag-and-drop not feasible |
| 3.11 | DnD from Files: file/group Move | SKIP | Drag-and-drop not feasible |
| 3.12 | DnD from Files: folder Link | SKIP | Drag-and-drop not feasible |
| 3.13 | DnD from Files: folder Copy | SKIP | Drag-and-drop not feasible |
| 3.14 | DnD from Files: folder Move | SKIP | Drag-and-drop not feasible |
| 3.15 | DnD from Packages | SKIP | Drag-and-drop not feasible |
| 3.16 | Alt + dragging (default action is link) | SKIP | Drag-and-drop not feasible |
| 4.1 | Preview entities (file, project, child, linked) | SKIP | No entities in space |
| 4.2 | Open file or project | SKIP | No entities in space |
| 4.3 | Search in the view | AMBIGUOUS | Search box present but filtering not verified |
| 5.1 | Share linked entity | SKIP | No linked entities |
| 5.2 | Share entity | SKIP | No entities in space |
| 5.3 | Share root space | PASS | Tested in Spaces.md scenario |
| 5.4 | Share child space | SKIP | Not tested separately |
| 6.1 | Delete linked space | SKIP | No linked spaces created |
| 7.1 | Moving functions from packages to spaces | SKIP | Drag-and-drop not feasible |
| 7.2 | Dragging root into its child | SKIP | Drag-and-drop not feasible |
| 7.3 | Dragging child to root Spaces level | SKIP | Drag-and-drop not feasible |

## Summary

1 step passed, 2 ambiguous, 27 skipped out of 30 total steps. This scenario focuses on UI-only interactions not covered by server autotests, primarily drag-and-drop operations. Browser automation via Chrome DevTools MCP cannot perform drag-and-drop between tree nodes reliably, making most steps untestable in this mode.

## Retrospective

### What worked well
- Share dialog opens correctly for spaces

### What did not work
- Almost all steps require drag-and-drop, which is not reliably supported by browser automation
- No entities were pre-loaded into spaces, so entity-level operations could not be tested

### Suggestions for the platform
- Provide JS API alternatives for drag-and-drop actions (e.g., `grok.dapi.projects.move(entityId, targetSpaceId)`)
- Add keyboard-driven alternatives for entity management in spaces

### Suggestions for the scenario
- Add setup instructions: "Before testing, create a space with at least one file, one project, and one linked entity"
- Consider providing a setup script that pre-populates test spaces with sample entities
- Mark drag-and-drop steps explicitly as "manual only" or provide API-based alternatives for automation
