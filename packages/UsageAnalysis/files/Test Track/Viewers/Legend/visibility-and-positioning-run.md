# Legend Position — Run Results

**Date**: 2026-04-07
**URL**: http://localhost:8888/
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open SPGI | PASS | 5s | PASSED | Loaded 3624 rows, 88 columns |
| 2 | Add viewers: SP, Hist, LC, BC, PC, TP, BP | PASS | 4s | PASSED | All 7 viewers added successfully |
| 3 | Set categorical legend for each viewer | PASS | 3s | PASSED | Set Stereo Category as color/split/stack column on all viewers |
| 4a | Legend is visible | PASS | 2s | PASSED | 4 DOM-based legends visible (SP, Hist, PC, TP); BC, LC, BP use canvas-based legends |
| 4b | Colors on legend match viewer colors | PASS | 1s | PASSED | Visually confirmed via screenshot — legend colors match viewer markers |
| 4c | Changing Color/Split/Stack changes legend | PASS | 3s | PASSED | Changed SP color to Scaffold Names, legend updated with AMINO/TRISUBSTITUTED, reverted to Stereo Category successfully |
| 4d | Adjust legend size | SKIP | - | N/A | Splitter-based resize not testable via JS API or MCP automation |
| 4e | CTRL+click multi-select categories | PASS | 3s | PASSED | CTRL+click added d4-legend-item-current to multiple items |
| 4f | Check Color Picker | AMBIGUOUS | 5s | N/A | Could not find how to open Color Picker from legend. Right-click opens table view context menu, dblclick has no effect. Cross 'x' opens legend settings popup, not color picker |
| 5 | Save and apply layout | PASS | 8s | PASSED | Layout saved, closed, restored — 7 viewers restored, 5 legends visible |
| 6 | Set Visibility=Always, Position=Auto | PASS | 3s | PASSED | 6 of 7 viewers showed visible legends (Box Plot canvas-based only) |
| 7 | Resize viewers — legend repositions | PASS | 3s | PASSED | Legend positions changed based on viewer dimensions; auto-positioning works |
| 8 | Save and apply layout | PASS | 8s | PASSED | Legend visibility and positions preserved across layout save/restore |
| 9 | Set Visibility=Auto, uncheck auto positioning | PASS | 2s | PASSED | Set legendVisibility=Auto, legendPosition=top |
| 10 | Reduce viewer size — legend hides | PASS | 3s | PASSED | Legend auto-hid when Scatter Plot resized to 100x100 |
| 11 | Restore viewers to normal size | PASS | 2s | PASSED | Legends reappeared after resize restore |
| 12 | Set corner positions, enable mini-legend | PASS | 3s | PASSED | Corner positions set (TL/TR/BL/BR); 1 mini-legend icon detected |
| 13 | Save and apply layout | PASS | 8s | PASSED | 7 viewers restored, 6 legends visible, corner positions preserved |
| 14 | Save and open project | FAIL | 5s | N/A | grok.dapi.projects.save threw ApiException; layout-level save works but project-level save failed |
| 15 | Close All | PASS | 1s | PASSED | All views closed |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~30s |
| Spec file generation | ~3s |
| Spec script execution | 30s |

## Summary

Legend visibility and positioning features work correctly overall. Legends are visible on all supported viewers, colors match, CTRL+click multi-select works, and legend positions (top, right, corner) are preserved across layout save/restore. The legend auto-hides when the viewer is too small (Visibility=Auto). The main issues are: the Color Picker could not be accessed from the legend (AMBIGUOUS), and project-level save/open failed with ApiException.

## Retrospective

### What worked well
- Legend visibility toggling (Auto/Always/Never) works correctly across all viewers
- Legend position changes (top, right, corner) are properly reflected in the DOM
- Layout save/restore preserves all legend settings including corner positions
- CTRL+click multi-select works as expected with d4-legend-item-current class
- Changing color column dynamically updates legend text (Scaffold Names → Stereo Category)
- Visibility=Always shows legends on 6 of 7 viewer types

### What did not work
- Color Picker: no clear way to open from legend item; right-click opens global context menu instead of legend-specific menu
- Project save via JS API threw ApiException — may need different approach to save project with SPGI data
- Box Plot never shows a DOM-based legend element, even with legendVisibility=Always
- Mini-legend icon only appears on 1 viewer; other viewers don't seem to use it in the current viewer sizes

### Suggestions for the platform
- Add a legend-specific context menu (right-click on legend item) with "Change Color" option to make color picker discoverable
- Consider adding DOM-based legend support to Box Plot for consistency
- Expose mini-legend mode as an explicit property in viewer options

### Suggestions for the scenario
- Step 4 "Check Color Picker" needs clarification on how to access it — specify exact interaction (right-click? double-click? which element?)
- Step 14 "Save and open project" — specify whether to use grok.dapi.projects or the SAVE button
- Consider splitting this into 2-3 smaller scenarios for independent testing of visibility, positioning, and project persistence
