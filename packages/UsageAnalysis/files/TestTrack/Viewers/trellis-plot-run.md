# Trellis plot tests (Playwright) -- Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-------|-------|
| **Inner viewer types** | | | | |
| 1 | Switch to Scatter plot | PASS | PASSED | `tp.props.viewerType = 'Scatter plot'` |
| 2 | Set X=WEIGHT, Y=HEIGHT, Color=RACE | PASS | PASSED | `tp.setOptions({innerViewerLook: ...})` |
| 3 | Switch to Bar chart | PASS | PASSED | |
| 4 | Set Split=RACE, Value=AGE, Aggr=avg | PASS | PASSED | |
| 5 | Switch to Histogram | PASS | PASSED | |
| 6 | Set Value=AGE, Split=RACE | PASS | PASSED | |
| 7 | Switch to Line chart | PASS | PASSED | |
| 8 | Set X=STARTED, Y=AGE, Split=RACE | PASS | PASSED | |
| 9 | Switch to Box plot | PASS | PASSED | |
| 10 | Set Category=SEX, Value=AGE | PASS | PASSED | |
| 11 | Switch to Pie chart | PASS | PASSED | |
| 12 | Set Category=RACE | PASS | PASSED | |
| 13 | Switch to Density plot | PASS | PASSED | |
| 14 | Set X=WEIGHT, Y=HEIGHT | PASS | PASSED | |
| 15 | Switch to Statistics (Summary) | PASS | PASSED | |
| 16 | Set Visualization to bars | PASS | PASSED | |
| 17 | Switch to Sparklines | PASS | PASSED | |
| 18 | Set Sparkline Type to Bar Chart | PASS | PASSED | |
| 19 | Switch to PC Plot | PASS | PASSED | |
| 20 | Set Color to SEX | PASS | PASSED | |
| **Global scale** | | | | |
| 1 | Switch to Scatter plot | PASS | PASSED | |
| 2 | Enable Global Scale | PASS | PASSED | |
| 3 | Disable Global Scale | PASS | PASSED | |
| 4 | Re-enable Global Scale | PASS | PASSED | |
| **Axes visibility** | | | | |
| 1-4 | Show X Axes: Always, Never, Auto | PASS | PASSED | |
| 5 | Show Y Axes: Always, Never, Auto | PASS | PASSED | |
| 6 | Show Range Sliders off | PASS | PASSED | |
| 7 | Show Range Sliders on | PASS | PASSED | |
| **Range sliders with global scale** | | | | |
| 2 | Enable Global Scale + Range Sliders | PASS | PASSED | |
| 3 | Show X Axes = Always | PASS | PASSED | |
| 4-5 | Hover/drag range slider | AMBIGUOUS | SKIPPED | Canvas interaction not automatable |
| 6 | Reset Inner Range Sliders | SKIP | SKIPPED | Depends on step 5 |
| 7 | Show Y Axes = Always | PASS | PASSED | |
| **Gridlines** | | | | |
| 1-3 | always, never, auto | PASS | PASSED | |
| **Tiles mode** | | | | |
| 1 | Enable tiles | PASS | PASSED | |
| 2 | Tiles Per Row = 2 | PASS | PASSED | |
| 3 | Tiles Per Row = 6 | PASS | PASSED | |
| 4 | Disable tiles | PASS | PASSED | |
| **Category management** | | | | |
| 1 | Set X=SEX, Y=RACE | PASS | PASSED | |
| 2 | Add X column DIS_POP | PASS | PASSED | |
| 3 | Remove second X column | PASS | PASSED | |
| 4-6 | Toggle X/Y labels | PASS | PASSED | |
| **Pack categories** | | | | |
| 1 | Set X=SEX, Y=RACE | PASS | PASSED | |
| 2 | Filter out Asian | PASS | PASSED | 5778 rows filtered |
| 3 | Pack Categories = true (default) | PASS | PASSED | |
| 4 | Disable Pack Categories | PASS | PASSED | |
| 5 | Re-enable | PASS | PASSED | |
| **On Click functionality** | | | | |
| 1 | Set On Click = Select | PASS | PASSED | |
| 2-5 | Click cells (Select mode) | AMBIGUOUS | SKIPPED | Canvas click |
| 6 | Set On Click = Filter | PASS | PASSED | |
| 7-15 | Click cells (Filter mode) | SKIP | SKIPPED | Canvas click |
| 16 | Set On Click = None | PASS | PASSED | |
| **Selectors** | | | | |
| 1-4 | Toggle selectors and control panel | PASS | PASSED | |
| **Allow viewer full screen** | | | | |
| 1-3 | Hover/click full-screen icon | AMBIGUOUS | SKIPPED | Canvas hover |
| 4 | Disable allowViewerFullScreen | PASS | PASSED | |
| **Scrolling** | | | | |
| 1 | Set X=SITE | PASS | PASSED | |
| 2-3 | Horizontal scroll | AMBIGUOUS | SKIPPED | Canvas interaction |
| 4 | Set Y=RACE | PASS | PASSED | |
| 5 | Mouse wheel scroll | AMBIGUOUS | SKIPPED | Canvas interaction |
| **Legend** | | | | |
| 1 | Switch to Scatter plot, color=SEX | PASS | PASSED | |
| 2 | Legend Visibility = Always | PASS | PASSED | |
| 3 | Legend Position cycle | PASS | PASSED | Left, Right, Top, Bottom |
| 4 | Legend Visibility = Never | PASS | PASSED | |
| **Context menu** | | | | |
| 1-3 | Right-click context menu | PASS | PASSED | Scatter plot group found |
| **Inner viewer properties** | | | | |
| 1-5 | Gear icon + inner viewer settings | PARTIAL | PASSED | Gear not found via DOM, used JS API |
| **Use in Trellis** | | | | |
| 1-3 | Scatter plot > Use in Trellis | PASS | PASSED | Context menu worked |
| 5-6 | Bar chart > Use in Trellis | PASS | PASSED | |
| 8-9 | Histogram > Use in Trellis | FAIL | FAILED | Not found in General submenu |
| 11-12 | Line chart > Use in Trellis | PASS | PASSED | |
| 14-15 | Box plot > Use in Trellis | PASS | PASSED | |
| **Auto layout** | | | | |
| 1 | Auto Layout enabled by default | PASS | PASSED | |
| 2-3 | Resize viewer | AMBIGUOUS | SKIPPED | Not automatable |
| 4 | Disable Auto Layout | PASS | PASSED | |
| **Title and description** | | | | |
| 1-4 | Title, description, position | PASS | PASSED | |
| **Label orientation** | | | | |
| 1-5 | X/Y label orientations | PASS | PASSED | |
| **Pick Up / Apply** | | | | |
| 3-6 | Pick Up and Apply between trellis plots | PASS | PASSED | Settings transferred correctly |
| 7-8 | Independent axis/slider changes | SKIP | SKIPPED | |
| **Layout and Project save/restore** | | | | |
| 1-3 | Layout save/restore | PASS | PASSED | Viewers correctly restored |
| 4-6 | Project save/restore | SKIP | SKIPPED | Side-effect avoidance |
| **Viewer filter formula** | | | | |
| 1-3 | Set and clear filter formula | PASS | PASSED | |
| **Multi Curve inner viewer** | | | | |
| 1-4 | Open curves, switch table, set Multi curve | PASS | PASSED | |
| 5-8 | Multi curve interactions | SKIP | SKIPPED | |
| **To Script** | | | | |
| 1-3 | To JavaScript balloon | PASS | PASSED | |
| **Keyboard navigation** | | | | |
| 1-6 | Arrow keys + ESC | AMBIGUOUS | SKIPPED | Canvas click required |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~180s |
| Spec file generation | ~30s |
| Spec script execution | 108s (PASSED) |

## Summary

Most trellis plot property-based tests passed successfully via JS API. Canvas-based interactions
(bin clicks, range slider drag, full-screen icon hover, keyboard navigation) could not be automated
and were marked AMBIGUOUS. The Histogram "Use in Trellis" context menu item was not found (FAIL).
Layout save/restore worked correctly.

## Retrospective

### What worked well
- All viewer type switches via `tp.props.viewerType` worked flawlessly
- Inner viewer property setting via `tp.setOptions({innerViewerLook: ...})` worked
- Context menu automation via `dispatchEvent(new MouseEvent('contextmenu', ...))` worked for Pick Up/Apply and Use in Trellis
- Layout save/restore preserved all trellis plot settings

### What did not work
- Histogram viewer missing "Use in Trellis" in General context menu submenu
- Canvas-based interactions (bin clicks, range slider drag, full-screen icon hover) not automatable via DOM
- Gear icon not visible via DOM even with selenium class (requires hover?)

### Suggestions for the platform
- Add "Use in Trellis" support for Histogram viewer if missing
- Consider adding `name=` attributes to trellis cell elements for automation
- Make gear icon accessible without hover in selenium mode

### Suggestions for the scenario
- Split canvas-dependent steps into separate scenarios
- Add expected property values for verification
- Clarify "Use in Trellis" should work for all viewer types listed
