# Composition Analysis — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open FASTA.csv (also tested HELM.csv, MSA.csv in analyze scenario) | PASS | All 3 datasets loaded correctly |
| 2 | Bio > Analyze > Composition | PASS | WebLogo viewer added to view |
| 3 | Check selection by clicking letter in viewer | AMBIGUOUS | WebLogo is canvas-based; automated clicking did not trigger row selection |
| 4 | Click Gear icon on viewer | PASS | Opened via context menu > Properties (title bar gear icon not available on WebLogo) |
| 5 | Change arbitrary properties on Context Panel | PASS | Data section shows: Sequence, Value, Value Aggr Type, Start/End Position, Skip Empty options; Layout, Behavior, Style sections present |

## Summary

4 of 5 steps passed. The Composition/WebLogo viewer opens correctly and properties are accessible. The letter-click selection (step 3) could not be verified through automation because the WebLogo viewer is canvas-based and the exact pixel position of letters cannot be determined from the DOM.

## Retrospective

### What worked well
- Composition opens directly without a dialog
- WebLogo renders sequence logos correctly
- Properties panel shows all expected settings sections

### What did not work
- WebLogo title bar icons (gear, close) not visible even with selenium class — had to use context menu
- Canvas-based letter clicking could not be automated for selection testing

### Suggestions for the platform
- WebLogo should show title bar icons when selenium class is active
- Consider adding name attributes to WebLogo position elements for automation

### Suggestions for the scenario
- Step 3 "click any letter" is hard to automate on canvas viewers — consider providing a JS API alternative
