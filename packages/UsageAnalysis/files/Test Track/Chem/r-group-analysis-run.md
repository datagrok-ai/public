# R Group Analysis — Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open sar_small.csv dataset | PASS | 8s | PASSED | 200 rows, 27 columns loaded |
| 2 | Open Chem > Analyze > R-Groups Analysis | PASS | 2s | PASSED | Dialog opened with sketcher, MCS, options |
| 3 | Click MCS button | PASS | 5s | PASSED | Core scaffold computed and displayed in sketcher |
| 4 | Check Visual analysis checkbox | PASS | 1s | PASSED | Already checked by default |
| 5 | Click OK — run with defaults | PASS | 10s | PASSED | Trellis plot appeared, columns increased to 34 (7 R-group cols added) |
| 6 | Run again, click MCS, uncheck Replace latest | PASS | 8s | PASSED | Columns increased to 41 — second set of R-group columns added |
| 7 | Run again, click MCS, check Replace latest | PASS | 8s | PASSED | Columns stay at 41 — latest columns replaced |
| 8 | Run without MCS — expect error | PASS | 3s | PASSED | Balloon "No core was provided" appeared |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~50s |
| Spec file generation | ~3s |
| Spec script execution | 1m 6s (PASSED) |

## Summary

All R-Group Analysis steps passed on the dev server. MCS computation, Visual analysis with trellis plot, "Replace latest" behavior (both checked and unchecked), and error handling for missing core all work correctly.

## Retrospective

### What worked well
- MCS computation is fast and produces correct core scaffold
- Trellis plot with color-coded R-groups displayed correctly
- "Replace latest" toggle works as expected — unchecked adds new columns, checked replaces
- Error handling shows clear balloon message when no core is provided

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- None

### Suggestions for the scenario
- Step numbering has gaps (no step 8, jumps from 7 to 9)
- Could clarify that Visual analysis checkbox is checked by default
