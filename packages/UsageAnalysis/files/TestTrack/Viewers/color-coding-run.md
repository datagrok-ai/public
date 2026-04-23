# Color Coding — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Notes |
|---|------|--------|------|-------|
| 1 | Open demog | PASS | 4s | 5850 rows, 11 cols |
| 2.1 | AGE: linear then conditional | PASS | 1s | `setLinear()` → "Linear"; `setConditional({'< 30': color, ...})` → "Conditional" |
| 2.2 | SEX: categorical with custom M/F colors | PASS | 1s | `setCategorical({'M': color, 'F': color})` → "Categorical" |
| 2.3 | CONTROL: categorical (default) | PASS | 1s | `setCategorical()` → "Categorical" |
| 2.4 | STARTED: linear with custom 3-stop scheme | PASS | 1s | `setLinear([c1, c2, c3])` → "Linear" |
| 3.1 | Disable AGE, SEX, STARTED | PASS | 1s | `setDisabled()` → "Off" for all three |
| 3.2 | Re-enable; custom colors preserved | PASS | 1s | Types restored; non-type tags unchanged (tagsPreserved=true) |
| 4.1 | Create Race_copy, apply categorical | PASS | 1s | 12 cols total; Race_copy type="Categorical" |
| 4.2 | Apply RACE coloring to Race_copy (Pick Up → Apply) | PASS | 1s | Tag copy via `.color-coding-*`; dstType === srcType |
| 4.3 | Apply STARTED coloring to HEIGHT | PASS | 1s | Tag copy; HEIGHT type matches STARTED |
| 5.1 | AGE: set custom 3-stop linear scheme | PASS | 1s | `setLinear([#1A237E, #F5F5F5, #B71C1C])` → "Linear" |
| 5.2 | AGE: invert the scheme | PASS | 1s | Tag `.color-coding-linear` reversed; first/last colors changed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec execution | 29s |

## API findings

| API call | Works | Notes |
|----------|-------|-------|
| `col.meta.colors.setLinear()` | ✅ | No args = default scheme |
| `col.meta.colors.setLinear([c1, c2, c3])` | ✅ | Color stop array accepted |
| `col.meta.colors.setConditional({'< 30': color, ...})` | ✅ | Object with range strings as keys |
| `col.meta.colors.setCategorical()` | ✅ | Default category colors |
| `col.meta.colors.setCategorical({'M': color, 'F': color})` | ✅ | Custom per-category color map |
| `col.meta.colors.setDisabled()` | ✅ | Sets type to "Off", preserves other tags |
| `col.meta.colors.getType()` | ✅ | Returns "Linear" / "Categorical" / "Conditional" / "Off" |
| Tag `.color-coding-linear` | ✅ | Stores linear color stops as JSON integer array |
| Tag `.color-coding-conditional` | ✅ | Stores conditional rules as JSON object with hex colors |
| Pick Up / Apply via tag copy | ✅ | Copy all `.color-coding-*` tags from src to dst column |

## Summary

All 11 steps passed. The entire test runs on the demog dataset (no SPGI_v2 needed).
UI-only steps (Grid Color Coding All/None/Auto, layout save/restore, Edit scheme dialog) are in `color-coding-ui.md`.
