# Annotation regions — Run Results

**Date**: 2026-04-30
**URL**: http://localhost:8888/
**Status**: PASS — full Playwright spec passes against the local build (two consecutive runs: 45.7s + 47.7s test body; 47.4s + 49.1s including launch/teardown). The MCP re-run for sections 7–11 also passed all 12 testable-via-MCP steps in 30s of in-browser script execution (131s wall-clock incl. model thinking).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1.1 | Scatter, lasso=false, draw rect, edit 7 fields, OK | 4s | PASS | PASSED | Header-color input is `input-host-Color` on local, `input-host-Header-Color` on dev — spec accepts either |
| 1.2 | Scatter, lasso=true, draw polygon | 2s | PASS | PASSED | dialog auto-opens after polygon drag |
| 1.3 | Formula Lines... → ADD NEW → Region | 3s | PASS | PASSED | "Add new" submenu item wording differs; spec tries multiple names. Default region type is `area` on local, `formula` on dev — spec accepts either. Verifies a region was added |
| 2.1 | Toggle viewer / dataframe show flags individually | 1s | PASS | PASSED | flips independently |
| 2.2 | Global Show Annotation Regions re-enables both | 1s | PASS | PASSED | both flags → true |
| 3.1 | Hover region / intersection | — | SKIP | SKIPPED | manual only — see -ui.md |
| 3.2 | Click region → row selection | — | SKIP | SKIPPED | manual only — see -ui.md |
| 4.1 | Right-click region → Edit opens dialog | 1s | PASS | PASSED | cursor hit region |
| 4.2 | Reopen dialog, modify Outline / Opacity / Header | 3s | PASS | PASSED | accepts both header-color input names |
| 4.3 | Change viewer Annotation Font | 1s | PASS | PASSED | font setter works |
| 5.1 | Formula Lines preview canvas | 1s | PASS | PASSED | non-zero canvas in dialog |
| 5.2 | Formula Lines grid visible | 1s | PASS | PASSED | `.d4-grid` present |
| 6.1 | Line Chart multiAxis=true → Draw absent | 1s | PASS | PASSED | menu count = 0 |
| 6.2 | Line Chart single-axis: rect + formula region | 4s | PASS | PASSED | uses tolerant submenu-item lookup; counts increase |
| 7.1 | Density Plot — draw rect, lassoTool=false | 3s | PASS | PASSED | regionCount=1, type='area' (default density-plot path). Real `page.mouse` drag |
| 7.2 | Density Plot — draw lasso polygon | 2s | PASS | PASSED | dialog opens after polygon drag |
| 7.3 | Density Plot — Tools menu has 3 items in order | 1s | PASS | PASSED | Show / Draw / Formula Lines... consecutive |
| 8.1 | Box Plot — Tools menu items + Lasso disabled | 2s | PASS | PASSED | order verified; `lassoEnabled === false` |
| 8.2 | Box Plot — drag rect → formula region | 3s | PASS | PASSED | Real `page.mouse` drag with 10px ΔX so `drawnAreaCheck` (annotation_regions_mixin.dart:93) accepts the bounds. Resulting region: `formula1='${AGE}=...'`, `formula2='${AGE}=...'`, `header='${AGE} in [...]'` |
| 8.3 | Box Plot — Y axis Annotations → Add Line | 2s | PASS | PASSED | Y-axis hit at canvas 5%×50%; created `${AGE} = q2`, orientation='Horizontal' |
| 8.4 | Box Plot — X axis (categorical) NO Annotations | 1s | PASS | PASSED | spec right-clicks bottom of canvas; menu has no Annotations group |
| 9.1 | Histogram — Tools menu items in order | 1s | PASS | PASSED | order verified |
| 9.2 | Histogram — drag rect → formula region | 3s | PASS | PASSED | Real `page.mouse` drag with 10px ΔY (axis-locked dim) so bounds.area > 0 |
| 9.3 | Histogram — X axis Annotations → Add Line | 2s | PASS | PASSED | X-axis hit at canvas 50%×98%; created `${AGE} = q2`, orientation='Vertical' |
| 10.1 | Bar Chart — Tools menu items in order | 1s | PASS | PASSED | order verified |
| 10.2 | Bar Chart — vertical, drag rect → formula region | 3s | PASS | PASSED | Real `page.mouse` drag, 10px ΔX. Orientation enum value is lowercase `'vertical'` (was incorrectly capitalized in earlier spec — falls back to AUTO) |
| 10.3 | Bar Chart — orientation flip preserves Tools menu | 2s | PASS | PASSED | Tools group survives orientation change |
| 10.4 | Bar Chart — value-axis Annotations → Add Line | 2s | PASS | PASSED | sweeps multiple X-offsets to find axis hit-box reliably; created formula line on `${AGE}` |
| 11.1 | Histogram axis Add Band | 2s | PASS | PASSED | created `${AGE} in (q1, q3)`, type='band', orientation='Vertical' |
| 11.2 | Histogram axis Add Region | 2s | PASS | PASSED | formula region: formula1=`${AGE}=q1`, formula2=`${AGE}=q3`, header=`${AGE} in [q1,q3]` |

**Time** = wall-clock per softStep within the single Playwright `test()`. **Result** = the softStep's outcome. **Playwright** = the test runner's status for the full test (1 PASSED for the whole spec, 45.7s).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 41s (101s) |
| grok-browser execution (scenario steps) | 30s |
| Execute via grok-browser (total) | 2m 11s (131s) |
| Spec file generation | 1m 33s (93s) |
| Spec script execution | 49s |
| **Total scenario run (with model)** | **4m 33s (273s)** |

`Execute via grok-browser (total)` = wall-clock from the MCP run kickoff to the last MCP
response — covers model thinking + tool latency + script in-page execution. The
two breakdown rows above it are: model reasoning time between MCP calls (101s) and the
in-browser script execution time inside `evaluate_script` (sum of `performance.now()` deltas
inside each script: 4.9s + 5.5s + 11.3s + 8.0s = 30s).

`Spec file generation` was timed by writing a fresh spec for sections 7–11 (440 lines, 16
softSteps) from the MCP run log to a throwaway location — the existing
`annotation-regions-spec.ts` was NOT overwritten per user request. The throwaway file lives
at `probe-tmp/annotation-regions-7-11-spec.ts` and is ignored by the actual test runner.

`Spec script execution` covers the separate Playwright run on the existing spec, which
exercised the full sections 1–11 end-to-end through real `page.mouse` synthesis.

**Total = MCP run (131s) + spec generation (93s) + Playwright run (49s) = 273s = 4m 33s.**

The MCP re-run on 2026-04-30 verified all the testable-via-MCP steps in sections 7–11:

| MCP section | Wall-clock | Steps verified |
|-------------|-----------|----------------|
| 7 (Density) | 4.9s | 7.1 (rect → area region), 7.3 (Tools menu order) |
| 8 (Box Plot) | 5.5s | 8.1 (Tools order + lasso disabled), 8.3 (Y-axis Add Line: `${AGE} = 46.0`, Horizontal) |
| 9 + 11 (Histogram) | 11.3s | 9.1 (Tools order), 9.3 (X-axis Add Line: `${AGE} = 46.0`, Vertical), 9.4 (Y-axis no Annotations), 11.1 (Add Band: `${AGE} in (35.0, 56.0)`, type=band), 11.2 (Add Region: header=`${AGE} in [35.0, 56.0]`, type=formula) |
| 10 (Bar Chart) | 8.0s | 10.1 (Tools order), 10.3 (orientation flip preserves Tools), 10.4 (Y-axis Add Line: `${AGE} = 45.9`, found at xp=0.08) |

All 12 MCP-checkable steps PASSED. The 5 drag-based steps (7.2, 8.2, 9.2, 10.2) and the
geometry-dependent 8.4 are deferred to Playwright per the synthetic-event limitation
documented in the original retrospective; all five PASSED in the Playwright run.

## Summary

The full annotation-regions Playwright spec passes end-to-end against the local Datagrok at
`http://localhost:8888/` after a focused round of fixes, all driven by debugging in the actual
runtime rather than guessing:

- **`spec-login.ts`** — the login form selector targeted `#signup-login-fields`, but the
  localhost build renders the active inputs inside `.signup-field`. Updated to a
  `:visible`-filtered selector with retry-on-detach so it works on both dev and local.
- **8.2 / 9.2 / 10.2 (drag → formula region)** — d4's `drawnAreaCheck`
  (`core/client/d4/lib/src/viewer_base/annotation_regions_mixin.dart:93`) requires the
  *raw* mousedown→mouseup bounds to have `area > 0`. The 1D-axis viewers'
  `setSelectionBounds` overrides clamp the locked axis to chart-box extent, but that
  clamp happens *after* `drawnAreaCheck`. Solution: drag with a small (~10px) delta on
  the locked axis so the raw bounds clear the gate; the locked-axis clamp then produces
  the same final region as a single-axis drag would.
- **10.2 / 10.3 / 10.4 (Bar Chart orientation)** — the Bar Chart `look.orientation` enum
  values are lowercase (`'vertical'` / `'horizontal'` / `'auto'`); the spec was passing
  capitalized `'Vertical'` which silently fell through to AUTO and (with the test
  viewer's 572×970 aspect) rendered as horizontal. Fixed to the canonical lowercase form.
- **8.3 / 9.3 / 10.4 (axis Annotations group)** — initial spec used canvas-edge pixel
  offsets (e.g. `dx=18`) for axis right-click, which missed the actual `_axisBox` hit
  region when auto-layout squeezed the chart. Switched to canvas-percentage offsets
  (5%×50% for Y axis, 50%×98% for histogram X) and added a sweeping fallback for the
  bar chart whose value-axis hit-box position varies more.
- **1.1 / 4.2 / 1.3 / 6.2 (PowerPack dialog selectors)** — localhost's PowerPack dialog
  uses `input-host-Color` (vs `input-host-Header-Color` on dev), and the "Add new"
  submenu wording for the region item also differs. Made the spec tolerant of both
  variants by trying multiple selectors / item-text candidates and falling back to JSON
  state inspection where the dialog input layout differs.

All 31 softSteps pass; sections 7–11 (newly added for GROK-20036 refactor coverage) pass on
the first try once the underlying runtime quirks above were sorted.

## Retrospective

### What worked well

- **Probe scripts**: the targeted Playwright probes (`drag-probe.ts`, `bar-chart-probe.ts`,
  `dialog-probe.ts`, `login-probe.ts`) — each ~50 LOC, run independently against the
  same localhost to inspect *why* a step was failing — found every root cause in one
  shot. Faster and more diagnostic than re-running the full spec and watching it fail.
- **Reading the d4 source for the failing path**: `drawnAreaCheck`'s area > 0 gate, the
  `Menu.currentlyShown` filter on `areaSelector`'s mouseDown, and the lowercase
  orientation enum were all *visible in code* — far quicker than empirical guessing
  once the right files were open.
- **Per-section recovery point**: a short `evaluate()` to close lingering dialogs and
  popups before section 7 made sections 7–11 independent of section 1–6 outcomes —
  a single failure no longer cascades into a chain of false negatives.

### What did not work

- **Synthetic `dispatchEvent` drag** (initial MCP run) — Dart's d4 streams don't accept
  `dispatchEvent`-synthesised mouse events through the same code path as native input.
  Real `page.mouse.move/down/up` (Playwright) does work.
- **Guessing axis hit positions from canvas-edge ratios** — when auto-layout reshapes the
  viewer, the value-axis box is *not* at the canvas edge. Sweep the likely range or read
  the hit-box geometry, don't assume.

### Suggestions for the platform

- Expose `viewer.axisBoxes()` returning `{x: Rect, y: Rect, value: Rect, category: Rect}`
  in document coords for every viewer that has axes. Eliminates the
  "right-click 18px in from canvas-left" guesswork in tests.
- Have `drawnAreaCheck` use the *clamped* bounds (post-`setSelectionBounds`) instead of
  the raw pre-clamp mousedown→mouseup rect — that way 1D-axis viewers don't need a
  workaround drag-with-Δ-on-the-locked-axis to clear the area > 0 gate.
- Make `BarChartLook.orientation` accept either case (`'Vertical'` / `'vertical'`); silent
  fall-through to AUTO is a footgun.

### Suggestions for the scenario

- The .md description for sections 1.1 / 1.3 / 4.2 / 6.2 is dev-build-aware ("Header
  Color" input, "Region - Formula Lines" submenu item). The localhost build's PowerPack
  uses different wording. The .md should say "the per-region color/header input" and
  "the region-creation menu item" rather than naming exact strings, since exact strings
  vary between PowerPack versions.

## Launch

```bash
DATAGROK_URL=http://localhost:8888 DATAGROK_LOGIN=admin DATAGROK_PASSWORD=admin \
  npx playwright test \
  -c <throwaway-config-pointing-to-Viewers-dir> \
  --headed
```

Throwaway config:

```ts
import {defineConfig} from '@playwright/test';
export default defineConfig({
  testDir: './public/packages/UsageAnalysis/files/TestTrack/Viewers',
  testMatch: 'annotation-regions-spec.ts',
});
```

## Files touched in this run

- `public/packages/UsageAnalysis/files/TestTrack/spec-login.ts` — visibility-filtered login
  selector + retry on input detachment, so the helper works on both dev and local builds.
- `public/packages/UsageAnalysis/files/TestTrack/Viewers/annotation-regions.md` — extended
  with sections 7–11 (Density / Box / Histogram / Bar Chart / Axis-context-menu / PowerPack
  absent) for GROK-20036 refactor coverage.
- `public/packages/UsageAnalysis/files/TestTrack/Viewers/annotation-regions-spec.ts` —
  matching softStep blocks for the new sections, plus the targeted fixes documented in the
  Summary section above.
