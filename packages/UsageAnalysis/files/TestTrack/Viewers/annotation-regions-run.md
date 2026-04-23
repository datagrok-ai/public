# Annotation regions — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai/
**Status**: PASS (45.5s)

## Steps

| # | Step | Result | Time | Notes |
|---|------|--------|------|-------|
| 1.1 | Scatter Plot, lassoTool=false, Draw Annotation Region, drag rect, edit 7 fields (Title, Description, Fill, Outline, OutlineWidth, Opacity, HeaderColor), OK | PASS | — | outlineWidth=3 and opacity=60 verified via `getOptions`; OutlineWidth/Opacity are raw `<input type="range">` set via evaluate |
| 1.2 | lassoTool=true, Draw Annotation Region, drag polygon | PASS | — | Dialog auto-opens after polygon drag; cancelled |
| 1.3 | Formula Lines... → ADD NEW → Region - Formula Lines | PASS | — | Formula pattern matches `${HEIGHT}.*${WEIGHT}` |
| 2.1 | Toggle Show Viewer / Show Dataframe Annotation Regions individually | PASS | — | Each flag flips independently |
| 2.2 | Global Show Annotation Regions re-enables both | PASS | — | Both flags → true |
| 3.1 | Hover region / intersection | SKIP | — | Manual only — see annotation-regions-ui.md (no worldToScreen API) |
| 3.2 | Click region → rows selected in grid | SKIP | — | Manual only — see annotation-regions-ui.md (canvas hit-test unreliable without worldToScreen) |
| 4.1 | Right-click region → Edit opens dialog | PASS* | — | "Edit" menu item found: true — cursor hit the region |
| 4.2 | Reopen Formula Lines dialog, click grid row, edit OutlineWidth/Opacity/HeaderColor, OK | PASS | — | outlineWidth=5, opacity=40 verified via `getOptions` |
| 4.3 | Change viewer-level Annotation Font | PASS | — | `annotationFont` → "18px bold Arial" |
| 5.1 | Formula Lines dialog: preview canvas renders | PASS | — | Any dialog canvas has non-zero width/height |
| 5.2 | Formula Lines dialog: grid visible | PASS | — | `.d4-grid` toBeVisible — column names not checkable (canvas-based) |
| 6.1 | Line Chart multiAxis=true — Draw Annotation Region absent | PASS | — | menu item count = 0 |
| 6.2 | Line Chart single-axis — draw rectangle + add formula region | PASS | — | annotationRegions.length increases after rect, then after formula |

## Launch

```bash
npx playwright test --config=playwright.config.ts --project=dev "annotation-regions-spec"
```

## Known limitations

- **3.2 / 4.1**: Canvas hit-test relies on approximate screen coordinates. Wide region in 3.2 maximizes hit probability; 4.1 is soft and logs AMBIGUOUS if cursor misses.
- **5.2**: Grid column names (Title / Formula / Show) are canvas-rendered — not accessible via DOM.
- **3.1**: Hover interaction requires manual verification — no reliable way to assert `mouseOverRowGroup` without `worldToScreen`.

## Startup pattern

Uses Playwright `{page}` fixture with standard login form detection:
```ts
await page.goto(baseUrl);
const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) { ... }
await page.locator('[name="Browse"]').waitFor({timeout: 120_000});
```

## Suggestions for the platform

- Expose `viewer.worldToScreen(dataX, dataY)` returning reliable screen coords — would unlock deterministic hover/click region tests.
- Give Annotation Regions visibility toggles stable `[name="input-host-…"]` selectors in the property panel so automation can reach them without the context menu fallback.
- Make dialog preview canvases report non-zero size synchronously to allow `toBeVisible()` without custom polling.
- **5.2**: Reconcile scenario description — dialog grid shows columns "Title / Formula / Show", not "type, header, fill/outline, formulas, description".
- **6.1**: Clarify expectation — `Draw Annotation Region` is **absent** from the menu (not just disabled), while `Formula Lines...` remains available.
