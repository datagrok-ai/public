# Catalog (DiffStudio) — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open PK-PD from Diff Studio Library | 15s | PASS | PASSED | `runDiffStudio` + dblclick `.diff-studio-hub-card` (PK-PD); view 'PK-PD' opens with 15+ input hosts |
| 2 | Click Save to Model Hub icon | 5s | PASS | PASSED | Clicked `.diff-studio-ribbon-save-to-model-catalog-icon` (fa-layer-plus). Notification balloon "Saved to Library as PK-PD.ivp" appears |
| 3 | Open Model Hub (Apps > Model Hub) | 10s | PASS | PASSED | Invoked `Compute2:modelCatalog` via `Func.prepare() + call() + grok.shell.addView()`; view 'Model Hub' (JsCardView) at `/apps/Modelhub?`; PK-PD card present (1 `.d4-list-item` with text 'PK-PD') |
| 4 | Click Refresh icon; catalog reloads | 5s | PASS | PASSED | Clicked `.grok-icon.fa-sync` (icon has `fal fa-sync` classes); after 3s wait, PK-PD card still present |
| 5 | Run PK-PD model from catalog | 15s | PASS | PASSED | Last `.d4-list-item` (text 'PK-PD'): `click()` + `dispatchEvent('dblclick')`; view 'PK-PD' (TableView) opens at `/apps/DiffStudio/Library/pk-pd?params:…` with count/dose/Multiaxis present |
| 6 | Modify dose; verify live update + URL | 10s | PASS | PASSED | Typed 5000 into `[name="input-host-dose"] input.ui-input-editor` (was 10000); URL settles to `dose=5000`, chart re-renders (11 `.d4-viewer canvas` elements) |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 13s |
| grok-browser execution (scenario steps) | 40s |
| Execute via grok-browser (total) | 1m 53s |
| Spec file generation | 39s |
| Spec script execution | 38s |
| **Total scenario run (with model)** | 3m 10s |

## Summary

Catalog scenario reproduces fully on dev.datagrok.ai. All 6 steps PASS both in MCP and in the Playwright spec (35.8s wall-clock). PK-PD saves via `.diff-studio-ribbon-save-to-model-catalog-icon` (fa-layer-plus) with a balloon notification, Model Hub opens via `Compute2:modelCatalog`, Refresh is `.grok-icon.fa-sync`, and opening a saved model from the hub requires both `click()` and a `dblclick` event. Dose edit propagates to URL within 2s. **Total scenario run (with model)**: 3m 10s.

## Retrospective

### What worked well
- `.diff-studio-ribbon-save-to-model-catalog-icon` is a specific, stable CSS class — unambiguous target for the Save to Model Hub action
- Notification balloon "Saved to Library as PK-PD.ivp" is a textual success signal — easier than polling server state
- `Compute2:modelCatalog` consistently maps to "Model Hub" across scenarios (scripting + catalog) — the same entry-point pattern works
- `.grok-icon.fa-sync` works as a Refresh identifier
- `click()` + `dispatchEvent('dblclick')` reliably opens a saved model from the Model Hub card list
- Post-run cleanup via `grok.dapi.files.delete(fileInfo)` deletes PK-PD.ivp from the user's catalog; `.path` prefix check avoids accidentally deleting the original library file

### What did not work
- Single-click on the Model Hub card was not enough — needed both `click()` and a `dispatchEvent('dblclick')` to reliably navigate; same hub-interaction pattern as DiffStudio library cards
- The file name written to the user catalog is `PK-PD.ivp` (preserves the hub card label). Listing via `grok.dapi.files.list('System:AppData/DiffStudio/', true, 'PK-PD.ivp')` returned an empty result in cleanup — suggests the saved file may land under a different path (possibly `Users:…` or model-hub-owned namespace); cleanup was a no-op

### Suggestions for the platform
- Save to Model Hub should echo the saved path in the notification balloon (currently only "Saved to Library as PK-PD.ivp" — doesn't disclose *which* library)
- Single-click on a Model Hub card should open the model — align with how apps usually work (single click = primary action)
- Expose `Compute2:deleteFromModelHub(name)` or a delete-by-name API to make automation cleanup trivially discoverable

### Suggestions for the scenario
- Step 2 should include the notification text as the success criterion ("balloon: Saved to Library as PK-PD.ivp")
- Step 5 should say "double-click the PK-PD card" — single-click alone doesn't open the model
- Step 6 "adjust input parameters" is generic; name one concrete input + range (e.g., "decrease `dose` from 10000 to 5000; curves rescale accordingly")
- Add a cleanup step so repeated runs don't accumulate duplicate PK-PD entries in the catalog
