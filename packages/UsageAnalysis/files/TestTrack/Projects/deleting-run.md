# Deleting — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Find the projects from the previous steps | 6s | PASS | PASSED | Setup: created two projects via `DG.Project.create()` + `dapi.projects.save()`. Browse → Dashboards opened gallery; both tiles found via `[name="div-{ProjectName}"]` (two AutoTestDelete tiles visible) |
| 2 | Right-click the project, select **Delete project** | 4s | PASS | FAILED | Real `MouseEvent('contextmenu')` dispatched on the tile opens the menu; `[name="div-Delete-Project"]` present. Confirm dialog opens with title "Are you sure?" and body `Delete project "{name}"?` — `[name="button-DELETE"]`/`[name="button-CANCEL"]` both rendered. Spec replay clicks the menu item but the subsequent DELETE click does not commit (see notes) |
| 3 | Click **DELETE** in the confirmation dialog | 3s | PASS | FAILED | In MCP, `[name="button-DELETE"].click()` closed the dialog and removed the project from `dapi.projects`. In the Playwright spec, the confirm dialog closes but `dapi.projects.filter('name = "X"').first()` continues to return the project (~30s wait); the OK action appears to not fire under token-injected sessions |
| 4 | Check that the project has been deleted | 4s | PASS | FAILED | After both DELETE clicks: tiles disappeared from the Dashboards gallery and `dapi.projects.filter(...).first()` returned null for both projects. Cleanup confirmed no AutoTestDelete projects remain server-side |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~3m |
| grok-browser execution (scenario steps) | ~17s |
| Execute via grok-browser (total) | ~3m 17s |
| Spec file generation | ~2m |
| Spec script execution | ~7m (4 attempts × ~1m–1m 45s) |
| **Total scenario run (with model)** | ~13m |

## Summary

All four scenario steps PASS in the live MCP run on `dev.datagrok.ai`: two test projects were
created via API, located in the Browse → Dashboards gallery, deleted through the right-click
context menu and confirm dialog, and verified gone from both UI and API. **Total scenario run
(with model)** ≈ 13 minutes including spec generation and four Playwright replay attempts.

The Playwright spec replay is partial: it reliably handles login, project setup, gallery
navigation, right-click → context menu, and dialog open. The specific `[name="button-DELETE"]`
click inside the confirm dialog does not propagate the delete to the server when the spec
runs under a token-injected session, even with native `locator.click()`, programmatic
`element.click()`, and full `mousedown`/`mouseup`/`click` event dispatch — the dialog closes
but `dapi.projects` still returns the entity. The same clicks worked instantly in the live
MCP browser session (cookie-based auth), so the platform itself is fine; the gap appears to
be specific to the spec replay environment. Spec retained as `deleting-spec.ts` for future
debugging; cleanup `finally` ensures no test projects are left behind.

## Retrospective

### What worked well
- Right-click context menu: `MouseEvent('contextmenu', {button: 2})` on the tile reliably opens
  the project context menu, with all expected `[name="div-..."]` items present
  (`Open`, `Rename...`, `Delete-Project`, `Save-as-Zip`, `Share...`, `Copy` submenu, `Add-to-favorites`)
- The confirm-dialog DOM contract from `references/projects.md` matches reality:
  title `Are you sure?`, body `Delete project "{friendlyName}"?` plus the warning hint, and
  both `[name="button-DELETE"]` and `[name="button-CANCEL"]` are present
- API path (`grok.dapi.projects.delete`) is reliable and idempotent — used both for setup
  cleanup and as the cleanup-after-spec safety net
- Gallery tile lookup by `[name="div-{ProjectName}"]` works once the gallery view is
  navigated (`/projects` URL or `Dashboards` tree click)

### What did not work
- **Playwright replay of the DELETE click does not commit the deletion.** Tested clicks:
  `locator.click()` (with `expect(...).toBeEnabled()` precheck), inline
  `element.click()` via `page.evaluate`, and a full `mousedown`/`mouseup`/`click`
  `MouseEvent` sequence. In all three: the confirm dialog closes, the
  topmost dialog disappears, but `dapi.projects.filter('name = "X"').first()` keeps
  returning the project (waited 30s after dialog close)
- `playwright.waitForFunction` with an inline `(async () => ...)()` predicate returns the
  outer Promise (truthy) and exits early without awaiting — switched to a manual polling
  loop (`waitProjectGone`) instead
- The Dashboards gallery view caches between renders. Re-clicking the `Dashboards` tree
  node alone is not enough — direct navigation to `/projects` is more reliable for the
  spec's first paint
- Setup using `grok.shell.project` (the scratchpad) and `grok.dapi.projects.save(p)` fails
  with `Unable to add entity ... to the project` because the scratchpad's table is already
  attached to a non-saved state. `DG.Project.create()` + `addChild(tableInfo)` +
  `tables.uploadDataFrame` + `tables.save` + `projects.save` is the working pattern (matches
  `ApiSamples/scripts/dapi/projects.js`)

### Suggestions for the platform
- Investigate why the destructive `[name="button-DELETE"]` OK action in `Modal.confirmDelete`
  closes the dialog without firing the delete handler when the click originates from a
  Playwright-driven page in a token-injected session. The MCP/cookie-session click works;
  the difference is worth tracing
- Consider exposing a stable `data-testid` (or keeping the `[name="button-DELETE"]`) on the
  confirm button and a `data-state="busy|idle|done"` so spec authors can wait for the
  command to finish instead of guessing on dialog-detached
- The `Dashboards` gallery would benefit from a programmatic refresh hook (e.g. a
  `dashboards.refresh()` JS API or a deterministic event after `dapi.projects.save`)

### Suggestions for the scenario
- The opening note `make sure this is the last test case` is fine, but the scenario itself
  has no setup section. Add an explicit "Pre-condition: at least two projects from
  `uploading.md` and `projects-copy-clone.md` must be present" so a fresh spec can decide
  whether to call setup or rely on prior state
- Step 4 (`Check that the project has been deleted and is no longer present`) is ambiguous
  about "where" — clarify it as "tile no longer appears in `Browse → Dashboards` AND
  `await grok.dapi.projects.filter('name = "{name}"').first()` returns null" so the
  verification surface is testable
