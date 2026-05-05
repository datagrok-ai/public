# Sparql Connection — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse > Databases, click ellipsis to list providers | 6s | PASS | PASSED | Page has 3 `[name="icon-ellipsis-h"]` icons; the relevant one is scoped inside the Databases group. After click, 78 providers listed including Sparql. |
| 2 | Right-click Sparql, choose "Add new connection" | 4s | AMBIGUOUS | PASSED | Actual menu label is `New connection...` (not literally "Add new connection"); same intent. Dialog title is "Add new connection". |
| 3 | Enter `test_sparql` in Name | 3s | PASS | PASSED | Focus + Ctrl+A + type_text — standard Dart-input pattern. |
| 4 | Endpoint, Requires Server=true, Prefixes empty | 5s | PASS | PASSED | Endpoint typed; Requires Server checkbox flipped false→true; Prefixes left blank. The Prefixes input is a `<textarea name="Prefixes">` (no `input-` prefix on the field itself). |
| 5 | Click TEST | 1s | FAIL | FAILED | Server returned `failed to connect: SocketException: Failed host lookup: 'data.ontotext.com' (OS Error: Name or service not known, errno = -2)`. Real platform/infra issue — dev cannot resolve `data.ontotext.com`. UI flow itself worked. |
| 6 | Click OK to save | 2s | PASS | PASSED | Connection persisted (id `bf8afaa0-21a0-11f1-8bf2-5f6e65ae0088`, dataSource `Sparql`). |
| 7 | Right-click test_sparql, Delete, confirm YES | 8s | AMBIGUOUS | PASSED | Confirm button is `DELETE`, not `YES`. Tree did not auto-refresh after save — required clicking the refresh icon and triggering the Sparql node's triangle/dblclick before the new connection appeared. |

**Time** = step 2b wall-clock per step (incl. thinking). **Result** = step 2b outcome. **Playwright** = step 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~1m 10s |
| grok-browser execution (scenario steps) | ~30s |
| Execute via grok-browser (total) | ~1m 40s |
| Spec file generation | ~2m |
| Spec script execution | 29s |
| **Total scenario run (with model)** | ~5m |

The two `scenario steps` rows sum to `Execute via grok-browser (total)`.

## Summary

The scenario reproduced cleanly through the UI on dev — the Add-new-connection dialog opened,
all fields accepted input, the connection saved, and the right-click delete cycle worked
end-to-end. **Step 5 (Test) is a real platform failure**: dev's connector backend cannot
DNS-resolve `data.ontotext.com`, so the connection test never reaches the SPARQL endpoint.
Steps 2 and 7 both required minor relabeling (menu item is "New connection..."; confirm
button is "DELETE"), so the scenario text is slightly off from what the UI actually shows.
Total scenario run (with model): ~5m.

## Retrospective

### What worked well
- `[name="..."]` selectors held up everywhere — every dialog input, button, menu item, and tree expander resolved on the first try.
- The shared `loginToDatagrok` + `specTestOptions` helper made the spec preamble trivial; token derived via `users/login/dev/<key>` from `~/.grok/config.yaml`.
- Synthetic `contextmenu` events on `.d4-tree-view-node` worked reliably for the right-click steps.
- Scoping the Step 1 ellipsis click *inside* the Databases group (not just the first ellipsis on the page) made the spec robust against pages with multiple section ellipses.

### What did not work
- **Step 5 (Test)** — `dev.datagrok.ai` server-side cannot resolve `data.ontotext.com`. Likely a deliberate egress/DNS restriction on the dev environment, but it makes any scenario hitting external SPARQL endpoints unrunnable on dev.
- **Tree did not auto-refresh after `OK`** — the just-saved `test_sparql` node was not visible under Sparql until the refresh icon was clicked and the provider node re-expanded (triangle click + dblclick fallback).
- The scenario assumed Tabs mode, but the spec navigates with `?browse=connections`, so we operate in the Browse view's tree. That worked, but the scenario doesn't specify which surface (left tree vs. center connection view).

### Suggestions for the platform
- Either allow the dev environment to reach common public test endpoints (data.ontotext.com, etc.) or document an internally-reachable SPARQL endpoint for QA scenarios so this case can pass on dev.
- Auto-refresh the Browse tree provider node when `grok.dapi.connections.save()` returns. Today the new connection only appears after a manual refresh + re-expand, which surprises both human testers and automated specs.
- Consider giving destructive-confirm dialogs a stable `name="button-CONFIRM"` (or matching the scenario's "YES" wording) so test code doesn't have to know whether it's `DELETE`, `YES`, `REMOVE`, etc.
- Give the Prefixes textarea a `name="input-Prefixes"` to match the convention used by every other field in the dialog (Name, Endpoint, Requires-Server, Secret-Name).

### Suggestions for the scenario
- Step 2: change "Add new connection" to "New connection..." to match the actual menu item.
- Step 7: change "press YES" to "click DELETE" (or "press the destructive confirm button").
- Step 1: clarify "three dots in the end of Databases list" → "click the **…** (ellipsis) icon at the bottom of the Databases section to reveal providers without connections" — readers may otherwise look for a context-menu "More" item.
- Add an explicit precondition that the test endpoint must be reachable from the target server, and call out which environments this is known to pass/fail on.
- After Step 6, add an instruction to refresh the tree (or note that it auto-refreshes if/when the platform bug is fixed) so Step 7 is reproducible.
