# Create metadata schema & entity type — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Log in with admin permissions | 0s | PASS | PASSED | Already logged in from prior MCP session — spec runs login from scratch. Login form is plain HTML, so `.fill()` is more reliable than `keyboard.type()` which drops characters on dev's first paint |
| 2 | Navigate to Browse → Platform → Sticky Meta → Types | 12s | PASS | PASSED | Tree click via Browse → Platform → Sticky Meta → Types. Direct `page.goto('/meta/types')` is unreliable: Datagrok's SPA ignores sibling-route goto() and reloads to the most-recently-visited Sticky-Meta tab |
| 3 | Click "NEW ENTITY TYPE..." button | 7s | PASS | PASSED | `[name="button-New-Entity-Type..."]` opens the "Create a new entity type" dialog |
| 4 | Enter name TestEntity1, matcher `semtype=molecule`, Save | 33s | PASS | PASSED | Focus input + Ctrl+A + `keyboard.type(..., {delay: 10})` for both fields. In the spec, OK is confirmed via `page.mouse.click()` at its bounding box — JS `.click()` and `keyboard.press('Control+Enter')` both silently no-op. Verified via search (list is paginated, 18 of 57) |
| 5 | Navigate to Browse → Platform → Sticky Meta → Schemas | 10s | PASS | PASSED | Tree click on "Schemas" (sibling of Types in the already-expanded tree) |
| 6 | Click "NEW SCHEMA..." button | 8s | PASS | PASSED | `[name="button-New-Schema..."]` opens the "Create a new schema" dialog |
| 7 | Enter TestSchema1 and associate with TestEntity1 | 34s | PASS | PASSED | Outer schema-name input scoped as `.d4-dialog .ui-form > [name="input-host-Name"] > [name="input-Name"]` (the same `name=` is duplicated inside the Properties table — direct-child scoping disambiguates). "Associated with" opens a nested type-picker dialog; checkbox at `[name="prop-view-testentity1"]` |
| 8 | Add 4 properties (rating/int, notes/string, verified/bool, review_date/datetime) and save | 59s | PASS | PASSED | Per-row: click Name input via `nth(i)` → Ctrl+A → `keyboard.type(name, {delay: 10})`; select type via `select.value = type; select.dispatchEvent('input'/'change')`. OK confirmed via `page.mouse.click()`. Verified by right-clicking TestSchema1 → Edit and reading the rows back: associated=TestEntity1, rows=[rating/int, notes/string, verified/bool, review_date/datetime] |

## Timing

All values below are for a single fresh run (dev pre-cleaned of TestEntity1 / TestSchema1, then the scenario executed end-to-end once). Iterative debugging time from earlier sessions is not included.

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~1m |
| grok-browser execution (scenario steps) | ~1m 48s |
| Execute via grok-browser (total) | 2m 48s |
| Spec file generation | 1m 24s |
| Spec script execution | 59s |
| **Total scenario run (with model)** | 5m 11s |

## Summary

All 8 scenario steps PASS in both the MCP browser reproduction and the generated Playwright spec. TestEntity1 is created with matcher `semtype=molecule`; TestSchema1 is associated with TestEntity1 and carries the four expected properties (rating/int, notes/string, verified/bool, review_date/datetime). Spec execution is idempotent: back-to-back runs each pass in ~58s — a search-based pre-cleanup + a duplicate-name recovery branch keep the test self-healing even when prior runs leave entities behind. **Total scenario run (with model)**: 5m 11s (2m 48s MCP reproduction + 1m 24s spec file generation + 59s spec execution).

## Retrospective

### What worked well
- Browse-tree navigation (Platform → Sticky Meta → Types / Schemas) is the only reliable way to switch between Sticky-Meta sibling routes; once you commit to it, the flow is mechanical
- `[name="button-New-Entity-Type..."]` / `[name="button-New-Schema..."]` / `[name="button-OK"]` / `[name="button-CANCEL"]` / `[name="button-DELETE"]` are stable selectors
- The "Associated with" type picker uses `.property-grid-item-editor-checkbox` rows with names like `[name="prop-view-<lowercase-name>"]` — trivially targetable
- Post-save verification via right-click → Edit reveals the schema's full structure (outer `input-Name`, `div-Associated-with-`, `table.d4-item-table` with per-row `input-Name` + `select[name="input-Property-Type"]`) and is perfect for an assertion
- Search-box + right-click → Delete as a pre-cleanup strategy handles the paginated list issue transparently

### What did not work
- `keyboard.type()` without `{delay: 10}` silently dropped characters on dev's login form (the first paint seems to eat fast keystrokes)
- `[name="input-Name"]` alone is ambiguous inside the schema dialog — it matches the outer schema-name input AND every row in the Properties table. Scoping via `.ui-form > [name="input-host-Name"]` (direct child) was necessary
- `page.goto('/meta/types')` after already having visited `/meta/schemas` in the same session silently redirects back to `/meta/schemas` — the Datagrok SPA router cached the last Sticky-Meta URL
- JS `.click()` and `keyboard.press('Control+Enter')` both failed to submit the "Create a new entity type" dialog even when OK was enabled and inputs were valid. Only `page.mouse.click()` at the button's bounding box worked
- The Types list is paginated (shows ~18 of 57). A plain label-by-text lookup misses any entry past the first page, which masked duplicates and caused misleading false-negatives

### Suggestions for the platform
- Give dialog-scoped inputs unique `name=` attributes: `input-host-Schema-Name` for the outer name, `input-host-Property-Name` for table rows. Without this, every Sticky-Meta dialog automation pays the strict-mode ambiguity tax
- Make `/meta/types` and `/meta/schemas` treat URL navigation as authoritative — the SPA currently redirects to whichever tab was opened first in the session, which is surprising
- Login form: debounce keystrokes properly so fast `keyboard.type()` input isn't swallowed
- Add a `name=` shortcut for navigating directly to Sticky-Meta tabs (e.g., `[name="div-Sticky-Meta---Types"]`) so tests don't need tree expansion logic
- The Types/Schemas lists should paginate in the DOM with a "load more" element, not silently truncate. Currently any automation that enumerates labels is brittle without using the search box

### Suggestions for the scenario
- Add an explicit "Precondition: TestEntity1 / TestSchema1 must not exist (delete if present via right-click → Delete from the list)"
- Note in the scenario body that the example names are consumed by the `add-and-edit.md` downstream scenario, so cleanup at end is intentional (leave them)
- Clarify that "Associated with:" is a link-action (`select entities`) that opens a nested type-picker dialog with checkboxes — not a combobox
