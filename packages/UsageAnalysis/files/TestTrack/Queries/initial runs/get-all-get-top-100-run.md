# Queries — Get All / Get Top 100 — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

### Part 1 — PostgresDart > NorthwindTest > Schemas > public > orders
| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1.1 | Browse → PostgresDart → NorthwindTest → Schemas → public → orders | 12s | PASS | PASSED | Tree node fully expanded; container `div-PostgresDart-PostgreSQLTest-Schemas` (uses internal name `PostgreSQLTest`). |
| 1.2 | Right-click orders → Get All | 8s | PASS | PASSED | New view "orders" opened with 830 rows × 14 cols. Initial render showed 50 rows then settled at 830. |
| 1.3 | Right-click orders → Get Top 100 | 6s | PASS | PASSED | New view "orders (2)" opened with 100 rows × 14 cols. |

### Part 2 — Postgres > Northwind > Schemas > public > orders
| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 2.1 | Browse → Postgres → Northwind → Schemas → public → orders | 6s | PASS | PASSED | Container `div-Postgres-PostgresNorthwind-Schemas`. |
| 2.2 | Right-click orders → Get All | 7s | PASS | PASSED | New view "orders (3)" opened with 830 rows × 14 cols. |
| 2.3 | Right-click orders → Get Top 100 | 6s | PASS | PASSED | New view "orders (4)" opened with 100 rows × 14 cols. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 20s |
| grok-browser execution (scenario steps) | 45s |
| Execute via grok-browser (total) | 2m 5s |
| Spec file generation | 4m 40s |
| Spec script execution | 31s |
| **Total scenario run (with model)** | 7m 16s |

All rows are full-phase wall-clock (incl. model thinking and retries).
The two `scenario steps` rows sum to `Execute via grok-browser (total)`.

## Summary

Both halves of the scenario PASS end-to-end on dev. Get All on the orders table
returns 830 rows × 14 cols and Get Top 100 returns 100 rows × 14 cols, for both
PostgresDart > NorthwindTest and Postgres > Northwind. The Playwright replay
went green on the third attempt (31s) after fixing two automation pitfalls. The
prior run (2026-04-24) had marked Part 1 as FAIL because NorthwindTest under
PostgresDart wasn't seeded — that connection now exists on dev (internal name
`PostgreSQLTest`, friendlyName `NorthwindTest`). Total scenario run (with model):
**7m 16s**.

## Retrospective

### What worked well
- Right-click on a Schemas table node IS automatable via synthetic
  `contextmenu` dispatch on the table's `.d4-tree-view-node` ancestor — both
  `Get All` and `Get Top 100` items appear and trigger correctly. This
  contradicts a prior run's finding that table-node context menus weren't
  automation-friendly; the path now reliably reproduces in MCP and replays
  in Playwright.
- The `name=` attribute on the Schemas group container
  (`div-{Provider}-{InternalConnName}-Schemas`) is stable enough to anchor
  Playwright lookups, but a top-down label-based traversal is more
  resilient against the connection's internal vs. friendly name mismatch.
- Get All initially renders ~50 rows asynchronously and then settles to the
  full row count — a stable-value polling loop (3 consecutive identical
  reads above the expected minimum) avoids both flakes and premature passes.

### What did not work
- **First spec attempt — fixed selector** assumed the Schemas container exists
  on initial mount. It doesn't: connection-level expansion triggers a remote
  schema load, so `[name="div-PostgresDart-PostgreSQLTest-Schemas"]` only
  appears after PostgresDart > NorthwindTest is clicked.
- **Second spec attempt — Browse toggle** unconditionally clicked
  `[name="Browse"]` to "open" the panel. After login the Browse panel is
  already mounted, so the click *toggled it off* — `treeRoot` count went from
  1 to 0 and the tree never rendered. Probe spec confirmed: 84 labels visible
  immediately after `loginToDatagrok`.

### Suggestions for the platform
- Add a stable `name=` attribute on `.d4-tree-view-node` for table-level rows
  (e.g. `name="tree-Databases---{Provider}---{Connection}---Schemas---{schema}---{table}"`)
  so right-click targets don't have to be located by label-walking.
- Make the Browse-panel toggle button reflect the panel's current open/closed
  state in its `aria-pressed` (or analogous) attribute, so automation can
  decide whether to click without probing the tree's children count.

### Suggestions for the scenario
- Replace plain "Northwind" / "NorthwindTest" with a hint about which provider
  exposes which (`Postgres → Northwind`, `PostgresDart → NorthwindTest`) — the
  current wording reads as if both providers expose both names.
- Mention that the Schemas group expansion can take a few seconds (remote
  schema load) so manual testers don't mistake the spinner for a hang.
