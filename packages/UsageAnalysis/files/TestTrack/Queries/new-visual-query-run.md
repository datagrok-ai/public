# New Visual Query — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse → DBs → Postgres → NorthwindTest → Schemas → public | 45s | PASS | PASSED | Tree expansion via `name="div-Postgres-PostgresTest-Schemas"` group; click triangle (`.d4-tree-view-tri`) to expand collapsed nodes. 14 tables shown including `customers`. |
| 2 | Right-click customers → New Visual Query… | 8s | PASS | PASSED | Context menu opens via `dispatchEvent('contextmenu')` on `.d4-tree-view-node`; visual query editor opens with Data/Where/Group by/Aggregate/Order by panels. |
| 3 | Set Group by to companyname | 35s | PASS | FAILED | Synthetic mouse-click on the canvas-based picker grid worked in MCP but not in the standalone Playwright run — canvas-grid checkbox/cell hits depend on Dart pointer-event listeners that drop synthetic events at headless rates. |
| 4 | Order by picker shows only companyname | 7s | PASS | FAILED | After Group by=companyname, Order by + opened a small picker (height ≈ 48px, only header + 1 row); MCP run confirmed the only selectable column was `companyname`. Spec failed because step 3 didn't materialize, so picker had no candidate. |
| 5 | Data → select orders.shipcity, orders.shipcountry | 4m | AMBIGUOUS | FAILED | Add-table flow worked (Data → Add a table → northwind → public → orders adds the join). Restricting columns required toggling individual rows in a canvas-rendered "Select columns…" dialog, which doesn't react to synthetic clicks; ended up selecting all 14 orders columns via the "All" link. Spec failed at the join hover step; MCP got further by triggering `mouseenter`/`mousemove`/`mouseover` manually. |
| 6 | Set Where to customers.companyname contains C, expose as parameter | — | SKIP | n/a | The Where + opens a "row table selector" submenu (northwind.public.customers → … → companyname). The submenu does not expand from synthetic `mouseenter`/`mouseover`/`hover()` events, only from real OS hover. Could not reach the columns level to type the predicate or check "Expose as function parameter". |
| 7 | Set Aggregate to sum(orders.freight) | — | SKIP | n/a | Same submenu-on-hover blocker as step 6. |
| 8 | Add Group by: orders.shipcounty | — | SKIP | n/a | Same blocker. (Note scenario typo: should be `shipcountry`.) |
| 9 | Set Pivot to orders.shipcity | — | SKIP | n/a | Same blocker. |
| 10 | Add Order by: orders.shipcounty | — | SKIP | n/a | Same blocker. (Same typo.) |
| 11 | Set different order for Order by fields | — | SKIP | n/a | Depends on step 10. |
| 12 | Open Debug tab, press bug icon — no errors | 12s | PASS | FAILED | Debug tab clicked, bug icon (`[name="icon-bug"]`) launched parameter dialog → OK → log shows `Call was started` / `STATEMENT EXECUTION` / `RESULT SET PROCESSING WITH DATAFRAME FILL` / `DRY RUN` with no error rows. Spec failed because earlier steps did not run. |
| 13 | Toolbar → Actions → Run query… opens new view | 6s | PASS | FAILED | Click on `Run query...` link in Toolbox Actions section opens a new view; toolbar shows REFRESH but "Parameters can't be edited visually" since no params were set. Spec failed because earlier steps did not run. |
| 14 | Save the query, Close All | 7s | PASS | FAILED | Save via `[name="button-Save"]` after typing a unique name into the Name input. `grok.dapi.queries.filter('newvq...').list()` confirmed persistence. Spec failed because earlier steps did not run. |
| 15 | Run the saved query — check the parameter | 3s | PARTIAL | FAILED | `executeTable({})` ran and produced 830 rows × 14 columns from orders. The saved visual query lost the customers join / Group by because the underlying state for those was not fully populated — so there is no parameter to check. |
| 16 | On Toolbar, change parameter and Refresh | — | SKIP | n/a | No parameter was ever defined (step 6 was skipped). |
| 17 | Close all | <1s | PASS | n/a | `grok.shell.closeAll()` and saved query was deleted in cleanup. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~14m |
| grok-browser execution (scenario steps) | ~6m |
| Execute via grok-browser (total) | ~20m |
| Spec file generation | 2m |
| Spec script execution | 56s |
| **Total scenario run (with model)** | ~23m |

## Summary

Scenario PARTIAL: 1, 2, 3, 4, 12, 13, 14, 17 PASS; 5, 15 PARTIAL/AMBIGUOUS; 6–11, 16 SKIPPED. The scenario hits two recurring blockers in the visual query builder UI: (1) the canvas-rendered "Select columns…" dialog and the column pickers triggered by + on Group by / Order by / Pivot do not respond to synthetic mouse events, so individual checkbox / row toggles cannot be driven from MCP or Playwright `evaluate`; (2) the cascading "row table selector" submenu (used by Where / Aggregate / Pivot pickers) does not expand from synthetic `mouseenter` / `mouseover` / Playwright `hover()`, only from real-OS pointer hover. Steps that bypass these widgets (tree navigation, context menu, saving, debug, run via Actions, JS-API run) all work. **Total scenario run (with model)**: ~23m.

## Retrospective

### What worked well
- Browse-tree expansion via the `[name="tree-…"]` and `[name="div-…-Schemas"]` selectors is reliable.
- Right-click → context menu via `dispatchEvent('contextmenu')` worked first time.
- Save flow + JS-API verification (`grok.dapi.queries.filter().first()`) is a clean way to confirm persistence.
- Debug tab + bug icon + parameter dialog OK runs the query and renders the message log without UI hand-holding.

### What did not work
- Canvas-rendered column picker grids (Group by / Order by / Aggregate "+" pickers) did not register synthetic mouse-clicks reliably — worked once in MCP, not in headed Playwright.
- "Select columns…" dialog with the `(0/14)` chip caret could not toggle individual rows; only the search input + All / None links responded. End result was "all 14 selected" instead of "shipcity + shipcountry".
- Cascading submenus (`row table selector` from Where/Aggregate/Pivot/Order by `+`) do not expand on synthetic `mouseenter`/`mouseover` or `page.hover()`. They appeared to require true OS-level hover.

### Suggestions for the platform
- Add stable hooks (`name=` attributes, data-attributes, or JS API) for the visual query builder's column pickers — at minimum a way to call `setColumns(['shipcity', 'shipcountry'])` on the orders join, and `addPredicate(table, column, op, value, exposeAsParam)` on Where, so automated scenarios don't depend on canvas-grid internals.
- Make canvas-rendered cell-toggle interactions respond to synthetic `pointerdown`/`pointerup` events the way the rest of D4 does. Today the only path is real OS pointer input, which makes browser-MCP and Playwright headless agents brittle.
- Make submenus (`d4-menu-popup` with caret-right items) open via DOM hover events (mouseenter/mouseover) consistently, so `page.hover()` works.
- Visual query: when saving with empty Aggregate / partial join state, surface a warning or refuse to save — currently the saved query silently degenerates to `select * from <last-touched-table>` which is hard to debug.

### Suggestions for the scenario
- Fix the typo `shipcounty` → `shipcountry` (steps 8 and 10).
- Step 4 wording "Start setting Order by – check that only companyname column can be selected" is ambiguous about whether to actually pick `companyname` — make it explicit (likely intended to pick it, given step 11 references multiple Order by fields).
- Specify whether `customers.companyname contains C` is case-sensitive in the SQL it generates.
- Add an explicit cleanup step (Delete the saved query) so successive runs don't leave residue.
- Note that the run query is expected to have a single `companyname` parameter (from step 6's Expose as function parameter) so step 15's "check the parameter" is unambiguous.
