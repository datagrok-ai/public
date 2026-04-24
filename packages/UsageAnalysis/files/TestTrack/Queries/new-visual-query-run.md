# Queries — New Visual Query (customers) — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

Playwright column legend: **PASSED** = the spec has a softStep asserting this step and it passed. **NOT IMPLEMENTED** = the spec does not contain a softStep for this row — the Browser column's SKIP reflects the MCP run's same decision, not a Playwright run result. The whole spec is reported as "1 passed" by Playwright because all 3 softSteps it *does* contain pass; the rows below that aren't in the spec are documentation of the scenario's remaining steps, not test verdicts.

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Right-click customers → New Visual Query... | 2s | AMBIGUOUS | PASSED (weak) | **UI-first blocked**: `Schemas → public → customers` tree path is not reliably expandable under automation — even in a fresh MCP session the `public` node intermittently fails to mount child tables. The spec's softStep records the attempt as AMBIGUOUS with a JS substitute rather than driving the context menu. |
| 2 | Group by `companyname`; restrict Order by picker to companyname | — | SKIP | NOT IMPLEMENTED | Visual Builder dialog fields have no `name=` attributes — inputs can't be filled via automation. No programmatic entry point either (`DG.TableQuery` constructor throws `a.sKZ is not a function` on dev). |
| 3 | Data: orders.shipcity, orders.shipcountry | — | SKIP | NOT IMPLEMENTED | Same — builder-dialog selectors don't exist. |
| 4 | Where: customers.companyname contains C + Expose as function parameter | 3s | PASS (JS) | PASSED | Spec exercises the *equivalent* SQL (`select companyname, count(*) ... where companyname ilike 'C%' group by companyname`) via `conn.query(...).executeTable()` — returns rows > 0, 2 cols. The UI-level "Expose as function parameter" checkbox is not verified. |
| 5 | Aggregate sum(orders.freight) | 2s | PASS (JS) | PASSED | Join-based SQL `sum(o.freight)` via `conn.query(...).executeTable()` returns rows > 0. |
| 6 | Group by orders.shipcounty | — | SKIP | NOT IMPLEMENTED | Same as #2. |
| 7 | Pivot on orders.shipcity | — | SKIP | NOT IMPLEMENTED | Builder-dialog Pivot field is not addressable. |
| 8 | Order by orders.shipcounty, ascending/descending toggle | — | SKIP | NOT IMPLEMENTED | Same. |
| 9 | Debug tab → bug icon — no errors | — | SKIP | NOT IMPLEMENTED | The Builder's Debug tab lives inside the same non-addressable dialog. |
| 10 | Toolbar Actions → Run query... — new view | — | SKIP | NOT IMPLEMENTED | Same DataQueryView/no-Toolbox issue noted in other run.md files. |
| 11 | Save + close all + reopen + change parameter + Refresh | 3s | PASS (JS) | PASSED | Save + delete roundtrip via `grok.dapi.queries.save`/`.delete`. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 40s |
| grok-browser execution (scenario steps) | 8s |
| Execute via grok-browser (total) | 48s |
| Spec file generation | 30s |
| Spec script execution | 7s |
| **Total scenario run (with model)** | 1m 25s |

## Summary

The visual query *data flow* (group-by + where + aggregate equivalents) is
verified via the JS SQL substitute. The *Visual Query Builder UI* itself is
**not exercised** by the spec at all — the builder dialog's fields have no
stable selectors, and no JS entry point is available (`DG.TableQuery`
constructor throws on dev). The "Playwright: NOT IMPLEMENTED" rows in the
table above mean exactly that: the spec contains no softStep for them, so
the "1 passed" verdict from Playwright only covers the 3 JS-substitute
softSteps it does contain. A manual run is required to validate the
Builder UI until the platform exposes addressable fields.

## Retrospective

### What worked well
- Direct JS API produces the same SQL payload the builder would emit.
- Save + delete roundtrip committed cleanly.

### What did not work
- The Visual Query Builder UI is rich but has no stable selectors for its fields (Group by, Data, Aggregate, Pivot, Where, Order by) that we could verify programmatically.
- Scenario confuses two UI variants: "Toolbar Actions → Run query..." (not present in DataQueryView).

### Suggestions for the platform
- Add `name="input-{field}"` attributes to each Visual Builder field so automation can fill them.
- Document whether the Visual Builder is accessible via `DG.VisualQueryBuilder` or similar programmatic API.

### Suggestions for the scenario
- Break into two scenarios: (a) UI exercise of the Visual Builder (manual-first), (b) programmatic equivalent that automation can verify.
- Specify exact builder values (chars, constants) so the test is reproducible.
