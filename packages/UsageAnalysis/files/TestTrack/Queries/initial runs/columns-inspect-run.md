# Queries — Column Inspection on PostgresDart / Postgres — Run Results

**Date**: 2026-05-04
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

### Part 1 — PostgresDart → NorthwindTest
| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1.1 | Browse → Databases → PostgresDart → NorthwindTest → Schemas → public | 1m 12s | PASS | PASSED | Tree label `NorthwindTest` is the friendly name for connection `Dbtests:PostgreSQLTest` (PostgresDart). 14 tables under public. Most of the 1m 12s went into discovering that `Schemas` and `public` are unnamed group nodes — must be located by label scoped under the parent's children host, not by `[name=...]`. |
| 1.2 | Expand each DB table to the column level | 12s | PASS | PASSED | Click on `.d4-tree-view-tri` of each table group reliably expanded to columns. 14 tables expanded. |
| 1.3 | Select each DB column for each DB table (click) | 37s | PASS | PASSED | 90 columns clicked in total (categories=4, customercustomerdemo=2, customerdemographics=2, customers=11, employees=18, employeeterritories=2, order_details=5, orders=14, products=10, region=2, shippers=3, suppliers=12, territories=3, usstates=4). |
| 1.4 | No errors on Context Panel | <1s | PASS | PASSED | Zero visible `.d4-balloon`/`.grok-balloon` and zero accordion panes with `.d4-error`. Last current object: column `stateregion` — accordion panes: General, Inspect, Database meta. |

### Part 2 — Postgres → Northwind
| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 2.1 | Browse → Databases → Postgres → Northwind → Schemas → public | 29s | PASS | PASSED | Tree label `Northwind` resolves under Postgres. 14 tables under public (same Northwind schema). |
| 2.2 | Expand each DB table to the column level | 8s | PASS | PASSED | Same chevron-click pattern as Part 1; all 14 tables expanded. |
| 2.3 | Select each DB column for each DB table (click) | 32s | PASS | PASSED | 90 columns clicked. |
| 2.4 | No errors on Context Panel | <1s | PASS | PASSED | Zero balloons, zero accordion errors. Same accordion shape as Part 1. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 30s |
| grok-browser execution (scenario steps) | 3m 13s |
| Execute via grok-browser (total) | 4m 43s |
| Spec file generation | 1m 30s |
| Spec script execution | 2m 45s |
| **Total scenario run (with model)** | 8m 58s |

`Spec script execution` is end-to-end including the one retry (attempt 1 failed because the connection-child node hadn't lazy-loaded after a fixed 1.5s wait; replaced with a polling waitFor and attempt 2 passed in 1m 30s).

## Summary

Both parts pass cleanly on dev. Walking the Browse tree from Databases →
Provider → Connection → Schemas → public and expanding every table works
reliably with chevron clicks; clicking each of the 90 columns under public
produces no error balloons and no error indicators on the Context Panel
accordion. The one Playwright retry was due to a fixed wait being too short
for first-load lazy mounting of provider children — replaced with a polling
`waitForSelector`. **Total scenario run (with model)**: 8m 58s.

This contradicts the previous (2026-04-24) run which marked column-level
expansion as "not automatable" and `NorthwindTest` as missing from
PostgresDart. The connection is present (it's the friendly name of
`Dbtests:PostgreSQLTest`), and `.d4-tree-view-tri` chevron clicks are a
stable expand affordance.

## Retrospective

### What worked well
- `[name="tree-Databases---{Provider}---{Connection}"]` selector resolves
  connection nodes uniquely, dodging label-collision issues (e.g. multiple
  `NorthwindTest` labels under different providers).
- `.d4-tree-view-tri` chevron click is a reliable, side-effect-free expand
  trigger for every level of the tree (provider, connection, table). No
  accidental Get-All / Navigate fires.
- Polling for child-node availability (vs fixed sleep) is the right pattern
  for tree-driven flows — fixed waits flake on a fresh Playwright context.

### What did not work
- `grok.dapi.connections.filter('name = "Northwind"')` and
  `filter('name = "NorthwindTest"')` both returned smart-search substring
  matches rather than exact-name matches, masking the fact that the
  connections actually exist on dev. Had to inspect tree DOM to confirm.
- Initial spec used a fixed 1500ms post-expand wait; on a fresh Playwright
  session, the Datagrok server's first DB-driven schema fetch can take
  longer, so the connection child wasn't in the DOM yet. Replaced with
  polling.

### Suggestions for the platform
- Consider giving `Schemas` / `public` / table group nodes their own
  `name=` attributes (e.g. `tree-Databases---Postgres---Northwind---Schemas`,
  `...---public`, `...---public---orders`), as already done for the
  provider/connection levels. Currently the absence forces automation to
  fall back to label-text scanning.
- `grok.dapi.connections.filter('name = "X"')` should respect `=` as
  exact match rather than substring, or expose a clearly-named exact-match
  variant. Alternatively, document the smart-search behavior in the
  reference (`connections.md`).

### Suggestions for the scenario
- Rename the scenario step labels to use `friendlyName`s exactly as they
  appear in the tree (`NorthwindTest`, `Northwind`) and note that the
  underlying connection `name` may differ (e.g. `PostgreSQLTest`).
- Tighten step 5 wording: "no errors on Context Panel" — clarify that this
  means (a) no visible error balloons platform-wide and (b) no error
  indicators on accordion panes for the selected column. Both are checked
  by this run.
