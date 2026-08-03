# Queries — query-postprocessing — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse → Databases → Postgres → NorthwindTest | 8s | PASS | PASSED | Tree expansion via `.d4-tree-view-tri`; Postgres `NorthwindTest` resolved to connection `PostgresTest`. |
| 2 | Right-click NorthwindTest → New Query…, set name + body | 12s | PASS | PASSED | Context-menu items found by `.d4-menu-item-label` text (no `[name="div-New-Query..."]` was attached on the dispatched contextmenu). Name set via `execCommand`; body set via `.CodeMirror.setValue`. |
| 3 | Run query — inline grid appears | 6s | PASS | PASSED | Ribbon Play button clicked; grid canvas appeared on first attempt. |
| 4 | Switch to Post-Process — replace line 7 with grok.shell.info | 8s | PASS | PASSED | Template seeds line 7 (`console.log(result.rowCount)`); replaced with `grok.shell.info(result.rowCount);` via CodeMirror `replaceRange`. |
| 5 | Switch to Layout — add scatter + correlation plots | 11s | PASS | PASSED | Layout tab inherited the run result; `[name="icon-scatter-plot"]` and `[name="icon-correlation-plot"]` added both viewers. |
| 6 | Save the query | 5s | PASS | PASSED | `[name="button-Save"]` clicked; query persisted as `TestPostprocessing` (server normalized) with friendlyName `Test_Postprocessing`. |
| 7 | Close all + re-navigate to Browse | 4s | PASS | PASSED | `grok.shell.closeAll()` + `/browse` URL navigation. |
| 8 | Preview Test_Postprocessing — viewers + balloon 77 | 14s | PASS | PASSED | Click on the query opened a TableView with the saved layout; both viewers rendered; `grok.shell.info` hook captured `"77"`. |
| 9 | Right-click Test_Postprocessing → Edit… | 7s | PASS | PASSED | Editor reopened with persisted body; layout/post-process saved with the query. |
| 10 | Run from Post-Process tab — balloon 77 | 6s | PASS | PASSED | Single Play click captured `"77"` from `grok.shell.info`. |
| 11 | Run from Layout tab — balloon 77 + viewers | 7s | PASS | PASSED | Both `viewer-Scatter-plot` and `viewer-Correlation-plot` present in the live view; `"77"` captured. |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 3m 46s |
| grok-browser execution (scenario steps) | 1m 30s |
| Execute via grok-browser (total) | 5m 16s |
| Spec file generation | 3m 12s |
| Spec script execution | 44s |
| **Total scenario run (with model)** | 9m 12s |

## Summary

All 10 scenario steps reproduced cleanly on dev. The ribbon Play button, post-process JS hook,
saved layout, and `Edit…` round-trip all worked on the first try; the `Test_Postprocessing` query
persisted correctly and re-applied its scatter + correlation layout and
`grok.shell.info(rowCount)` post-process on every preview/run path. The generated Playwright spec
passed end-to-end in 44s. **Total scenario run (with model)**: 9m 12s.

## Retrospective

### What worked well
- The `grok.shell.info` monkey-patch turned the visual "green balloon 77" check into a deterministic
  assertion (`__capturedInfo` array) — far more reliable than scraping balloon DOM.
- Layout tab automatically inherited the result DataFrame from the prior Run, so adding viewers via
  Toolbox icons "just worked" without re-running.
- `name = "Test_Postprocessing"` server-normalizes to `TestPostprocessing`; filtering with
  `name in ("Test_Postprocessing", "TestPostprocessing")` covers both for pre-cleanup.

### What did not work
- `New Query…` context-menu items did not have `[name="div-New-Query..."]` attributes on the
  dispatched-contextmenu menu, so had to fall back to `.d4-menu-item-label` text matching. (Worked,
  but worth noting that the `name=` selectors documented in `connections.md` aren't always present.)
- Running the spec required deriving `DATAGROK_AUTH_TOKEN` from `~/.grok/config.yaml`'s dev key via
  `POST /api/users/login/dev/<key>` — `loginToDatagrok` errors out otherwise. Not a failure, just
  an obscure prerequisite for ad-hoc spec runs outside `grok test`.

### Suggestions for the platform
- Surface `name=` attributes consistently on connection / query context-menu items
  (e.g. `New Query...`) — the `connections.md` reference promises them but they aren't applied
  uniformly.
- `q.options` / `q.script` / `q.layoutId` returned by `grok.dapi.queries.find()` came back empty
  even after a successful save; if those fields were exposed on the `Query` entity, specs could
  verify persistence directly without round-tripping through the UI.

### Suggestions for the scenario
- Step 4 says "On line 7, add: `grok.shell.info(result.rowCount);`" — but the post-process
  template already has `console.log(result.rowCount)` on line 7. Clarify whether to **replace** or
  **insert** (replacement makes the balloon assertion clean; insertion still works but leaves a
  stray `console.log`).
- Step 9 says "Preview and run" — on dev, single-clicking the query in Browse opens a full
  TableView (effectively running it); there is no distinct preview step in the Context Panel.
  Tighten the wording to match the observed behavior, or call out the Context Panel vs. full
  TableView difference.
- Step 10 sub-bullets verify the balloon and viewers separately for Post-Process and Layout tabs —
  but both run the same query, so the balloon will fire identically; the meaningful difference is
  that Layout tab also renders the viewers. Worth making that distinction explicit.
