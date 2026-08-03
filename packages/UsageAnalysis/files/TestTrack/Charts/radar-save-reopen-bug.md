---
feature: charts
target_layer: playwright
coverage_type: edge
priority: p1
realizes_atlas: [radar-table-rebind-project-roundtrip]
realizes: [charts.radar]
realized_as:
  - radar-save-reopen-bug-spec.ts
pyramid_layer: bug-focused
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/radar-save-reopen-bug.md
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
related_bugs:
  - GROK-18085
---

# Radar viewer — table-rebind on project save/reopen (GROK-18085)

Bug-focused regression scenario for GROK-18085: a Radar viewer's
table-rebind state serializes incorrectly across project save →
close-all → reopen, causing a null error on deserialization. The
canonical `radar.md` smoke scenario does not exercise project
save/reopen, so this scenario exists specifically to lock in the
GROK-18085 reproduction: rebind a Radar viewer's bound table
mid-session, save the project, reopen it, and verify Radar deserializes
correctly with the new table binding.

## Setup

A clean Datagrok session. Single scenario; uses `grok.dapi.projects` for
the save/reopen flow programmatically (UI Save Project dialog requires
selectors deferred per cycle charts-migrate-2026-05-07).

## Scenarios

### Scenario 1: Radar table-rebind survives project save → reopen (GROK-18085 invariant)

Steps:

1. Open `System:DemoFiles/demog.csv` and `System:DemoFiles/chem/SPGI.csv` in the
   same session (both as TableViews via `grok.shell.addTableView`).
2. On the SPGI view, add a Radar viewer (`tv.addViewer('Radar')`).
   **Expected:** Radar attached to SPGI without console error.
3. Rebind the Radar's bound table from SPGI to demog via
   `radar.setOptions({table: 'demog'})`. Wait for re-render.
   **Expected:** Radar's `props.get('table')` returns `'demog'` (or the
   read-back races to null — race-tolerant per cycle lessons).
4. Save the project: save the CURRENT WORKSPACE project
   (`grok.shell.project`) — it holds the open demog/SPGI views, the Radar,
   and the rebind layout — via `grok.dapi.projects.save`, after setting a
   unique name like `radar-rebind-${Date.now()}`. Do NOT construct a fresh
   `DG.Project.create()`: an empty project carries none of the rebind state,
   so Steps 6-7 (verify the rebind survived) would have nothing to verify.
   **Expected:** save resolves; project carries non-empty `id`.
5. Close all views (`grok.shell.closeAll()`).
6. Reopen the project via `grok.dapi.projects.find(<id>)` + open API.
   **Expected (GROK-18085 invariant):** project reopens WITHOUT a null
   error / deserialization exception. Console errors during reopen are
   either empty or filtered to benign noise (Failed to load resource etc.).
7. Verify the Radar viewer is present in the reopened table view and
   `radar.props.get('table')` returns `'demog'` (or null on race).
8. Cleanup: delete the saved project via `grok.dapi.projects.delete`.

## Notes

- GROK-18085 invariant: Step 6 — the project must reopen without a null
  error / deserialization exception. Console-error capture filters out
  benign network noise (404s, favicon, "Failed to load resource").
- Cleanup: Step 8 deletes the saved project; the spec wraps this in
  try/finally so cleanup runs even if an assertion fails.

## Dataset metadata

```json
{
  "order": 32,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/chem/SPGI.csv"]
}
```
