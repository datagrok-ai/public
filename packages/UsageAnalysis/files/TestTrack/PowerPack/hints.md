---
feature: powerpack
target_layer: playwright
coverage_type: smoke
priority: p1
realizes_atlas: []
realizes: []
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-function-hint-tooltip
ui_coverage_delegated_to: add-new-column.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/hints.md
migration_date: 2026-05-20
source_text_fixes:
  - fucntion-typo-corrected-to-function
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
realized_as:
  - hints-spec.ts
gate_verdicts:
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T09:49:00Z
    spec_runs:
      - spec: hints-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 23
        failure_keys: []
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-hints
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-hints
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
    review_round: 1
---

# Add New Column dialog — function-hint tooltip

Focused UI scenario verifying that hovering over an inserted function
name in the Add New Column dialog's formula editor surfaces a tooltip
containing the function's signature.

The dialog's basic open-and-add flow is covered by `add-new-column.md`.

## Setup

A clean Datagrok session is the only shared setup. The scenario opens
its own dataset; no fixture chaining is required.

## Scenarios

### Hover over inserted function name → tooltip with signature appears

1. Open the demog dataset (`System:DemoFiles/demog.csv`) — for example
   via **Browse** > **Files**, **File** > **Open**, or the equivalent
   JS-API loader. **Verify:** the demog grid renders.
2. Open the **Add New Column** dialog — click the `Add New Column`
   icon on the table view toolbar (or **Edit** > **Add New Column**
   from the top menu). **Verify:** the Add New Column dialog opens
   with its formula editor (CodeMirror), functions panel, and preview
   grid visible.
3. Insert any function into the formula text field using any
   supported approach (typing the function name and confirming via
   autocomplete, clicking a function entry in the functions panel, or
   dragging a function onto the formula field). For example, type
   `a` and select `Abs` from the autocomplete tooltip; the function
   is inserted as `Abs(num)`.
4. Hover the mouse pointer over the inserted function name (e.g.
   over `Abs` in `Abs(num)`). **Verify:** a tooltip appears
   containing the function's signature (function name plus its input
   parameter list, e.g. `Abs(num)`).

---
{
  "order": 3
}
