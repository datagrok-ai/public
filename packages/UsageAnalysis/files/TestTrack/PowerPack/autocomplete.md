---
feature: powerpack
target_layer: playwright
coverage_type: smoke
priority: p1
realizes: []
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-autocomplete
  - add-new-column-ctrl-space
  - add-new-column-dollar-column-suggestions
ui_coverage_delegated_to: add-new-column.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/autocomplete.md
migration_date: 2026-05-20
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
realized_as:
  - autocomplete-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-autocomplete
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T09:34:58Z
    spec_runs:
      - spec: autocomplete-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 51
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-autocomplete
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
    review_round: 1
---

# Add New Column — autocomplete (demog)

Focused UI scenario for the autocomplete mechanics of the
**Add New Column** dialog: the type-triggered function tooltip,
explicit invocation via **Ctrl+Space**, and the **$**
column-suggestion tooltip. The dialog's basic open/close and
OK/Cancel behavior is covered by `add-new-column.md`.

## Setup

1. Open `demog.csv` from `System:DemoFiles` (e.g. via
   `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` + a
   table view).
2. Open the **Add new column** dialog (toolbar icon
   `[name="icon-add-new-column"]` on the grid ribbon, or
   **Edit | Add New Column** top-menu).

## Scenarios

### Autocomplete: type-triggered function tooltip

1. Click the formula editor (CodeMirror `.d4-dialog .cm-content`)
   to give it focus.
2. Type the letter `a`.
3. **Verify:** an autocomplete tooltip appears listing functions
   whose names start with `a` (e.g. `Abs`, `Acos`, `Avg`).
4. Select a function from the list — try **both** completion
   actions:
   - Press **Enter** on the highlighted entry, AND
   - Use mouse click on an entry (re-do the type-`a` step first
     to reopen the tooltip).
5. **Verify:** the function is inserted into the editor in the
   form `Abs(num)` — function name followed by parenthesised
   parameter-type placeholders.
6. Remove the inserted function from the editor (select all +
   delete, or backspace).

### Autocomplete: Ctrl+Space explicit invocation

1. With the editor focused and empty, press **Ctrl+Space**.
2. **Verify:** an autocomplete tooltip appears with the full
   function list (same widget as the type-triggered case, just
   triggered explicitly).

### Autocomplete: `$` column-suggestion tooltip

1. With the editor focused (clear any open tooltip with
   **Escape** first if needed), type the **`$`** character.
2. **Verify:** the autocomplete tooltip appears listing the
   **columns** of the current dataset (for `demog.csv`,
   entries such as `HEIGHT`, `WEIGHT`, `AGE` appear), distinct
   from the function-name list of the previous two scenarios.

## Notes

- **Related bug GROK-17004** (paste-handler crash on complex
  formulas) touches the dialog's autocomplete / paste surface
  but is NOT directly reproduced here — this scenario tests
  autocomplete on **typing**, not paste. The paste-crash regression
  is covered by `highlight.md`.
- **Source dataset.** `demog.csv` is platform-provided
  (`System:DemoFiles`); no other scenario in this section depends
  on it.

---
{
  "order": 2
}
