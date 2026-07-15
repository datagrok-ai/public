---
feature: powerpack
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [powerpack.cp.add-new-column-persists]
realizes: []
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - add-new-column-dialog
  - context-panel-formula-edit
  - add-new-column-formula-recalc-dependency
  - save-project-with-formula-columns-persistence
ui_coverage_delegated_to: add-new-column.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/formula-refreshing.md
migration_date: 2026-05-20
source_text_fixes:
  - step-1-demog-open-spelled-out-as-system-demofiles-demog-csv
  - step-5-context-panel-formula-edit-clarified-as-formula-info-panel-widget
  - additional-notes-persistence-check-promoted-to-explicit-scenario-step
candidate_helpers:
  - helpers.powerpack.openDemog
  - helpers.powerpack.addCalculatedColumn
  - helpers.powerpack.editFormulaViaContextPanel
  - helpers.powerpack.saveProjectWithDatasync
  - helpers.powerpack.reopenProject
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-17109
realized_as:
  - formula-refreshing-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-formula-refreshing
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T11:13:49Z
    spec_runs:
      - spec: formula-refreshing-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 46
        failure_keys: []
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-formula-refreshing
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
---

# Add New Column — formula dependency chain + recalc + save/reopen persistence

Bug-focused regression scenario for chained calculated columns on the
Demog dataset: build a 3-step dependency chain (`Weight2` → `Weight3` →
`Weight4`), modify the upstream formulas via the `Formula` info panel
widget, verify recalculation propagates through the dependency chain
without errors, then save the project and reopen it to verify the
chained calculated columns persist and recompute correctly.

This is one of several scenarios covering GROK-17109 (calculated
columns must persist across save, data sync, and reopen with correct
recalculation) — this one focuses on the dependency-chain
recalculation and save/reopen persistence; the basic add-column
dialog flow is covered by `add-new-column.md`.

## Setup

1. Open `System:DemoFiles/demog.csv` so the Demog table view is the
   active view. (Source-text fix: original step 1 said only "Open the
   demog dataset" — explicit source per the trailing JSON metadata
   block on the original scenario.)

## Scenarios

### Build a 3-step calculated-column dependency chain

1. **Add column `Weight2 = ${WEIGHT} + 100`.**
   - Open the Add New Column dialog (toolbar "Add new column" icon on
     the Demog table view ribbon).
   - Enter formula `${WEIGHT} + 100` in the formula editor.
   - Enter `Weight2` in the name input.
   - Click OK.
   - A new column named `Weight2` is added to the Demog dataset with
     values correctly computed as `WEIGHT + 100` per row.
2. **Add column `Weight3 = ${Weight2} + 100`.**
   - Open the Add New Column dialog.
   - Enter formula `${Weight2} + 100` (referencing the column added
     in step 1).
   - Enter `Weight3` in the name input.
   - Click OK.
   - A new column named `Weight3` is added; values are correctly
     computed based on `Weight2` per row (i.e. equivalent to
     `WEIGHT + 200`).
3. **Add column `Weight4 = Log10(${Weight3}) - 0.2`.**
   - Open the Add New Column dialog.
   - Enter formula `Log10(${Weight3}) - 0.2` (referencing the column
     added in step 2).
   - Enter `Weight4` in the name input.
   - Click OK.
   - A new column named `Weight4` is added; values are correctly
     computed based on `Weight3` per row.

### Verify formula dependency recalculation via the Formula info panel

1. **Open the `Formula` info panel for `Weight2`.**
   - Select the `Weight2` column header.
   - Open the Context Panel (right-side panel) for the column.
   - Locate the `Formula` info panel widget (registered for columns
     tagged with a formula via `powerpack.formula.is-formula-column`;
     the widget opens an `AddNewColumnDialog` pre-bound to the
     column's existing formula per `powerpack.formula.widget`).
2. **Modify the `Weight2` formula and apply.**
   - Change the formula to a different valid expression that still
     yields numeric values (for example `${WEIGHT} + 200`).
   - Apply the change.
   - `Weight2` recomputes with the new formula.
   - `Weight3` recomputes automatically (it depends on `Weight2`).
   - `Weight4` recomputes automatically (it depends on `Weight3`).
   - No calculation errors or unexpected values appear in the
     dependent columns; the dependency hierarchy `Weight2 → Weight3
     → Weight4` recalculates in the correct order.
3. **Modify the `Weight3` formula and apply.**
   - Open the `Formula` info panel for `Weight3`.
   - Change the formula to a different valid expression (for example
     `${Weight2} + 50`).
   - Apply the change.
   - `Weight3` recomputes with the new formula.
   - `Weight4` recomputes automatically based on the new `Weight3`.
   - `Weight2` is unaffected (it is upstream of `Weight3`).
   - No errors are surfaced; downstream-only propagation honored.
4. **Modify the `Weight4` formula and apply.**
   - Open the `Formula` info panel for `Weight4`.
   - Change the formula to a different valid expression (for example
     `Log10(${Weight3}) - 0.1`).
   - Apply the change.
   - `Weight4` recomputes with the new formula.
   - `Weight2` and `Weight3` are unaffected (terminal node in the
     dependency chain).
   - No errors are surfaced.

### Save the project, reopen it, verify chained columns persist (GROK-17109 regression surface)

1. **Save the current view as a project.**
   - Save the Demog view as a new project containing the
     `Weight2` / `Weight3` / `Weight4` calculated columns; enable
     data sync where the project save dialog offers it.
2. **Close the project.** Close the table view / project so the
   working state is cleared.
3. **Reopen the saved project.** Open the project just saved.
4. **Verify all three chained calculated columns persist.**
   - `Weight2`, `Weight3`, and `Weight4` are present in the dataset
     upon reopen.
   - Each column carries its formula tag (a column flagged by
     `powerpack.formula.is-formula-column` exposes the `Formula`
     info panel widget when selected).
   - Values match the last-applied formula state (i.e. the values
     reflect the most recent formula edits from the previous
     scenario, not the original formulas from the build step).
   - No errors are surfaced on reopen; the dependency chain
     `Weight2 → Weight3 → Weight4` is intact and recomputes
     correctly.

## Notes

- **Related bug.** GROK-17109: calculated columns were not saved to
  the project when using data sync (fixed in 1.23.0). The Save +
  Reopen scenario block above is the regression check for this bug.
