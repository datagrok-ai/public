---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column-func
  - powerpack.dialogs.add-new-column
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-dialog
  - add-new-column-toolbar-icon
  - add-new-column-recent-activities
  - add-new-column-autocomplete
  - add-new-column-drag-n-drop-columns
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/add-new-column.md
migration_date: 2026-05-23
source_text_fixes:
  - testtrack-star-icon-replaced-with-explicit-demog-open
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-17109
  - GROK-17004
realized_as:
  - add-new-column-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T00:00:00Z
    review_round: 1
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T09:42:13Z
    spec_runs:
      - spec: add-new-column-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 31
        failure_keys: []
---

# Add New Column — Demog smoke (toolbar + UI sanity + Recent Activities autofill)

Single-source UI smoke for the Add New Column dialog on the Demog
dataset. Exercises the toolbar entry path (`addNewColumnDialog`
editor-for `AddNewColumn`), the dialog's built-in UI sanity (no
overlapping content, resizable, tooltips), the autocomplete hint
mechanic, the column drag-n-drop into the formula field, and the
Recent Activities autofill after the dialog is reopened.

Chain witness role: this scenario is the PowerPack chain's smoke
witness (`smoke_scenario: add-new-column.md`) for the dialog UI
surface; the sibling AddNewColumn/ scenarios delegate their basic
dialog-open-and-add flows here per chain
`ui_coverage_plan.delegated_scenarios`.

## Setup

1. Open `System:DemoFiles/demog.csv` so the Demog table view is the
   active view (replaces the TestTrack "star" / "Open test data"
   icon — same dataset, explicit source per the JSON metadata block
   on the original scenario).

## Scenarios

### Add a column via the dialog and autofill from Recent Activities

1. **Open the Add New Column dialog.** Click the "Add new column"
   toolbar icon on the Demog table view ribbon. A dialog opens.
2. **Verify dialog UI sanity.** Hover across the entire dialog and
   verify:
   - No overlapping text anywhere in the dialog.
   - No unnecessary scrollbars, icons, or list contents extending
     beyond the dialog boundaries.
   - Tooltips on dialog controls render with clear descriptions.
3. **Verify the dialog resizes.** Resize the dialog larger and then
   smaller; the layout adjusts appropriately to both directions
   without clipping or overlap.
4. **Add a new column named "New" with formula
   `Round(${HEIGHT} + ${WEIGHT})`.**
   - Type the column name `New` in the name input.
   - Compose the formula in the editor using BOTH mechanics:
     - **Autocomplete hint.** Start typing `Rou` in the formula
       editor and accept the `Round` function from the
       autocomplete tooltip suggestions.
     - **Column drag-n-drop.** Drag the `HEIGHT` column header
       from the columns panel into the formula editor at the
       argument position; do the same for `WEIGHT` so the
       formula becomes `Round(${HEIGHT} + ${WEIGHT})`.
5. **Click OK.** A new column named `New` is added to the Demog
   dataset.
6. **Reopen the Add New Column dialog.** Click the "Add new
   column" toolbar icon again. The dialog reopens.
7. **Open Recent Activities and select the most recent entry.**
   Locate the Recent Activities icon in the dialog and click it;
   select the most recent activity (the entry corresponding to the
   `New` column added at Step 5).
8. **Verify the form autofills.** The dialog's name input is
   prefilled with `New` and the formula editor is prefilled with
   the formula composed at Step 4 (`Round(${HEIGHT} + ${WEIGHT})`).

## Notes

- **Chain context.** This is the PowerPack chain's
  `smoke_scenario` per
  `scenario-chains/powerpack.yaml :: ui_coverage_plan`. Owns all
  five flows in `ui_coverage_responsibility:` directly; no
  delegation.
- **Related bugs.** `GROK-17109` (calculated columns persist
  across save+datasync+reopen) and `GROK-17004` (paste handler
  crash on complex formulas) both touch
  `powerpack.dialogs.add-new-column-func` +
  `powerpack.dialogs.add-new-column`. This smoke exercises the
  dialog's headline path but does NOT walk the full
  save+datasync+reopen invariant (covered by
  `AddNewColumn/add-new-column.md` + `AddNewColumn/formula-refreshing.md`)
  nor the complex-paste invariant (covered by
  `AddNewColumn/highlight.md`). Cross-cutting candidates are
  emitted at the chain level (`bug_focused_candidates`).
- **Sibling spec.** A Playwright spec already exists at
  `public/packages/PowerPack/src/tests/add-new-column.ts` (see
  `existing-test-index.yaml`); house-style anchor for Automator
  when authoring the migrated scenario's `-spec.ts`.
- **Source-text fix.** The original referenced opening the Demog
  dataset via the TestTrack "star" / "Open test data" icon
  (TestTrack-runner-specific UI). The migrated scenario references
  the dataset explicitly as `System:DemoFiles/demog.csv` per the
  original's trailing JSON metadata block. Same dataset, explicit
  source.
- **Original trailing JSON metadata.** The original scenario ended
  with `{"order": 2, "datasets": ["System:DemoFiles/demog.csv"]}`.
  The `order` field is captured in chain
  `order_from_files`; the `datasets` field is captured by the
  explicit `System:DemoFiles/demog.csv` reference in Setup.
