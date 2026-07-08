---
feature: powerpack
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [powerpack.cp.add-new-column-persists]
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - add-new-column-dialog
  - column-rename-context-action
  - save-project-with-datasync
  - project-reopen-with-formula-recalc
ui_coverage_delegated_to: add-new-column.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/add-new-column.md
migration_date: 2026-05-20
source_text_fixes:
  - step-3-todo-resolved-formula-weight-plus-100-column-weight2
  - step-4-resolved-formula-weight2-plus-100-column-weight3
  - close-all-spelled-out-as-close-all-views-via-shell
candidate_helpers:
  - helpers.powerpack.openTableFromLocalStorage
  - helpers.powerpack.openTableFromHomeDir
  - helpers.powerpack.openQueryResult
  - helpers.powerpack.saveProjectWithDatasync
  - helpers.powerpack.reopenProject
  - helpers.powerpack.addCalculatedColumn
  - helpers.powerpack.renameColumnViaContextAction
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      Save-with-datasync provenance is not wired in the test environment, so
      save-with-datasync degrades to snapshot-only and GROK-17109 cannot be
      exercised here. The add-new-column advanced flow is asserted; the
      datasync-provenance branch is deferred.
    verdict_status: SCOPE_REDUCTION
related_bugs:
  - GROK-17109
realized_as:
  - add-new-column-advanced-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-add-new-column-subdir
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-add-new-column-subdir
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T11:00:00Z
    spec_runs:
      - spec: add-new-column-advanced-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 37
        failure_keys: []
---

# Add New Column — multi-source datasync persistence + formula recalc on column rename

Multi-source regression test for GROK-17109: open data from five
different table sources, add two calculated columns that reference
each other, rename and edit the source column, save the project with
data sync, close and reopen, and verify the calculated columns persist
with their formulas intact and recompute correctly after the rename.

This scenario focuses on the persistence and rename-recalc behavior;
the basic dialog open-and-add flow is covered by `add-new-column.md`.

## Setup

The scenario walks five table sources in sequence; for each source,
the same Scenarios block is executed end-to-end. Setup establishes
which source the cycle is currently using:

1. Pick one source from the matrix below and open the resulting
   table view:
   - **Source 1 — Local storage.** Open a table previously saved to
     local storage from the Browse tree.
   - **Source 2 — Home dir.** Open a table from the user's Home
     directory.
   - **Source 3 — Query result.** Run `Postgres:Northwind:OrdersByEmployee`
     from the Browse tree and open the result table.
   - **Source 4 — DB GetTop100 result.** Run the auto-generated
     `GetTop100` query against the `products` table in
     `Postgres:Northwind` and open the result table.
   - **Source 5 — DB GetAll result.** Run the auto-generated `GetAll`
     query against the `products` table in `Postgres:Northwind` and
     open the result table.
2. Pick a numeric source column on the open table that will be the
   subject of the chained formulas; call it `WEIGHT` in the
   instructions below. For Northwind `products`, `unitprice` is the
   natural choice; substitute the column name as needed per source.

## Scenarios

### Add chained calculated columns, mutate, save with datasync, reopen, verify

1. **Open the Add New Column dialog (first time).** Click the "Add
   new column" toolbar icon on the open table view ribbon. The
   `AddNewColumnDialog` opens.
2. **Add the first calculated column.**
   - Enter the formula `${WEIGHT} + 100` in the formula editor
     (substituting the source column name picked at Setup step 2;
     the example matches the `Weight2`/`Weight3` chain in the
     sibling `formula-refreshing.md` scenario for cross-scenario
     consistency).
   - Set the new column name to `Weight2`.
   - Click OK. A new column `Weight2` is added to the table; values
     are computed as the chosen source column plus 100.
3. **Open the Add New Column dialog (second time).** Click the
   "Add new column" toolbar icon again. The dialog reopens.
4. **Add the second calculated column referencing the first.**
   - Enter the formula `${Weight2} + 100` in the formula editor.
   - Set the new column name to `Weight3`.
   - Click OK. A new column `Weight3` is added; values are computed
     as `Weight2 + 100` (and therefore the chosen source column
     plus 200 transitively).
5. **Mutate the source column in the grid.**
   - Rename the source column header from `WEIGHT` to a new name —
     e.g. `BaseWeight` — using the context menu / Rename column
     action on the column header in the grid.
   - Change one or more cell values in the renamed source column.
   - **Expected result.** `Weight2` updates its formula text to
     reference the new source column name (`${BaseWeight} + 100`)
     and its values recompute accordingly. `Weight3`'s values
     recompute transitively from the updated `Weight2`.
6. **Save the project with datasync (where the source supports it).**
   - Save the current view as a project from the toolbar / Save
     dialog. In the Save Project dialog, enable the "Data sync"
     option for the open table.
   - **Note for matrix.** Some sources may not support datasync
     (e.g. tables opened from local storage have no upstream source
     to sync against). For sources where datasync is unavailable,
     save the project without datasync and note the source in the
     run log; the persistence-on-reopen invariant (Step 8) still
     applies to the calculated columns themselves regardless of
     datasync availability.
7. **Close all views.** Close the open table view (and any other
   views opened during this cycle) so the workspace is clean before
   reopen.
8. **Reopen the saved project.** Open the project saved at Step 6
   from Recent Projects / Dashboards / context menu.
   - **Expected result (GROK-17109 invariant).** Both calculated
     columns `Weight2` and `Weight3` are present in the dataset
     upon reopen, with formula tags preserved and values intact.
     Columns must NOT be missing — this is the canonical GROK-17109
     regression invariant.
9. **Rename the source column inside the formula (post-reopen).**
   - In the grid of the reopened project, rename the source column
     used by the formulas (the column carrying the new name from
     Step 5; if datasync rewrote it back on reopen, the current
     name is whatever the persisted project carries) to a different
     name — e.g. `BaseWeight2`.
   - **Expected result.** The formula on `Weight2` updates to
     reference the new source column name; the formula on `Weight3`
     remains referencing `${Weight2}` and remains valid.
10. **Change values in the source column (post-reopen).**
    - Edit one or more cell values in the source column on the
      reopened project's grid.
    - **Expected result.** Values in `Weight2` recompute according
      to its formula. Values in `Weight3` recompute transitively
      from the updated `Weight2`. The calculated columns track the
      source-column changes accordingly.

**Overall expected result.** The calculated columns are added,
their values change according to changes in the table source
columns AND persist across save-project-with-datasync + close-all
+ reopen with formula tags intact. Renaming the source column at
any point (pre-save or post-reopen) updates the formula text
automatically. This is the GROK-17109 reproduction path; before
the fix in 1.23.0, columns disappeared on reopen.

## Notes

- **Formula naming.** The `Weight2 = ${WEIGHT} + 100`, `Weight3 =
  ${Weight2} + 100` chain mirrors the formula pattern used by
  `formula-refreshing.md`, which extends the chain with a third
  column (`Weight4`) and additional Context Panel formula-edit
  coverage.
- **Datasync availability.** Sources 1 (local storage) and 2 (Home
  dir) typically do not support data sync, since there is no
  upstream source to sync against; Sources 3-5 (Northwind query
  results) do support it. The persistence-on-reopen check (Step 8)
  applies to all five sources regardless of data-sync availability.
