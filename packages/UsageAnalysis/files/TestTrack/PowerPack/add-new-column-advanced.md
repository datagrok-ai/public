---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column
  - powerpack.dialogs.add-new-column-func
  - powerpack.dialogs.prepare-add-column-call
  - powerpack.formula.is-formula-column
  - powerpack.dialogs
target_layer: playwright
coverage_type: regression
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
scope_reductions: []
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
    scope_reduction_proposal: |
      Independent re-derivation for cycle 2026-05-28-powerpack-automate-02
      reaches SCOPE_REDUCTION on the same axis prior cycles found: the
      scenario "## Setup" step 1 enumerates five distinct source classes
      walked in sequence (scenario body: Source 1 local storage, Source 2
      Home dir, Source 3 Postgres:Northwind:OrdersByEmployee, Source 4
      Northwind products GetTop100, Source 5 Northwind products GetAll)
      and prescribes "the same Scenarios block is executed end-to-end"
      for each source. The spec walks only the file source
      (System:DemoFiles/demog.csv) for the full 10-step chain (softStep
      "Setup: open System:DemoFiles/demog.csv with datasync provenance"
      at spec L167 via openTableFromFile) and documents the four
      un-walked sources in its leading block-comment "Scope reductions"
      section (spec L74-105) with rationale: the 5x10 step matrix (50
      combinations) is impractical in a single ~600s Playwright test,
      and the GROK-17109 invariant (calculated columns persist across
      save+datasync+reopen with formula tags intact) is
      source-class-independent at the calc-column-persistence layer
      (the bug surfaces on the add-new-column dialog -> datasync save
      path, not on a specific source-class binding).

      Scenario frontmatter still records `scope_reductions: []` (L34) —
      the narrowing is NOT declared at the scenario layer; it lives only
      in the spec's leading block comment. Independent application of
      the Gate E checklist routes this to SCOPE_REDUCTION rather than
      FAIL: E-TRACE-02 PASSES because all 10 numbered scenario steps
      ARE covered in the spec (Step 1 L187, Step 2 L202, Step 3 L310,
      Step 4 L327, Step 5 L427, Step 6 L490, Step 7 L507, Step 8 L525,
      Step 9 L577, Step 10 L616). The narrowing lives in Setup step 1's
      source multiplicity (5 -> 1), not in any numbered-step omission.
      Per the spec-mode rubric "SCOPE_REDUCTION: acceptable if ... the
      spec covers fewer assertions than the scenario explicitly because
      of a documented technical limitation" — the limitation (5x10 =
      50 combinations exceeding the 600s playwright runtime bound) is
      documented in the spec's local notes.

      Proposed resolution (preferred 1 + 2 together):

      1. Backfill scenario frontmatter `scope_reductions:` with an SR
         entry `id: SR-01, check: setup-multi-source-walk-narrowed,
         rationale: "Walk 1 representative file source
         (System:DemoFiles/demog.csv) for the full 10-step chain in
         this -spec.ts; un-walked sources (local-storage, Home-dir,
         OrdersByEmployee, GetTop100, GetAll) move to per-source matrix
         specs at a future expansion. GROK-17109 invariant is source-
         class-independent at the calc-column persistence layer.",
         verdict_status: SCOPE_REDUCTION`. Re-emit Gate D and Gate E so
         the SR is declared at the layer where it is materialized.

      2. Surface the per-source matrix at the chain level as a
         `bug_focused_candidates` entry in
         `scenario-chains/powerpack.yaml` (proposed spec name
         `powerpack-add-new-column-multi-source-matrix-spec.ts`) so the
         un-walked sources are tracked as future work, not silently
         dropped. Scenario Notes already gesture at this with the
         GROK-17109 cross-cutting-candidate language.

      3. (NOT recommended.) Extend this spec to walk all five sources
         in a per-source loop. Runtime impact ~3000s — needs
         long-running-test budget allowance and `test.slow()`. Reject
         on runtime grounds.

      Retry-context check (E-RETRY-IGNORES-GATE-B): the scenario's
      `gate_verdicts.b` block (cycle_id 2026-05-26-powerpack-automate-01,
      the most-recent-prior-cycle of this scenario) records
      `verdict: PASS` (3/3 attempts, failure_keys: []). The
      retry-context predicate requires `gate_verdicts.b.verdict ==
      "FAIL"`; the current Gate B verdict is PASS, so retry-context
      detection does NOT fire and E-RETRY-IGNORES-GATE-B never fires.
      (The prior-cycle SR prose narrating a -03 B-STAB FAIL is stale —
      the actual recorded Gate B verdict is now PASS, validated against
      the per-attempt Playwright JSON reporter output.) Independent of
      the predicate: the spec source is unchanged from prior reads and
      retains no diagnosed-failing code path, so even under a
      hypothetical FAIL predicate the ignored-evidence pattern would
      not be present.

      Gate E evaluates spec-vs-scenario traceability and discipline on
      the authored source. The SCOPE_REDUCTION verdict stands on the
      source-matrix narrowing axis (5 -> 1 in the spec leading block
      comment vs empty scenario frontmatter scope_reductions[]); all
      other E-* checks PASS independently.
    claims:
      - check_id: E-STRUCT-MECH-01
        status: PASS
      - check_id: E-STRUCT-MECH-02
        status: PASS
      - check_id: E-STRUCT-MECH-03
        status: PASS
      - check_id: E-STRUCT-MECH-04
        status: PASS
      - check_id: E-STRUCT-MECH-05
        status: PASS
      - check_id: E-STRUCT-MECH-06
        status: PASS
      - check_id: E-TRACE-01
        status: PASS
      - check_id: E-TRACE-02
        status: PASS
      - check_id: E-TRACE-03
        status: PASS
      - check_id: E-SEL-01
        status: PASS
      - check_id: E-SEL-02
        status: PASS
      - check_id: E-SEL-03
        status: PASS
      - check_id: E-HELP-01
        status: PASS
      - check_id: E-HELP-02
        status: PASS
      - check_id: E-LAYER-01
        status: PASS
      - check_id: E-LAYER-02
        status: NA
      - check_id: E-LAYER-COMPLIANCE-01
        status: PASS
      - check_id: E-BOUND-01
        status: PASS
      - check_id: E-BOUND-02
        status: PASS
      - check_id: E-RETRY-IGNORES-GATE-B
        status: NA
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-add-new-column-subdir
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
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

Multi-source matrix-like sequence walking the canonical GROK-17109
reproduction surface: open data from five distinct sources, add two
chained calculated columns, mutate the source column (rename + edit
values), save the project with datasync where the source supports it,
close all, reopen, and verify the calculated columns persist with
formulas intact AND recompute correctly when the source column is
renamed.

Chain witness role: this scenario is bug-focused per chain
`pyramid_layer: bug-focused` and owns the specialty persistence flows
(`save-project-with-datasync`, `project-reopen-with-formula-recalc`,
`column-rename-context-action`) for the PowerPack chain. Basic
dialog-open-and-add flow is delegated to the chain's smoke
witness at `add-new-column.md` (top-level PowerPack scenario,
`ui_coverage_delegated_to: add-new-column.md`).

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

- **Bug focus and chain role.** This scenario is bug-focused per
  chain `pyramid_layer: bug-focused`; it walks the canonical
  GROK-17109 reproduction (Step 6 save-with-datasync + Step 8
  reopen-and-verify) plus the formula-recalc-on-rename invariant
  (Step 9). Cross-cutting bug candidate `GROK-17109` is emitted at
  the chain level (`bug_focused_candidates` in
  `scenario-chains/powerpack.yaml`); proposed spec
  `powerpack-grok-17109-spec.ts` spans this scenario plus
  `add-new-column.md` (top-level smoke) and
  `AddNewColumn/formula-refreshing.md` (dependency-chain recalc).
- **UI delegation.** Basic dialog-open-and-add flow (toolbar icon,
  formula editor, OK button) is owned by the chain's smoke witness
  at `add-new-column.md`. This scenario owns the specialty
  persistence flows: `save-project-with-datasync`,
  `project-reopen-with-formula-recalc`, `column-rename-context-action`.
  See chain `ui_coverage_plan.delegated_scenarios` entry.
- **Sibling spec.** A Playwright spec already exists at
  `public/packages/PowerPack/src/tests/add-new-column.ts` (see
  `existing-test-index.yaml`); house-style anchor for Automator
  when authoring the migrated scenario's `-spec.ts`.
- **Sibling scenario alignment.** The `Weight2 = ${WEIGHT} + 100`,
  `Weight3 = ${Weight2} + 100` chain mirrors the formula pattern
  used by `AddNewColumn/formula-refreshing.md` (Weight2 → Weight3 →
  Weight4) for cross-scenario consistency; that sibling extends the
  chain with `Weight4 = Log10(${Weight3}) - 0.2` and adds Context
  Panel formula-edit coverage which this scenario does not duplicate.
- **Source-text fixes.** The original scenario contained the TODO
  marker `(TODO: specify which formula to use)` on Step 3; the
  migration resolves it to `${WEIGHT} + 100` for Step 3 and
  `${Weight2} + 100` for Step 4 per the chain-level recommendation
  (a) in `unresolved_ambiguities :: formula-not-specified` — the
  resolution mirrors the Weight2/Weight3 chain in the sibling
  `formula-refreshing.md` scenario. The original "Close All"
  shorthand is spelled out as closing all open views via the shell
  (Step 7) so the action is unambiguous when this scenario is
  automated.
- **Candidate helpers.** Several patterns recur and warrant
  registry candidates: opening tables from local storage / Home dir
  / query results (`openTableFromLocalStorage`, `openTableFromHomeDir`,
  `openQueryResult`), saving a project with datasync
  (`saveProjectWithDatasync`), reopening a saved project
  (`reopenProject`), adding a calculated column via the dialog
  (`addCalculatedColumn`), and renaming a column via the grid
  context action (`renameColumnViaContextAction`). None of these
  exist in `helpers-registry.yaml` yet; surfaced here as candidates
  for the registry per the Migrator candidate-helper convention.
- **Datasync availability matrix.** Sources 1 (local storage) and
  2 (Home dir) typically do not support datasync since there is no
  upstream source to sync against; Sources 3-5 (Northwind query
  results) do support datasync. The scenario body's Step 6
  "(where available)" qualifier preserves the original's allowance
  for variable datasync support across sources; the persistence-
  on-reopen invariant (Step 8) applies to all five sources.
- **Original trailing JSON metadata.** The original scenario ended
  with `{"order": 1}`. The `order` field is captured in chain
  `order_from_files` under the path-relative key
  `AddNewColumn/add-new-column.md`.
- **Name-collision awareness.** There is a sibling top-level
  scenario at `PowerPack/add-new-column.md` (the chain's smoke
  witness for the Demog flow). The two scenarios share the same
  basename but are DISTINCT and not interchangeable; the chain
  encodes them with path-relative keys. See chain
  `unresolved_ambiguities :: naming-collision-add-new-column` for
  future rename suggestions.
