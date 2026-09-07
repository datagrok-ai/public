---
feature: chem
realizes_atlas:
  - chem.int.empty-input-analyses
realizes:
  - chem.analyze.r-groups-analysis
  - chem.analyze.chemical-space
priority: p1
target_layer: playwright
coverage_type: edge
pyramid_layer: bug-focused
produced_from: atlas-driven
realized_as:
  - chem-empty-input-analyses-spec.ts
related_bugs:
  - id: GROK-16329
    status: fixed
  - id: GROK-18407
    status: open
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    check: E-EXPECT-COVERAGE-01
    rationale: |
      The empty-RESULT branch of R-Groups Analysis is NOT covered. There are two
      distinct outcomes and this scenario, being the empty-INPUT boundary,
      exercises only the first: (a) empty input — no molecules at all — where MCS
      derives no core and rGroupDecomp short-circuits with "No core was provided"
      before decomposing; (b) real molecules that share no common substructure,
      where the decomposition genuinely runs and returns an empty result, which
      is the branch that emits "No R-Groups were found"
      (analysis/r-group-analysis.ts#L181-189, reachable only when the call
      returned a result object). Branch (b) needs a fixture of molecules with no
      shared core plus the Visual analysis checkbox set, and is a separate
      scenario. Measured on dev 2026-08-21: 4/4 attempts on an all-empty column
      produced "No core was provided" with OK enabled, confirming (b) is
      unreachable from this input.
  - id: SR-02
    check: E-EXPECT-COVERAGE-01
    rationale: |
      No cell-tooltip behaviour is covered by this scenario at all. S1.5 drives
      the cell click and asserts only what the click does to product state: the
      cell the grid's own hit test resolved becomes current (current row = the
      clicked row, current column = structure), from a pre-click discriminator
      confirming the current cell was a different one. That is the scenario's
      subject — the app stays usable after a rejected analysis — and a frozen
      grid cannot produce it.
      The tooltip is not merely unasserted, it is unreachable on this fixture by
      design: the cell-tooltip column filter (grid_core.dart:357-359) keeps a
      column only when it is NOT already visible in the grid, unless
      showVisibleColumnsInTooltip is set, and that defaults to false
      (grid_look.dart:299). Both of this fixture's columns are visible, so every
      column is filtered out, the row tooltip has no rows, and tooltip.dart:449
      returns before _show — TOOLTIP_SHOWN (:365) never fires. Gate B confirmed
      this 3/3 on 2026-08-22: host present but empty, no row table, no text.
      An earlier revision of this SR planned to let a live run supply a measured
      .d4-tooltip-empty-value count and then decide whether the empty-cell text
      could be asserted. That plan is WITHDRAWN, not deferred: the count it would
      have consumed was read out of a tooltip that never rendered, so it measured
      nothing about the empty-cell rendering branch. Any future coverage of that
      branch needs a fixture built for it — a column hidden from the grid, or
      showVisibleColumnsInTooltip set, plus the isOnlyQualityColumn preconditions
      (quality tag, a registered HtmlCellRenderer for it, no second column sharing
      the tag; utils.dart:1178-1181) — and belongs to a tooltip scenario, not to
      this empty-input boundary.
gate_verdicts:
  b:
    verdict: KNOWN_BUG_REPRODUCED
    cycle_id: direct-gate-b-2026-08-22-chem-empty-input-analyses-r4
    timestamp: 2026-08-22T16:15:10Z
    spec_runs:
      - spec: chem-empty-input-analyses-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 15
        failure_keys: []
        run_mode: headless-cold
  a:
    verdict: PASS
    cycle_id: 2026-08-20-chem-new-15
    timestamp: 2026-08-20T00:00:00Z
    failure_keys: []
    review_round: 1
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-empty-input-analyses-r3
    timestamp: 2026-08-22T16:20:04Z
    failure_keys: []
expected_results:
  - anchor: "S1.4"
    expectation: >-
      Run against a molecular column whose every cell is empty, so the input is
      empty by construction rather than by a guess about which core a populated
      dataset contains. Empty input must be REJECTED VISIBLY: a "No core was
      provided" balloon appears, no null-reference console error fires, and the
      dataset stays fully unfiltered (filter.trueCount === rowCount) with no
      R-group columns and no Trellis plot. The failure this guards against is a
      silent dialog close or a null crash, not a particular wording. Measured on
      dev 2026-08-21, 4/4 attempts: pressing MCS on an all-empty column derives
      no core, OK stays ENABLED, and confirming raises "No core was provided".
      The sibling message "No R-Groups were found" belongs to a DIFFERENT branch
      — a decomposition that ran over real molecules and returned nothing — and
      is unreachable from empty input, because rGroupDecomp short-circuits on an
      empty core before decomposing. That empty-RESULT branch is not covered by
      this scenario; see scope_reductions.
  - anchor: "S1.5"
    expectation: >-
      No R-group result columns are appended, no Trellis plot is created, and
      the grid remains visible and interactive. Interactivity is read from the
      grid's response to a real click, not from visibility alone: a `structure`
      cell several rows down is located through the grid's own hit test and
      clicked, after which that cell is current (current row = the clicked row,
      current column = `structure`). No cell tooltip is expected — with both of
      this fixture's columns visible in the grid, the product deliberately shows
      none (grid_core.dart:357-359, grid_look.dart:299, tooltip.dart:449).
  - anchor: "S2.3-5"
    expectation: >-
      Empty-column Chemical Space is rejected — no Embed_X embedding column is
      appended (GROK-18407, currently open) — and no JavaScript crash fires.
---

# Chem — Empty Input Boundary: R-Groups Analysis and Chemical Space

## Setup

1. Open Datagrok in a clean browser session and sign in.
2. Build a small table whose molecular column is **entirely empty**: ten rows with
   an `id` column and a `structure` column in which every cell is blank, with
   `structure` marked as the Molecule semantic type. Confirm the Grid viewer is
   visible and that every `structure` cell is empty.

   This input is what makes the decomposition empty **by construction**. Running
   the same steps against a populated file such as `smiles.csv` leaves the outcome
   dependent on whether a core happens to match that data — on `smiles.csv` one
   does, so the decomposition succeeds, no "no R-groups" message appears, and the
   check fails on correct product behaviour.

## Scenarios

### Scenario 1: R-Groups Analysis with MCS strategy on an all-empty molecular column

Steps:

1. With the all-empty table open, open the top menu and navigate to
   Chem > Analyze > R-Groups Analysis. The R-Groups Analysis dialog opens.
2. In the dialog, locate the strategy selector and choose **MCS** (Maximum Common
   Substructure). Leave the molecule column set to the detected Molecule column.
   Do not enter a core structure manually — the MCS strategy derives the core
   from the data.
3. Click **OK** to start the decomposition.
4. Observe the UI response. Empty input must be rejected visibly rather than
   closing the dialog silently: a balloon reading **"No core was provided"**
   appears, because with no molecules the MCS strategy derives no core and the
   decomposition short-circuits before it runs. Verify that no JavaScript console
   error is logged in the browser developer tools, and that the grid row count is
   unchanged from the pre-analysis state (all rows still visible, no filtering
   applied by the rejected run).
5. Dismiss the message or balloon. Verify that the table's column list does not
   contain any newly appended R-group columns, and that the application remains
   interactive — click a `structure` cell a few rows down and confirm the grid
   still handles input: that cell becomes the current cell (the current row moves
   to the clicked row and the current column becomes `structure`).

   Do not expect a cell tooltip on this table. Datagrok deliberately shows none
   when every column is already visible in the grid: the cell-tooltip column
   filter keeps a column only when it is NOT visible, unless the grid's
   `showVisibleColumnsInTooltip` property is set, and that defaults to false
   (`grid_core.dart:357-359`, `grid_look.dart:299`). This fixture has two columns
   and both are visible, so no tooltip row is left to render and the tooltip is
   never shown (`tooltip.dart:449`). A tooltip here would be a defect, not a
   requirement.

### Scenario 2: Chemical Space dialog rejects empty-column input without silently closing

Steps:

1. Open a second dataset: navigate to Browse > Demo Files > Chem > smiles.csv
   (or reuse the already-open table). Alternatively open a table that contains
   a column where all molecule cells are empty strings.
2. Navigate to top menu Chem > Analyze > Chemical Space. The Chemical Space
   editor dialog opens. Verify the dialog title is present (contains "Chemical
   Space" or "Chem Space") and the molecule-column selector is visible.
3. In the column selector, choose a column whose cells are all empty (if the
   test dataset has one), or clear the column selection entirely so the
   selected column field is blank or shows a placeholder. If the dialog pre-
   selects a valid column, use the selector dropdown to switch to a column with
   no molecule data.
4. Click **Run** (or **OK**) without entering any valid column selection.
5. Observe the dialog response. The dialog must NOT close silently with no
   feedback. One of the following must be true: (a) the Run button is disabled
   and cannot be clicked when the column is blank or invalid; (b) an inline
   validation message appears within the open dialog identifying the problem;
   (c) a balloon notification appears immediately after the click explaining the
   rejection. Confirm no Embed_X_N scatter-plot or embedding column is appended
   to the table, and no JavaScript console error fires.

## Notes

- target_layer rationale: both scenarios require multi-step UI interaction
  with Chem top-menu dialogs and DOM-level assertion of balloon/notification
  content. Playwright handles dialog open/close detection and balloon DOM
  queries natively.
- GROK-16329 (fixed in 1.20.0): R-Group MCS null crash when decomposition
  returns empty — Scenario 1 locks in the regression guard. The bug is marked
  fixed; this scenario detects any re-introduction.
- GROK-18407 (regression-risk, no fixed_in): Chemical Space silent-close on
  empty-column input — Scenario 2 locks in the input-validation requirement.
  Run button disabled or inline error are both acceptable implementations;
  the assertion fails only on a silent close with zero UI feedback.
- The chain's existing `r-group-analysis.md` covers the GROK-16329 repro path on
  `smiles.csv`. Note that the two scenarios exercise DIFFERENT branches: this one
  is empty INPUT ("No core was provided", measured on dev 2026-08-21), while the
  empty-RESULT branch that emits "No R Groups were found" requires real molecules
  sharing no common core and is not reached from here — see scope_reductions SR-01.
  This scenario is authored for the interaction `chem.int.empty-input-analyses`,
  which names BOTH the R-Groups empty-result AND the Chemical Space empty-column
  invariant together as a shared empty-input boundary class — the two scenarios
  here close that combined id. The existing `r-group-analysis.md` does not carry
  `realizes_atlas: [chem.int.empty-input-analyses]` in its frontmatter (it was
  migrated before the interaction was added to the atlas in rev 7); authoring a
  new scenario for the full interaction id is therefore the correct coverage
  move rather than extending the migrated file.
- The existing `chemical-space.md` explicitly does not exercise the empty-column
  edge case (chain analysis note: "Body does NOT test empty-column-input edge
  case").
- `derived_from` citation: atlas entry `chem.int.empty-input-analyses` was added
  in atlas rev 7 (2026-08-19) from bug-library GROK-16329 and GROK-18407
  patterns.

## Automation notes

- For Scenario 1 balloon detection, query for the Datagrok notification/balloon
  container (class pattern `d4-balloon` or `ui-notification`) and assert its
  text content matches a "no R groups" pattern before dismissal.
- For Scenario 2 Run-button-disabled branch: `expect(runButton).toBeDisabled()`
  is a valid assertion path. For the inline-message branch: assert a DOM element
  with a non-empty text describing the column validation error is visible within
  the dialog before clicking away.
- Do not fabricate an "empty" column by injecting data via JS API — the test
  must drive the column selector through the UI dropdown so the dialog's own
  validation gate is exercised.
