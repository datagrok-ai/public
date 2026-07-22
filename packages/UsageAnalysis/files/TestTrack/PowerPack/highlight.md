---
feature: powerpack
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: [add-new-column-paste-complex-formula]
realizes: []
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/highlight.md
migration_date: 2026-05-20
source_text_fixes:
  - text-filed-to-text-field-typo
  - step-renumbering-duplicate-step-3
  - drag-n-drop-phrase-normalization
candidate_helpers:
  - addNewColumn.openDialog
  - addNewColumn.pasteFormula
  - addNewColumn.assertColumnNameHighlighted
unresolved_ambiguities:
  - blue-highlight-css-token-not-specified
  - autocomplete-vs-drag-n-drop-step-5-alternation
scope_reductions: []
related_bugs:
  - GROK-17004
pyramid_layer: bug-focused
realized_as:
  - highlight-spec.ts
ui_coverage_responsibility:
  - add-new-column-column-name-highlight
  - add-new-column-paste-formula
  - add-new-column-drag-n-drop-columns
ui_coverage_delegated_to: add-new-column.md
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-add-new-column-highlight
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T11:26:03Z
    spec_runs:
      - spec: highlight-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 55
        failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-add-new-column-highlight
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
    review_round: 1
---

# Add New Column — Column-Name Highlighting in Formulas (GROK-17004 regression)

Regression test for GROK-17004: when a formula in the Add New Column
dialog references a column (via `${col}` or `$[col]`), the referenced
column name should be highlighted in blue in the formula editor —
whether the formula was pasted, typed with autocomplete, or built by
dragging a column into the editor. Before the fix, pasting a complex
multi-column formula could crash the dialog instead of highlighting
the column references.

## Setup

1. Open the Demog dataset (`System:DemoFiles/demog.csv`) — a TableView with the Demog grid should be the active view.
2. Open the `Add New Column` dialog (toolbar `+`-column icon or `Edit | Add New Column`).

## Scenarios

### Column-name highlight on pasted formula referencing one column

Verifies that pasting a formula referencing a column (`${col}` form) highlights
the column name in blue — this is the code path that crashed under GROK-17004
for complex formulas.

1. Paste the literal text `Abs(${age})` into the expression editor (`text field` in the
   dialog body).
2. The referenced column name `age` inside `${...}` is highlighted in blue.

### Column-name highlight on pasted formula using bracket-reference form

Verifies that the highlight also works for the `$[col]` bracket-reference form,
which is parsed independently from `${col}`.

1. Paste the literal text `Avg($[age])` into the expression editor.
2. The referenced column name `age` inside `$[...]` is highlighted in blue.

### Column-name highlight on column added via autocomplete

Verifies that a column reference inserted via autocomplete (typing `${` or `$[`
triggers the column-suggestion menu) is also highlighted.

1. Type a different function into the expression editor (for example `Round(`).
2. Trigger column autocomplete (type `${` inside the function arguments) and select a
   column from the suggestion menu — the editor inserts a `${<column>}` reference at the
   caret.
3. The inserted column name is highlighted in blue.

### Column-name highlight on column added via drag-and-drop

Verifies the highlight invariant for column references inserted by dragging a column
header from the Demog grid into the expression editor.

1. Type a different function into the expression editor (for example `Sin(`).
2. Drag a column header from the Demog grid into the expression editor — the editor
   inserts a `${<column>}` reference at the drop position.
3. The inserted column name is highlighted in blue.

### Complex multi-column nested-conditional paste (GROK-17004 verbatim repro)

Verifies the highlight against the exact formula that triggered GROK-17004: a
complex, multi-column pharmacokinetic formula with nested conditionals. This is
the paste that crashed the dialog with `TypeError: Cannot read properties of
undefined (reading 'to')` on builds before 1.23.0.

Bug-invariant assertions:
- The paste handler does NOT throw (no `TypeError` / uncaught exception in console).
- `cm-column-name` spans appear for each `${col}` reference in the formula — five
  distinct referenced columns (`Whole blood assay 1`, `Route Admin`, `Chemical Space X`,
  `Average Mass`, `Species`).
- At least one rendered span has a blue computed color (the
  `getColumnNamesAndSelections → addColHighlight` pipeline completed successfully).

1. Switch the active dataset to SPGI by opening `System:DemoFiles/chem/SPGI.csv` — a
   TableView with the SPGI grid should become the active view (SPGI carries the
   `Whole blood assay 1`, `Route Admin`, `Chemical Space X`, `Average Mass`, and
   `Species` columns referenced by the GROK-17004 repro formula).
2. Open the `Add New Column` dialog against the SPGI grid (toolbar `+`-column icon
   or `Edit | Add New Column`).
3. Paste the GROK-17004 repro formula verbatim into the expression editor as a
   single paste action:

   ```
   if(${Whole blood assay 1} != null, ${Whole blood assay 1}, if(${Route Admin}=="PO", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}=="Dog", 30.9, if(${Species}=="Monkey", 43.6, if(${Species}=="Minipig", 39, null)))))*100
   ```

4. **Verify** the paste handler completes without throwing — no `TypeError: Cannot
   read properties of undefined (reading 'to')` or other uncaught exception appears
   in the browser console during or after the paste.
5. **Verify** the editor renders at least one `.cm-column-name` highlight span. The
   five referenced column names (`Whole blood assay 1`, `Route Admin`, `Chemical
   Space X`, `Average Mass`, `Species`) are each highlighted in blue at every
   `${col}` reference site.

## Notes

- This scenario covers the GROK-17004 regression: scenarios 1-4 exercise the
  highlight-on-paste / autocomplete / drag-and-drop code path with simple
  single-column formulas on Demog; scenario 5 replays the exact multi-column
  formula from the original bug report on SPGI.
