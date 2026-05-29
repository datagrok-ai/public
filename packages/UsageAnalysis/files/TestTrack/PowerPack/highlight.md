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
---

## Setup

1. Open the Demog dataset (`System:DemoFiles/demog.csv`) — a TableView with the Demog grid should be the active view.
2. Open the `Add New Column` dialog (toolbar `+`-column icon or `Edit | Add New Column`).

## Scenarios

### Column-name highlight on pasted formula referencing one column

Verifies the paste-handler highlight invariant (GROK-17004 regression surface) for the
canonical `${col}` reference form. The pasted-formula path is the entry point at which
GROK-17004's `getColumnNamesAndSelections` crash manifested for complex pIC50-style
formulas; the simple-paste case here exercises the same highlight-on-paste codepath.

1. Paste the literal text `Abs(${age})` into the expression editor (`text field` in the
   dialog body).
2. The referenced column name `age` inside `${...}` is highlighted in blue.

### Column-name highlight on pasted formula using bracket-reference form

Verifies that the highlight invariant also covers the `$[col]` bracket-reference form.
The two reference forms are independently parsed by the column-resolution path, and the
bug surface includes both.

1. Paste the literal text `Avg($[age])` into the expression editor.
2. The referenced column name `age` inside `$[...]` is highlighted in blue.

### Column-name highlight on column added via autocomplete

Verifies the highlight invariant for column references inserted through the in-editor
autocomplete path (typed `${` or `$[` triggers the column-suggestion menu).

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

Verifies the highlight invariant against the exact GROK-17004 bug-repro paste —
a multi-column pIC50-style pharmacokinetic formula with nested `if` / `Contains` /
species discrimination / dose-normalization arithmetic. This is the canonical
NIBR-workflow paste that triggered `TypeError: Cannot read properties of undefined
(reading 'to')` on `add-new-column.ts:611:67` (in `getColumnNamesAndSelections`)
on pre-1.23.0 builds (`bug-library/powerpack.yaml:115-156`).

Bug-invariant assertions:
- The paste handler does NOT throw (no `TypeError` / uncaught exception in console).
- `cm-column-name` spans appear for each `${col}` reference in the formula — five
  distinct referenced columns (`Whole blood assay 1`, `Route Admin`, `Chemical Space X`,
  `Average Mass`, `Species`).
- At least one rendered span has a blue computed color (the
  `getColumnNamesAndSelections → addColHighlight` pipeline completed successfully).

1. Switch the active dataset to SPGI by opening `System:DemoFiles/SPGI.csv` — a
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

- Bug coverage: this scenario walks the GROK-17004 regression surface (column-name
  highlight on paste). Scenarios 1-4 exercise the highlight-on-paste / autocomplete /
  drag-n-drop invariant on Demog with simple single-column references — sufficient
  to detect any regression in the `updateListener → getColumnNamesAndSelections →
  addColHighlight` pipeline that GROK-17004 broke. Scenario 5 (added 2026-05-28)
  exercises the verbatim GROK-17004 bug-repro paste on SPGI: a multi-column nested-
  conditional pharmacokinetic formula that triggered the original `TypeError` on
  pre-1.23.0 builds. Together the two halves close the `test_coverage: needed` slot
  on `bug-library/powerpack.yaml :: GROK-17004` once Validator confirms Scenario 5
  PASS — Gate F (coverage adjudication) is the authoritative gate for that flip.
- UI coverage delegation: basic dialog-open and toolbar-icon flows are owned by the
  section smoke `AddNewColumn/add-new-column.md` per the chain
  `ui_coverage_plan.smoke_covers`. This scenario owns the paste-handler column-name-
  highlight specialty flow and the per-reference-form coverage on its
  `ui_coverage_responsibility[]`.
- Highlight-color token: the original TestTrack scenario specifies "blue" without a
  CSS token. The migrated assertions remain at semantic "highlighted in blue" level;
  downstream Automator may resolve this to a concrete CSS selector / computed-style
  check against the editor's highlight class. Surfaced in `unresolved_ambiguities`.
- Helpers used by name in this scenario do not yet exist in
  `helpers-registry.yaml`; surfaced as `candidate_helpers` for the registry curator.
