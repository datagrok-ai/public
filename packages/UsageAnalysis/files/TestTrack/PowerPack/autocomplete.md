---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column-func
  - powerpack.dialogs.add-new-column
target_layer: playwright
coverage_type: smoke
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
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML and carries all required fields:
          feature=powerpack, sub_features_covered=[powerpack.dialogs.add-new-column-func,
          powerpack.dialogs.add-new-column], target_layer=playwright,
          coverage_type=smoke. Both sub_features_covered entries resolve to
          atlas sub_features (lines 677 and 709 of powerpack.yaml).
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains "## Setup" and "## Scenarios" level-2 headings plus
          three "### Autocomplete: ..." scenario sub-headings under
          ## Scenarios. Scenario headings are present.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Each of the three scenario sub-sections has numbered steps
          starting with "1." (type-triggered tooltip — 6 steps;
          Ctrl+Space — 2 steps; dollar column-suggestion — 2 steps).
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          No empty scenarios. Each of the three scenarios contains at
          least two numbered steps with concrete actions and a Verify
          assertion.
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer=playwright is one of the allowed values
          {playwright, apitest, manual-only}.
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type=smoke is one of the allowed values
          {smoke, regression, edge, perf}.
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          File-level coverage_type=smoke label is present in frontmatter
          and applies to all three scenarios in the body. No severity-axis
          (p0..p3) value present.
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          Common setup (opening demog.csv and the Add new column dialog)
          is factored into a single "## Setup" section. The three scenario
          bodies do not duplicate the setup; they begin from the focused
          editor state established by Setup.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          pyramid_layer=ui-smoke and coverage_type=smoke are aligned per
          the hard rule (ui-smoke requires smoke). Orchestrator pre-emptive
          correction from regression to smoke restored alignment before
          Gate A entry.
      - check: A-CONT-01
        status: PASS
        evidence: |
          Scenario body uses real names throughout: dataset demog.csv
          from System:DemoFiles; columns HEIGHT, WEIGHT, AGE; functions
          Abs, Acos, Avg; selectors [name="icon-add-new-column"],
          .d4-dialog .cm-content, .cm-tooltip-autocomplete. No
          angle-bracket placeholders, square-bracket placeholders, or
          generic stand-ins detected.
      - check: A-BUG-01
        status: PASS
        evidence: |
          Atlas known_issues for powerpack is [] (line 1224 of
          powerpack.yaml — "Imported from bug-library/powerpack.yaml when
          curated. No bug-library file exists for powerpack yet").
          PASS-by-vacuity. Body Notes references GROK-17004 as a
          cross-cutting bug-focused candidate emitted at chain level —
          this is informational and does not affect A-BUG-01 since atlas
          has no needed entries.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          No scenarios opted out for effort or complexity. The scenario
          set covers three distinct autocomplete mechanics (type-triggered
          function tooltip, Ctrl+Space invocation, dollar column-suggestion)
          with explicit completion-action variants (Enter and mouse click).
      - check: A-MERIT-02
        status: PASS
        evidence: |
          No "TODO: add later" or deferred-to-next-phase markers in the
          scenario body or frontmatter. The Notes section explains current
          scope-split decisions (delegation to add-new-column.md for the
          basic dialog surface; GROK-17004 at chain level) with explicit
          rationale citing scenario-chains/powerpack.yaml rev 1.
---

# Add New Column — autocomplete (demog)

Short focused UI scenario for the autocomplete mechanics of the
**Add New Column** dialog: type-triggered function tooltip,
**Ctrl+Space** explicit invocation, and **$** column-suggestion
tooltip. Delegates the basic dialog UI surface (open / close,
preview grid, OK/Cancel) to the chain smoke `add-new-column.md`
per `scenario-chains/powerpack.yaml` rev 1
`ui_coverage_plan.delegated_scenarios`.

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

- **Chain context.** `scenario-chains/powerpack.yaml` rev 1
  classifies this scenario as `pyramid_layer: ui-smoke` /
  `classification: simple`. The chain's elected smoke is the
  top-level `add-new-column.md` (Rule 1 edge case — shortest
  step-count + broader sub-feature coverage); this scenario
  owns the specialty autocomplete mechanics on its
  `ui_coverage_responsibility:` list
  (`add-new-column-autocomplete`,
  `add-new-column-ctrl-space`,
  `add-new-column-dollar-column-suggestions`).
- **Related bug GROK-17004** (paste-handler crash on complex
  formulas) touches the dialog's autocomplete / paste surface
  but is NOT directly reproduced here — this scenario tests
  autocomplete on **typing**, not paste. The cross-cutting
  bug-focused candidate for GROK-17004 is emitted at chain
  level (`bug_focused_candidates`) spanning `highlight.md`
  + the top-level `add-new-column.md`; frontmatter
  `related_bugs: []` for this scenario reflects that
  scope-split.
- **Sibling spec.** A Playwright spec already exists at
  `autocomplete-spec.ts` (same directory). Its test name is
  `AddNewColumn: autocomplete (demog)`. It drives the dialog
  via the toolbar icon (`[name="icon-add-new-column"]`),
  uses CodeMirror selectors (`.d4-dialog .cm-content`,
  `.cm-tooltip-autocomplete`), and asserts columns
  `HEIGHT`, `WEIGHT`, `AGE` for the `$` case — patterns to
  reuse if the Automator regenerates the spec.
- **Source dataset** — `demog.csv` is platform-provided
  (`System:DemoFiles`); no fixture is produced or consumed by
  any other scenario in this section.

---
{
  "order": 2
}
