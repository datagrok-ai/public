# Scatter plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `scatterplot/CLAUDE.md` "Automation
surface", the TestTrack `scatterplot-*.md` and `scatter-plot-ui.md`. The "(no value)" legend click is
in `scatter-plot-legend-no-value.feature`.

## Needs a library change

- **GROK-17282 — an unsupported line formula in the Formula Lines dialog**
  (`scatterplot-regression-and-formula-lines.md`, Scenario 3). Typing `${HEIGHT} - ${WEIGHT} = 1`
  into the Line pane puts `d4-forced-invalid` on its textarea and the class `disabled` on OK
  (`PowerPack/src/dialogs/formula-lines.ts`, `setFormulaValidationResult` / `onFormulaValidation`),
  and the dialog stays open (measured). Two things stop the claim:
  1. `{element} should be disabled` does not see it: `assertions.ts` treats `aria-disabled`,
     `u2-input-disabled`, `d4-disabled`, `d4-filter-disabled`, `d4-menu-item-disabled` and the native
     `disabled` property as disabled, not a `.ui-btn` with the plain `disabled` class (the button
     keeps `enabled` too). A click on it is then swallowed by the command bar under it. Add
     `ui-btn` + `disabled` to the disabled check.
  2. The Line pane's formula editor has no label (placeholder "Formula", host without a name), so
     no `{element}` resolves it. Either a name on the input in PowerPack, or a library step
     `user types the formula {string} into the Formula Lines dialog` beside the existing
     `user adds the formula line …` (`bindings/tiers/viewers/formula-lines.ts`, which already
     reaches `[name="pane-Line"] textarea`).

  ```gherkin
  Scenario: A formula with two columns on its left side is refused, and the fixed one is kept (GROK-17282)
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line" from the open menu
    And user types the formula "${HEIGHT} - ${WEIGHT} = 1" into the Formula Lines dialog
    Then the formula editor in "Formula Lines" dialog should be invalid
    And OK button in "Formula Lines" dialog should be disabled
    And "Formula Lines" dialog should be visible
    When user types the formula "${HEIGHT} = ${WEIGHT} + 1" into the Formula Lines dialog
    Then the formula editor in "Formula Lines" dialog should be valid
    And OK button in "Formula Lines" dialog should be enabled
    When user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then the formula editor in "Formula Lines" dialog should be valid
    And no errors should have been logged
  ```

## Needs a reading

- **A legend entry dimmed by a legend click** (`scatterplot-legend.md`, "the other color entries are
  dimmed"; "no entry dimmed" after the release). `should be selected` reads the clicked entry; nothing
  reports the others' dimmed state. Needs a `dimmed` state on legend items (the class the legend puts
  on an entry it filters out, read through `item {string} of {element} should be {state}`), then
  `"Caucasian" legend item in legend of scatter plot viewer should be dimmed` after the "Asian" click
  and `should not be dimmed` after the second click.
- **Structures on the X axis** (`scatter-plot-ui.md`, "Structures as markers and as an axis", items
  7-8, spgi-100). With a molecule column on X the axis draws structures instead of tick labels;
  the status reports the axis column and range, not what the ticks are drawn as. Needs an
  `x axis labels drawn as` reading (`text` / `structure`, from the axis painter's renderer), claimed
  after the pick and again after a saved layout is applied.
