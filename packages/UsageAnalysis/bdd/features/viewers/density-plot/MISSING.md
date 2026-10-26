# Density plot — what the features cannot say yet

Written during the viewer gap review of 23 Oct 2026 against `density_plot/CLAUDE.md` "Automation
surface". What could be said is in `density-plot-bin-selector-escape-tooltip.feature`.

## Needs a step

- **The bin selector is on screen, and gone when Show Bin Selector is off** (`density-plot-ui.md`,
  "Show Bin Selector"). The `bin selector shown` reading is `showBinSelector`, the auto-layout's
  computed result (the root must be over 300×300), not a check of the screen — and the slider is
  hover-revealed, so `true` means the layout keeps it, not that a person sees it this instant
  (`density_plot/CLAUDE.md`, readings table). The viewer already publishes the slider as the
  `bin selector` **part** (the bins range `<input>`, `features/density_plot_status.dart`), but no step
  claims a part's visibility. Needs `the {string} part of {widget} should be visible` / `… hidden`
  in the library (an on-screen check of the part's element: in the document, non-zero box, not
  `display:none` / `visibility:hidden` / `opacity:0`):

  ```gherkin
  Scenario: Show Bin Selector takes the bin selector off the screen and puts it back
    When user hovers over the "densest bin" area of density plot viewer
    Then the "bin selector" part of density plot viewer should be visible
    When user moves the pointer away from density plot viewer
    And user sets "showBinSelector" property of density plot viewer to "false"
    And user hovers over the "densest bin" area of density plot viewer
    Then the "bin selector" part of density plot viewer should be hidden
    When user moves the pointer away from density plot viewer
    And user sets "showBinSelector" property of density plot viewer to "true"
    Then no errors should have been logged
  ```

## Edge cases (`density-plot-edge-cases.feature`)

Added with `density-plot-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **ApiTests to write (core#1048)**: Bins set to 0, -5 or null through `setOptions` (the slider
   cannot produce them) and a viewer root of 0 × 0 must not throw on every frame.
2. **Steps requested**: the bin tooltip's range on a logarithmic axis (core#1048 "Wrong 4") needs the
   tooltip's range as numbers, not text.
3. **Dart test to write**: listener count across table swaps (core#1048 "Leaks").
