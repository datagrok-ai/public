# Histogram — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against the histogram's "Automation surface" (`core/client/d4/lib/src/viewers/histogram/CLAUDE.md`)
and the TestTrack checklists (`files/TestTrack/Viewers/Histogram/`). What could be said went into
`histogram-selectors-and-mouse-over.feature` (Show Split Selector, Show Mouse Over Row), the bins
slider against the Bins property (github-2296) in `histogram.feature`, and exact counts in
`histogram-range-filter.feature`.

## Still open

- **X Axis Height (core#1086).** With Auto Layout off the layout takes the value
  (`histogram_core.dart` cuts `look._calcXAxisHeight` off the bottom): the bars get shorter, which
  `histogram-selectors-and-mouse-over.feature` claims. The drawn axis strip itself looks the same for
  every value, with Auto Layout on and off, and nothing reports it: that half needs a reading of the
  drawn axis, e.g. `x axis strip height` (the height of the band the x axis labels and ticks were last painted in, from the axis
  render in `histogram_core.dart`). Viewer gaps round, 23 Oct 2026. Then:

  ```gherkin
  Scenario: With Auto Layout off, X Axis Height sets the height of the drawn x axis strip
    When user sets "Auto Layout" property of histogram viewer to "false"
    And user sets "X Axis Height" property of histogram viewer to "60"
    Then the "x axis strip height" reading of histogram viewer should be 60
    When user sets "X Axis Height" property of histogram viewer to "30"
    Then the "x axis strip height" reading of histogram viewer should be 30
    When user sets "Auto Layout" property of histogram viewer to "true"
    Then no errors should have been logged
  ```

- **The Show Values counts follow the value column (GROK-19759).** `histogram-setup-tune-select.md`
  step 8: with Show Values on, setting Value to WEIGHT must redraw the per-bin counts for the WEIGHT
  distribution. The status reports a `bin labels` area (the strip Show Values reserves) but not the
  text drawn in it, so the only claim today is "the strip repainted" — which any repaint satisfies,
  stale counts included. Needs a reading of the drawn labels, e.g. `bin label <n>` (the count text the
  last frame drew above bin n, from the Show Values loop in `histogram_core.dart`), absent while Show
  Values is off. Then:

  ```gherkin
  Scenario: Show Values redraws the counts for the new value column (GROK-19759)
    Given user sets "Show Values" property of histogram viewer to "true"
    And the "bin label 8" reading of histogram viewer should be "99"
    When user sets "Value" property of histogram viewer to "WEIGHT"
    Then the "bin label 8" reading of histogram viewer should not be "99"
    And the "bin label 8" reading of histogram viewer should be "<count of WEIGHT bin 8, measured>"
    And no errors should have been logged
  ```

## Edge cases

The edge-case round of 2026-10-26 (localhost). Its checks are now steps of the journeys: `histogram-range-filter.feature`. Only the features were written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Steps requested**: the range inputs' own validation. Min equal to max (50 and 50) is flagged by
   the input validator ("Max should be greater than min", `histogram_core.dart`), while the viewer's
   `error` reading stays empty (measured on localhost; the viewer reports only max below min). No step
   reads the inputs' invalid state.
2. **Decided**: Bins typed as 201 is not taken (the out-of-bounds rule): a known failure (core#1200)
   in `histogram-range-filter.feature`, its last scenario.
3. **Decided (2026-10-30)**, claimed in `histogram-range-filter.feature`: `abc` in a range input is not
   taken and the filter stays; switching Value while a range filters releases the range. Bin labels
   dropped above 100 bins are right as they are (no scenario: no reading of the axis labels).
4. **Not written yet**: the selectors hidden at 200 × 200 and shown at 201 × 201. The range checks
   for a date value column (`STARTED`, date inputs) are in `histogram-range-filter.feature`: a maximum
   below the minimum is reported, a bound past the dates is clamped, a range normalized to the
   filter logs no error with every row filtered out (green, 2026-10-29).

### API or Dart tests to write

- A layout whose Value column (HEIGHT) the table no longer has, loaded back, leaves the histogram
  with no bin and the message `value: Column "HEIGHT" does not exist.`, and logs no error.
- With every row filtered out the histogram shows 0 rows and reports no error; once the filter is
  reset it shows 1000 rows again.
- Bins set to 201 or 0 keeps 20: a known failure, core#1200, still open.

Moved out of BDD (no user gesture); requested from ApiTests.

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed histogram issues into scenarios, now in the histogram's own journeys.

- **The distribution lines follow a filter made elsewhere (core#998).** The lines are computed in
  `render` over the rows shown and cached by `combinedFilter` (`histogram/CLAUDE.md` "Split"); the bug
  left them on the old rows after any filter but the histogram's own. A repaint claim cannot tell the
  lines from the bars, which move with the same filter. Needs a reading of what a line was computed
  from — `distribution line <category> mean` (and `… stdev`) — then:

  ```gherkin
  Scenario: The distribution lines follow a filter made in another place (core#998)
    Given user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value                   | AGE  |
      | Split                   | SEX  |
      | Show Distribution Lines | true |
    When user remembers the "distribution line F mean" reading of histogram viewer
    And user filters rows where "AGE" is between 18 and 40
    Then the "distribution line F mean" reading of histogram viewer should be lower than remembered
    When user resets the filter
    Then the "distribution line F mean" reading of histogram viewer should be as remembered
    And no errors should have been logged
  ```
