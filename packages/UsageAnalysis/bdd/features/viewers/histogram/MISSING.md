# Histogram — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against the histogram's "Automation surface" (`core/client/d4/lib/src/viewers/histogram/CLAUDE.md`)
and the TestTrack checklists (`files/TestTrack/Viewers/Histogram/`). What could be said went into
`histogram-selectors-and-mouse-over.feature` (Show Split Selector, the bins slider against the Bins
property — github-2296, X Axis Height, Show Mouse Over Row) and exact counts in
`histogram-range-filter.feature`.

## Still open

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
