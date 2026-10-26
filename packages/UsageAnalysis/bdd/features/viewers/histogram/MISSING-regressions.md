# Histogram — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed histogram issues into `histogram-regressions.feature`.

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
