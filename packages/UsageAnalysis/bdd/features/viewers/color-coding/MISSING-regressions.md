# Color coding — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed colour coding issues into `color-coding-regressions.feature`.

### Regressions found 26 Oct 2026

- **A colour-coded category column paints the bar chart's category names (public#3621).** Measured twice
  on localhost: on demog-1000 with RACE coloured categorically (Caucasian #FF0000, Other #00FF00) and
  a bar chart split by RACE, the `y axis` area (the category names) is painted #FF0000 — 7009 px of it
  in the 75×200 area, every one of the four names on a red background, the green "Other" included,
  while the grid cells are coloured as set. `_renderLabel` in `bar_chart_core.dart` sets the style's
  `backColor` to transparent, but the category renderer it calls still paints the cell. Expected (the
  issue): no colouring on the category names.

  ```gherkin
  @realizes:public#3621
  Scenario: A colour-coded category column does not paint the bar chart's category names (public#3621)
    Given user is logged in
    And user opens demog-1000 dataset
    When user colors "RACE" column categorically:
      | Caucasian | #FF0000 |
      | Other     | #00FF00 |
    Then the "cell 1 of RACE" area of grid should contain the color "#FF0000"
    When user adds a bar chart viewer with:
      | Split | RACE |
    Then bar chart viewer should have a "y axis" area
    And the "y axis" area of bar chart viewer should not contain the color "#FF0000"
    And the "y axis" area of bar chart viewer should not contain the color "#00FF00"
    And no errors should have been logged
  ```

### Questions for owners

- **Color Coding > Apply Coloring in another table (public#1887's note).** Pick Up / Apply > Pick Up on
  demog-1000's RACE, then Pick Up / Apply > Apply on demog-100's RACE, carries the colouring (claimed).
  Color Coding > Pick Up Coloring on demog-1000's RACE, then the header menu of demog-100's RACE: Apply
  Coloring is disabled with "Pick up coloring first" (seen 26 Oct 2026). The issue's note says applying
  the colour coding should work from a style picked up in a different dataframe: is the Color Coding
  pair meant to be per table?
