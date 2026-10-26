# Box plot — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed box plot issues into `box-plot-regressions.feature`.

- **The violin style over a value column holding a real NaN (core#960).** The crash needed NaN, not a
  missing value: a column built with `DG.Column.fromList`, where NaN is not none, threw `NaN.floor()` in
  `_defineCatKdeValues` and on every hover. Every data step writes NaN as a missing value: `user sets
  "WEIGHT" column in row 2 to "NaN"` stores the float null (`isNone` true, measured 26 Oct 2026), and so
  do the calculated columns `Sqrt(-1)`, `0/0` and `${WEIGHT} * NaN`. Needs a step that writes a raw NaN
  into a float column's buffer (`col.getRawData()[i] = NaN` followed by `fireValuesChanged`), then:

  ```gherkin
  Scenario: The violin style draws and hovers over a value column holding a NaN (core#960)
    Given user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | WEIGHT |
      | Category 1 | RACE   |
      | Plot Style | violin |
    When user writes a raw NaN into "WEIGHT" column in row 2
    And user hovers over the "Caucasian values" area of box plot viewer
    Then box plot viewer should report no error
    And box plot viewer should have a "category Caucasian" area
    And no errors should have been logged
  ```

- **The statistics strip lines up under each box and uses the column's format (public#3693,
  public#3773).** The status reports one `stats` area for the whole strip and no per-category cell or
  text. Needs `stats <category>` areas (one per category column of the strip) and a `stats <category>`
  reading of the text drawn there:

  ```gherkin
  Scenario: Each category's statistics sit under its own box in the column's format (public#3693)
    Given user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Category 1      | RACE |
      | Value           | AGE  |
      | Show Statistics | true |
    Then the "stats Caucasian" area of box plot viewer should lie below the "Caucasian values" area
    And the "stats Caucasian" area of box plot viewer should lie to the right of the "stats Black" area
    And no errors should have been logged
  ```
