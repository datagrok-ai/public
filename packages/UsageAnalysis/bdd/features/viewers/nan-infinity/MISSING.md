# NaN and Infinity — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026, against
`files/TestTrack/Viewers/working-with-nan-infinity.md` and its spec. Added in that round:
`working-with-nan-infinity-histogram.feature` (an Infinity in WEIGHT and a NaN in HEIGHT leave the
histogram's value axis on its finite ends, its 20 bins and no error).

## Still open

- **The SPGI_v2_infinity table** (the md's step 1, "SPGI_v2_infinity.csv — add viewers + layouts").
  The old spec opens `System:AppData/ApiTests/datasets/SPGI_v2_infinity.csv`, whose Chemical Space
  X / Y columns hold Infinity values as the file stores them — parsed from the CSV, not written into
  a cell afterwards. The library has no dataset for it (`bindings/platform/datasets.ts` lists
  `spgi-linked1` / `spgi-linked2` from the same ApiTests folder, not this one). Needs a
  `dataset('spgi-infinity', {path: 'System:AppData/ApiTests/datasets/SPGI_v2_infinity.csv', …})`
  entry (the ApiTests package published, as for the linked tables), and then:

  ```gherkin
  Scenario: Infinity values read from the file leave the scatter plot and the histogram finite
    Given user opens spgi-infinity dataset
    And user adds a scatter plot viewer with:
      | xColumnName | Chemical Space X |
      | yColumnName | Chemical Space Y |
    And user adds a histogram viewer with:
      | valueColumnName | Chemical Space X |
    Then the "x axis max" reading of scatter plot viewer should be a finite number
    And the "axis max" reading of histogram viewer should be a finite number
    And histogram viewer should report no error
    And no errors should have been logged
  ```
