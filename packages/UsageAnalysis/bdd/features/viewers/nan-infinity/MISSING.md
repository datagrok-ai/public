# NaN and Infinity — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026, against
`files/TestTrack/Viewers/working-with-nan-infinity.md` and its spec. Added in that round:
`working-with-nan-infinity-histogram.feature` (an Infinity in WEIGHT and a NaN in HEIGHT leave the
histogram's value axis on its finite ends, its 20 bins and no error).

## Still open

- **The SPGI_v2_infinity table** (the md's step 1, "SPGI_v2_infinity.csv — add viewers + layouts").
  Reachable now (`System:AppData/ApiTests/datasets/SPGI_v2_infinity.csv` is on the local stand and
  the Browse tree opens it), but not a BDD check: the claim is about the data, not a gesture (the
  Chemical Space X / Y columns hold Infinity values parsed from the CSV, and the viewers' axes must
  stay finite), and the file is 12 MB, 3624 rows with structures, against the smallest-dataset rule.
  Requested from ApiTests: read the file, add a scatter plot on Chemical Space X / Y and a histogram
  on Chemical Space X, and check that the scatter plot's x axis max and the histogram's axis max are
  finite and neither viewer reports an error. A three-row CSV with `Infinity` / `-Infinity` cells
  shows the same parse.
