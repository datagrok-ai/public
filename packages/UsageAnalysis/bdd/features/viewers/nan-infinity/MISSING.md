# NaN and Infinity — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026, against
`files/TestTrack/Viewers/working-with-nan-infinity.md` and its spec. Added in that round:
`working-with-nan-infinity-histogram.feature` (an Infinity in WEIGHT and a NaN in HEIGHT leave the
histogram's value axis on its finite ends, its 20 bins and no error).

## Won't do

- The SPGI_v2_infinity table (the md's step 1): infinities in every viewer are xamgle `d4 | Viewers | Non-finite values
  | <viewer>`, in the column statistics ddt `non_finite_values_test.dart`; a 12 MB file for one parse.
