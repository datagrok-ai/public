# Timelines — what the features cannot say yet

## Edge cases

The edge cases of the Timelines viewer (26 Oct 2026) that no feature claims, and why. Nothing here
asks for a library, binding or code change.

### Still open

#### 1. Steps requested

1. **An interval as a hit area.** The status reports `view` and `lane <value>` only, so no step can
   hover or click one interval. The mixed-types warning ("The columns of different types cannot be
   used for representing dates.", `timelines-viewer.ts:514-519`) is raised only from a hover or a
   click on an interval (`isSameDate`), so it has no scenario. Needs `interval <lane> <n>` areas
   from the series data; then: End = a numeric column, hover the first interval → that warning
   balloon.

#### 2. Decided

1. **An event whose End is earlier than its Start** is drawn as a point: right for now (2026-10-30),
   API test 4.
2. **Date Format set by the user** goes back to the default on any column change: public#4459, no
   scenario yet (no reading of the axis labels).
3. **An event with a missing start or end date while Show Open Intervals is off** is drawn as a
   point at its known date, not as an interval: what the help means by the option (2026-10-30).

#### 3. Questions for owners

None.

#### 4. API or Dart tests to write

1. Color set to grade gives a legend of 2 items; Color cleared hides the legend, the chart stays
   painted, no error (public#617).
2. Every row filtered out (dose 50–60): 0 rows shown and 0 intervals; the filter reset gives 4 rows;
   no error.
3. Two lanes with 200-character names: 2 lanes, 2 intervals, the chart painted, no error
   (public#1973).
4. An event whose End is earlier than its Start: 2 lanes, 2 intervals, the chart painted, no error.

Moved out of BDD (no user gesture); requested from ApiTests.

### Won't do

Nothing.
