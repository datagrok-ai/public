# Timelines edge cases — what the features cannot say yet

Written for `timelines-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

## Still open

### 1. Steps requested

1. **An interval as a hit area.** The status reports `view` and `lane <value>` only, so no step can
   hover or click one interval. The mixed-types warning ("The columns of different types cannot be
   used for representing dates.", `timelines-viewer.ts:514-519`) is raised only from a hover or a
   click on an interval (`isSameDate`), so it has no scenario. Needs `interval <lane> <n>` areas
   from the series data; then: End = a numeric column, hover the first interval → that warning
   balloon.

### 2. Decided

1. **An event whose End is earlier than its Start** is drawn as a point: right for now (2026-10-30),
   claimed in `timelines-edge-cases.feature` (the interval count).
2. **Date Format set by the user** goes back to the default on any column change: public#4459, no
   scenario yet (no reading of the axis labels).
3. **An event with a missing start or end date while Show Open Intervals is off** is drawn as a
   point at its known date, not as an interval: what the help means by the option (2026-10-30).

### 3. Questions for owners

None.

### 4. API or Dart tests to write

None.

## Won't do

Nothing.
