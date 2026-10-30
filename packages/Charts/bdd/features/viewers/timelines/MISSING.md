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

### 3. Questions for owners

1. **Marker Size and Line Width 0 or negative** (plain ints with no bounds,
   `timelines-viewer.ts:81`, `:84`): declare bounds, or clamp?
2. **An invalid Date Format** (free text, `:87`): what is shown?
3. **The mixed date types warning**: is a warning on hover the intended answer to a Start and an
   End of different types?
4. **Events with a missing start or end date while Show Open Intervals is off.** Help
   (`timelines.md`) says "Include events with missing dates by enabling Show Open Intervals", yet
   with the option off such an event is still drawn, as a point at its known date; with the option
   on it becomes an interval to the edge of the data. Should an event with a missing date be left
   out while the option is off, or is the help to say that it is drawn as a point?

### 4. API or Dart tests to write

None.

## Won't do

Nothing.
