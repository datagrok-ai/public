<!-- TITLE: Tests: Shape map -->
<!-- SUBTITLE: -->

# Tests: shape map

Shows a map that is applicable for the specified dataset. Typically, it would represent a geographical area (countries,
states, counties, etc), but also it can show an arbitrary shapes (such as a store floor plan, brain regions, or EEG
electrodes).

## Testing scenario

1. Open "ua_population" dataset

2. Add [Shape Map](../viewers/shape-map.md) viewer

    * Viewer added
    * Map of Ukraine regions is drawn on viewer

3. Change aggregation in selector on viewer

    * Aggregations changes affect colors on the map

4. Click on region on map

    * Row (s) selected

5. Select not matching rows (from viewer context menu or from "hamburger" menu)

    * One row ("Kiev City 1") selected

6. Change display color (from viewer context menu or from "hamburger" menu)

    * Color changes according to the selected palette

7. Test common features for viewers ([viewers-test](../viewers/viewers-test.md))

See also:

* [Shape map](../viewers/shape-map.md)
