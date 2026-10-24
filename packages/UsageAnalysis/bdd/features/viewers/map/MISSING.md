# Map — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against `files/TestTrack/Viewers/Map/map.md` and
the GIS package's `CLAUDE.md` "Automation surface". `map.feature` and `map-navigation.feature` (wheel
zoom, drag pan) hold what could be said. Every claim on this viewer has to be a reading, never the
picture: the base maps are internet tiles.

## Still open

- **Marker Min Size / Size column.** The md's "Marker Min Size 12 — the markers grow" and "Size by
  Depth" have no reading: `map.feature` claims only that no feature is lost. Needs a `marker size min`
  / `marker size max` reading (the radii the WebGL marker style was built with, `gis-openlayer.ts`
  marker style) so the scenario can be:

  ```gherkin
  Scenario: Marker Min Size grows the smallest marker
    Given user remembers the "marker size min" reading of map viewer
    When user sets "markerMinSize" property of map viewer to "12"
    Then the "marker size min" reading of map viewer should be 12
    And the "marker size min" reading of map viewer should be higher than remembered
    And the "markers" reading of map viewer should be 2426
    And no errors should have been logged
  ```

- **Selected markers.** The Ctrl-drag in `map.feature` claims "some rows should be selected" — the
  map's own count of highlighted features is not reported, so a selection made elsewhere (the grid)
  reaching the map cannot be claimed. Needs a `selected markers` reading (the features the selection
  layer draws, `gis-viewer.ts` selection handling):

  ```gherkin
  Scenario: Rows selected in the grid are the markers the map highlights
    When user selects rows where "MagType" is "Mb"
    Then 72 rows should be selected
    And the "selected markers" reading of map viewer should be 72
    When user selects no rows
    Then the "selected markers" reading of map viewer should be 0
    And no errors should have been logged
  ```

  Also the md's "Selected Color" property (step 20) needs the same reading plus the colour it uses.

- **Custom geo layers (doc.kml, ecoregions2.geojson, boothbay2018.kmz, uk.topojson).** The md drops the
  files onto the map, unchecks the new layer and deletes it. Missing: (1) the fixtures in a place the
  suite can read (the GIS package's `files/` or DemoFiles); (2) a library step that drops a local file
  onto an element and lets the drop complete (`a file dragged from disk enters {element} and the drag's
  end never reaches the page` is the only file-drag step, and it never drops); (3) a
  `layer "<name>" delete` area — the `del` cell of the layers grid, next to the `vis` cell
  `layer "<name>" visibility` already reports (`gis-viewer-status.ts`).

  ```gherkin
  Scenario: A KML file dropped on the map becomes a layer, which can be hidden and deleted
    When user drops the file "doc.kml" onto map viewer
    Then the "layers" reading of map viewer should be 6
    And the "layer \"doc.kml\" visible" reading of map viewer should be "true"
    When user clicks on "Map layers" button
    And user clicks on the "layer \"doc.kml\" visibility" area of map viewer
    Then the "layer \"doc.kml\" visible" reading of map viewer should be "false"
    When user clicks on the "layer \"doc.kml\" delete" area of map viewer
    Then the "layers" reading of map viewer should be 5
    And map viewer should not have a "layer \"doc.kml\"" area
    And no errors should have been logged
  ```
