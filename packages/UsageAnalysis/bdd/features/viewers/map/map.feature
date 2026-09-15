@journey @viewers @realizes:viewers.map-viewer
Feature: Map viewer layers, zoom, selection and the point tooltip
  The GIS package's Map — an OpenLayers map that draws a table's Latitude / Longitude columns as
  WebGL markers or a heat map over a base map. Not the Dart Shape Map, which is an SVG choropleth.

  NOT TRANSLATED, and the reason: four of the eleven scenarios of the spec this replaces asserted
  only that a screenshot of the viewer hashed differently after a change — "Color by Magnitude",
  "Size by Depth", "Marker Min Size redraws the markers" and "Render Type cycles markers / heatmap
  / both" each ended in a pixel comparison. `initMap` builds a `BingMaps` aerial layer and an
  `OSM()` base layer before anything of ours is drawn, and both fetch tiles from the public
  internet, so that hash is a claim about a tile server's availability and latency, not about
  Datagrok: offline it fails, online a tile landing one frame later changes it. The tile layers and
  the WebGL marker layer also render onto SEPARATE canvases inside `.ol-viewport`, so "the canvas"
  of this viewer is whichever comes first in the DOM — the picture is not even well defined here.
  What those four scenarios could still say honestly is kept: the properties round-trip, Render
  Type moves the two product layers (`layer "Heatmap" visible` / `layer "Markers GL" visible`), and
  no feature is lost. Nothing below compares pictures.

  Two corrections to the old spec. The layers panel it clicked at `panel.x + 11, panel.y + 11 +
  row * 24` is an OpenLayers control whose body is a real `DG.Grid` with `rowHeight: 24` — that is
  where the hardcoded 24-pixel row pitch came from; the viewer now reports the rows itself, and
  there are two areas per layer: `layer "<name>"` is the row (clicking it makes the layer current)
  and `layer "<name>" visibility` is the checkbox cell (clicking it toggles the layer). And the
  panel starts hidden, so those areas are absent until the layers button is clicked — the
  `layer "<name>" visible` READINGS need no panel and are what the claims below rest on.

  Fixture: earthquakes, 2426 rows; MagType is Mw 2312, Mb 72, Me 25, ML 7, Unk 6, Ms 4.

  Background:
    Given user is logged in
    And user opens earthquakes dataset
    And user adds a map viewer
    Then map viewer should be visible
    And 2426 rows should pass the filter
    And the "markers" reading of map viewer should be 2426
    And the "rows shown" reading of map viewer should be 2426

  Scenario: The geo columns are detected and every filtered row becomes a marker
    Then "latitudeColumnName" property of map viewer should be "Latitude"
    And "longitudeColumnName" property of map viewer should be "Longitude"
    And the "markers" and "rows shown" readings of map viewer should be the same
    And the "render type" reading of map viewer should be "markers"
    And the "layer \"Markers GL\" visible" reading of map viewer should be "true"
    And the "layer \"Heatmap\" visible" reading of map viewer should be "false"
    And the "layers" reading of map viewer should be 5
    And map viewer should report no error
    And no errors should have been logged

  Scenario: A colour column and a size column keep every feature on the map
    When user sets properties of map viewer:
      | colorColumnName | Magnitude |
      | sizeColumnName  | Depth     |
    Then "colorColumnName" property of map viewer should be "Magnitude"
    And "sizeColumnName" property of map viewer should be "Depth"
    And the "markers" reading of map viewer should be 2426
    And the "rows shown" reading of map viewer should be 2426
    And the "layer \"Markers GL\" visible" reading of map viewer should be "true"
    When user sets properties of map viewer:
      | colorColumnName |  |
      | sizeColumnName  |  |
    Then the "markers" reading of map viewer should be 2426
    And no errors should have been logged

  Scenario: The layers button reveals a row per layer, and closes them again
    Then map viewer should not have a "layers panel" area
    And map viewer should not have a "layer \"Heatmap\"" area
    When user clicks on "Map layers" button
    Then map viewer should have a "layers panel" area
    And map viewer should have a "layer \"Bing sat\"" area
    And map viewer should have a "layer \"BaseLayer\"" area
    And map viewer should have a "layer \"Heatmap\"" area
    And map viewer should have a "layer \"Markers GL\"" area
    And map viewer should have a "layer \"Heatmap\" visibility" area
    When user clicks on "Map layers" button
    Then map viewer should not have a "layers panel" area
    And map viewer should not have a "layer \"Heatmap\"" area
    And no errors should have been logged

  Scenario: The checkbox cell of a row toggles that layer, the row itself makes it current
    When user clicks on "Map layers" button
    Then the "layer \"Heatmap\" visible" reading of map viewer should be "false"
    When user clicks on the "layer \"Heatmap\" visibility" area of map viewer
    Then the "layer \"Heatmap\" visible" reading of map viewer should be "true"
    And the "visible layers" reading of map viewer should be 5
    When user clicks on the "layer \"Markers GL\" visibility" area of map viewer
    Then the "layer \"Markers GL\" visible" reading of map viewer should be "false"
    And the "visible layers" reading of map viewer should be 4
    When user clicks on the "layer \"BaseLayer\"" area of map viewer
    Then "currentLayer" property of map viewer should be "BaseLayer"
    And the "layer \"BaseLayer\" visible" reading of map viewer should be "true"
    When user clicks on the "layer \"Markers GL\" visibility" area of map viewer
    And user clicks on the "layer \"Heatmap\" visibility" area of map viewer
    Then the "layer \"Markers GL\" visible" reading of map viewer should be "true"
    And the "layer \"Heatmap\" visible" reading of map viewer should be "false"
    When user clicks on "Map layers" button
    Then map viewer should not have a "layers panel" area
    And no errors should have been logged

  Scenario: Render Type moves the two product layers and leaves the base maps alone
    When user sets "renderType" property of map viewer to "heatmap"
    Then the "render type" reading of map viewer should be "heatmap"
    And the "layer \"Heatmap\" visible" reading of map viewer should be "true"
    And the "layer \"Markers GL\" visible" reading of map viewer should be "false"
    And the "layer \"BaseLayer\" visible" reading of map viewer should be "true"
    When user sets "renderType" property of map viewer to "both"
    Then the "render type" reading of map viewer should be "both"
    And the "layer \"Heatmap\" visible" reading of map viewer should be "true"
    And the "layer \"Markers GL\" visible" reading of map viewer should be "true"
    And the "visible layers" reading of map viewer should be 5
    When user sets "renderType" property of map viewer to "markers"
    Then the "render type" reading of map viewer should be "markers"
    And the "layer \"Heatmap\" visible" reading of map viewer should be "false"
    And the "layer \"Markers GL\" visible" reading of map viewer should be "true"
    And the "visible layers" reading of map viewer should be 4
    And no errors should have been logged

  Scenario: The zoom buttons move the view by one level each way
    Given user remembers the "zoom" reading of map viewer
    When user clicks on the "zoom in" area of map viewer
    Then the "zoom" reading of map viewer should be higher than before
    And the "zoom" reading of map viewer should not be as remembered
    When user clicks on the "zoom out" area of map viewer
    Then the "zoom" reading of map viewer should be lower than before
    And the "zoom" reading of map viewer should be as remembered
    And no errors should have been logged

  Scenario: Ctrl and a drag select the points inside the rectangle, Escape clears them
    Given user clears the row selection
    When user drags a box over the "view" area of map viewer holding Control
    Then some rows should be selected
    And every selected row should pass the filter
    When user presses Escape
    Then no rows should be selected
    And no errors should have been logged

  Scenario: A table filter narrows what the map holds
    When user adds a categorical filter on "MagType" keeping "Mw"
    Then 2312 rows should pass the filter
    And the "markers" reading of map viewer should be 2312
    And the "rows shown" reading of map viewer should be 2312
    And the "markers" and "rows shown" readings of map viewer should be the same
    When user hovers over "MagType" filter card
    And user clicks on close of "MagType" filter card
    Then 2426 rows should pass the filter
    And the "markers" reading of map viewer should be 2426
    And no errors should have been logged

  Scenario: Show Tooltip decides whether a point under the pointer says anything
    Then "showTooltip" property of map viewer should be "false"
    When user hovers over the "point 5" area of map viewer
    Then tooltip should be hidden
    When user sets "showTooltip" property of map viewer to "true"
    And user moves the pointer away from map viewer
    And user hovers over the "point 5" area of map viewer
    Then exactly one tooltip should be shown
    And tooltip should contain text "Latitude"
    And tooltip should contain text "Longitude"
    And tooltip should contain text "Magnitude"
    And tooltip should contain text "MagType"
    When user sets "showTooltip" property of map viewer to "false"
    And user moves the pointer away from map viewer
    And user hovers over the "point 5" area of map viewer
    Then tooltip should be hidden
    And no errors should have been logged

  Scenario: Closing the map under the pointer disposes it without an error
    Given user remembers the place of map viewer
    When user clicks on close icon of map viewer
    Then map viewer should be absent
    When user moves the pointer across the remembered place
    Then the open tableview should have 0 map viewers
    And no errors should have been logged
