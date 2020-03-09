<!-- TITLE: Tests: Shape Map -->
<!-- SUBTITLE: -->

# Tests: Shape Map

Shows a map that is applicable for the specified dataset. Typically, it would 
represent a geographical area (countries, states, counties, etc), but also
it can show an arbitrary shapes (such as a store floor plan, brain regions,
or EEG electrodes).

## Testing scenario

1. Open "ua_population" dataset

1. Add [Shape Map](../viewers/shape-map.md) viewer 
   * Viewer added
   * Map of Ukraine regions is drawn on viewer
     
1. Change aggregation in selector on viewer
   * Aggregations changes affect colors on the map
   
1. Click on region on map 
   * Row (s) selected 
   
1. Select not matching rows (from viewer context menu or from "hamburger" menu)   
   * One row ("Kiev City 1") selected
   
1. Change display color (from viewer context menu or from "hamburger" menu) 
   * Color changes according to the selected palette
   
1. Test common features for viewers ([viewers-test](../viewers/viewers-test.md))   
   
See also:
 * [Shape Map](../viewers/shape-map.md)  
