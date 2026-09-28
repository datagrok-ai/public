@apps @demos @realizes:views.browse
Feature: The Visualization demos open from Browse > Apps > Demo with their viewer
  Each viewer demo is opened by its node in the Browse tree; the Tutorials demo app says when it has
  run and its view is named (`demo-loaded`, with the demo's path), and the row claims the viewer the
  demo is about and the table it shows. A row starts with the capability gate on the package the
  demo comes from (Tutorials, Charts, Dendrogram). Translated from the TestTrack case Apps/apps.md 2
  and playwright-public/browse/demo_apps.test.ts; the case's "select points: the highlight follows"
  is the last scenario (the demo's histogram and bar chart paint the selection in the demo's own
  colour, which the pixel reading of a selection does not know, so the scatter plot is the witness).

  The @full-stand rows read a demo file a smaller stand's System:DemoFiles may lack (Line Chart,
  Correlation Plot, Chord and Sankey read energy_uk.csv-like files, Word Cloud word_cloud.csv, Network
  Diagram got-s1-edges.csv and fetches its node pictures from an outside image service), or draw
  base-map tiles from an outside tile server (Map); they claim the viewer and that it reports no error.
  Data Annotations is a server project only a dev stand carries, not a demo function.

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Demo tree node inside browse tree is expanded
    And Apps---Demo---Visualization tree node inside browse tree is expanded

  Scenario Outline: The <demo> demo opens with its <viewer> viewer
    Given the "<package>" package is installed
    And Apps---Demo---Visualization---<group> tree node inside browse tree is expanded
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Visualization---<group>---<node> tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Visualization | <section> | <demo>"
    And the "<demo>" view should be current
    And <viewer> viewer should be visible
    And the table should have <rows> rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | package    | group                   | section                 | demo            | node            | viewer          | rows |
      | Tutorials  | General                 | General                 | Scatter Plot    | Scatter-Plot    | scatter plot    | 5850 |
      | Tutorials  | General                 | General                 | Bar Chart       | Bar-Chart       | bar chart       | 5850 |
      | Tutorials  | General                 | General                 | Histogram       | Histogram       | histogram       | 5850 |
      | Tutorials  | General                 | General                 | Pie Chart       | Pie-Chart       | pie chart       | 5850 |
      | Tutorials  | General                 | General                 | 3D Scatter Plot | 3D-Scatter-Plot | 3d scatter plot | 5850 |
      | Tutorials  | General                 | General                 | Density Plot    | Density-Plot    | density plot    | 5850 |
      | Tutorials  | General                 | General                 | Filters         | Filters         | filters         | 200  |
      | Tutorials  | General                 | General                 | Markup          | Markup          | markup          | 5850 |
      | Tutorials  | General                 | General                 | Tile Viewer     | Tile-Viewer     | Tile Viewer     | 200  |
      | Tutorials  | Data-Flow-and-Hierarchy | Data Flow and Hierarchy | Tree Map        | Tree-Map        | tree map        | 5850 |
      | Tutorials  | Data-Separation         | Data Separation         | Trellis Plot    | Trellis-Plot    | trellis plot    | 5850 |
      | Tutorials  | Data-Separation         | Data Separation         | Matrix Plot     | Matrix-Plot     | matrix plot     | 5850 |
      | Tutorials  | Input-and-Edit          | Input and Edit          | Grid            | Grid            | grid            | 100  |
      | Tutorials  | Input-and-Edit          | Input and Edit          | Form            | Form            | form            | 200  |
      | Tutorials  | Statistical             | Statistical             | Box Plot        | Box-Plot        | box plot        | 5850 |
      | Tutorials  | Statistical             | Statistical             | PC Plot         | PC-Plot         | PC Plot         | 5850 |
      | Tutorials  | Statistical             | Statistical             | Pivot Table     | Pivot-Table     | pivot table     | 5850 |
      | Tutorials  | Statistical             | Statistical             | Statistics      | Statistics      | statistics      | 5850 |
      | Tutorials  | Time-and-Date           | Time and Date           | Calendar        | Calendar        | calendar        | 5850 |
      | Charts     | General                 | General                 | Radar           | Radar           | radar           | 5850 |
      | Charts     | General                 | General                 | Sunburst        | Sunburst        | sunburst        | 5850 |
      | Charts     | General                 | General                 | Surface Plot    | Surface-Plot    | surface plot    | 121  |
      | Charts     | General                 | General                 | Timelines       | Timelines       | timelines       | 143  |
      | Charts     | Data-Flow-and-Hierarchy | Data Flow and Hierarchy | Tree            | Tree            | tree            | 5850 |
      | Dendrogram | General                 | General                 | Heatmap         | Heatmap         | grid            | 31   |

  @full-stand
  Scenario Outline: The <demo> demo opens with its <viewer> viewer on a full stand
    Given the "<package>" package is installed
    And Apps---Demo---Visualization---<group> tree node inside browse tree is expanded
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Visualization---<group>---<node> tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Visualization | <section> | <demo>"
    And the "<demo>" view should be current
    And <viewer> viewer should be visible
    And <viewer> viewer should report no error
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | package   | group                   | section                 | demo             | node             | viewer           |
      | Tutorials | General                 | General                 | Line Chart       | Line-Chart       | line chart       |
      | Tutorials | Statistical             | Statistical             | Correlation Plot | Correlation-Plot | correlation plot |
      | Tutorials | Data-Flow-and-Hierarchy | Data Flow and Hierarchy | Network Diagram  | Network-Diagram  | network diagram  |
      | Charts    | General                 | General                 | Chord            | Chord            | chord            |
      | Charts    | General                 | General                 | Sankey           | Sankey           | sankey           |
      | Charts    | General                 | General                 | Word Cloud       | Word-Cloud       | word cloud       |
      | Gis       | Geographical            | Geographical            | Map              | Map              | map              |

  # apps.md 2: "select points on charts — selected cells or points highlight properly"
  Scenario: Rows selected in the Scatter Plot demo reach its scatter plot
    Given the "Tutorials" package is installed
    And Apps---Demo---Visualization---General tree node inside browse tree is expanded
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Visualization---General---Scatter-Plot tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Visualization | General | Scatter Plot"
    When user takes a snapshot of scatter plot viewer
    And user selects rows where "SEX" is "F"
    Then only rows where "SEX" is "F" should be selected
    And the "rows selected" reading of scatter plot viewer should be 3243
    And scatter plot viewer should show more selection highlight than before
    And no errors should have been logged
