@viewers @realizes:viewers.chrome
Feature: Viewer chrome — title and description, Pick Up / Apply, the context menu, docking and help
  The title and the description every viewer inherits from the base viewer: shown when set, placed
  where Description Position says, gone under Never, and cleared again. One outline over the
  viewers whose own features used to carry the same scenario; a viewer with a chrome of its own
  (the tile viewer's above/below content, the trellis plot's description slot) keeps its own.
  Then what the `*-ui.md` checklists ask of several viewers alike: Pick Up on one viewer and Apply
  on a second of the same type copies a setting, and a later change of the first leaves the second
  as it was; Tooltip > Edit... opens the tooltip editor ("Edit Aggregated Tooltip" for the line
  chart, which aggregates its points) and CANCEL closes it; the context menu carries the General and Tooltip groups, and General's Clone and Close act; the viewer
  added to the right of the grid, docks along the left edge of the view by its title bar, and its "?"
  icon opens its help page (the page's title is searched in the help panel).
  Not translated: the General items that leave the page (Save to Gallery, Save as PNG, Embed...)
  beyond being offered, editing the tooltip's columns inside the dialog, and undocking the viewer
  into a floating window.
  Then closed issues every viewer inherits, claimed on the scatter plot they were reported on: the
  gear switches the context panel to its viewer after a property was edited, in the same view or in
  another (public#667, public#719); a viewer stacked over the grid and closed does not come back from
  a layout (public#1745; a stacked viewer has no title bar, so it is closed with General > Close); a
  title typed into the header becomes the Title property (public#2925; the title is set first, as an
  empty title shows no text to double-click); a layout with a viewer bound to another table applies
  over the tables opened anew with no warning and keeps the binding (public#338). The layouts saved
  to the server go at feature end.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario Outline: <viewer> shows and clears its title and description
    Given user adds a <viewer> viewer
    When user sets properties of <viewer> viewer:
      | Show Title | true         |
      | Title      | Demographics |
    Then title of <viewer> viewer should have text "Demographics"
    When user sets properties of <viewer> viewer:
      | Description                 | By race |
      | Description Visibility Mode | Always  |
    Then description of <viewer> viewer should have text "By race"
    And the description of <viewer> viewer should be above its content
    When user sets "Description Position" property of <viewer> viewer to "Bottom"
    Then description of <viewer> viewer should be visible
    And the description of <viewer> viewer should be below its content
    When user sets "Description Position" property of <viewer> viewer to "Left"
    Then description of <viewer> viewer should be visible
    When user sets "Description Position" property of <viewer> viewer to "Right"
    Then description of <viewer> viewer should be visible
    When user sets "Description Visibility Mode" property of <viewer> viewer to "Never"
    Then description of <viewer> viewer should be absent
    When user sets properties of <viewer> viewer:
      | Show Title                  | false |
      | Title                       |       |
      | Description                 |       |
      | Description Visibility Mode | Auto  |
      | Description Position        | Top   |
    Then description of <viewer> viewer should be absent
    And no errors should have been logged

    Examples:
      | viewer           |
      | bar chart        |
      | box plot         |
      | histogram        |
      | pc plot          |
      | pie chart        |
      | line chart       |
      | density plot     |
      | correlation plot |

  Scenario Outline: Pick Up / Apply copies the <viewer>'s <property> onto a second <viewer>, and the two stay apart
    Given user adds a <viewer> viewer
    And user adds a <viewer> viewer
    Then the open tableview should have 2 <viewer> viewers
    When user sets "<property>" property of first <viewer> viewer to "<value>"
    Then "<property>" property of second <viewer> viewer should not be "<value>"
    When user picks "Pick Up / Apply > Pick Up" from the viewer menu of first <viewer> viewer
    And user picks "Pick Up / Apply > Apply" from the viewer menu of second <viewer> viewer
    Then "<property>" property of second <viewer> viewer should be "<value>"
    When user sets "<property>" property of first <viewer> viewer to "<other>"
    Then "<property>" property of first <viewer> viewer should be "<other>"
    And "<property>" property of second <viewer> viewer should be "<value>"
    And no errors should have been logged

    Examples:
      | viewer       | property    | value     | other     |
      | scatter plot | Color       | RACE      | SEX       |
      | box plot     | Category 1  | RACE      | DIS_POP   |
      | trellis plot | Viewer Type | Bar chart | Histogram |
      | line chart   | lineWidth   | 3         | 5         |
      | bar chart    | Stack       | SEX       | RACE      |
      | histogram    | Bins        | 15        | 30        |

  Scenario Outline: Tooltip > Edit... opens the tooltip editor of the <viewer>, and CANCEL closes it
    Given user adds a <viewer> viewer
    Then <viewer> viewer should be painted
    When user picks "Tooltip > Edit..." from the viewer menu of <viewer> viewer
    Then "<dialog>" dialog should be visible
    When user clicks on CANCEL button in "<dialog>" dialog
    Then "<dialog>" dialog should be absent
    And no errors should have been logged

    Examples:
      | viewer     | dialog                  |
      | line chart | Edit Aggregated Tooltip |
      | box plot   | Edit Tooltip            |
      | pc plot    | Edit Tooltip            |

  Scenario Outline: The <viewer>'s context menu offers the General and Tooltip groups, and General clones and closes it
    Given user adds a <viewer> viewer
    When user opens the viewer menu of <viewer> viewer
    Then the open menu should list "General > Clone"
    And the open menu should list "General > Full Screen"
    And the open menu should list "General > Save to Gallery"
    And the open menu should list "General > Save as PNG"
    And the open menu should list "General > Embed..."
    And the open menu should list "General > Close"
    And the open menu should list "Tooltip > Edit..."
    And the open menu should list "<own item>"
    When user closes the context menu
    And user picks "General > Clone" from the viewer menu of <viewer> viewer
    Then the open tableview should have 2 <viewer> viewers
    When user picks "General > Close" from the viewer menu of second <viewer> viewer
    Then the open tableview should have 1 <viewer> viewer
    And no errors should have been logged

    Examples:
      | viewer       | own item                  |
      | pc plot      | Columns > Columns         |
      | trellis plot | Pick Up / Apply > Pick Up |

  Scenario Outline: The <viewer> docks along the left edge by its title bar
    Given user adds a <viewer> viewer
    Then <viewer> viewer should not be docked along the left edge of the view
    When user docks <viewer> viewer to the left edge of the view
    Then <viewer> viewer should be docked along the left edge of the view
    And no errors should have been logged

    Examples:
      | viewer       |
      | pc plot      |
      | trellis plot |

  Scenario Outline: The <viewer> opens its help from the title bar
    Given the stand serves the help pages
    And user adds a <viewer> viewer
    When user opens the help of <viewer> viewer
    Then help panel should contain text "<help>"
    And no errors should have been logged

    Examples:
      | viewer       | help                      |
      | pc plot      | Parallel coordinates plot |
      | trellis plot | Trellis                   |

  @realizes:public#667 @realizes:public#719
  Scenario: After a property is edited in the context panel, the gear of another viewer, in this view or another, switches the panel to it (public#667, public#719)
    Given the context panel is open
    And user adds a scatter plot viewer
    And user adds a histogram viewer
    And user opens cars dataset
    And user adds a histogram viewer
    When user clicks on settings icon of histogram viewer
    Then the context panel should show "Histogram"
    Given "X Axis" category in context panel is expanded
    When user enters "30" in "Bins" property in context panel
    Then "Bins" property of histogram viewer should be "30"
    Given user switches to the "demog-1000" table view
    When user clicks on settings icon of scatter plot viewer
    Then the context panel should show "Scatter plot"
    And "Bins" property in context panel should be absent
    Given "Color" category in context panel is expanded
    When user selects "RACE" in "Color" property in context panel
    Then "Color" property of scatter plot viewer should be "RACE"
    When user clicks on settings icon of histogram viewer
    Then the context panel should show "Histogram"
    And "Marker Type" property in context panel should be absent
    And no errors should have been logged

  @realizes:public#1745
  Scenario: A viewer stacked over the grid and closed does not come back as a floating window from the layout (public#1745)
    Given user adds a scatter plot viewer
    And user adds a histogram viewer
    When user stacks scatter plot viewer onto grid as a tab
    Then the viewers "demog-1000, Scatter plot" should be the tabs of one panel
    When user picks "General > Close" from the viewer menu of scatter plot viewer
    Then the open tableview should have 0 scatter plot viewers
    And the open tableview should have 1 histogram viewer
    When user saves the layout of the current table view to the server
    And user picks "General > Close" from the viewer menu of histogram viewer
    Then the open tableview should have 0 histogram viewers
    When user loads the saved layout
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 0 scatter plot viewers
    And grid should be visible
    And no errors should have been logged

  @realizes:public#2925
  Scenario: A title typed into the viewer's header becomes its Title (public#2925)
    Given user adds a scatter plot viewer
    And user sets "Title" property of scatter plot viewer to "Plot"
    Then title of scatter plot viewer should have text "Plot"
    When user double-clicks on title of scatter plot viewer
    And user types "Heights" at the caret
    And user presses Enter
    Then title of scatter plot viewer should have text "Heights"
    And "Title" property of scatter plot viewer should be "Heights"
    And no errors should have been logged

  @realizes:public#338
  Scenario: A layout with a viewer bound to another table applies over the tables opened anew with no warning (public#338)
    Given the context panel is open
    And user opens cars dataset
    And user adds a scatter plot viewer
    When user clicks on settings icon of scatter plot viewer
    And user selects "demog-1000" in "Table" property in context panel
    Then scatter plot viewer should be bound to table "demog-1000"
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user opens cars dataset
    And user loads the saved layout
    Then scatter plot viewer should be bound to table "demog-1000"
    And no error or warning balloon should have been shown
    And no errors should have been logged
