@guide @tutorial
Feature: How-to · End to end
  An end-to-end how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films
  it as a demo clip that outlines what each step presses). One task from start to finish, every move on camera:
  Browse > Files > Demo is expanded and demog.csv opened with a double click; the Context Panel is switched on from
  the status bar and the AGE header's Stats pane sums the column up; the toolbox's bar chart, scatter plot (WEIGHT by
  HEIGHT, colored by SEX) and histogram of AGE chart it; the filter panel keeps the RA patients; the ribbon's SAVE
  saves it all as the dashboard "Rheumatoid arthritis overview", and the Share dialog the first save opens shares it
  with the fixture user bddsecond (the library's sharing user) at View and use, with Send notifications off. Then
  the Browse tab, Dashboards and a search show the dashboard's card.

  The project, with its table, view, picture and grant, is deleted when the feature starts and ends; the sharing
  user's notifications are cleared now and at feature end.

  Scenario: From a file to a shared dashboard
    Open a file, look at it, chart it, narrow it to one group, then save it as a dashboard and share it.

    Given user is logged in
    And simple mode is off
    And the sharing user has no notifications, now and when the feature ends
    And no project named "Rheumatoid arthritis overview" is on the server
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    When the video says "The task: a chart overview of one patient group, shared with a colleague"
    # caption: Start from the file: Browse > Files > Demo holds sample data
    And user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    # caption: Double-click demog.csv to open it
    And user double-clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And the table should have 5850 rows
    # caption: First, look: switch on the Context Panel and click a column
    When user clicks on "Context Panel" status bar toggle
    And user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    # caption: Stats sums the column up
    When user expands Stats pane in context panel
    Then Stats pane in context panel should contain text "avg"
    When user clicks on "Context Panel" status bar toggle
    Then toolbox should be visible
    # caption: Chart it from the toolbox: a bar chart counts the patients per group
    When user clicks on "bar chart" icon in toolbox
    Then bar chart viewer should be visible
    When user picks "DIS_POP" in the "split" column selector of bar chart viewer
    # caption: A scatter plot: weight against height, colored by sex
    And user clicks on "scatter plot" icon in toolbox
    Then scatter plot viewer should be visible
    When user picks "HEIGHT" in the "x" column selector of scatter plot viewer
    And user picks "WEIGHT" in the "y" column selector of scatter plot viewer
    And user picks "SEX" in the "color" column selector of scatter plot viewer
    # caption: A histogram shows the ages
    And user clicks on "histogram" icon in toolbox
    Then histogram viewer should be visible
    When user picks "AGE" in the "value" column selector of histogram viewer
    # caption: Narrow it to the group that matters: the filter panel
    And user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user clicks on the "category RA of DIS_POP" area of filter panel
    Then the filter should pass exactly the rows where "DIS_POP" is "RA"
    When the video says "Every chart now shows only the RA patients"
    And user hovers over scatter plot viewer
    # caption: SAVE keeps the table, the charts and the filter together
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "Rheumatoid arthritis overview" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "Rheumatoid arthritis overview" should be on the server
    And "Share Rheumatoid arthritis overview" dialog should be visible
    # caption: Share it: type a colleague's name, keep View and use
    When user types "bddsecond" into "User, group, or email" input in "Share Rheumatoid arthritis overview" dialog
    And user clicks on "bddsecond" table row in "Share Rheumatoid arthritis overview" dialog
    # caption: Uncheck Send notifications to share without a message
    And user unchecks "Send notifications" input in "Share Rheumatoid arthritis overview" dialog
    And user clicks on OK button in "Share Rheumatoid arthritis overview" dialog
    Then the "Share Rheumatoid arthritis overview" dialog should close
    And an info balloon containing "Shared" should have been shown
    # caption: Done. To find it later: Browse > Dashboards
    When user clicks on browse tab
    Then browse panel should be visible
    When user collapses Files tree node inside browse tree
    And user clicks on Dashboards tree node inside browse tree
    And user enters "Rheumatoid arthritis overview" into gallery search
    Then "Rheumatoid arthritis overview" gallery card should be visible
    When the video says "There it is, for you and for them"
    And user hovers over "Rheumatoid arthritis overview" gallery card
    Then no errors should have been logged
