@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Four ways to open a table: from Browse (Files > Demo > cars.csv:
  a click previews the file, a double-click opens it), from the computer through the Open local file icon of the
  Browse toolbar (files/assay-results.csv of this folder), by dropping a file on the window
  (files/compound-yields.csv), and by pasting text into the Import text view the Open text icon opens. A file opened
  from the computer stays in the browser: nothing is uploaded, nothing is made on the server. My stuff > Recent is
  left out: it lists entities (dashboards, connections), not the files opened here.

  Scenario: Open your data
    Open a file from Browse or from your computer, drop it on the window, or paste it as text.

    Given user is logged in
    And simple mode is off
    # caption: Browse: the files and databases on the server
    When user clicks on browse tab
    # caption: Files > Demo: sample data that comes with Datagrok
    And user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    # caption: One click previews a file
    And user clicks on Files---Demo---cars.csv tree node inside browse tree
    # caption: A double-click opens it as a table
    And user double-clicks on Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    And the table should have 30 rows
    # caption: A file on your computer: the open-file icon of Browse
    When user clicks on browse tab
    And user uploads "features/howto/01-getting-started/files/assay-results.csv" through "Open local file" icon inside browse toolbar
    Then the "assay-results" view should be current
    And the table should have 24 rows
    # caption: Or drag a file from your computer onto the window
    When user drops the "features/howto/01-getting-started/files/compound-yields.csv" file of the project onto grid
    Then the "compound-yields" view should be current
    And the table should have 12 rows
    # caption: Copied cells or text: the Open text icon
    When user clicks on browse tab
    And user clicks on "Open text" icon inside browse toolbar
    Then the "Import text" view should be current
    # caption: Paste them in; the preview below follows
    When user clicks on code editor
    And user pastes "compound,yield\nCMP-201,72.5\nCMP-202,65.0\nCMP-203,81.2\nCMP-204,58.9" into code editor
    # caption: DONE opens it as a table
    And user clicks on DONE button
    Then the table should have 4 rows
    # caption: Each table has its own tab: click one to go back
    When user clicks on "cars" view tab
    Then the "cars" view should be current
    And no errors should have been logged
