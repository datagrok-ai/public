@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The toolbar's Export icon offers the formats; As CSV
  downloads the whole table (demog-1000, 1000 rows, Psoriasis among them). With the rows filtered to DIS_POP = RA
  (434 of 1000), As CSV (options) opens Save as CSV, whose Filtered rows only downloads just those rows (the
  dialog also has the delimiter, selected rows and columns only, and more). demog-1000 is opened off camera; the
  files go to the browser's downloads, nothing is made on the server.

  Scenario: Export data
    Download a table as a file, or only the rows that pass the filter.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Export: download the table as a file
    When user clicks on Export icon in toolbar
    # caption: Pick a format; the file goes to your downloads
    And user downloads a file through "As CSV" text in toolbar
    Then the downloaded file should contain 434 occurrences of ",RA,"
    And the downloaded file should contain "Psoriasis"
    # caption: Filter the rows you want to keep
    When user clicks on filter icon in toolbar
    And user clicks on the "category RA of DIS_POP" area of filter panel
    Then 434 rows should pass the filter
    # caption: CSV with options: what goes into the file
    When user clicks on Export icon in toolbar
    And user clicks on "As CSV (options)..." text in toolbar
    Then "Save as CSV" dialog should be visible
    # caption: Filtered rows only: just the 434 rows on screen
    When user checks "Filtered rows only" input in "Save as CSV" dialog
    And user downloads a file through OK button in "Save as CSV" dialog
    Then the downloaded file should contain 434 occurrences of ",RA,"
    And the downloaded file should not contain "Psoriasis"
    And no errors should have been logged
