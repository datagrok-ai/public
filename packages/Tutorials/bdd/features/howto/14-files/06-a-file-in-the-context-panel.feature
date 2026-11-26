@guide @tutorial
Feature: How-to · Files
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). What the context panel shows for a file of Browse: cars.csv of
  Files > Demo clicked (its preview opens), the context panel opened from the status bar, its Details (the share,
  the folder, the size and the dates) and its Actions (Open, Download, Rename..., Delete... and more,
  `file_info_meta.dart`), and Download used: the browser gets the file. Nothing is made on the server; the
  downloaded copy stays in the test's download folder.

  Scenario: A file's details and actions
    Click a file in Browse, and the context panel shows its details and what you can do with it.

    Given user is logged in
    And simple mode is off
    And user watches downloads
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    And user clicks on Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on "Context Panel" status bar toggle
    Then context panel should be visible
    And the context panel should show "cars.csv"
    When the video says "The context panel shows the file you clicked"
    And user expands Details pane in context panel
    And user hovers over the text "Size" in Details pane in context panel
    And the video says "Details: where it is, its size, when it changed"
    And user expands Actions pane in context panel
    And the video says "Actions: what you can do with it"
    And user downloads a file through "Download" label in Actions pane in context panel
    Then a file "cars.csv" should have been downloaded
    And the downloaded file "cars.csv" should contain text "alfaromeo"
    And no errors should have been logged
    When the video says "Download saves it to your computer"
