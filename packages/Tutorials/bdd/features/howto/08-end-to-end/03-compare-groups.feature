@guide @tutorial
Feature: How-to · End to end
  An end-to-end how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films
  it as a demo clip that outlines what each step presses). One task from start to finish: do men and women differ in
  weight, and does it hold in every race. demog is opened off camera as "Demographics" (the stand holds many layouts
  named "demog", and a new one would be saved as "demog (2)"), with the toolbox open. The filter panel keeps the RA
  patients; a box plot of WEIGHT by SEX shows the difference and its p-value; the box plot's Category 1 set to RACE
  and its Category2 (in its settings) to SEX split every race by sex on one axis; and Save in the toolbox's Layouts
  section, expanded on camera, keeps the arrangement under the table's name, renamed "Weight by sex and race" in the
  context panel (the save shows it there although the Layouts header was clicked just before). A trellis plot of box plots was tried and dropped: its cells
  scale apart and its control panel offers only the inner Value, so the comparison did not read. Every move after
  the start is on camera.

  The layout, under either name, with its picture, is deleted when the feature ends; one an earlier run of this
  account left is swept at the start (bindings/howto-dashboards.ts).

  Scenario: Compare groups
    Answer "do the groups differ?" with a filter and a box plot, then keep the arrangement as a layout.

    Given user is logged in
    And simple mode is off
    And the layouts named "Demographics, Weight by sex and race" are deleted when the feature ends
    And the user's own layouts named "Demographics, Weight by sex and race" are removed now and at feature end
    And user opens demog dataset keeping the first 5850 rows as "Demographics"
    And the toolbox pane is shown
    And the context panel is closed
    And the browse panel is closed
    When the video says "The question: do men and women differ in weight, in every race"
    # caption: First choose who to compare: the RA patients
    And user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user clicks on the "category RA of DIS_POP" area of filter panel
    Then the filter should pass exactly the rows where "DIS_POP" is "RA"
    # caption: A box plot compares weight between the sexes
    When user clicks on box plot icon on toolbox
    Then box plot viewer should be visible
    When user picks "WEIGHT" in the "Value" column selector of box plot viewer
    And user picks "SEX" in the "Category 1" column selector of box plot viewer
    # caption: Men weigh more; the tiny p-value says it is unlikely to be chance
    And user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "Welch's t-test"
    # caption: Now check every race: make race the category
    When user picks "RACE" in the "Category 1" column selector of box plot viewer
    # caption: and split each race by sex, in the settings
    And user hovers over box plot viewer
    And user clicks on settings icon of box plot viewer
    Then the context panel should show "Box plot"
    When user selects "SEX" in "Category2" property in context panel
    Then box plot viewer should have a "category Asian | F" area
    When the video says "Within every race, the men's box sits higher"
    And user hovers over the "category Black | M" area of box plot viewer
    # caption: Keep this arrangement: Toolbox > Layouts > Save
    And user expands Layouts accordion header in toolbox
    And user clicks on Save button in layouts pane
    Then "Demographics" layout card should be visible
    And the context panel should show "Demographics"
    # caption: Give it a name you will find again
    When user types "Weight by sex and race" into Name details field in context panel
    And user presses Enter
    Then "Weight by sex and race" layout card should be visible
    When the video says "Its card applies it to any table with these columns"
    And user hovers over "Weight by sex and race" layout card
    Then no errors should have been logged
