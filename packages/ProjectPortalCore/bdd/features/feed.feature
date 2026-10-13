@portal
Feature: Activity feed
  A program's activity is read from the portal's own records: gate decisions and their evidence,
  runs saved and selected, documents, snapshots, compound role changes, what the sync changed, new
  data. The program page lists it newest first, by kind; My feed puts together everything the user
  follows.

  Background:
    Given user opens the Project Portal app

  Scenario: A gate decision shows in the project's activity and in My feed
    Given a scratch project "BDD-FEED" at stage "TV" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on feed project card
    And user clicks on open program button
    And user clicks on follow button
    And user clicks on assay criterion met button
    And user clicks on advance button
    Then activity heading should contain text "Activity"
    And activity entries should contain text "Advance: Target Validation →"
    When user selects "Gate decisions (1)" in activity kind input
    Then activity entries should contain text "with 1 piece of evidence"
    When user clicks on my feed button
    Then my feed dialog should be visible
    And my feed entries should contain text "Advance: Target Validation →"
    When user clicks on first my feed program link
    Then my feed dialog should not be visible
    And program title should contain text "BDD-FEED"
