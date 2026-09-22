@serial @realizes:views.scripts
Feature: A script in every language runs its template
  Every language of the New menu opens a template that takes the sample table and answers the
  number of its cells; run with cars (30 rows by 17 columns) it answers 510, and the JavaScript
  template raises a browser alert instead. Nothing is saved. Translated from
  files/TestTrack/Scripts/create.md (13-20) and the language loop of
  playwright-public/scripts/scripts-create-debugged.test.ts.

  R, Python, Octave, NodeJS and Julia run in containers, so their rows are @full-stand (dev has
  them all); Grok and Pyodide run in the browser. Julia is not in the md but is in the New menu and
  runs, so it is here. Octave is scripts-octave.feature: it answers 527 for the same table.

  Not translated, and why: Flow Script... opens a flow designer, not the script editor — another
  feature's subject. The md's "F5 or the button": the button, since F5 in a browser that the
  platform has not focused reloads the page.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And user opens the Scripts view

  Scenario Outline: The <language> template runs with cars
    When user clicks on New button
    And user picks "<language> Script..." from the open menu
    Then the "Template" view should be current
    And code editor should contain the text "sample: cars.csv"
    When user clicks on "Open script sample table" icon
    Then table "cars" should be open
    When user clicks on "Run script (F5)" icon
    Then "Template" dialog should be visible
    When user selects "cars" in Table input in "Template" dialog
    And user clicks on OK button in "Template" dialog
    Then the "Template" dialog should close
    And the script results should show "count" as "510"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    @full-stand
    Examples: container languages
      | language |
      | R        |
      | Python   |
      | NodeJS   |
      | Julia    |

    Examples: languages that run in the browser
      | language |
      | Grok     |
      | Pyodide  |

  Scenario: The JavaScript template raises its alert
    Given browser alerts are recorded
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    And code editor should contain the text "alert('Hello World!')"
    When user clicks on "Run script (F5)" icon
    Then the browser should have shown the alert "Hello World!"
    And no errors should have been logged
    And no error or warning balloon should have been shown
