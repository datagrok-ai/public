@browse @realizes:views.browse
Feature: Opening what the address names
  The address the platform writes for what is open leads back to it, and an address that names
  nothing leaves the shell working. Translated from the TestTrack case Browse/browse.md 5 ("copy the
  URL and open it — the same entity opens") and the manual cases Browse-Route-01..04
  (browse_manual_tests2.md section 12, playwright-public/browse/route.test.ts).

  The address is followed inside the running app (`grok.shell.route`), as a pasted link is by the
  shell's router: a feature has one page, so the case's "open it in a new tab" becomes the same
  address taken again after every view was closed. A query's string parameter goes in quoted, as the
  platform writes it itself (`?shipCountry=%22Germany%22`); unquoted, the value is read as a
  variable. The parameterized query is the Dbtests package's PostgresByStringChoices on the
  NorthwindTest connection (PostgresTest), which the queries features run on too: Germany has 122
  orders; a stand that does not reach its database skips the scenario.

  Background:
    Given user is logged in

  # the project's name has no dash: a project named through the JS API with dashes gets an address
  # (/p/admin.bdd-x-y/...) whose path parse stops at the first dash ("Unable to get project asset
  # "bdd"") — a candidate finding, not claimed until it is walked by hand
  Scenario: The address of an open project opens it again after everything was closed
    Given no project named "bddbrowseroute" is on the server
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    When user saves the current view as project "bddbrowseroute"
    And user closes all views
    And user opens the "bddbrowseroute" project
    Then the page address should contain "bddbrowseroute"
    When user remembers the page address
    And user closes all views
    Then the page address should not contain "bddbrowseroute"
    When user opens the remembered address
    Then the page address should contain "bddbrowseroute"
    And the current view should hold at least 2 viewers
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A folder address opens that folder
    When user opens the address "/files/System.DemoFiles/chem"
    Then the "Demo/chem" view should be current
    And the gallery counter should show as many items as the "System:DemoFiles/chem/" folder holds on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A query address runs the query with the parameter it carries
    # the query runs on db.datagrok.ai, outside the stand
    Given the stand has a reachable "PostgresTest" connection
    When user opens the address "/func/Dbtests.PostgresByStringChoices?shipCountry=%22Germany%22"
    Then the "PostgresByStringChoices" view should be current
    And the table should have 122 rows
    And every value of "shipcountry" column should match "^Germany$"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An address that names nothing says so and leaves the shell as it was
    Given user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    When user opens the address "/p/no.such_project/none"
    Then an error balloon containing "Unable to get project asset" should have been shown
    And the "demog-1000" view should be current
    And no errors should have been logged
