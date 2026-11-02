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

  An address that names nothing else leaves the shell working too. A space address, and the address
  of a space deleted since, lead to the Spaces list with no error; a function address says "Unable to
  find function" and keeps the open view; an app address (`/apps/<package>/<app>`, `routing.dart`
  `_appsHandler`) says "Application not found" and opens the Apps view; the address of a project
  deleted since says "Unable to get project asset", as a never-saved name does. A folder or a script
  address that names nothing is expected to say so; today the folder opens an empty folder view with
  an "error" label and logs "File or Directory ... not exists", and the script opens an empty Script
  view and logs the server's 404 (core#1152, a known failure). The space is this feature's own,
  named with the run's time and removed through the API before its address is followed; the project
  is removed the same way, and both are swept again at feature end.

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

  Scenario: The address of a deleted project says so and keeps the open view
    Given no project named "bddedgeaddrproj{time}" is on the server
    And user opens demog-1000 dataset
    When user saves the current view as project "bddedgeaddrproj{time}"
    And user closes all views
    And user opens the "bddedgeaddrproj{time}" project
    Then the page address should contain "bddedgeaddrproj{time}"
    When user remembers the page address
    And user closes all views
    Given the user's own project "bddedgeaddrproj{time}" is removed now and at feature end
    And user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    When user opens the remembered address
    Then an error balloon containing "Unable to get project asset" should have been shown
    And the "demog-1000" view should be current
    And no errors should have been logged

  Scenario: A folder address opens that folder, and a space address that names nothing or a deleted space leads to the Spaces list
    When user opens the address "/files/System.DemoFiles/chem"
    Then the "Demo/chem" view should be current
    And the gallery counter should show as many items as the "System:DemoFiles/chem/" folder holds on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user opens the address "/s/BDDEdgeNoSpace{time}"
    Then the "Spaces" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown
    Given the browse panel is open
    And a space named "BDDEdgeAddrSpace{time}" is on the server
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDEdgeAddrSpace{time} tree node inside browse tree
    Then the "BDDEdgeAddrSpace{time}" view should be current
    And the page address should contain "/s/BDDEdgeAddrSpace{time}"
    When user remembers the page address
    And user closes all views
    Given no space named "BDDEdgeAddrSpace{time}" is on the server
    When user opens the remembered address
    Then the "Spaces" view should be current
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

  Scenario: An address that names no project or function says so and keeps the open view, and one that names no app opens the Apps view
    Given user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    When user opens the address "/p/no.such_project/none"
    Then an error balloon containing "Unable to get project asset" should have been shown
    And the "demog-1000" view should be current
    And no errors should have been logged
    When user opens the address "/func/BddNoSuchPackage.BddNoSuchFunc{time}"
    Then an error balloon containing "Unable to find function" should have been shown
    And the "demog-1000" view should be current
    When user opens the address "/apps/BddNoSuchApp{time}"
    Then an error balloon containing "Application not found" should have been shown
    And the "Apps" view should be current

  Scenario Outline: A <kind> address that names nothing says the <kind> is not there
    The address of a folder the file share does not have opens an empty folder view named after it,
    with an "error" label in the context panel and no message, and logs "File or Directory ... not
    exists"; the address of a script the server does not have opens an empty, unnamed Script view
    and logs the server's 404, and the next address followed throws "NullError: method not found"
    in the router (core#1152).
    Given user opens demog-1000 dataset
    When user opens the address "<address>"
    Then an error or warning balloon matching "not found|not exist|n't exist" should have been shown

    @known-failure:123 @realizes:core#1152
    Examples:
      | kind   | address                                                       |
      | folder | /files/System.AppData/UsageAnalysis/bdd-edge-no-folder-{time} |
      | script | /script/BddNoSuchScript{time}                                 |
