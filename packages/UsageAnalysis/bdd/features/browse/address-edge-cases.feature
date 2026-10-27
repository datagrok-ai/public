@browse @edge-cases @realizes:views.browse
Feature: Address edge cases: a space, a folder, a function, a script and a project that are not there
  Guards for an address that names an entity the server does not have, or no longer has. The
  contract is the manual case Browse-Route-04 (browse_manual_tests2.md section 12): an address that
  names nothing shows a "not found" message or leads to a safe screen, the application does not
  break, and no error is shown in the page. browse-routing claims it for a project address with a
  name nobody has; here it is claimed for the other kinds of address, and for a project and a space
  that existed and were deleted.

  The address is followed inside the running app (`grok.shell.route`), as a pasted link is by the
  shell's router. A space address that names nothing leads to the Spaces list; a function address
  that names nothing says "Unable to find function" and keeps the view that was open; a deleted
  project's address says "Unable to get project asset" as a never-saved name does.
  Two scenarios describe defects. A folder address that names nothing opens an empty folder view
  with an "error" label and logs "File or Directory ... not exists"; a script address that names
  nothing opens an empty Script view with an "error" label and logs the server's 404. Both are
  core#1152.

  The project and the space are this feature's own, named with the run's time, and removed through
  the API before their address is followed; they are swept again at feature end.

  Background:
    Given user is logged in

  Scenario: A space address that names nothing leads to the Spaces list
    Given user opens demog-1000 dataset
    When user opens the address "/s/BDDEdgeNoSpace{time}"
    Then the "Spaces" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The address of a deleted space leads to the Spaces list
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

  Scenario: A function address that names nothing says so and keeps the open view
    Given user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    When user opens the address "/func/BddNoSuchPackage.BddNoSuchFunc{time}"
    Then an error balloon containing "Unable to find function" should have been shown
    And the "demog-1000" view should be current

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

  @known-failure:79 @realizes:core#1152
  Scenario: A folder address that names nothing says the folder is not there
    The defect: the address of a folder the file share does not have opens an empty folder view named
    after it, with an "error" label in the context panel and no message, and logs "File or Directory
    ... not exists" with a stack trace (core#1152). Remove the tag when such an address says the folder is not
    there without logging an error (if the fix leads to a safe screen instead, this claim changes).
    Given user opens demog-1000 dataset
    When user opens the address "/files/System.AppData/UsageAnalysis/bdd-edge-no-folder-{time}"
    Then an error or warning balloon matching "not found|not exist|n't exist" should have been shown
    And no errors should have been logged

  @known-failure:91 @realizes:core#1152
  Scenario: A script address that names nothing says the script is not there
    The defect: the address of a script the server does not have opens an empty, unnamed Script view
    with an "error" label and no message, and logs the server's 404 with a stack trace; the next
    address followed while that view is open throws "NullError: method not found" in the router (core#1152).
    Remove the tag when such an address says the script is not there without logging an error (if
    the fix leads to a safe screen instead, this claim changes).
    Given user opens demog-1000 dataset
    When user opens the address "/script/BddNoSuchScript{time}"
    Then an error or warning balloon matching "not found|not exist|n't exist" should have been shown
    And no errors should have been logged
