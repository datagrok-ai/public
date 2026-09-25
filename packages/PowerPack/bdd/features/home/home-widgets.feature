@journey @serial @realizes:powerpack.view.welcome @realizes:powerpack.dashboard.spotlight @realizes:powerpack.dashboard.community
Feature: The widgets of the Home page
  The Home page of an administrator: the search box, the Spotlight, Community, Usage and Reports
  widgets, how a widget is closed and brought back, what the page keeps across a reload, the search
  that replaces the widgets, and what each widget holds. Translated from TestTrack
  PowerPack/Widgets/home_widgets_manual_tests.md (cases Display-01, Controls-01/02, Customize-01/02,
  Search-01, Nav-01, Spotlight-01/02, Community-01, Usage-01, Reports-01, Perm-02); Perm-01 and
  Spotlight-03 need the second account and are in home-as-sharing-user.feature.

  The page lays the widgets out by their `order`, not in the order they are built: Spotlight first,
  then Reports, Usage and Community. The widget settings are read from the server, which is what a
  reload starts from, and every scenario that hides a widget brings it back through the Customize form;
  the settings the account had come back when the feature ends in any case.

  The tip at the bottom of Spotlight changes with the weekday — a demo on Monday and Friday, a
  tutorial on Wednesday, a plain tip on Tuesday and Thursday — so what the tip is claimed to open
  depends on the day the feature runs. The Community widget's links are checked by their addresses:
  the community site is not opened.

  What a Spotlight tab lists is the account's own history (pins, favorites, activity): each tab is
  claimed to show its own page, not particular items: the Spotlight page shows the account's recent
  items, or, for an account with none, interactive tutorials and demo apps to start with. The
  Workspace hint of the md ("Select a pinned item...") shows only when something is pinned. The
  System block of Usage lists the services by name (Health Scheduler, Garbage Collector, Credentials
  Server, Core, Jupyter, Grok Spawner, Grok Connect), not as links, and no Datlas entry.

  Not translated: the tooltip "Remove" of the close icon (the icon carries it as its aria label, which
  is how the scenario finds it); "the widget scrolls through recent reports" (the Reports list is
  claimed to be there, its scrolling is not).

  A page keeps the widget settings it read when it loaded, and writes them all back to the server once
  its widgets are built: a page of the same account that loads while this feature has a widget hidden
  or shown again writes the older setting back. So the feature signs in, on its page, as an
  administrator account of its own ("bddhomeadmin", a member of Administrators, made once per stand),
  whose settings no other page reads or writes, and every claim about what is stored is read from the
  server after the gesture that changed it. The running account comes back when the feature ends. It
  runs @serial with the other features that sign in on their page.

  Background:
    Given user is logged in
    And an administrator account "bddhomeadmin" is on the server
    And user is signed in as "bddhomeadmin" on this page
    And the widget settings of the Home page come back when the feature ends

  Scenario: The Home page shows the search box and four loaded widgets after a reload
    When user reloads the page
    Then home search should be visible
    And the placeholder of home search should start with "Search everywhere"
    And the Home page should show the widgets "Spotlight, Reports, Usage, Community"
    And every widget of the Home page should show content
    And title of Community home widget should have text "Community"
    And title of Spotlight home widget should be hidden
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The close icon shows on hover and removes the widget
    Then close icon of Community home widget should be hidden
    When user hovers over Community home widget
    Then close icon of Community home widget should be visible
    When user clicks on close icon of Community home widget
    Then Community home widget should be absent
    And the Home page should show the widgets "Spotlight, Reports, Usage"
    And the Community widget should be stored as hidden
    When user clicks on "Customize widgets..." link
    And user checks "Community" input in context panel
    Then Community home widget should be visible
    And the Community widget should be stored as shown
    And no errors should have been logged

  Scenario: A closed widget stays closed after a reload
    When user hovers over Community home widget
    And user clicks on close icon of Community home widget
    Then the Community widget should be stored as hidden
    When user reloads the page
    Then Community home widget should be absent
    And the Home page should show the widgets "Spotlight, Reports, Usage"
    When user clicks on "Customize widgets..." link
    Then "Community" input in context panel should not be checked
    When user checks "Community" input in context panel
    Then the Community widget should be stored as shown
    When user reloads the page
    Then the Home page should show the widgets "Spotlight, Reports, Usage, Community"
    And no errors should have been logged

  Scenario: The Customize form hides and shows a widget
    When user clicks on "Customize widgets..." link
    Then the following elements should be visible:
      | "Spotlight" input in context panel |
      | "Community" input in context panel |
      | "Usage" input in context panel     |
      | "Reports" input in context panel   |
    And "Community" input in context panel should be checked
    When user unchecks "Community" input in context panel
    Then Community home widget should be absent
    And the Community widget should be stored as hidden
    When user checks "Community" input in context panel
    Then Community home widget should be visible
    And the Community widget should be stored as shown
    And no errors should have been logged

  Scenario: A widget hidden in the Customize form stays hidden after a reload
    When user clicks on "Customize widgets..." link
    And user unchecks "Community" input in context panel
    Then the Community widget should be stored as hidden
    When user reloads the page
    Then Community home widget should be absent
    When user clicks on "Customize widgets..." link
    Then "Community" input in context panel should not be checked
    When user checks "Community" input in context panel
    Then Community home widget should be visible
    And the Community widget should be stored as shown
    And no errors should have been logged

  Scenario: A search replaces the widgets, and clearing it brings them back
    When user types "aspirin" into home search
    Then home widgets panel should be hidden
    And home search results should be visible
    And the page address should contain "search?q=aspirin"
    And the search should have finished
    When user clears home search
    Then home widgets panel should be visible
    And home search results should be hidden
    And the page address should not contain "?q="
    And no errors should have been logged

  Scenario: The Home icon brings back the same widgets, before and after a reload
    Given user opens demog dataset
    And the browse panel is open
    Then the "demog" view should be current
    When user clicks on "Home" icon in browse toolbar
    Then the "Home" view should be current
    And the Home page should show the widgets "Spotlight, Reports, Usage, Community"
    When user reloads the page
    Then the Home page should show the widgets "Spotlight, Reports, Usage, Community"
    And no errors should have been logged

  Scenario: Spotlight has six tabs, and each shows its own page
    Then the following elements should be visible:
      | Workspace tab in Spotlight home widget     |
      | Spotlight tab in Spotlight home widget     |
      | Favorites tab in Spotlight home widget     |
      | Notifications tab in Spotlight home widget |
      | "My Activity" tab in Spotlight home widget |
      | Learn tab in Spotlight home widget         |
    When user clicks on Spotlight tab in Spotlight home widget
    Then the "Spotlight" tab of Spotlight home widget should be showing
    And spotlight page of Spotlight home widget should contain one of the texts "Recent | Interactive Tutorials"
    When user clicks on Favorites tab in Spotlight home widget
    Then the "Favorites" tab of Spotlight home widget should be showing
    When user clicks on Notifications tab in Spotlight home widget
    Then the "Notifications" tab of Spotlight home widget should be showing
    When user clicks on "My Activity" tab in Spotlight home widget
    Then the "My Activity" tab of Spotlight home widget should be showing
    When user clicks on Learn tab in Spotlight home widget
    Then the "Learn" tab of Spotlight home widget should be showing
    And the following elements should be visible:
      | VIDEO tab in Spotlight home widget     |
      | WIKI tab in Spotlight home widget      |
      | DEMO tab in Spotlight home widget      |
      | TUTORIALS tab in Spotlight home widget |
    And Spotlight home widget should contain text "Cheminformatics"
    When user clicks on WIKI tab in Spotlight home widget
    Then the "WIKI" tab of Spotlight home widget should be showing
    When user clicks on DEMO tab in Spotlight home widget
    Then the "DEMO" tab of Spotlight home widget should be showing
    When user clicks on TUTORIALS tab in Spotlight home widget
    Then the "TUTORIALS" tab of Spotlight home widget should be showing
    When user clicks on VIDEO tab in Spotlight home widget
    Then the "VIDEO" tab of Spotlight home widget should be showing
    When user clicks on Workspace tab in Spotlight home widget
    Then the "Workspace" tab of Spotlight home widget should be showing
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The tip of the day at the bottom of Spotlight opens what it names
    Then the tip of the day of Spotlight home widget should open what it names
    And no errors should have been logged

  Scenario: Community lists links to the community site
    Then every link of Community home widget should point to "https://community.datagrok.ai/"
    And Community home widget should contain text "Platform Releases"

  Scenario: Usage shows its user and error charts and the state of the services
    Then Usage home widget should contain text "Users"
    And Usage home widget should contain text "Errors"
    And Usage home widget should contain text "System"
    And Usage home widget should contain text "Jupyter"
    And Usage home widget should contain text "Grok Spawner"
    And Usage home widget should contain text "Grok Connect"
    And there should be 2 visible line chart viewer in Usage home widget
    And the "lines" reading of first line chart viewer in Usage home widget should be 1
    And the "rows shown" reading of first line chart viewer in Usage home widget should be at least 1
    And the "lines" reading of second line chart viewer in Usage home widget should be 1
    And the "rows shown" reading of second line chart viewer in Usage home widget should be at least 1
    And "Open Usage Analysis" link in Usage home widget should be visible
    And no errors should have been logged

  Scenario: Open Usage Analysis in the Usage widget opens the app's Overview view, with nothing logged
    When user clicks on "Open Usage Analysis" link in Usage home widget
    Then the "Overview" view should be open
    And the "Home" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # The click, bubbling up to the Home view's root, makes Home current again: the Reports widget
  # stops that propagation, the Usage widget does not.
  @known-failure
  Scenario: The Overview view Open Usage Analysis opened is brought to the front
    Then the "Overview" view should be current

  Scenario: Reports lists recent reports and opens the reports view
    Then "Open Reports" link in Reports home widget should be visible
    When user clicks on "Open Reports" link in Reports home widget
    Then the "Reports" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown
