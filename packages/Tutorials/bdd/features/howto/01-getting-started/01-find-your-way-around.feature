@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The long tour of the window, every move filmed (the lead's
  "no jumps" round): Home and its search box, the sidebar's icons, then Browse with each section opened to show
  what is in it (My stuff's buckets, Spaces, the apps, the Dashboards gallery searched for the two activity-cliffs demo dashboards, the
  Biologics connection with its queries and tables, Platform and its Plugins view, Files > Demo), demog opened from
  Files > Demo, then the table view around it: the ribbon's menus and icons, the toolbox and its sections, the
  context panel on a column, the status bar's table name, column and row counts and the toggles at its bottom
  right, each named by its tooltip; last, the ribbon follows the data: smiles.csv (molecules) opened from
  Files > Demo > chem brings a Chem menu and peptides.csv (sequences) from Files > Demo > bio a Bio menu. The
  Biologics connection is the stand's own Postgres (no outside database). Nothing is made on the server.

  Scenario: Find your way around Datagrok
    Home, the sidebar and Browse, then an open table with its menus, toolbox, context panel and status bar.

    Given user is logged in
    And simple mode is off
    When the video says "Home: one search box for everything on the server"
    And user hovers over home search
    # caption: The sidebar opens the main panels
    And user hovers over browse tab
    # caption: Toolbox: tools for the table you work on
    And user hovers over toolbox tab
    # caption: Dashboards: the open dashboards and their tables
    And user hovers over "Dashboards" tab
    # caption: Favorites: what you starred
    And user hovers over favorites sidebar tab
    # caption: Help: documentation, community, feedback
    And user hovers over "Help" tab
    # caption: Settings: your preferences and the server's
    And user hovers over "Settings" tab
    # caption: Browse: everything on the server, in one tree
    And user clicks on browse tab
    Then browse panel should be visible
    # caption: My stuff: recent items, favorites, what others shared with you
    When user expands "My stuff" tree node inside browse tree
    And user collapses "My stuff" tree node inside browse tree
    # caption: Spaces: folders a team keeps together
    And user expands Spaces tree node inside browse tree
    And user collapses Spaces tree node inside browse tree
    # caption: Apps: the applications installed on the server
    And user expands Apps tree node inside browse tree
    And user collapses Apps tree node inside browse tree
    # caption: Dashboards: every saved dashboard, with a preview
    And user clicks on Dashboards tree node inside browse tree
    # caption: Find one by name
    And user enters "cliffs" into gallery search
    Then "demo_activity_cliffs" gallery card should be visible
    # caption: Databases: connections, grouped by database type
    When user expands Databases tree node inside browse tree
    And user expands Databases---PostgresDart tree node inside browse tree
    # caption: A connection holds its queries and its tables
    And user expands Databases---PostgresDart---Biologics tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics---Schemas tree node inside browse tree
    And user expands Databases---PostgresDart---Biologics---Schemas---biologics tree node inside browse tree
    And user collapses Databases tree node inside browse tree
    # caption: Platform: users, groups, plugins, functions, settings
    And user expands Platform tree node inside browse tree
    # caption: Plugins: what is installed on this server
    And user clicks on Platform---Plugins tree node inside browse tree
    And user collapses Platform tree node inside browse tree
    # caption: Files: file shares; Demo holds sample data
    And user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    # caption: Double-click a file to open it as a table
    And user double-clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    # caption: Edit: add rows and columns, find and replace
    When user hovers over "Edit" menu item
    # caption: View: layout, tooltips, presentation mode
    And user hovers over "View" menu item
    # caption: Select: all, random rows, duplicates, missing values
    And user hovers over "Select" menu item
    # caption: Data: aggregate, join, link, split
    And user hovers over "Data" menu item
    # caption: ML: clustering, dimensionality reduction, imputation
    And user hovers over "ML" menu item
    # caption: Save: keep the table and its viewers as a dashboard
    And user hovers over Save button in toolbar
    # caption: Export: download the table as a file
    And user hovers over "arrow-to-bottom" icon in toolbar
    # caption: Add viewer: the gallery of charts
    And user hovers over "Add viewer" icon in toolbar
    # caption: Filters: a filter for every column
    And user hovers over "Toggle filters" icon in toolbar
    # caption: Select all rows, none, or the others
    And user hovers over "All" icon in toolbar
    # caption: Add a column computed by a formula
    And user hovers over "Add New Column..." icon in toolbar
    Then toolbox should be visible
    # caption: On the left, the toolbox took Browse's place
    When user hovers over toolbox tab
    # caption: Search: find rows in this table
    And user hovers over Search accordion header inside toolbox
    # caption: Viewers: click one to add it to the view
    And user hovers over Viewers accordion header inside toolbox
    # caption: Layouts: saved arrangements of viewers
    And user hovers over Layouts accordion header inside toolbox
    # caption: Context panel: details of what you click
    And user clicks on "Context Panel" status bar toggle
    And user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    # caption: Status bar: the table's name, columns and rows
    When user hovers over "demog" text in status bar
    And user hovers over "Columns: 11" text in status bar
    And user hovers over "Rows: 5,850" text in status bar
    # caption: The toggles at the right show and hide panels
    And user hovers over "AI" status bar toggle
    # caption: Tabs: the tabs of the open views
    And user hovers over "Tabs" status bar toggle
    # caption: Toolbox: the panel on the left
    And user hovers over "Toolbox" status bar toggle
    # caption: Context Panel: the details on the right (F4)
    And user hovers over "Context Panel" status bar toggle
    # caption: Context Help: the help page of what you click (F1)
    And user hovers over "Context Help" status bar toggle
    # caption: Variables: values you set in the console
    And user hovers over "Variables" status bar toggle
    # caption: Console: type commands, read the log
    And user hovers over "Console" status bar toggle
    # caption: Presentation mode: only the view, no panels (F7)
    And user hovers over "Presentation mode" status bar toggle
    # caption: The menu follows the data: open a table of molecules
    And user clicks on browse tab
    And user expands Files---Demo---chem tree node inside browse tree
    And user double-clicks on Files---Demo---chem---smiles.csv tree node inside browse tree
    Then the "smiles" view should be current
    # caption: Molecules bring the Chem menu
    When user hovers over "Chem" menu item
    Then "Chem" menu item should be visible
    # caption: Now a table of peptide sequences
    When user clicks on browse tab
    And user expands Files---Demo---bio tree node inside browse tree
    And user double-clicks on Files---Demo---bio---peptides.csv tree node inside browse tree
    Then the "peptides" view should be current
    # caption: Sequences bring the Bio menu
    When user hovers over "Bio" menu item
    Then "Bio" menu item should be visible
    And no errors should have been logged
