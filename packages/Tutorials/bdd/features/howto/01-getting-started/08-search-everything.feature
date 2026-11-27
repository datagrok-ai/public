@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Home page's search box looks through the whole server, in one
  clip (it absorbed the former "Search everything from one box", "Ask a question in the search box" and "Look up an
  ID" clips): "scatter plot" brings functions, demos, help pages and samples, each kind under its own header (a help
  link's tooltip says how it opens); "Threshold" (in no query's search pattern, so no suggestion list covers the
  results) brings dashboards, files and queries with that word in the name, a dashboard previews on hover and opens on
  a double-click, and the house tab goes back to Home; a query with a search pattern (`--meta.searchPattern`) answers a
  typed sentence: "biologics" makes the box suggest the sentences (one is pointed at so the list stays in view; its
  items are the Dart popup menu PowerPack shows under the box), "biologics assays for E. coli" runs the query and
  shows its table under the box; an identifier a plugin registers (GROKADC- for ADCs, GROKMOL- for drugs) shows its
  record as a card, several IDs a card each, and "Open table" puts them in a new table whose IDs stay links (a click
  shows the record in the context panel). The dashboards, queries and IDs are the Biologics package's, on the
  stand's own Postgres (no outside service; a four-letter query would also ask a public protein database for a PDB
  ID, so none is typed). Nothing is made on the server: the dashboard is only viewed, the new table never saved.

  Scenario: Search everything from Home
    One search box on Home finds functions, help, dashboards, files and queries, answers questions and looks up IDs.

    Given user is logged in
    And simple mode is off
    When the video says "The Home search box looks through the whole server"
    And user clicks on home search
    # caption: Type what you look for: results come by kind
    And user types "scatter plot" into home search
    Then the Home search should have finished
    And the Home search results should list the kinds "Functions, Demos, Help, Samples"
    # caption: Functions, demos, code samples and help pages
    When user hovers over "Scatterplot" link inside home search results
    # caption: A word from a name finds dashboards, files and queries
    And user types "Threshold" key by key into home search
    Then the Home search should have finished
    And the Home search results should list the kinds "Projects, Files, Queries"
    # caption: Point at a dashboard to preview it
    When user hovers over second "ADCs by IC50 Threshold" link inside home search results
    # caption: Double-click to open it
    And user double-clicks on second "ADCs by IC50 Threshold" link inside home search results
    Then the "ADCs by IC50 Threshold" view should be current
    # caption: The house tab takes you back to Home
    When user clicks on the tab of the "Home" view
    Then the "Home" view should be current
    # caption: Some queries answer a sentence. Start typing one
    When user types "biologics" key by key into home search
    # caption: The box suggests the sentences it knows
    And user hovers over "biologics assays for ${organism}" menu item inside context menu
    # caption: Finish the sentence and the query runs
    And user types "biologics assays for E. coli" key by key into home search
    Then the Home search should have finished
    And the table the Home search answered with should have loaded
    # caption: Its answer: the assays run in E. coli
    When user hovers over home search results
    # caption: Type an ID to see its record
    And user types "GROKADC-000001" into home search
    Then the Home search should have finished
    When user hovers over home search results
    # caption: Several IDs, separated by commas, show a card each
    And user types "GROKMOL-000001, GROKMOL-000002, GROKMOL-000003" into home search
    Then the Home search should have finished
    When user hovers over home search results
    # caption: Open table puts them in a table
    And user clicks on "Open table" link inside home search results
    Then the table should have 3 rows
    # caption: The IDs stay links: click one to see its record
    When user clicks on "Context Panel" status bar toggle
    And user clicks on the "cell 2 of DG_BIOLOGICS_DRUG_ID" area of grid
    Then context panel should be visible
    And no errors should have been logged
