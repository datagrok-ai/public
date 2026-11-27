@guide @tutorial
Feature: How-to · Databases
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). A query's parameter is an `--input:` line in its header (type, name and a
  default value) used in the SQL as `@name`. The clip opens New Query... on the Biologics connection, names the query,
  types the parameter line and a select over biologics.assay_types that keeps one category, and runs it: the play icon
  asks for the category in a dialog (Developability by default); the clip types Legacy (9 assay types; Developability
  has 4) and the result shows under the editor. Nothing is saved: the editor's view is closed by the shell reset.

  Scenario: Add a parameter to a query
    Turn a value in your SQL into a parameter that is asked for each time the query runs.

    Given user is logged in
    And simple mode is off
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    When the video says "A parameter makes one query answer many questions"
    And user expands Databases tree node inside browse tree
    And user expands Databases---PostgresDart tree node inside browse tree
    # caption: A new query on the connection
    And user picks "New Query..." from the context menu of Databases---PostgresDart---Biologics tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "Assay types in a category" into Name input
    And user replaces the code of code editor with '--input: string category = "Developability"'
    # caption: An --input line declares the parameter and its default
    And user hovers over the text "--input" in code editor
    And user appends "select name, description from biologics.assay_types where category = @category" to code editor
    # caption: The SQL uses it as @category
    And user hovers over the text "@category" in code editor
    # caption: Run it: it asks for the value
    And user runs the query with the play icon
    Then "Assay types in a category" dialog should be visible
    And Category input in "Assay types in a category" dialog should have value "Developability"
    # caption: It starts with the default, Developability
    When user hovers over Category input in "Assay types in a category" dialog
    # caption: Type another category
    And user enters "Legacy" into Category input in "Assay types in a category" dialog
    And user clicks on OK button in "Assay types in a category" dialog
    Then the "Assay types in a category" dialog should close
    And grid should be visible
    And the "rows" reading of grid should be 9
    When the video says "Only the Legacy assay types come back" at the bottom
    And user hovers over grid
    Then no errors should have been logged
