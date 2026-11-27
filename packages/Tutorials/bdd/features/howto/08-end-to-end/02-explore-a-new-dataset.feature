@guide @tutorial
Feature: How-to · End to end
  An end-to-end how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films
  it as a demo clip that outlines what each step presses). One task from start to finish: getting to know cars (30
  rows), opened off camera and cut to 8 of its 17 columns, price among them (with all 17 the correlation plot ran off
  the right edge of its dock, hiding the cells the clip points at), the toolbox and the Context Panel open. The
  status bar gives the size; the price header and its Stats pane sum the column up; a histogram of price from the
  toolbox shows a gap and three cars beyond it (its 20 bins: bin 17 holds the Mercedes at 31600, bin 20 the Jaguar at
  35550 and the Porsche at 37028), and a selection box over bins 17 to 20 selects those three. A correlation plot
  points at engine size (r 0.91 with price; curb weight and horsepower 0.83), a double click on that cell adds a
  scatter plot of the pair with the three selected cars, the three biggest engines (258, 194, 183), at its top right,
  and Select > Selection to Filter shows only them. The cars table has no text column to group the prices by (model
  is one per row), so there is no box plot. Every move after the start is on camera. Nothing is made on the server.

  Scenario: Explore a new dataset
    Get to know a table in a few minutes: its numbers, its spread, the rows that stand out, and why.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And user removes "diesel" column
    And user removes "two.doors" column
    And user removes "hatchback" column
    And user removes "wheel.base" column
    And user removes "length" column
    And user removes "width" column
    And user removes "height" column
    And user removes "peak.rpm" column
    And user removes "symbol" column
    And the toolbox pane is shown
    And the context panel is open
    And the browse panel is closed
    When the video says "A table you have not seen before. Get to know it"
    # caption: The status bar: 30 rows, 8 columns
    And user hovers over the text "Columns: 8" in status bar
    # caption: Click a column: Stats gives its range and average
    And user clicks on the "header price" area of grid
    Then the context panel should show "price"
    When user expands Stats pane in context panel
    Then Stats pane in context panel should contain text "avg"
    When user clicks on "Context Panel" status bar toggle
    # caption: A histogram shows the spread: most cars are cheap
    And user clicks on "histogram" icon in toolbox
    Then histogram viewer should be visible
    When user picks "price" in the "value" column selector of histogram viewer
    # caption: Three cost far more. Shift+drag over their bars to select them
    And user drags a selection box from the "bin 17" area to the "bin 20" area of histogram viewer
    Then 3 rows should be selected
    And only rows where "model" is one of "jaguar, mercedes, porsche" should be selected
    When the video says "The grid marks the selected rows too"
    And user hovers over the "cell 10 of model" area of grid
    # caption: Why so expensive: a correlation plot shows what goes with price
    And user clicks on correlation plot icon on toolbox
    Then correlation plot viewer should be visible
    # caption: Blue: one falls as the other rises, as mileage does with price
    When user hovers over the "cell city.mpg x price" area of correlation plot viewer
    # caption: Red: they rise together. Engine size goes with price most
    And user hovers over the "cell eng.size x price" area of correlation plot viewer
    # caption: Double-click the cell to plot the two
    And user double-clicks on the "cell eng.size x price" area of correlation plot viewer
    Then the open tableview should have 1 scatter plot viewer
    When the video says "The three selected cars have the biggest engines: top right"
    And user hovers over the "marker of row 19" area of scatter plot viewer
    # caption: Selection to Filter shows only those rows, everywhere
    And user picks "Select > Selection to Filter" from the top menu
    Then 3 rows should pass the filter
    And every selected row should pass the filter
    When the video says "Jaguar, Mercedes and Porsche: the priciest, with the biggest engines"
    And user hovers over the "cell 10 of model" area of grid
    Then no errors should have been logged
