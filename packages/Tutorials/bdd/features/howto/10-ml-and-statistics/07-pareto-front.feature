@guide @tutorial
Feature: How-to · ML
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). ML > Pareto Front... (EDA, in the browser, no dialog) docks an
  Optimize panel on the left (one choice per numeric column: minimize, maximize or nothing), a Pareto front viewer on
  the right and a PC plot under the grid. On cars it starts by minimizing the last two numeric columns, highway.mpg
  and price; the clip clears highway.mpg and maximizes horsepower, so the front is the cars for which no other car is
  both cheaper and more powerful (green, "optimal": alfaromeo, chevrolet, dodge1, isuzu, jaguar, mercury, mitsubishi,
  porsche, subaru, toyota1, toyota3, computed from data/demo/cars.csv). Then what the front means, on three cars the
  pointer rests on (rows 12, 1 and 19 of cars, read from the file): mercedes, grey, 123 hp for 31,600, beaten by
  alfaromeo, 154 hp for 16,500 (52% of the price); porsche, the most powerful (207 hp) and the dearest (37,028), green
  as nothing beats its power. The PC plot draws every car as a line across the objectives and the first numeric
  columns (EDA shows at most 10 axes), colored by the same result.
  Setting the objectives in the viewer's own properties instead (Objectives > Minimize / Maximize) leaves the Optimize
  panel and the PC plot on the old ones (in the report). cars is opened off camera; nothing is made on the server.

  Scenario: Find the best trade-offs with a Pareto front
    A Pareto front marks the rows you cannot improve on one goal without giving up another.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Cheap or powerful: the Pareto front finds the best compromises" at the bottom
    # caption: Pareto Front is in the ML menu
    And user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    Then "price" input should have value "minimize"
    # caption: Optimize: minimize or maximize each column; leave mileage out
    When user selects "" in "highway.mpg" input
    # caption: Horsepower: higher is better. Price stays at minimize
    And user selects "maximize" in horsepower input
    Then "Minimize" property of pareto front viewer should be "price"
    And "Maximize" property of pareto front viewer should be "horsepower"
    And scatter plot viewer in pareto front viewer should be painted
    When the video says "Green: no other car is both cheaper and more powerful" at the bottom
    # caption: Grey Mercedes: 123 hp for 31,600
    And user hovers over the "marker of row 12" area of scatter plot viewer in pareto front viewer
    Then the "hovered row" reading of scatter plot viewer in pareto front viewer should be 12
    And the value of "model" column in row 12 should be "mercedes"
    And the value of "Pareto optimality" column in row 12 should be "non-optimal"
    # caption: Alfa Romeo beats it: more power for about half the price
    When user hovers over the "marker of row 1" area of scatter plot viewer in pareto front viewer
    Then the "hovered row" reading of scatter plot viewer in pareto front viewer should be 1
    And the value of "model" column in row 1 should be "alfaromeo"
    And the value of "Pareto optimality" column in row 1 should be "optimal"
    # caption: Porsche: the dearest, but no car has more power
    When user hovers over the "marker of row 19" area of scatter plot viewer in pareto front viewer
    Then the "hovered row" reading of scatter plot viewer in pareto front viewer should be 19
    And the value of "model" column in row 19 should be "porsche"
    And the value of "Pareto optimality" column in row 19 should be "optimal"
    And "horsepower" column should have its maximum in row 19
    And "price" column should have its maximum in row 19
    # caption: Below, one line per car across its columns, green on the front
    When user hovers over pc plot viewer
    Then no error or warning balloon should have been shown
    And no errors should have been logged
