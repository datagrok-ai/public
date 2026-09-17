@journey @viewers @realizes:viewers.matrix-plot
Feature: Matrix plot — the cells, the plots inside them and the rows they keep
  One cell per pair of numerical columns, a density or scatter plot in each off-diagonal one and a
  histogram on the diagonal whatever Cell Plot Type says; each cell holds its own copy of the frame,
  so a blank in either of its columns keeps a row out of that cell alone.
  The four specs this replaces shared `matrix-helpers.ts`, 106 lines whose whole job was to decide
  when a cell had stopped painting: `settledCellInk` polled `getImageData` over twelve rounds 300 ms
  apart until two readings agreed within 40, because the plot fired VIEWER_RENDERED before its inner
  viewers had drawn anything. That is fixed in the core, and a cell's picture is now
  `cell signature HEIGHT x AGE` — a hash of the frame it actually drew — so "the cell repainted" is
  an inequality of two numbers and "Density to Scatter and back returns to the same picture" is an
  equality, where the old spec could only demand `|delta| > 500` and `|back - dens| < 40`.
  The cells are also addressed by name rather than by `querySelectorAll(...)[1]`, which is what made
  the old wheel-zoom scenario's "the neighbour did not move" a claim about index 2.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a matrix plot viewer
    Then 1000 rows should pass the filter
    And the "rows shown" reading of matrix plot viewer should be 1000
    And the "cells" reading of matrix plot viewer should be 16
    And the cells of matrix plot viewer should be 4 wide and 4 tall
    And the "cell viewer type" reading of matrix plot viewer should be "Density plot"
    And the "cells drawn" reading of matrix plot viewer should be 16

  Scenario: One cell per pair of numerical columns, each with a picture of its own
    Then the "x columns" reading of matrix plot viewer should be 4
    And the "y columns" reading of matrix plot viewer should be 4
    And matrix plot viewer should have a "cell HEIGHT x AGE" area
    And matrix plot viewer should have a "cell body HEIGHT x AGE" area
    And matrix plot viewer should have a "cell 1,1" area
    And matrix plot viewer should not have a "cell SEX x AGE" area
    And the "blank cells" reading of matrix plot viewer should be 0
    And the "distinct cell signatures" reading of matrix plot viewer should be 16
    And the "error" reading of matrix plot viewer should be ""
    And no errors should have been logged

  Scenario: The diagonal hosts a histogram whatever Cell Plot Type says
    Then the "cell viewer type of AGE x AGE" reading of matrix plot viewer should be "Histogram"
    And the "cell viewer type of WEIGHT x WEIGHT" reading of matrix plot viewer should be "Histogram"
    And the "cell viewer type of HEIGHT x AGE" reading of matrix plot viewer should be "Density plot"
    When user sets "cellPlotType" property of matrix plot viewer to "Scatter plot"
    Then the "cell viewer type of HEIGHT x AGE" reading of matrix plot viewer should be "Scatter plot"
    And the "cell viewer type of AGE x AGE" reading of matrix plot viewer should be "Histogram"
    And the "cell viewer type of WEIGHT x WEIGHT" reading of matrix plot viewer should be "Histogram"
    When user sets "cellPlotType" property of matrix plot viewer to "Density plot"
    Then the "cell viewer type of HEIGHT x AGE" reading of matrix plot viewer should be "Density plot"
    And no errors should have been logged

  Scenario: Cell Plot Type redraws every off-diagonal cell and the round trip returns to the picture it started from
    Then the "distinct cell signatures" reading of matrix plot viewer should be 16
    When user remembers the "cell signature HEIGHT x AGE" reading of matrix plot viewer
    And user sets "cellPlotType" property of matrix plot viewer to "Scatter plot"
    Then the "cell signature HEIGHT x AGE" reading of matrix plot viewer should not be as remembered
    And the "cell signature WEIGHT x HEIGHT" reading of matrix plot viewer should differ from before
    And the "cell signature AGE x AGE" reading of matrix plot viewer should be the same as before
    And the "cells" reading of matrix plot viewer should be 16
    And the "blank cells" reading of matrix plot viewer should be 0
    When user sets "cellPlotType" property of matrix plot viewer to "Density plot"
    Then the "cell signature HEIGHT x AGE" reading of matrix plot viewer should be as remembered
    And no errors should have been logged

  Scenario: Each cell keeps its own rows, so a blank in either column drops the row from that cell only
    Then the "rows shown" reading of matrix plot viewer should be 1000
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 872
    And the "cell rows shown of WEIGHT x AGE" reading of matrix plot viewer should be 1000
    And the "cell rows shown of AGE x AGE" reading of matrix plot viewer should be 1000
    And the "cell rows shown of STARTED x HEIGHT" reading of matrix plot viewer should be 872
    And no errors should have been logged

  Scenario: Narrowing X re-tiles the grid and takes the labels with it
    When user sets "xColumnNames" property of matrix plot viewer to "AGE, HEIGHT"
    Then the "cells" reading of matrix plot viewer should be 8
    And the cells of matrix plot viewer should be 2 wide and 4 tall
    And the "cells drawn" reading of matrix plot viewer should be 8
    And matrix plot viewer should not have a "cell WEIGHT x AGE" area
    And matrix plot viewer should have a "x label HEIGHT" area
    And matrix plot viewer should not have a "x label WEIGHT" area
    And matrix plot viewer should have a "y label WEIGHT" area
    When user sets "xColumnNames" property of matrix plot viewer to "AGE, HEIGHT, WEIGHT, STARTED"
    Then the "cells" reading of matrix plot viewer should be 16
    And matrix plot viewer should have a "cell WEIGHT x AGE" area
    And no errors should have been logged

  Scenario: Cycling the column sets ends where it started, with no error (GROK-16473)
    When user sets "xColumnNames" property of matrix plot viewer to "AGE, HEIGHT, WEIGHT"
    Then the "cells" reading of matrix plot viewer should be 12
    When user sets "yColumnNames" property of matrix plot viewer to "AGE, HEIGHT"
    Then the "cells" reading of matrix plot viewer should be 6
    And the cells of matrix plot viewer should be 3 wide and 2 tall
    When user sets properties of matrix plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of matrix plot viewer should be 16
    And the "cells drawn" reading of matrix plot viewer should be 16
    And the "blank cells" reading of matrix plot viewer should be 0
    And the "error" reading of matrix plot viewer should be ""
    And no errors should have been logged

  Scenario: The viewer's own filter narrows every cell and leaves the table alone
    When user sets "filter" property of matrix plot viewer to "${AGE} > 30"
    Then the "rows shown" reading of matrix plot viewer should be 844
    And 1000 rows should pass the filter
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 745
    And the "cell rows shown of AGE x AGE" reading of matrix plot viewer should be 844
    And the "cell signature HEIGHT x AGE" reading of matrix plot viewer should differ from before
    When user sets "filter" property of matrix plot viewer to ""
    Then the "rows shown" reading of matrix plot viewer should be 1000
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 872
    And no errors should have been logged

  Scenario: Row Source Selected redraws every cell over the selection alone
    Then the "distinct cell signatures" reading of matrix plot viewer should be 16
    When user remembers the "cell signature HEIGHT x AGE" reading of matrix plot viewer
    And user selects the first 50 rows
    Then 50 rows should be selected
    When user sets "rowSource" property of matrix plot viewer to "Selected"
    Then the "rows shown" reading of matrix plot viewer should be 50
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 50
    And the "cell rows shown of AGE x AGE" reading of matrix plot viewer should be 50
    And the "cell signature HEIGHT x AGE" reading of matrix plot viewer should not be as remembered
    When user sets "rowSource" property of matrix plot viewer to "Filtered"
    And user clears the row selection
    Then the "rows shown" reading of matrix plot viewer should be 1000
    And the "cell signature HEIGHT x AGE" reading of matrix plot viewer should be as remembered
    And no errors should have been logged

  Scenario: The cell tooltip names the pair the cell stands for
    When user hovers over the "cell HEIGHT x AGE" area of matrix plot viewer
    Then exactly one tooltip should be shown
    And tooltip should contain the text "X: HEIGHT"
    And tooltip should contain the text "Y: AGE"
    When user moves the pointer away from matrix plot viewer
    Then no errors should have been logged

  Scenario: The expand icon shows on hover and opens the cell as the viewer it hosts
    When user moves the pointer away from matrix plot viewer
    Then matrix plot viewer should not have an "expand icon of HEIGHT x AGE" area
    And the open tableview should have 0 density plot viewers
    When user hovers over the "cell HEIGHT x AGE" area of matrix plot viewer
    Then matrix plot viewer should have an "expand icon of HEIGHT x AGE" area
    When user clicks on the "expand icon of HEIGHT x AGE" area of matrix plot viewer
    Then the open tableview should have 1 density plot viewer
    When user clicks on close icon of density plot viewer
    Then the open tableview should have 0 density plot viewers
    When user hovers over the "cell AGE x AGE" area of matrix plot viewer
    And user clicks on the "expand icon of AGE x AGE" area of matrix plot viewer
    Then the open tableview should have 1 histogram viewer
    When user clicks on close icon of histogram viewer
    Then the open tableview should have 0 histogram viewers
    And user moves the pointer away from matrix plot viewer
    Then no errors should have been logged

  Scenario: The wheel zooms the cell under the pointer and leaves its neighbour alone
    # A cell that still owes a frame reports NO signature, so the grid is asked for all sixteen
    # first: the previous scenario opened and closed a standalone viewer, and the re-tile it caused
    # left eleven of the cells unpainted at the moment the next step read them (seen once in three
    # runs). This line polls until every cell has a picture again; it is not a wait for a signal
    # nobody sends, it is the signal.
    Then the "distinct cell signatures" reading of matrix plot viewer should be 16
    When user remembers the "cell signature HEIGHT x AGE" reading of matrix plot viewer
    And user scrolls the mouse wheel up over the "cell body HEIGHT x AGE" area of matrix plot viewer
    Then the "cell signature HEIGHT x AGE" reading of matrix plot viewer should not be as remembered
    And the "cell signature WEIGHT x AGE" reading of matrix plot viewer should be the same as before
    When user scrolls the mouse wheel down over the "cell body HEIGHT x AGE" area of matrix plot viewer
    Then the "cell signature HEIGHT x AGE" reading of matrix plot viewer should be as remembered
    When user moves the pointer away from matrix plot viewer
    Then no errors should have been logged
