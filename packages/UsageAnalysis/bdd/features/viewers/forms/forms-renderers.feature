@journey @viewers @realizes:viewers.forms
Feature: Forms viewer renderers, renderer size and the twenty-field cap
  A column with a semantic type and a cell renderer is drawn on a canvas inside the card instead of
  in an input, and the canvas is sized by the Renderer Size property. On spgi-100 `Structure` and
  `Core` carry Molecule; `Primary Series Name` is a plain string next to them, which is how "two
  molecular canvases per card" is claimed without counting DOM nodes.
  The Renderer Size ladder is claimed as MONOTONE GROWTH. The spec this replaces re-implemented the
  product's own sizing formula in the test — `w * 0.66`, `w`, `floor(w * 1.5)`, floored and
  multiplied by the device pixel ratio — and compared it with the canvas it had just derived it
  from, which can only fail when the arithmetic is retyped wrongly. What a user can see is that
  large is bigger than normal and normal bigger than small, and that is what is asserted.
  The twenty-field cap is silent: spgi-100 has far more than twenty columns and the viewer takes
  the first twenty in table order with no message and no balloon. `COLS_LIMIT_EXCEEDED_WARNING` is
  still declared in the source but nothing reads it, the message having been removed deliberately.
  The last scenario is the fit renderer on curves.csv, where the semantic type promotes Renderer
  Size to `normal` without anyone setting it — the same claim the old Step 4a made by re-deriving
  the ladder, made here by reading the property and the canvas box.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a forms viewer
    And user makes row 1 current
    Then forms viewer should be visible
    And "Structure" column should have semantic type "Molecule"
    And "Core" column should have semantic type "Molecule"

  Scenario: The field set stops at twenty columns, in table order, with no message
    Then the "fields shown" reading of forms viewer should be 20
    And the "fields" reading of forms viewer should be "Id, Structure, CAST Idea ID, Last Published Date, Chemist, Lab Notebook, Stereo Category, Series, Scaffold Names, Primary Series Name, Primary Scaffold Name, Has Unlabeled R-Groups, Core, R1, R2, R3, R100, R101, Chemical Space X, Chemical Space Y"
    And the "fields" and "header labels" readings of forms viewer should be the same
    And forms viewer should report no error
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A molecule column is drawn on a canvas and a plain string is not
    When user sets "fieldsColumnNames" property of forms viewer to "Structure, Core, Primary Series Name"
    Then the "field kind of Structure" reading of forms viewer should be "canvas"
    And the "field kind of Core" reading of forms viewer should be "canvas"
    And the "field kind of Primary Series Name" reading of forms viewer should be "input"
    And the "Structure of card 1" reading of forms viewer should be "canvas"
    And the "Primary Series Name of card 1" reading of forms viewer should be "Pyrrolidines"
    And forms viewer should have a "field Structure of card 1" area
    And no errors should have been logged

  Scenario: Renderer Size grows the canvas from small through normal to large
    Then "rendererSize" property of forms viewer should be "small"
    And the "width of Structure of card 1" reading of forms viewer should be at least 1
    When user sets "rendererSize" property of forms viewer to "normal"
    Then the "width of Structure of card 1" reading of forms viewer should be higher than before
    And the "height of Structure of card 1" reading of forms viewer should be higher than before
    When user sets "rendererSize" property of forms viewer to "large"
    Then the "width of Structure of card 1" reading of forms viewer should be higher than before
    And the "height of Structure of card 1" reading of forms viewer should be higher than before
    And the "field kind of Structure" reading of forms viewer should be "canvas"
    When user sets "rendererSize" property of forms viewer to "small"
    Then the "width of Structure of card 1" reading of forms viewer should be lower than before
    And the "height of Structure of card 1" reading of forms viewer should be lower than before
    And no errors should have been logged

  Scenario: Every card of a selection carries its own molecule canvas, and a filter keeps them
    The selection is kept to three rows on purpose: the card list is virtual, so a selection larger
    than the strip can hold leaves the rest of the cards out of the DOM, and "the record cards are
    exactly the selected rows that pass the filter" would then be false for a reason that is not
    the viewer's fault.
    When user sets "fieldsColumnNames" property of forms viewer to "Structure, Core, Primary Series Name"
    And user sets "rendererSize" property of forms viewer to "normal"
    And user selects the first 3 rows
    Then 3 rows should be selected
    And the record cards of forms viewer should show rows "1, 2, 3"
    And the record cards of forms viewer should be exactly the selected rows that pass the filter
    And every record card of forms viewer should show "canvas" in "Structure"
    And every record card of forms viewer should show "canvas" in "Core"
    When user adds a categorical filter on "Primary Series Name" keeping "Pyrrolidines"
    Then fewer than 100 rows should pass the filter
    And the record cards of forms viewer should be exactly the selected rows that pass the filter
    And every record card of forms viewer should show "canvas" in "Structure"
    And every record card of forms viewer should show "Pyrrolidines" in "Primary Series Name"
    When user hovers over "Primary Series Name" filter card
    And user clicks on close of "Primary Series Name" filter card
    And user clears the row selection
    Then 100 rows should pass the filter
    And no errors should have been logged

  Scenario: A fit column promotes Renderer Size to normal without anyone setting it
    Given user closes all views
    And user opens curves dataset
    And user adds a forms viewer with:
      | fieldsColumnNames | smiles, multiple prefit |
    And user makes row 1 current
    Then forms viewer should be visible
    And "multiple prefit" column should have semantic type "fit"
    And "rendererSize" property of forms viewer should be "normal"
    And the "field kind of multiple prefit" reading of forms viewer should be "canvas"
    And the "field kind of smiles" reading of forms viewer should be "canvas"
    And the "multiple prefit of card 1" reading of forms viewer should be "canvas"
    And the "width of multiple prefit of card 1" reading of forms viewer should be 200
    And the "height of multiple prefit of card 1" reading of forms viewer should be 100
    And forms viewer should report no error
    And no errors should have been logged
