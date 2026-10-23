@viewers @realizes:viewers.forms
Feature: Forms viewer — every selected row gets a card
  Step 4c of `forms.md` (Scenario 4), which `forms-renderers.feature` stops short of: with three
  rows selected (as the md has it — a fourth card would sit past the viewer's right edge), every
  selected row gets a card of its own. Cards are read by position: card 1 is the current row's,
  card 2 the mouse-over card (blank while nothing is hovered), and the selected rows follow from
  card 3. The md's molecule and curve fields are on the cards; what their renderers draw is not a
  feature (the Chem and Curves package tests), so only the cards and their records are claimed.

  Not translated: `forms-ui.md` "Molecules occupy the left part of the card". On spgi-100 with
  Renderer Size normal the card stacks its fields in one column (Structure at 182,5 200x100,
  Primary Series Name at 182,105 below it), so the md's "text fields to its right" is not the
  layout the viewer has; whether it should be is a question for the viewer's owner, not a claim.

  Scenario: With three rows selected, every selected row gets a card
    Given user is logged in
    And the "Chem" package is installed
    And the "Curves" package is installed
    And user opens curves dataset
    And user adds a forms viewer with:
      | fieldsColumnNames | smiles, multiple prefit |
    And user makes row 1 current
    When user selects the first 3 rows
    Then 3 rows should be selected
    And the "card kind of card 1" reading of forms viewer should be "current"
    And the "card kind of card 2" reading of forms viewer should be "mouse-over"
    And the "record of card 2" reading of forms viewer should be ""
    And the "record of card 3" reading of forms viewer should be "1"
    And the "record of card 4" reading of forms viewer should be "2"
    And the "record of card 5" reading of forms viewer should be "3"
    And no errors should have been logged
