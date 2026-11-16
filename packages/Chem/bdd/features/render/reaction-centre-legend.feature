@journey @realizes:chem.reactions.centre-legend @eval-b7
Feature: Every drawing that colours a reaction's centre names its colours in a legend under it
  The reaction-centre colours (green made, red broken, blue order changed) are kept, with a legend under every drawing;
  the app's colours once had no legend. The legend is one line under a picture: "Reaction centre:", then a swatch and
  words for each colour the picture has, right under it and centred on it. A picture that colours nothing has none. A
  route's steps are coloured as the route is and share its legend, once, right under the route's picture (it sat at the
  end, under the steps). The app's default route makes and breaks bonds and changes no bond's order, so its legend
  names green "bond made" and red "bond broken" only. Which
  bonds each step colours (step 1's new amide C-N bond green, no ring bond coloured) is Chem's package tests
  ('Reactions: centre'), from the table of every step's centre.

  Background:
    Given user is logged in
    And the Chemical Engine is "Crux"

  Scenario: The Reactions app's preview and its full-screen view name the default route's colours
    Given user opens Chem's Reactions app on its default example
    Then the Reactions app should name the centre's colours once, right under its route: "bond made, bond broken"
    When user clicks the Full screen icon of the Reactions app's preview
    Then the full-screen Reaction dialog should name the centre's colours once, right under its route: "bond made, bond broken"
    And no errors should have been logged

  Scenario: A PNG exported with the centre's colours carries their legend
    Given user opens Chem's Reactions app on its default example
    And user watches downloads
    When user clicks Export... in the Reactions app's ribbon
    And user chooses "Screen", "double column (17.8 cm)" and "PNG 300 dpi" in the Export panel
    Then the Export panel's preview should carry the centre's legend
    When user clicks Save in the Export panel
    Then a file "reaction.png" should have been downloaded
    And the saved PNG should have the centre's green and red in it
    And no errors should have been logged
