@sketcher-controls @serial
Feature: Apps | Chem | Sketcher: Crux filling its view, the drawing the current object, the address its link, and Export
  The owner (2026-11-19): a regular chem sketcher app under Apps | Chem, the whole view the sketcher's (no scrollbars,
  no other inputs), the query options on, an example preloaded and the context panel open, so that a change of the
  drawing makes it the current object and the panel renders its panes; opened from a URL whose argument is the
  drawing's CXSMILES without coordinates, which drawing rewrites through the platform's own mechanism (the app's view
  path, which core writes into the address), and an Export dialog like the reactions' (SVG, PNG and CDXML in the ACS
  1996 and as-drawn styles; MOL V2000 and V3000, SDF, SMILES, CXSMILES, SMARTS, InChI, InChIKey), copied or downloaded.

  The owner, later that day: "the sketcher app by default should be in normal mode, not query mode": the app opens in
  molecule mode, the query options a click away (the gear's Query mode); a link holding a query opens in query mode
  (main's call: a shared query reopens as one).

  The session's sketcher is pinned to OpenChemLib: the app's sketcher host has its own, Crux (js-api's
  `{sketcherType: 'Crux', showInput: false}`: no input row), and leaves the session's as it was. The app opens the
  context panel and shuts the help panel, as the Reaction Sketcher does (the owner, 2026-11-19). The feature drives
  Crux's own controls (@sketcher-controls): a run that pins another sketcher for the whole run skips it.
  What each export holds exactly, the link's rules (a query as CXSMARTS or SMARTS, an empty drawing, a malformed or too
  long link) and the address changed in place are Chem's package tests ('Sketcher: app').

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the molecule sketcher is "OpenChemLib"

  Scenario: Opened from Apps: the example in Crux, molecule mode, the context panel on it, the help panel shut, the address its link
    Given the browse panel is open
    And the help panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Chem tree node inside browse tree is expanded
    When user clicks on Apps---Chem---Sketcher tree node inside browse tree
    Then the "Sketcher" view should be current
    And the "ready" reading of crux sketcher widget should be "true"
    And the "mode" reading of crux sketcher widget should be "molecule"
    And the "smiles" reading of crux sketcher widget should be the molecule "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
    And the current object should be the molecule "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
    And the context panel should show the "Chemistry" pane
    And help panel should be hidden
    And the address should be the Sketcher app's
    And the address's link should be "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
    And the Sketcher app should fill its view with Crux, no scrollbar and no field
    And the session's sketcher should be "OpenChemLib"
    And no errors should have been logged

  Scenario: A bond drawn makes the drawing the current object and rewrites the link; a cleared canvas leaves it empty
    Given user opens Chem's Sketcher app
    And the "ready" reading of crux sketcher widget should be "true"
    And the molecules made the current object are counted
    When user clicks on crux single bond tool
    And user clicks on Crux canvas 6% across and 8% down
    Then a molecule should have become the current object 1 time, the molecule "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O.CC"
    And the address's link should be the molecule "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O.CC"
    When user clicks on crux clear button
    Then the address's link should be ""
    And no errors should have been logged

  Scenario: A link with a molecule, followed in a new page load, opens on that molecule, its AND group kept
    When user follows the link "/apps/Chem/Sketcher?smiles=C%5BC%40H%5D(O)%5BC%40%40H%5D(C)Cl%20%7C%261%3A1%2C3%7C" in a new page load
    Then the "ready" reading of crux sketcher widget should be "true"
    And the "smiles" reading of crux sketcher widget should be the molecule "C[C@H](O)[C@@H](C)Cl"
    And the address's link should be "C[C@H](O)[C@@H](C)Cl |&1:1,3|"
    And the current object should be the molecule "C[C@H](O)[C@@H](C)Cl"
    And the "mode" reading of crux sketcher widget should be "molecule"
    And no errors should have been logged

  # The sketcher gives back the query as the link gave it until the first edit (the host's caller's string); the address
  # holds Crux's own SMARTS of the drawing.
  Scenario: A link with a query opens on it in query mode, and the address keeps it as SMARTS
    When user opens the address "/apps/Chem/Sketcher?smiles=%5B%236%3BR%5D-%5B%237%2C%238%5D"
    Then the "Sketcher" view should be current
    And the "ready" reading of crux sketcher widget should be "true"
    And the Crux sketcher should hold the query "[#6;R]-[#7,#8]"
    And the "query" reading of crux sketcher widget should be "true"
    And the "mode" reading of crux sketcher widget should be "query"
    And the address's link should be "[#6&R]-[#7,#8]"
    And no errors should have been logged

  # Main's F1 (run19): a query's link is its CXSMARTS where the SMARTS in it is Crux's own, so what SMARTS cannot hold
  # travels with the query; the link is written from the drawing, so it holds what Crux drew.
  Scenario: A link with a query holding an AND group opens on it in query mode, and the address keeps the group
    When user opens the address "/apps/Chem/Sketcher?smiles=%5B%236%5D-%5B%236%40H%5D(-%5B%238%5D)-%5B%237%2C%238%5D%20%7C%261%3A1%7C"
    Then the "Sketcher" view should be current
    And the "ready" reading of crux sketcher widget should be "true"
    And the "mode" reading of crux sketcher widget should be "query"
    And the address's link should be "[#6]-[#6@&H1](-[#8])-[#7,#8] |&1:1|"
    And no errors should have been logged

  Scenario: A link with a query holding an attachment point opens on it in query mode, and the address keeps the point
    When user opens the address "/apps/Chem/Sketcher?smiles=%5B%236%5D-%5B%236%5D-%5B%237%2C%238%5D%20%7CatomProp%3A2.molAttchpt.1%7C"
    Then the "Sketcher" view should be current
    And the "ready" reading of crux sketcher widget should be "true"
    And the "mode" reading of crux sketcher widget should be "query"
    And the address's link should be "[#6]-[#6]-[#7,#8] |atomProp:2.molAttchpt.1|"
    And no errors should have been logged

  Scenario: The gear's Query mode turns the query tools on and off, the drawing and its link as they were
    Given user opens Chem's Sketcher app
    And the "ready" reading of crux sketcher widget should be "true"
    And the "mode" reading of crux sketcher widget should be "molecule"
    When user clicks on Crux settings button
    And user clicks on Crux query mode item
    Then the "mode" reading of crux sketcher widget should be "query"
    When user clicks on Crux settings button
    And user clicks on Crux query mode item
    Then the "mode" reading of crux sketcher widget should be "molecule"
    And the address's link should be the molecule "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
    And no errors should have been logged

  Scenario Outline: Export saves <format> in the <style> style
    Given user opens Chem's Sketcher app
    And the "ready" reading of crux sketcher widget should be "true"
    When user clicks Export... in the Sketcher app's ribbon
    And user chooses the "<style>" style and the "<format>" format in the molecule Export dialog
    And user clicks Save in the molecule Export dialog
    Then the saved file should be "<file>", holding "<holds>"
    And no errors should have been logged

    Examples:
      | style             | format        | file              | holds                       |
      | ACS Document 1996 | SVG           | molecule.svg      | pt'                         |
      | As drawn          | SVG           | molecule.svg      | px'                         |
      | ACS Document 1996 | PNG (300 dpi) | molecule.png      | PNG                         |
      | As drawn          | PNG (300 dpi) | molecule.png      | PNG                         |
      | ACS Document 1996 | CDXML         | molecule.cdxml    | <CDXML                      |
      | As drawn          | CDXML         | molecule.cdxml    | <CDXML                      |
      | As drawn          | MOL V2000     | molecule.mol      | V2000                       |
      | As drawn          | MOL V3000     | molecule.mol      | V3000                       |
      | As drawn          | SDF           | molecule.sdf      | $$$$                        |
      | As drawn          | SMILES        | molecule.smi      | CC(C)c1c(C(=O)Nc2ccccc2)    |
      | As drawn          | CXSMILES      | molecule.cxsmiles | [C@@H](O)C[C@@H](O)CC(=O)O  |
      | As drawn          | SMARTS        | molecule.sma      | [#6]-[#6](-[#6])            |
      | As drawn          | InChI         | molecule.inchi    | InChI=1S/C33H35FN2O5/       |
      | As drawn          | InChIKey      | molecule.inchikey | XUKUURHRXDUEBC-KAYWLYCHSA-N |

  Scenario: Export copies the SMILES it shows; the style is off for a text format
    Given user opens Chem's Sketcher app
    And the "ready" reading of crux sketcher widget should be "true"
    When user clicks Export... in the Sketcher app's ribbon
    And user chooses the "As drawn" style and the "SMILES" format in the molecule Export dialog
    Then the molecule Export dialog's style should be off
    And the molecule Export dialog should show the text "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
    When user clicks Copy in the molecule Export dialog
    Then the clipboard should hold the molecule "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
    And no errors should have been logged

  Scenario Outline: No scrollbar at <width> by <height>
    Given the browser window is <width> by <height>
    And user opens Chem's Sketcher app
    Then the "ready" reading of crux sketcher widget should be "true"
    And the Sketcher app should fill its view with Crux, no scrollbar and no field

    Examples:
      | width | height |
      | 1280  | 720    |
      | 1440  | 900    |
      | 1920  | 1080   |
