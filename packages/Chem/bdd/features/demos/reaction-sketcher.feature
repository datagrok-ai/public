@guide @demo
Feature: Chemistry · The Reaction Sketcher
  Demo clips of the Reaction Sketcher (Apps > Chem > Reactions > Reaction Sketcher; features/demos, filmed with
  `grok-bdd guide`: the @demo tag leaves the set-up and the checks off camera). The app edits a reaction or a route in
  Crux Sketch on the left and draws it for publication on the right. On the canvas an arrow drawn after the last product
  continues the route, an arrow's "+" handle adds a reactant, a reagent or a product to its step, a condition typed over
  an arrow can be dragged and sized, a catalytic cycle's slots fill with the molecules clicked, and the Map tool numbers
  an atom through every step; the right pane maps each step by itself and says how sure the mapper is. Export previews
  the route for a journal and saves it; Add to table puts it into a table. Each clip opens the app on a reaction of its
  own (kept in the browser as the app keeps its document; what was kept comes back at feature end). What the canvas does
  is crux-sketch's own BDD; what the app keeps is Chem's package tests ('Reactions: app', 'Reactions: export',
  'Reactions: add panel').

  Scenario: The Reaction Sketcher, on its example route
    Apps, Chem, Reactions: the Reaction Sketcher opens on a five-step route, drawn to edit on the left and for print on the right.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the Reaction Sketcher opens on its example
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Chem tree node inside browse tree is expanded
    When user expands Apps---Chem---Reactions tree node inside browse tree
    And user clicks on Apps---Chem---Reactions---Reaction-Sketcher tree node inside browse tree
    Then the Reaction Sketcher should show its example route
    When user hovers over reactions app drawing
    Given the browse panel is closed
    When the video says "A five-step route: edit it on the left, read it on the right"
    And user scrolls the Reaction Sketcher's preview to its steps
    Then no errors should have been logged

  Scenario: Building a route on the canvas
    A new arrow continues the route from the last product; the arrow's "+" handle adds its product and a reagent pasted beside.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the Reaction Sketcher is open on "Brc1ccc(C=O)cc1.OB(O)c1ccccc1>>O=Cc1ccc(-c2ccccc2)cc1", titled "4-Phenylbenzyl alcohol in two steps"
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "A new arrow continues from the last product: step 2"
    And user clicks on crux reaction arrow tool
    And user draws an arrow 6 bond lengths long under the drawing
    And the video says "The arrow's + handle adds to its step"
    And user hovers over the middle of arrow 2
    And user clicks the "+" handle of the arrow
    And user chooses "Add Product" in the sketcher's menu
    And user pastes "OCc1ccc(-c2ccccc2)cc1"
    And user presses Escape
    And the video says "A reagent pasted, then added to step 2"
    And user pastes "CC(C)C[AlH]CC(C)C" at 75% across and 25% down the canvas
    And user hovers over the middle of arrow 2
    And user clicks the "+" handle of the arrow
    And user chooses "Add Reagent" in the sketcher's menu
    And user clicks on the last atom of the sketcher
    And user presses Escape
    And user presses Escape
    Then the Reaction Sketcher's Text should be the route "Brc1ccc(C=O)cc1.OB(O)c1ccccc1>>O=Cc1ccc(-c2ccccc2)cc1--**--O=Cc1ccc(-c2ccccc2)cc1>CC(C)C[AlH]CC(C)C>OCc1ccc(-c2ccccc2)cc1"
    And no errors should have been logged

  Scenario: Conditions over an arrow
    A condition typed over the arrow; dragged, it shows where it will land; its menu sizes it.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the Reaction Sketcher is open on "Brc1ccc(C=O)cc1.OB(O)c1ccccc1>>O=Cc1ccc(-c2ccccc2)cc1", titled "Suzuki coupling"
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "Conditions typed over the arrow"
    And user clicks the "reaction.text" tool of the sketcher
    And user writes "Pd(PPh3)4, K2CO3" over arrow 1
    Then the sketcher should draw the text "Pd(PPh3)4, K2CO3"
    When the video says "Dragged, it shows where it will land"
    And user clicks the "select" tool of the sketcher
    And user drags the text "Pd(PPh3)4, K2CO3" under arrow 1
    And the video says "Its menu sizes it"
    And user right-clicks the text "Pd(PPh3)4, K2CO3"
    And user chooses "Size" in the sketcher's menu
    And user chooses "14" in the sketcher's menu
    Then the sketcher should draw the text "Pd(PPh3)4, K2CO3" in 14 pt
    And no errors should have been logged

  Scenario: A catalytic cycle
    The Cycle tool places an empty cycle; the molecules clicked fill its slots in turn, joined by its arcs.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the Reaction Sketcher is open on "[Pd].Br[Pd]c1ccccc1.c1ccc([Pd]c2ccccc2)cc1", titled "A catalytic cycle"
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "A cycle placed empty"
    And user clicks the "reaction.cycles" tool of the sketcher
    And user clicks the canvas 3 bond lengths under the drawing
    And the video says "Each molecule clicked fills the next slot"
    And user clicks on the "atom 0" area of crux sketcher widget
    And user clicks on the "atom 1" area of crux sketcher widget
    And user clicks on the "atom 12" area of crux sketcher widget
    Then the Reaction Sketcher's Text should be the route "[Pd]>>Br[Pd]c1ccccc1--**--Br[Pd]c1ccccc1>>c1ccc([Pd]c2ccccc2)cc1--**--c1ccc([Pd]c2ccccc2)cc1>>[Pd]"
    And no errors should have been logged

  Scenario: Map numbers through a route
    The right pane maps each step by itself and says how sure; the Map tool numbers one atom through three steps.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the Reaction Sketcher is open on "CC(=O)O.CCN>>CCNC(C)=O--**--CCNC(C)=O>>CCNCC--**--CCNCC>>CCN(CC)C(C)=O", titled "Acylation, reduction, acylation"
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "Each step mapped by itself: how sure, and what to check"
    And user scrolls the Reaction Sketcher's preview to its mapping
    Then the Reaction Sketcher's preview should say how sure the mapper is of step 1
    And the Reaction Sketcher's preview should warn of step 3 "Check this mapping"
    When the video says "The Map tool: one nitrogen through three steps"
    And user clicks the "reaction.map" tool of the sketcher
    And user drags the "atom 6" area of crux sketcher widget to the "atom 9" area
    And user drags the "atom 9" area of crux sketcher widget to the "atom 15" area
    And user drags the "atom 15" area of crux sketcher widget to the "atom 20" area
    Then the Reaction Sketcher's Text should be the route "CC(=O)O.CC[NH2:1]>>CC[NH:1]C(C)=O--**--CC[NH:1]C(C)=O>>CC[NH:1]CC--**--CC[NH:1]CC>>CC[N:1](CC)C(C)=O"
    When the video says "An atom's menu clears its number"
    And user right-clicks atom 20 of the sketcher
    And user chooses "Clear Map Number" in the sketcher's menu
    And the video says "Auto-map numbers every step"
    And user clicks the "reaction.map.open" tool of the sketcher
    And user chooses "Auto-map" in the sketcher's menu
    Then every step of the Reaction Sketcher's Text should carry map numbers
    And no errors should have been logged

  Scenario: Export for a journal
    Export previews the route for a journal: ACS 1996, a line a step, map numbers, or as drawn; Save writes the file.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the Reaction Sketcher is open on its example
    And the toolbox pane is hidden
    And the context panel is closed
    And user watches downloads
    When the video says "Export: the route in a journal's style, previewed"
    And user clicks Export... in the Reactions app's ribbon
    And user selects "ACS Document 1996" in "Style" input in "Export reaction" dialog
    And user checks "Fit each step on one line" input in "Export reaction" dialog
    And user checks "Map numbers" input in "Export reaction" dialog
    And the video says "As drawn: the canvas as you drew it"
    And user checks "As drawn" input in "Export reaction" dialog
    And the video says "The preview follows the dialog"
    And user resizes "Export reaction" dialog by 320 and 0 pixels from its corner
    And user selects "SVG" in "Format" input in "Export reaction" dialog
    And user clicks on SAVE button in "Export reaction" dialog
    Then a file "reaction.svg" should have been downloaded
    And no errors should have been logged

  Scenario: The route into a table
    Add to table puts the route into a table of reactions, drawn in its cell as one chain.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the Reaction Sketcher is open on its example
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "Add to table: the route into a table"
    And user clicks on "Add to table" button
    Then the "Reactions" table view should open with 1 rows
    Given the "reaction" column of the grid is 1300 px wide
    And the rows of the grid are 330 px tall
    When user hovers over the "cell 1 of reaction" area of grid
    Then no errors should have been logged
