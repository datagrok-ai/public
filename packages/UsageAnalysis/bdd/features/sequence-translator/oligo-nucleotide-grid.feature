@journey @realizes:sequencetranslator.oligo-renderer.convert-helm-to-oligo @realizes:sequencetranslator.cell.oligo-nucleotide
@realizes:sequencetranslator.panel.oligo-nucleotide @realizes:sequencetranslator.panel.oligo-structures
Feature: OligoNucleotide duplex column: conversion, context panels and cell actions
  SequenceTranslator builds an OligoNucleotide duplex column from a HELM column of sirna-demo
  (44 siRNA and ASO rows). Oligo > Convert HELM to Oligo in the cell's context menu adds
  "oligo_helm (oligo)"; its Oligo-Nucleotide pane summarises the current duplex (lengths,
  modifications, conjugates and the Duplex alignment line), its Oligo Structures pane builds the
  strand structures, and its cell menu copies the HELM or a picture, opens the HELM editor and the
  oligo enumerator; a double click shows the duplex full screen. Translated from the TestTrack case
  SequenceTranslator/oligo-nucleotide-grid (Blocks A, C, D, E, F, H, I; Block G, on cyclized.csv, is
  oligo-polytool.feature). Oligo > Combine Sense+Antisense to Oligo (Block H) comes last: it adds
  a second duplex column that the earlier scenarios' column counts do not expect.

  Rows used: 1 siR-0001 (19/19 duplex with a GalNAc-L3 conjugate), 2 siR-0002, 34 aso-0034
  (single strand), 38 siR-0038 (3' overhangs), 40 siR-0040 (explicit HELM base pairs). The grid
  rows of a HELM column are tall, so a row past the first screen is reached by the mouse wheel.
  The Duplex texts of rows 38 and 40 were read on localhost (SequenceTranslator 1.11.5).

  Not here: Block B and the "renders as a duplex / single strand" parts of Blocks A and H (the
  duplex drawing and the monomer tooltip are canvas-only — no reading); the structure pictures of the Oligo Structures pane and of the full-screen dialog (canvas); the HELM
  editor's OK path and a real enumeration (both need a gesture inside the HELM Web Editor canvas).
  The cell actions sit under Current Value in the grid's context menu. The Enumerate dialog draws
  the loaded HELM as numbered monomer labels (m1 G2 sp3 …), so row 1's sense start is read as
  "m1G2sp3m4A5sp6m7C8p9m10U11p12" (row 2 starts with a Chol conjugate).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded
    When user double-clicks Files---App-Data---SequenceTranslator---samples---sirna-demo.csv tree node inside browse tree
    Then the "sirna-demo" view should be current
    And the table should have 44 rows
    And "oligo_helm" column should have semantic type "Macromolecule"
    When user scrolls the grid to the "oligo_helm" column
    And user picks "Oligo > Convert HELM to Oligo" from the context menu of the "cell 1 of oligo_helm" area of grid
    Then the table should have a column "oligo_helm (oligo)"

  Scenario: Convert HELM to Oligo appends an OligoNucleotide column and leaves the HELM column as it was
    Then "oligo_helm (oligo)" column should have semantic type "OligoNucleotide"
    And "oligo_helm (oligo)" column should hold the same values as "oligo_helm" column
    And the value of "oligo_helm" column in row 1 should be "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"
    And "oligo_helm" column should have semantic type "Macromolecule"
    And the table should have 13 columns
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The Oligo-Nucleotide pane summarises the current duplex and follows the current cell
    Given the context panel is open
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user clicks on the "cell 1 of oligo_helm (oligo)" area of grid
    And user expands "Oligo-Nucleotide" accordion header in context panel
    Then "Sense length" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 nt"
    And "Antisense length" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 nt"
    And "Modifications used" table row in "Oligo-Nucleotide" pane in context panel should contain the text "2'-OMe ×38, PS ×8"
    And "Conjugates" table row in "Oligo-Nucleotide" pane in context panel should contain the text "GalNAc-L3 linker ×1"
    And "Duplex" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 bp, blunt (auto-aligned)"
    And "Oligo-Nucleotide" pane in context panel should contain the text "2'-O-Methyl"
    And "Oligo-Nucleotide" pane in context panel should contain the text "Phosphorothioate (linkage)"
    When user clicks on the "cell 2 of oligo_helm (oligo)" area of grid
    Then "Modifications used" table row in "Oligo-Nucleotide" pane in context panel should contain the text "2'-OMe ×20, 2'-F ×18, PS ×8"
    And "Conjugates" table row in "Oligo-Nucleotide" pane in context panel should contain the text "Cholesterol ×1, GalNAc-L3 linker ×1"
    And no errors should have been logged

  Scenario: A single-strand ASO has no Duplex line; the Duplex line reports overhangs and HELM base pairs
    Given the context panel is open
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user clicks on the "cell 1 of oligo_helm (oligo)" area of grid
    And user expands "Oligo-Nucleotide" accordion header in context panel
    Then "Duplex" table row in "Oligo-Nucleotide" pane in context panel should be visible
    When user scrolls the mouse wheel down 40 times over the "cell 1 of oligo_helm (oligo)" area of grid
    And user clicks on the "cell 38 of oligo_helm (oligo)" area of grid
    Then row 38 should be current
    And "Duplex" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 bp, overhangs: 3' antisense +2, 3' sense +2 (auto-aligned)"
    When user clicks on the "cell 40 of oligo_helm (oligo)" area of grid
    Then row 40 should be current
    And "Duplex" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 bp, blunt (from HELM pairs)"
    When user scrolls the mouse wheel up 1 times over the "cell 40 of oligo_helm (oligo)" area of grid
    And user clicks on the "cell 34 of oligo_helm (oligo)" area of grid
    Then row 34 should be current
    And "Antisense length" table row in "Oligo-Nucleotide" pane in context panel should contain the text "single-strand"
    And "Modifications used" table row in "Oligo-Nucleotide" pane in context panel should contain the text "LNA ×6, PS ×18"
    And "Duplex" table row in "Oligo-Nucleotide" pane in context panel should be absent
    When user scrolls the mouse wheel up 40 times over the "cell 34 of oligo_helm (oligo)" area of grid
    Then grid should have a "cell 1 of oligo_helm (oligo)" area
    When user clicks on the "cell 1 of oligo_helm (oligo)" area of grid
    Then row 1 should be current
    And "Antisense length" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 nt"
    And no errors should have been logged

  Scenario: The Oligo Structures pane builds the sense and antisense structures
    Given the context panel is open
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user clicks on the "cell 1 of oligo_helm (oligo)" area of grid
    Then row 1 should be current
    And "Duplex" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 bp, blunt (auto-aligned)"
    When user expands "Oligo Structures" accordion header in context panel
    Then "Sense" accordion header in "Oligo Structures" pane in context panel should be visible
    And "Antisense" accordion header in "Oligo Structures" pane in context panel should be visible
    When user expands "Sense" accordion header in "Oligo Structures" pane in context panel
    Then "Explore" accordion header in "Sense" pane in context panel should be visible
    When user expands "Antisense" accordion header in "Oligo Structures" pane in context panel
    Then "Explore" accordion header in "Antisense" pane in context panel should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Copy as HELM copies the duplex HELM of the right-clicked cell
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user right-clicks on the "cell 1 of oligo_helm (oligo)" area of grid
    Then the open menu should list "Current Value > Copy as HELM"
    And the open menu should list "Current Value > Copy as Image"
    And the open menu should list "Current Value > Edit HELM"
    And the open menu should list "Enumerate Oligos"
    When user closes the context menu
    And user picks "Current Value > Copy as HELM" from the context menu of the "cell 1 of oligo_helm (oligo)" area of grid
    Then an info balloon containing "HELM copied to clipboard" should have been shown
    And the clipboard should have the text "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"
    And no errors should have been logged

  Scenario: Copy as Image puts a PNG picture of the cell on the clipboard
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user picks "Current Value > Copy as Image" from the context menu of the "cell 1 of oligo_helm (oligo)" area of grid
    Then an info balloon containing "Image copied to clipboard" should have been shown
    And the clipboard should hold a PNG image of at least 10000 bytes
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Edit HELM opens the HELM editor on the duplex, and CANCEL leaves the cell as it was
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user picks "Current Value > Edit HELM" from the context menu of the "cell 1 of oligo_helm (oligo)" area of grid
    Then HELM notation tab should be visible
    When user clicks on HELM notation tab
    Then HELM notation should contain the text "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p"
    And HELM notation should contain the text "|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p"
    When user clicks on CANCEL button
    Then HELM notation tab should be hidden
    And the value of "oligo_helm (oligo)" column in row 1 should be "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"
    And no errors should have been logged

  Scenario: A double click shows the duplex full screen in a dialog without buttons
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user double-clicks on the "cell 1 of oligo_helm (oligo)" area of grid
    Then "Oligonucleotide" dialog should be visible
    And CANCEL button in "Oligonucleotide" dialog should be hidden
    And OK button in "Oligonucleotide" dialog should be hidden
    When user clicks on Close icon in "Oligonucleotide" dialog
    Then the "Oligonucleotide" dialog should close
    And the "sirna-demo" view should be current
    And the value of "oligo_helm (oligo)" column in row 1 should be "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Enumerate Oligos opens the HELM enumeration on the cell, and CANCEL adds nothing
    When user scrolls the grid to the "oligo_helm (oligo)" column
    And user picks "Enumerate Oligos" from the context menu of the "cell 1 of oligo_helm (oligo)" area of grid
    Then "PolyTool Helm Enumeration" dialog should be visible
    And "PolyTool Helm Enumeration" dialog should contain the text "m1G2sp3m4A5sp6m7C8p9m10U11p12"
    And "PolyTool Helm Enumeration" dialog should contain the text "L3"
    When user clicks on CANCEL button in "PolyTool Helm Enumeration" dialog
    Then the "PolyTool Helm Enumeration" dialog should close
    And the table should have 13 columns
    And the table should have 44 rows
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Combine Sense+Antisense to Oligo pairs the sense and antisense columns into a duplex column
    When user scrolls the grid to the "sense_helm" column
    And user picks "Oligo > Combine Sense+Antisense to Oligo..." from the context menu of the "cell 1 of sense_helm" area of grid
    Then "Combine Sense + Antisense to Oligonucleotide" dialog should be visible
    And Sense input in "Combine Sense + Antisense to Oligonucleotide" dialog should contain the text "sense_helm"
    When user selects "antisense_helm" in Antisense input in "Combine Sense + Antisense to Oligonucleotide" dialog
    And user clicks on OK button in "Combine Sense + Antisense to Oligonucleotide" dialog
    Then the "Combine Sense + Antisense to Oligonucleotide" dialog should close
    And the table should have a column "sense_helm+antisense_helm (oligo)"
    And "sense_helm+antisense_helm (oligo)" column should have semantic type "OligoNucleotide"
    And the table should have 14 columns
    And the value of "sense_helm+antisense_helm (oligo)" column in row 1 should be "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"
    And the value of "sense_helm+antisense_helm (oligo)" column in row 34 should be "RNA1{[lna](C)[sp].[lna](A)[sp].[lna](G)[sp].d(T)[sp].d(G)[sp].d(T)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(T)[sp].d(G)[sp].d(C)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(A)[sp].[lna](T)[sp].[lna](A)[sp].[lna](A)}$$$$"
    Given the context panel is open
    When user scrolls the grid to the "sense_helm+antisense_helm (oligo)" column
    And user clicks on the "cell 1 of sense_helm+antisense_helm (oligo)" area of grid
    Then row 1 should be current
    And "Sense length" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 nt"
    And "Antisense length" table row in "Oligo-Nucleotide" pane in context panel should contain the text "19 nt"
    When user scrolls the mouse wheel down 40 times over the "cell 1 of sense_helm+antisense_helm (oligo)" area of grid
    And user scrolls the mouse wheel up 1 times over the "cell 40 of sense_helm+antisense_helm (oligo)" area of grid
    And user clicks on the "cell 34 of sense_helm+antisense_helm (oligo)" area of grid
    Then row 34 should be current
    And "Antisense length" table row in "Oligo-Nucleotide" pane in context panel should contain the text "single-strand"
    And "Duplex" table row in "Oligo-Nucleotide" pane in context panel should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged
