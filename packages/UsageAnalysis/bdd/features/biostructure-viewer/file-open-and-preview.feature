@viewers @realizes:biostructureviewer.preview.biostructure @realizes:biostructureviewer.import.pdb @realizes:biostructureviewer.import.pdbqt
Feature: Opening and previewing structure files from the Files browser
  A click on a structure file in Browse > Files previews it with the Mol* engine, in a view named
  after the file; a double-click opens that view on its own. A .pdb file's menu in the folder view
  has Open table residues. Translated from the TestTrack case
  BiostructureViewer/biostructureviewer-file-open-and-preview (GROK-17654, GROK-18999, GROK-14442,
  GROK-16968, CLAUDE-33, GROK-21118).

  Facts of the stand the md did not have: the view a double-click opens is named after the file
  ("1bdq.pdb"), not "Mol*" — the open shows the file's previewer as a view rather than calling the
  package's Import PDB handler; a preview becomes current under the file's name, which is how the
  preview title of GROK-18999 is claimed; Open table residues is offered in the folder view's
  context menu, not in the Browse tree's. The residue table's rows and columns are tested in
  BiostructureViewer src/tests/pdb-helper-tests.ts 'pdbToDf'.

  The Mol* engine is claimed by its Reset Camera overlay button. A preview also gets a view tab
  under the file's name, so a double-clicked view is told from the preview by previewing another
  file afterwards: a preview's tab goes with the next preview (shown first in the same scenario),
  an opened view's tab stays.

  Help in the sidebar opens the datagrok.ai help site in a new browser tab; the help panel the md
  means by "open Help and close it" is toggled with F1, which is what the CLAUDE-33 scenario does.


  The NGL-only formats (Scenario 2 and Scenario 6) preview and open in an NGL host, claimed as the
  "NGL host" element (`data-u2-name="ngl-host"` on the `.d4-ngl-viewer` that holds the NGL canvas)
  with no Mol* engine beside it.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded
    And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded

  Scenario: PDB, mmCIF, CIF and PDBQT files preview with the Mol* engine under their own names (GROK-17654, GROK-18999)
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree
    Then the "1bdq.pdb" view should be current
    And "Reset Camera" button should be visible
    When user clicks on Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree
    Then the "1RQ9.mmcif" view should be current
    And "Reset Camera" button should be visible
    When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree
    Then the "1rq9-assembly1.cif" view should be current
    And "Reset Camera" button should be visible
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.autodock-gpu.pdbqt tree node inside browse tree
    Then the "1bdq.autodock-gpu.pdbqt" view should be current
    And "Reset Camera" button should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: NGL-only formats preview with the NGL engine (GROK-13650)
    Given Files---Demo tree node inside browse tree is expanded
    And Files---Demo---bio tree node inside browse tree is expanded
    And Files---Demo---bio---ngl-formats tree node inside browse tree is expanded
    When user clicks on Files---Demo---bio---ngl-formats---1blu.mmtf tree node inside browse tree
    Then the "1blu.mmtf" view should be current
    And "NGL host" element should be visible
    And "Reset Camera" button should be absent
    When user clicks on Files---Demo---bio---ngl-formats---1crn.ply tree node inside browse tree
    Then the "1crn.ply" view should be current
    And "NGL host" element should be visible
    When user clicks on Files---Demo---bio---ngl-formats---1crn.obj tree node inside browse tree
    Then the "1crn.obj" view should be current
    And "NGL host" element should be visible
    When user clicks on Files---Demo---bio---ngl-formats---1lee.ccp4 tree node inside browse tree
    Then the "1lee.ccp4" view should be current
    And "NGL host" element should be visible
    When user clicks on Files---Demo---bio---ngl-formats---3pqr.cns tree node inside browse tree
    Then the "3pqr.cns" view should be current
    And "NGL host" element should be visible
    And "Reset Camera" button should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A double-click on an NGL-only format opens it in a view of its own with the NGL engine
    Given simple mode is off
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---bio tree node inside browse tree is expanded
    And Files---Demo---bio---ngl-formats tree node inside browse tree is expanded
    When user double-clicks on Files---Demo---bio---ngl-formats---1blu.mmtf tree node inside browse tree
    Then "1blu.mmtf" view should be visible
    And the "1blu.mmtf" view should be current
    And "NGL host" element should be visible
    And "Reset Camera" button should be absent
    When user clicks on Files---Demo---bio---ngl-formats---1crn.obj tree node inside browse tree
    Then the "1crn.obj" view should be current
    And "1blu.mmtf" view should be visible
    When user double-clicks on Files---Demo---bio---ngl-formats---1lee.ccp4 tree node inside browse tree
    Then "1lee.ccp4" view should be visible
    And the "1lee.ccp4" view should be current
    And "NGL host" element should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A double-click on a PDB opens it in a view of its own with the Mol* engine (GROK-14442, GROK-16968)
    Given simple mode is off
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree
    Then the "1bdq.pdb" view should be current
    And "Reset Camera" button should be visible
    When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree
    Then the "1rq9-assembly1.cif" view should be current
    And "1bdq.pdb" view should be absent
    When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree
    Then "1bdq.pdb" view should be visible
    And the "1bdq.pdb" view should be current
    And "Reset Camera" button should be visible
    And grid should be absent
    And "Open file" dialog should be absent
    When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree
    Then the "1rq9-assembly1.cif" view should be current
    And "1bdq.pdb" view should be visible
    When user clicks on "1bdq.pdb" view
    Then the "1bdq.pdb" view should be current
    And "Reset Camera" button should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A double-click on an mmCIF opens it in a view of its own
    Given simple mode is off
    When user double-clicks on Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree
    Then "1RQ9.mmcif" view should be visible
    And the "1RQ9.mmcif" view should be current
    And "Reset Camera" button should be visible
    When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree
    Then the "1rq9-assembly1.cif" view should be current
    And "1RQ9.mmcif" view should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A double-click on a PDBQT opens it in a view of its own
    Given simple mode is off
    When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq.autodock-gpu.pdbqt tree node inside browse tree
    Then "1bdq.autodock-gpu.pdbqt" view should be visible
    And the "1bdq.autodock-gpu.pdbqt" view should be current
    And "Reset Camera" button should be visible
    When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree
    Then the "1rq9-assembly1.cif" view should be current
    And "1bdq.autodock-gpu.pdbqt" view should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Closing an unrelated view while a structure preview is shown raises nothing (CLAUDE-33)
    Given simple mode is off
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree
    Then the "1bdq.pdb" view should be current
    And "Reset Camera" button should be visible
    Given Files---Demo tree node inside browse tree is expanded
    When user double-clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user closes demog view
    Then "demog" view should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    Given the browse panel is open
    When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree
    Then the "1bdq.pdb" view should be current
    And "Reset Camera" button should be visible
    And help panel should be absent
    When user presses F1
    Then help panel should be visible
    When user presses F1
    Then help panel should be absent
    And the "1bdq.pdb" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Open table residues from a PDB file's menu opens a residue table with an NGL viewer
    When user clicks on Files---App-Data---BiostructureViewer---samples tree node inside browse tree
    Then 1bdq.pdb link in gallery should be visible
    When user picks "Open table residues" from the context menu of 1bdq.pdb link in gallery
    Then the "Table" view should be current
    And NGL viewer should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @known-failure
  Scenario: Open table residues gives the residue names in compId (GROK-21118)
    When user clicks on Files---App-Data---BiostructureViewer---samples tree node inside browse tree
    Then 1bdq.pdb link in gallery should be visible
    When user picks "Open table residues" from the context menu of 1bdq.pdb link in gallery
    Then the "Table" table view should open with 198 rows
    And the value of "compId" column in row 1 should be "PRO"

  @known-failure
  Scenario: A Biostructure viewer added to a structure table takes the structure column by itself (GROK-21119)
    When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree
    Then the "pdb_data" table view should open with 6 rows
    When user clicks on "Add viewer" icon in toolbar
    And user types "Biostructure" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Biostructure" button in "Add Viewer" dialog
    Then Biostructure viewer should be visible
    And "Biostructure Id" property of Biostructure viewer should be "pdb"
