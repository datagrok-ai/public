@guide @tutorial
Feature: How-to · Files
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A single click on a file in Browse > Files > Demo previews it in
  the view area without opening it: a CSV as a table (cars.csv), an image (images/cats/cat2.jpg, the platform's
  image viewer), a PDF page by page (texts/pdf/research-paper.pdf, FileEditors), an SDF as a table of structures
  (chem/sdf/ApprovedDrugs2015.sdf, Chem) and a PDB file in 3D (bio/ngl-formats/1crn.pdb, BiostructureViewer's
  Mol* viewer; the file is read from the server, nothing is fetched from outside). Each folder is opened on camera
  and closed again after its file, so the tree stays short. Opening a file with a double-click is
  01-getting-started/02's. Nothing is made on the server.

  Scenario: Preview any file from Browse
    One click on a file in Browse shows it in place: tables, images, PDFs, molecules and 3D structures.

    Given user is logged in
    And simple mode is off
    And the "FileEditors" package is installed
    And the "Chem" package is installed
    And the "BiostructureViewer" package is installed
    And the package autostarts have completed
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Browse > Files > Demo: sample files of many kinds"
    And user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    And user clicks on Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    And no table should be left in the workspace
    When the video says "One click shows a file here, without opening it"
    And user expands Files---Demo---images tree node inside browse tree
    And user expands Files---Demo---images---cats tree node inside browse tree
    And user clicks on Files---Demo---images---cats---cat2.jpg tree node inside browse tree
    Then the "cat2.jpg" view should be current
    When the video says "Images show as they are"
    And user collapses Files---Demo---images tree node inside browse tree
    And user expands Files---Demo---texts tree node inside browse tree
    And user expands Files---Demo---texts---pdf tree node inside browse tree
    And user clicks on Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree
    Then the "research-paper.pdf" view should be current
    When the video says "A PDF shows its pages"
    And user collapses Files---Demo---texts tree node inside browse tree
    And user expands Files---Demo---chem tree node inside browse tree
    And user expands Files---Demo---chem---sdf tree node inside browse tree
    And user clicks on Files---Demo---chem---sdf---ApprovedDrugs2015.sdf tree node inside browse tree
    Then the "ApprovedDrugs2015" view should be current
    When the video says "A molecule file (SDF) shows its structures"
    And user collapses Files---Demo---chem tree node inside browse tree
    And user expands Files---Demo---bio tree node inside browse tree
    And user expands Files---Demo---bio---ngl-formats tree node inside browse tree
    And user clicks on Files---Demo---bio---ngl-formats---1crn.pdb tree node inside browse tree
    Then the "1crn.pdb" view should be current
    And "Reset Camera" button should be visible
    And no table should be left in the workspace
    And no errors should have been logged
    When the video says "A protein structure (PDB) shows in 3D"
