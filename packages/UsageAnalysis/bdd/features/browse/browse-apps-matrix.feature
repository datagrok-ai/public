@browse @realizes:views.browse
Feature: Every application a stand runs in the browser opens from Browse > Apps with its own content
  The TestTrack case Apps/apps.md 1 ("open each application: no errors, no visual glitches") and the
  apps matrix of playwright-public/browse/apps_matrix.test.ts, whose check came down to "some node
  exists and nothing logged an unfiltered error" — it never asserted that a view opened. Here each
  row names its node, the view the application opens, and two things only that application draws:
  its own inputs and buttons, the list it fills, the form it builds. The manual note of
  browse_manual_tests2.md section 6 (a double click as well as a single one) is its own outline.

  Each row starts with the capability gate on its package: a stand without the package skips the
  row with the reason instead of failing it. Covered by their own packages' features: MPO profiles
  (Chem), Diff Studio, Monomer Libraries and Collections (Bio), U2 Demo (U2Demo), Model Hub
  (browse-apps-and-dashboards), the Tutorials app (Tutorials).

  Not rows: HELM Enumerator and Oligo Batch Calculator open a preview card whose only content is
  RUN — the application itself is behind it; Markush Enumerator's view shows no control of its own
  to claim; Excalidraw logs a React error on opening (a candidate finding, walked before anything is
  claimed or filed); Plates > Create writes plates. Left out by the rule (a container, a server
  script or an outside service behind the application): Boltz-1, Docking, Admetica, MolTrack,
  Preclinical Case, Benchling, CDD Vault, Chemspace, Signals, Revvity Signals, KNIME, and Clinical Case
  (a Python reader for its studies, and a file written into App Data).

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded

  Scenario Outline: <view> opens from the tree with its own content
    Given the "<package>" package is installed
    And <group> tree node inside browse tree is expanded
    And <parent> tree node inside browse tree is expanded
    When user clicks on <node> tree node inside browse tree
    Then the "<view>" view should be current
    And <content> should be visible
    And <more> should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | package            | group         | parent                     | node                                               | view                     | content                            | more                            |
      | Chem               | Apps---Chem   | Apps---Chem---Reactions    | Apps---Chem---Reactions---Reaction-Enumerator      | Reaction Enumerator      | "Number of steps" input            | "Next: Reactions" button        |
      | Chem               | Apps---Chem   | Apps---Chem---Reactions    | Apps---Chem---Reactions---Transformation-Reactions | Transformation Reactions | "Molecules" input                  | "Run Reaction" button           |
      | Chem               | Apps---Chem   | Apps---Chem---Reactions    | Apps---Chem---Reactions---Two-Component-Reactions  | Two-Component Reactions  | "Reactant 1" input                 | "Run Reaction" button           |
      | HitTriage          | Apps---Chem   | Apps---Chem                | Apps---Chem---Hit-Triage                           | Hit Triage               | "Number Of Molecules" input        | "Start" button                  |
      | HitTriage          | Apps---Chem   | Apps---Chem                | Apps---Chem---Hit-Design                           | Hit Design               | "Target molecules" input           | "Start" button                  |
      | HitTriage          | Apps---Peptides | Apps---Peptides          | Apps---Peptides---PeptiHit                         | PeptiHit                 | "Chemist" input                    | "Start" button                  |
      | HitTriage          | Apps---Peptides | Apps---Peptides          | Apps---Peptides---PepTriage                        | PepTriage                | "Peptide Count" input              | "Start" button                  |
      | SequenceTranslator | Apps---Peptides | Apps---Peptides---Oligo-Toolkit | Apps---Peptides---Oligo-Toolkit---Oligo-Translator | Oligo Translator    | "Input format" input               | "Convert" button                |
      | SequenceTranslator | Apps---Peptides | Apps---Peptides---Oligo-Toolkit | Apps---Peptides---Oligo-Toolkit---Oligo-Pattern    | Oligo Pattern       | "Sense strand length" input        | "Edit strands" button           |
      | SequenceTranslator | Apps---Peptides | Apps---Peptides---Oligo-Toolkit | Apps---Peptides---Oligo-Toolkit---Oligo-Structure  | Oligo Structure     | "AS direction" input               | "Save SDF" button               |
      | Metabolicgraph     | Apps---Misc   | Apps---Misc                | Apps---Misc---MetabolicGraph                       | Metabolic Graph          | "PGK" text                         | "ICDHyr" text                   |
      | Flow               | Apps          | Apps                       | Apps---Flow                                        | Flow                     | "Create your first flow" button    | "Workflow demo" text            |
      | Plates             | Apps---Plates | Apps---Plates              | Apps---Plates---Search-plates                      | Search Plates            | "Imaging device" input             | "Plate cell count" input        |
      | Plates             | Apps---Plates | Apps---Plates              | Apps---Plates---Search-analyses                    | Search Analyses          | "IC50" input                       | "Hill Slope" input              |
      | U2demo             | Apps---Dev    | Apps---Dev                 | Apps---Dev---Reports-Browser                       | Reports                  | first item in list                 | first row actions               |
      | U2demo             | Apps---Dev    | Apps---Dev                 | Apps---Dev---U2-Designer                           | U2 Designer              | "nameInput" element                | "Save" button                   |

  Scenario Outline: A double click keeps <view> open through the next click in the tree
    Given the "<package>" package is installed
    And <group> tree node inside browse tree is expanded
    And <parent> tree node inside browse tree is expanded
    When user double-clicks on <node> tree node inside browse tree
    Then the "<view>" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    Then Projects view should be visible
    And "<view>" view should be present
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | package   | group       | parent                  | node                                               | view                     |
      | Chem      | Apps---Chem | Apps---Chem---Reactions | Apps---Chem---Reactions---Transformation-Reactions | Transformation Reactions |
      | HitTriage | Apps---Chem | Apps---Chem             | Apps---Chem---Hit-Triage                           | Hit Triage               |
      | U2demo    | Apps---Dev  | Apps---Dev              | Apps---Dev---U2-Designer                           | U2 Designer              |
