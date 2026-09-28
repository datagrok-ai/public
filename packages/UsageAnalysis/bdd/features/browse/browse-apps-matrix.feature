@browse @realizes:views.browse
Feature: Every application a stand runs in the browser opens from Browse > Apps with its own content
  The TestTrack case Apps/apps.md 1 ("open each application: no errors, no visual glitches") and the
  apps matrix of playwright-public/browse/apps_matrix.test.ts, whose check came down to "some node
  exists and nothing logged an unfiltered error" — it never asserted that a view opened. Here each
  row names its node, the view the application opens, and something only that application draws:
  its own inputs and buttons, the list it fills, the form it builds. The manual note of
  browse_manual_tests2.md section 6 (a double click as well as a single one) is its own outline.

  The rows are the client-side applications of the packages a local stand carries (Chem,
  U2Demo). Others are covered by their own packages' features: MPO profiles (Chem), Diff Studio,
  Monomer Libraries and Collections (Bio), U2 Demo (U2Demo), Model Hub (browse-apps-and-dashboards),
  and the Tutorials app (Tutorials). The client-side applications of packages a local stand does not
  carry — Flow, Excalidraw, MetabolicGraph, Hit Triage and Hit Design, PeptiHit and PepTriage, the
  Oligo Toolkit family and Oligo Batch Calculator, Markush and HELM Enumerators, Plates search — are
  not rows yet: a row is written once it has been run against the application. Left out by the rule
  (a container, a server script or an outside service behind the application): Boltz-1, Docking,
  Admetica, MolTrack, Preclinical Case, Benchling, CDD Vault, Chemspace, Signals, Revvity Signals,
  KNIME, and Clinical Case (a Python reader for its studies, and a file written into App Data).

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded

  Scenario Outline: <node> opens from the tree with its own content
    Given Apps---<group> tree node inside browse tree is expanded
    And <parent> tree node inside browse tree is expanded
    When user clicks on <node> tree node inside browse tree
    Then the "<view>" view should be current
    And the page address should contain "<address>"
    And <content> should be visible
    And <more> should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | group | parent                | node                                             | view                     | address                          | content                                | more                         |
      | Chem  | Apps---Chem---Reactions | Apps---Chem---Reactions---Reaction-Enumerator      | Reaction Enumerator      | /apps/Chem/ReactionEnumerator    | "Number of steps" input                | "Next: Reactions" button     |
      | Chem  | Apps---Chem---Reactions | Apps---Chem---Reactions---Transformation-Reactions | Transformation Reactions | /apps/Chem/TransformationReactions | "Molecules" input                     | "Run Reaction" button        |
      | Chem  | Apps---Chem---Reactions | Apps---Chem---Reactions---Two-Component-Reactions  | Two-Component Reactions  | /apps/Chem/TwoComponentReactions | "Reactant 1" input                     | "Run Reaction" button        |
      | Dev   | Apps---Dev              | Apps---Dev---Reports-Browser                       | Reports                  | /apps/U2demo/ReportsBrowser      | first item in list                     | first row actions            |
      | Dev   | Apps---Dev              | Apps---Dev---U2-Designer                           | U2 Designer              | /apps/U2demo/U2Designer          | "nameInput" element                    | "Save" button                |

  Scenario Outline: A double click keeps <node> open through the next click in the tree
    Given Apps---<group> tree node inside browse tree is expanded
    And <parent> tree node inside browse tree is expanded
    When user double-clicks on <node> tree node inside browse tree
    Then the "<view>" view should be current
    When user clicks on Dashboards tree node inside browse tree
    Then Projects view should be visible
    And "<view>" view should be present
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | group | parent                  | node                                               | view                     |
      | Chem  | Apps---Chem---Reactions | Apps---Chem---Reactions---Transformation-Reactions | Transformation Reactions |
      | Dev   | Apps---Dev              | Apps---Dev---U2-Designer                           | U2 Designer              |
