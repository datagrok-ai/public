@apps @demos @realizes:views.browse
Feature: The Data Access, Compute and Curves demos open from Browse > Apps > Demo with their content
  Each demo is opened by its node in the Browse tree; the Tutorials demo app says when it has run
  and its view is named (`demo-loaded`, with the demo's path), and the row claims what the demo
  builds. A row starts with the capability gate on the package the demo comes from. Translated from
  the TestTrack case Apps/apps.md 2 and playwright-public/browse/demo_apps.test.ts; the case's
  interactivity is claimed on Table Linking, where a current row in one table filters the other.

  Domain Databases turns on the account's "Domain databases" beta setting; the feature remembers the
  setting and puts it back. The Databases demo's connection list is the stand's own content, so only
  the view is claimed there. Plates' Assay Plates demo is not in the tree (its category is not one of
  the demo app's), so it has no row.

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Demo tree node inside browse tree is expanded

  Scenario Outline: The <demo> demo opens with its <content>
    Given the "<package>" package is installed
    And Apps---Demo---<group> tree node inside browse tree is expanded
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---<group>---<node> tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "<section> | <demo>"
    And the "<demo>" view should be current
    And <content> should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | package    | group       | section     | demo                  | node                  | content                        |
      | Tutorials  | Data-Access | Data Access | Table Linking         | Table-Linking         | second grid viewer             |
      | DiffStudio | Compute     | Compute     | Diff Studio           | Diff-Studio           | line chart viewer              |
      | DiffStudio | Compute     | Compute     | PK-PD Modeling        | PK-PD-Modeling        | line chart viewer              |
      | DiffStudio | Compute     | Compute     | Bioreactor            | Bioreactor            | line chart viewer              |
      | Eda        | Compute     | Compute     | Multivariate Analysis | Multivariate-Analysis | bar chart viewer               |
      | Curves     | Curves      | Curves      | Curve Fitting         | Curve-Fitting         | grid                           |
      | Curves     | Curves      | Curves      | Assay Curves          | Assay-Curves          | MultiCurveViewer viewer        |

  Scenario Outline: The <demo> demo opens the platform's <type> view
    Given the "Tutorials" package is installed
    And Apps---Demo---Data-Access tree node inside browse tree is expanded
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Data-Access---<demo> tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Data Access | <demo>"
    And the "<demo>" view should be current
    And the current view should be a <type> view
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | demo      | type      |
      | Files     | files     |
      | Databases | databases |

  Scenario: Table Linking filters the demographics by the category made current
    Given the "Tutorials" package is installed
    And Apps---Demo---Data-Access tree node inside browse tree is expanded
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Data-Access---Table-Linking tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Data Access | Table Linking"
    And table "Categories" should have 8 rows
    And table "Demographics" should have 5850 rows
    When user makes row 1 of table "Categories" current
    Then 37 rows of table "Demographics" should pass the filter
    When user makes row 7 of table "Categories" current
    Then 2444 rows of table "Demographics" should pass the filter
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Domain Databases demo turns the beta setting on, and the feature puts it back
    Given the "Tutorials" package is installed
    And the "enableDomainDatabases" shell setting is put back at feature end
    And Apps---Demo---Data-Access tree node inside browse tree is expanded
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Data-Access---Domain-Databases tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Data Access | Domain Databases"
    And the "Domain Databases" view should be current
    When user clicks on "START" button
    Then the "enableDomainDatabases" shell setting should be true
    And no errors should have been logged
