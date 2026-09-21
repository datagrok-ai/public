@browse @realizes:views.browse
Feature: The My stuff section of the Browse tree
  The section that gathers what belongs to the current user. Translated from the manual case
  Browse-MyStuff-01 (playwright-public/browse/mystuff.test.ts, browse_manual_tests2.md section 4).

  Browse-MyStuff-02 (Recent holds what was opened recently), -03 and -05 (Add to Favorites from
  My Files, GROK-19848) are not translated and no other file holds them: the Favorites and Recent
  nodes do not announce that they have reloaded, so a claim made right after the change reads the
  old list. GROK-19848 therefore has no regression test here.
  -06 (a new script appears under My stuff) is writable — the old spec made its own fixture and
  deleted it — and is simply not written yet. -04 depends on what others have shared with this
  account.

  Recent, Favorites and Shared with me are the section's own nodes. The rest are the buckets of
  the user's personal project (project_meta.dart bucketOf: Connections, Files, Dashboards,
  Scripts, Tables, Flows, Spaces), and a bucket is there only while the account owns something
  of that kind, so only the three fixed nodes are claimed. A bucket carries the same name as the
  section elsewhere in the tree; the full path tells them apart.

  A node below the top level is named by its full tree path ("Files---Demo"), which is what the
  platform writes into its own `name` attribute. Several sections carry a node called Demo, Files
  or App Data, and a bare name matches whichever of them another feature happened to leave
  open: the tree remembers its expanded set per user, across features and across runs.

  Background:
    Given user is logged in
    And the browse panel is open
    And My stuff tree node inside browse tree is expanded

  Scenario: My stuff gathers the user's own things
    Then the following elements should be visible:
      | My-stuff---Recent tree node inside browse tree          |
      | My-stuff---Favorites tree node inside browse tree       |
      | My-stuff---Shared-with-me tree node inside browse tree  |
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The section closes on its own twistie
    Given My-stuff---Recent tree node inside browse tree should be visible
    When user collapses My stuff tree node inside browse tree
    Then My-stuff---Recent tree node inside browse tree should be hidden
    And My-stuff---Favorites tree node inside browse tree should be hidden
    And My stuff tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
