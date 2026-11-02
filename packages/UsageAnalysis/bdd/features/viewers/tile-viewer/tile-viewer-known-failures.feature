@viewers @realizes:viewers.tile-viewer
Feature: Tile viewer — the known failures
  Temporary: the known failures of the tile viewer while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into tile-viewer-lanes.feature as an ordinary @edge-cases scenario.

  On a table opened for it with lanes on a column holding an empty value, the lane of that value
  should be offered as "(empty)", so that a card dropped into it gets an empty value (core#1217).
  Today that lane has no name (the lane names read "x, y, "), so there is no "(empty)" lane to drop a
  card into; the scenario ends there, and the drop and the empty value it writes are claimed once the
  lane exists.

  Background:
    Given user is logged in

  @edge-cases @known-failure:26 @realizes:core#1217
  Scenario: The lane of an empty value is offered as (empty) (core#1217)
    Given user opens a table "with-empty" with:
      | site | arm |
      | A    | x   |
      | A    |     |
      | B    | y   |
    And user adds a tile viewer with:
      | Lanes Column Name | arm |
    Then the "lane of row 1" reading of tile viewer should be "x"
    And tile viewer should have a "lane (empty)" area
