@browse @realizes:views.browse
Feature: The quick filters of a Browse gallery
  The filter panel a gallery opens from its toolbar, with the quick filters on top. Translated from
  the manual cases Browse-Filter-02 and -03 (browse_manual_tests2.md section 11,
  playwright-public/browse/filter.test.ts), whose old specs clicked a tag and claimed nothing.

  A quick filter writes its query into the gallery's search ("created < 1w") and the gallery reloads
  on it; All clears the query and brings the counter back to what it was. Whether the filter
  narrows the list depends on the stand — a local stand whose demo files were copied this week keeps
  all of them, and the query is recursive, so it lists the files of the subfolders too — so the
  claim is the query and the way back, not a smaller count. Each scenario ends on All, the state the next one
  expects. Browse-Filter-01 (no nameless filter property, GROK-19691) and -04 (Users by group) are in
  users-view.feature as far as the Users filter panel goes; -06 (a search that finds nothing) there
  and in spaces-search.feature.

  Browse-Filter-05 (a filter and a search together, GROK-19690) is claimed on the Scripts gallery,
  whose quick filters include Favorites, over scripts the feature makes: two whose names hold the
  search ("BddGallery<time>Starred", starred, and "...Other", not) and, for the order where the
  filter comes first, a second starred one the search does not reach ("BddFavoriteKeptAside"). The
  search the box applies to the gallery is the user's text; a quick filter picked over it is added
  to it, and one picked first is written into the box until the user types over it, after which it
  still holds. All empties the box, and the same search typed again applies again (it did not until
  2026-10-13: the box still took the cleared text for the one it had last reported, so typing it
  back was no change). The search text is a fuzzy match, so these names share nothing with any other script
  but the search; the gallery then counts exactly what both let through. The Users gallery of the
  old spec lists every login on the stand, which no count can be stated over. Both scenarios run
  in the serial lane with the other features of the Scripts view.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: Created recently puts its query in the search and All takes it back
    Given Files tree node inside browse tree is expanded
    When user clicks on Files---Demo tree node inside browse tree
    Then the "Demo" view should be current
    When user clicks on "Toggle filters" icon
    And user remembers the gallery counter
    And user clicks on "Created recently" tag
    Then gallery search should have value "created < 1w"
    When user clicks on "All" tag
    Then gallery search should have value ""
    And the gallery counter should not be lower than remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # the Dockers view warns in a balloon on a stand whose Grok Spawner is off
  Scenario: The Dockers gallery offers Used by me and Created recently
    Given the stand runs the "Grok Spawner" service
    And Platform tree node inside browse tree is expanded
    When user clicks on Platform---Dockers tree node inside browse tree
    Then the "Dockers" view should be current
    When user clicks on "Toggle filters" icon
    Then "Used by me" tag should be visible
    And "Created recently" tag should be visible
    When user remembers the gallery counter
    And user clicks on "Used by me" tag
    Then gallery search should have value "usedBy = @current"
    When user clicks on "All" tag
    Then gallery search should have value ""
    And the gallery counter should not be lower than remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @serial @GROK-19690
  Scenario: A quick filter picked over a search keeps only what both let through
    Given a script "BddGallery{time}Starred" is on the server:
      """
      //language: javascript
      //output: int x
      x = 1;
      """
    And a script "BddGallery{time}Other" is on the server:
      """
      //language: javascript
      //output: int x
      x = 2;
      """
    And the script "BddGallery{time}Starred" is in favorites
    And user opens the Scripts view
    When user enters "BddGallery{time}" into gallery search
    Then gallery counter should have text "2"
    When user remembers the gallery counter
    And user clicks on "Toggle filters" icon
    And user clicks on "Favorites" tag
    Then gallery counter should have text "1"
    And "BddGallery{time}Starred" link in gallery should be visible
    And "BddGallery{time}Other" link in gallery should be absent
    And gallery search should have value "BddGallery{time}"
    When user clicks on "All" tag
    Then gallery search should have value ""
    And the gallery counter should be higher than remembered
    When user enters "BddGallery{time}" into gallery search
    Then gallery counter should have text "2"
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @serial @GROK-19690
  Scenario: A search typed over a quick filter narrows what the filter lets through
    Given a script "BddGallery{time}Starred" is on the server:
      """
      //language: javascript
      //output: int x
      x = 1;
      """
    And a script "BddGallery{time}Other" is on the server:
      """
      //language: javascript
      //output: int x
      x = 2;
      """
    And a script "BddFavoriteKeptAside" is on the server:
      """
      //language: javascript
      //output: int x
      x = 3;
      """
    And the script "BddGallery{time}Starred" is in favorites
    And the script "BddFavoriteKeptAside" is in favorites
    And user opens the Scripts view
    When user clicks on "Toggle filters" icon
    And user remembers the gallery counter
    And user clicks on "Favorites" tag
    Then gallery search should have value "starredBy = @current"
    And the gallery counter should be lower than remembered
    And "BddGallery{time}Starred" link in gallery should be visible
    And "BddFavoriteKeptAside" link in gallery should be visible
    And "BddGallery{time}Other" link in gallery should be absent
    When user enters "BddGallery{time}" into gallery search
    Then gallery counter should have text "1"
    And "BddGallery{time}Starred" link in gallery should be visible
    And "BddFavoriteKeptAside" link in gallery should be absent
    And "BddGallery{time}Other" link in gallery should be absent
    When user clicks on "All" tag
    Then gallery search should have value ""
    And the gallery counter should not be lower than remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown
