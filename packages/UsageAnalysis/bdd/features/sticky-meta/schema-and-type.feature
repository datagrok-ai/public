@journey @serial @sticky-meta @realizes:views.entity-type @realizes:views.entity-property-schemas
Feature: An entity type and a metadata schema, from creation to deletion
  The administrator's side of Sticky Meta under Browse > Platform > Sticky Meta: a new entity type
  (what objects the metadata attaches to), a schema associated with it holding four typed
  properties, the schema read back in its Edit dialog, and both deleted from their lists.
  Translated from TestTrack StickyMeta/01-schema-and-type.md (Tests 1.1-1.4), its primary
  create-schema-and-type case, and playwright-tests/e2e/stickymeta/01-schema-and-type.test.ts.

  Types and Schemas are switched by their tree nodes: the views share one route family, and
  /meta/types opened by address shows Schemas. Each list is paginated, so a card is looked for after
  a search for its name. The type matches molecules ("semtype=molecule"), as the case asks.

  Everything is made through the UI under names unique to the run, and the last two scenarios
  delete it through the UI as well; a journey runs them even when an earlier scenario failed. What
  the library does not have yet — reading the type and the schema back from the server, and sweeping
  them at the feature's start and end — is listed in sticky-meta/MISSING.md, as are the property
  types the Edit dialog shows (the case's "rating/int, notes/string…" is claimed for the names and
  the order only).

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Sticky Meta" tree node inside browse tree

  Scenario: A new entity type needs a name and a matching expression (1.1)
    When user clicks on "Platform > Sticky Meta > Types" tree node inside browse tree
    Then the "Entity types" view should be current
    When user clicks on "New Entity Type..." button
    Then "Create a new entity type" dialog should be visible
    And OK button in "Create a new entity type" dialog should be disabled
    When user enters "bdd-sm-type-{time}" into "Name" input in "Create a new entity type" dialog
    Then OK button in "Create a new entity type" dialog should be disabled
    When user enters "semtype=molecule" into "Matching expression" input in "Create a new entity type" dialog
    Then OK button in "Create a new entity type" dialog should be enabled
    When user clicks on OK button in "Create a new entity type" dialog
    Then the "Create a new entity type" dialog should close
    When user types "bdd-sm-type-{time}" into gallery search
    Then "bdd-sm-type-{time}" link in gallery should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A schema is associated with the type and given four typed properties (1.2)
    When user clicks on "Platform > Sticky Meta > Schemas" tree node inside browse tree
    Then the "Schemas" view should be current
    When user clicks on "New Schema..." button
    Then "Create a new schema" dialog should be visible
    And OK button in "Create a new schema" dialog should be disabled
    And "Property Type" input in "Create a new schema" dialog should offer "string, int, bool, double, datetime, string_list"
    When user enters "bdd-sm-schema-{time}" into "Name" input in "Create a new schema" dialog
    And user clicks on "select entities" action in "Create a new schema" dialog
    Then "Select types for bdd-sm-schema-{time}" dialog should be visible
    When user checks "bdd-sm-type-{time}" property in "Select types for bdd-sm-schema-{time}" dialog
    And user clicks on OK button in "Select types for bdd-sm-schema-{time}" dialog
    Then the "Select types for bdd-sm-schema-{time}" dialog should close
    And "Associated with" input in "Create a new schema" dialog should contain text "bdd-sm-type-{time}"
    When user enters "rating" into second "Name" input in "Create a new schema" dialog
    And user selects "int" in "Property Type" input in "Create a new schema" dialog
    Then "Property Type" input in "Create a new schema" dialog should have value "int"
    When user clicks on "Add new property to schema" button in "Create a new schema" dialog
    And user enters "notes" into third "Name" input in "Create a new schema" dialog
    And user selects "string" in second "Property Type" input in "Create a new schema" dialog
    Then second "Property Type" input in "Create a new schema" dialog should have value "string"
    When user clicks on "Add new property to schema" button in "Create a new schema" dialog
    And user enters "verified" into fourth "Name" input in "Create a new schema" dialog
    And user selects "bool" in third "Property Type" input in "Create a new schema" dialog
    Then third "Property Type" input in "Create a new schema" dialog should have value "bool"
    When user clicks on "Add new property to schema" button in "Create a new schema" dialog
    And user enters "review_date" into fifth "Name" input in "Create a new schema" dialog
    And user selects "datetime" in fourth "Property Type" input in "Create a new schema" dialog
    Then fourth "Property Type" input in "Create a new schema" dialog should have value "datetime"
    Then OK button in "Create a new schema" dialog should be enabled
    When user clicks on OK button in "Create a new schema" dialog
    Then the "Create a new schema" dialog should close
    When user types "bdd-sm-schema-{time}" into gallery search
    Then "bdd-sm-schema-{time}" link in gallery should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Edit dialog shows the schema as it was saved (1.3)
    When user picks "Edit" from the context menu of "bdd-sm-schema-{time}" link in gallery
    Then "Edit schema" dialog should be visible
    And first "Name" input in "Edit schema" dialog should have value "bdd-sm-schema-{time}"
    And "Associated with" input in "Edit schema" dialog should contain text "bdd-sm-type-{time}"
    And second "Name" input in "Edit schema" dialog should have value "rating"
    And third "Name" input in "Edit schema" dialog should have value "notes"
    And fourth "Name" input in "Edit schema" dialog should have value "verified"
    And fifth "Name" input in "Edit schema" dialog should have value "review_date"
    And sixth "Name" input in "Edit schema" dialog should be absent
    When user clicks on CANCEL button in "Edit schema" dialog
    Then the "Edit schema" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The schema is deleted from its list (1.4)
    When user remembers the gallery counter
    And user picks "Delete" from the context menu of "bdd-sm-schema-{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And the gallery counter should be lower than remembered
    And "bdd-sm-schema-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The entity type is deleted from its list (1.4)
    When user clicks on "Platform > Sticky Meta > Types" tree node inside browse tree
    Then the "Entity types" view should be current
    When user types "bdd-sm-type-{time}" into gallery search
    Then "bdd-sm-type-{time}" link in gallery should be visible
    When user remembers the gallery counter
    And user picks "Delete" from the context menu of "bdd-sm-type-{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And the gallery counter should be lower than remembered
    And "bdd-sm-type-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown
