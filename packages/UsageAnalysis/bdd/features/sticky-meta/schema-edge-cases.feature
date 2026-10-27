@serial @sticky-meta @edge-cases @realizes:views.entity-property-schemas @realizes:views.entity-type
Feature: Sticky Meta schema edge cases: a taken name, a malformed matching expression, a bad address, delete from the tree
  Guards for the schema and entity type editors under Browse > Platform > Sticky Meta. A new schema
  under the name of an existing one would drop that schema's properties, so OK asks first ("This
  action will delete properties …", `property_schema_meta.dart` `save`), and CANCEL there keeps the
  schema as it was and the New Schema dialog open. A matching expression not in "tag=value" form is
  refused ('Use "tag=value" format', `entity_type_meta.dart` `edit`). An address naming no schema says
  so in an error balloon ('Schema "…" not found', `routing.dart` `_metaSchemasHandler`).

  Deleting a schema from the Browse tree shows no balloon and leaves its node in the tree until the tree
  is refreshed (core#677): that scenario is a known failure, and the tag goes when core#677 is fixed.

  The schema and its entity type are made through the API, matching molecules, and removed at feature
  end with the type the malformed expression is typed for; serial with the other Sticky Meta features,
  since a schema that matches molecules shows in every molecule column's Sticky meta pane; only CANCEL
  in the confirmation is claimed.

  Background:
    Given user is logged in
    And a Sticky Meta schema "bdd-sme-schema-{time}" for molecules, with the entity type "bdd-sme-type-{time}", is on the server with the properties:
      | rating | int    |
      | notes  | string |
    And the Sticky Meta schema "bdd-sme-none-{time}" and entity type "bdd-sme-badtype-{time}" are removed now and at feature end
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Sticky Meta" tree node inside browse tree

  Scenario: A new schema under a taken name asks before dropping its properties, and CANCEL keeps them
    When user clicks on "Platform > Sticky Meta > Schemas" tree node inside browse tree
    Then the "Schemas" view should be current
    When user clicks on "New Schema..." button
    Then "Create a new schema" dialog should be visible
    When user enters "bdd-sme-schema-{time}" into "Name" input in "Create a new schema" dialog
    And user clicks on "select entities" action in "Create a new schema" dialog
    Then "Select types for bdd-sme-schema-{time}" dialog should be visible
    When user checks "bdd-sme-type-{time}" property in "Select types for bdd-sme-schema-{time}" dialog
    And user clicks on OK button in "Select types for bdd-sme-schema-{time}" dialog
    Then the "Select types for bdd-sme-schema-{time}" dialog should close
    When user enters "score" into second "Name" input in "Create a new schema" dialog
    And user selects "int" in "Property Type" input in "Create a new schema" dialog
    Then OK button in "Create a new schema" dialog should be enabled
    When user clicks on OK button in "Create a new schema" dialog
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "delete properties rating, notes"
    When user clicks on CANCEL button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And "Create a new schema" dialog should be visible
    And the Sticky Meta schema "bdd-sme-schema-{time}" should be on the server with the properties "rating: int, notes: string"
    When user clicks on CANCEL button in "Create a new schema" dialog
    Then the "Create a new schema" dialog should close
    And the Sticky Meta schema "bdd-sme-schema-{time}" should be on the server with the properties "rating: int, notes: string"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A matching expression not in tag=value form is refused
    When user clicks on "Platform > Sticky Meta > Types" tree node inside browse tree
    Then the "Entity types" view should be current
    When user clicks on "New Entity Type..." button
    Then "Create a new entity type" dialog should be visible
    When user enters "bdd-sme-badtype-{time}" into "Name" input in "Create a new entity type" dialog
    And user enters "semtype" into "Matching expression" input in "Create a new entity type" dialog
    Then "Matching expression" input in "Create a new entity type" dialog should be invalid
    And OK button in "Create a new entity type" dialog should be disabled
    When user enters "semtype=" into "Matching expression" input in "Create a new entity type" dialog
    Then "Matching expression" input in "Create a new entity type" dialog should be invalid
    And OK button in "Create a new entity type" dialog should be disabled
    When user enters "semtype=molecule" into "Matching expression" input in "Create a new entity type" dialog
    Then "Matching expression" input in "Create a new entity type" dialog should be valid
    And OK button in "Create a new entity type" dialog should be enabled
    When user clicks on CANCEL button in "Create a new entity type" dialog
    Then the "Create a new entity type" dialog should close
    And no entity type named "bdd-sme-badtype-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An address naming no schema says so and leaves the current view
    When user clicks on "Platform > Sticky Meta > Schemas" tree node inside browse tree
    Then the "Schemas" view should be current
    When user opens the address "/meta/schemas/bdd-sme-none-{time}"
    Then an error balloon containing "Schema \"bdd-sme-none-{time}\" not found" should have been shown
    And the "Schemas" view should be current

  @known-failure:95 @realizes:core#677
  Scenario: A schema deleted from the Browse tree leaves the tree at once
    core#677: deleting a schema from the Browse tree shows no balloon, and the schema stays in the tree
    until the tree is refreshed. Remove the tag when core#677 is fixed.
    When user expands "Platform > Sticky Meta > Schemas" tree node inside browse tree
    Then "Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}" tree node inside browse tree should be visible
    When user picks "Delete" from the context menu of "Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}" tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And no Sticky Meta schema named "bdd-sme-schema-{time}" should be on the server
    And "Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}" tree node inside browse tree should be absent
    And an info balloon should have been shown
    And no errors should have been logged
