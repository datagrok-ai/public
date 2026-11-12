@demo @realizes:u2.image-input @realizes:u2.filter-builder
Feature: Catalog looks
  Claims on how Catalog samples look, kept apart from catalog.feature, which gen-catalog.mjs writes
  from the manifest: the image preview beside its caption, and the filter builder's read-only and
  invalid states.

  Background:
    Given user opens the "Catalog" demo page

  Scenario: The image preview sits beside its caption
    Then editor of "image-input-default" image input should be right of label of "image-input-default" image input
    And editor of "image-input-default" image input should be top-aligned with label of "image-input-default" image input

  Scenario: A read-only filter builder offers nothing that edits it
    Then Advanced link in "filter-builder-readOnly" filter builder should be hidden
    And "Add condition" button in "filter-builder-readOnly" filter builder should be hidden
    And Advanced link in "filter-builder-default" filter builder should be visible
    And "Add condition" button in "filter-builder-default" filter builder should be visible

  Scenario: A blank required filter builder shows its error
    Then error of "filter-builder-invalid" filter builder should be visible
    And error of "filter-builder-invalid" filter builder should have text "Value can't be empty"
