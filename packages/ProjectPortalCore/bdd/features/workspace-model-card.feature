@portal
Feature: A catalogued model's card beside its preview in the Workspace on Home
  Power Pack's Workspace tab on Home shows, beside a model's own preview, what other packages say
  about the model. For a model in the portal's catalogue that is its card: approval, the portfolios
  it belongs to, the question it answers, its inputs with their validity range and its outputs.

  Background:
    Given user opens the Project Portal app

  Scenario: Pin the example model, select it in the Workspace, read its card
    Given the user has pinned the "ProjectPortalCore:FedBatchTiter" model
    When user selects the "ProjectPortalCore:FedBatchTiter" model in the Workspace on Home
    Then workspace model card should contain text "What final titer"
    And workspace model card should contain text "draft"
    And workspace model card should contain text "EX-BIOC"
