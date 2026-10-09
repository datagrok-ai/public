@portal
Feature: A catalogued model's card in the context panel
  Selecting a model anywhere, here in Model Hub, shows its card in the context panel when a
  catalogue of the user's portfolios lists it. The card is labelled with those portfolios.

  Scenario: Select the example model in Model Hub, read its card in the context panel
    When user selects "Fed-batch titer (example)" in Model Hub
    Then model card pane portfolios should contain text "In the catalogue of EX-BIOC"
    And model card pane should contain text "What final titer"
    And model card pane should contain text "draft"
