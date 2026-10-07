@guide
Feature: Several portfolios on one instance
  A guide for the team: one Datagrok runs small-molecule discovery, biologics discovery and
  biologics development side by side. Each portfolio has its own registry feed and its own
  stage-gate model, and says whether the registry or the portal owns the stage. Development mirrors
  the registry's gates and refuses to move them; biologics discovery reads the same registry but not
  its gates; small-molecule discovery decides its gates in the portal, one stage at a time, once the
  next stage's criteria are met. The portfolios run on synthetic example registries.

  Scenario: Three portfolios, three stage models, two kinds of gate authority
    Given user is logged in
    And simple mode is off
    And a scratch project "KINASE-X" in the "EX-SMD" portfolio
    And user opens the Project Portal app
    And the side panels are hidden
    Then portfolio list should contain text "Example Biologics Development"
    When user clicks on biologics development portfolio
    Then board header should contain text "stages from registry"
    When user clicks on EXB-201 card
    Then gate reasons should contain text "The registry owns this portfolio's stages"
    When user clicks on portfolios link
    And user clicks on biologics discovery portfolio
    Then board header should contain text "EX-BIOD-V1"
    When user clicks on EXB-101 card
    Then program detail should contain text "Registry stage: 3"
    When user clicks on portfolios link
    And user clicks on small molecule discovery portfolio
    Then board header should contain text "stages decided here"
    When user clicks on KINASE-X card
    And user clicks on advance button
    Then KINASE-X card in target validation column should be visible
    And gate reasons should contain text "required criteria of Hit Finding"
    When user clicks on assay criterion met button
    And user clicks on advance button
    Then transition history should contain text "Hit Finding"
