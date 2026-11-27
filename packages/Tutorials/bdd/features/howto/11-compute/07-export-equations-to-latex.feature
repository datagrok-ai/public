@guide @tutorial
Feature: How-to · Compute
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A model's equations exported as LaTeX: the Lotka-Volterra model of
  Diff Studio's library (Browse > Apps > Diff Studio), the ribbon's download icon (Save or export model), As
  LaTeX..., and the Export as LaTeX dialog: its live preview of the LaTeX source, Standalone turned off for a
  fragment to paste into a document, the copy icon (claimed on the clipboard), then CANCEL. Nothing is downloaded or
  made on the server; opening the model adds it to Diff Studio's recent models, a file in the account's home that is
  put back at feature end.

  Scenario: Export a model's equations to LaTeX
    Diff Studio writes a model's equations as LaTeX or Markdown, ready to paste into a paper or a report.

    Given user is logged in
    And simple mode is off
    And the "DiffStudio" package is installed
    And the package autostarts have completed
    And the file "diff-studio-recent.d42" of the user's home is put back at feature end
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Diff Studio is in Browse > Apps
    When user clicks on browse tab
    And user clicks on Apps tree node inside browse tree
    And user double-clicks on Diff-Studio gallery card
    Then the "Diff Studio" view should be current
    When user double-clicks on "Lotka-Volterra" model card
    Then "Prey" input should be visible
    # caption: This icon saves or exports the model
    When user clicks on "Save or export model" icon
    And user picks "As LaTeX..." from the open menu
    Then "Export as LaTeX" dialog should be visible
    When the video says "The preview is the LaTeX source, ready to use"
    # caption: Standalone off: only the part to paste into your document
    And user unchecks "Standalone" input in "Export as LaTeX" dialog
    # caption: Copy it, or download it as a .tex file
    And user clicks on "Copy the generated source to the clipboard." icon in "Export as LaTeX" dialog
    Then the clipboard should contain the text "alpha"
    When user clicks on CANCEL button in "Export as LaTeX" dialog
    Then the "Export as LaTeX" dialog should close
    And no errors should have been logged
