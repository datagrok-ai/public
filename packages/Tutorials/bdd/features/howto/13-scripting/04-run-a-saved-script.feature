@guide @tutorial @viewport:1440x810
Feature: How-to · Scripting
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A saved script run again without opening it: Browse > Platform >
  Functions > Scripts, searched for "Body mass", the card clicked; the context panel's Details say what it does and
  what it takes and returns, its Script pane shows the code (a plain text area), its Run pane is the same form as the
  run dialog (Height 160, Weight 55), RUN, and the
  output lands in the Variables panel (bmi 21.48 = 55 / 1.6²). Arranged off camera: the script "Body mass index"
  saved through the API (deleted at feature end, swept at the start), the toolbox hidden, the context panel and
  Browse open. The Scripts gallery's search (an account setting) is emptied again at feature end. Filmed in a
  1440x810 window, so the panel's text reads larger. The Activity pane is left out: it lists the creation, and the
  run only once the panel is drawn again.

  Scenario: Run a saved script without opening it
    Click a saved script to see what it takes, and run it from the context panel.

    Given user is logged in
    And simple mode is off
    And a script "Body mass index" is on the server:
      """
      //language: javascript
      //description: BMI from height and weight
      //input: double height = 175 {caption: Height, cm}
      //input: double weight = 70 {caption: Weight, kg}
      //output: double bmi
      bmi = weight / (height / 100) ** 2;
      """
    And the Scripts gallery search is emptied at feature end
    And the toolbox pane is hidden
    And the context panel is open
    And the browse panel is open
    # caption: Saved scripts are in Browse, under Platform > Functions > Scripts
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Functions" tree node inside browse tree
    And user clicks on "Platform > Functions > Scripts" tree node inside browse tree
    Then the "Scripts" view should be current
    When user types "Body mass" into gallery search
    # caption: Click a script: what it does, its inputs and outputs
    And user clicks on "Body mass index" link in gallery
    Then the context panel should show "Body mass index"
    # caption: Script: its code, without opening an editor
    When user clicks on "Script" pane header in context panel
    Then the Script pane of the context panel should show "bmi = weight / (height / 100) ** 2;"
    # caption: Run: the same form, right in the context panel
    When user clicks on "Run" pane header in context panel
    And user enters "160" into "Height, cm" input in context panel
    And user enters "55" into "Weight, kg" input in context panel
    And user clicks on context panel run button
    # 55 kg / 1.6² m
    Then the Variables panel should show "bmi" as "21.48"
    # caption: The result goes to the Variables panel
    When user hovers over the text "bmi" in "Variables" dock panel
    Then no errors should have been logged
