@guide @tutorial @viewport:1280x720
Feature: How-to · Scripting
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A JavaScript script written from the template: Browse >
  Platform > Functions > Scripts, New > JavaScript Script, the header lines typed in (name, description, two
  inputs with captions and defaults, an output) and the code line (BMI from height and weight); Run opens a form
  built from the header, Height set to 180, OK, and the Results under the editor read 21.60 (70 / 1.8²); Save
  stores it, and the Scripts view finds it by its name. Starts with Browse open and the toolbox and the context
  panel hidden. The script "Body mass index" is deleted at feature end and swept at its start, and the Scripts
  gallery's search (an account setting) is emptied again; nothing else is made on the server. The gallery shows
  the stand's other scripts for a few seconds behind the New menu. Filmed in a 1280x720 window, so the code reads larger in the video.

  Scenario: Write a JavaScript script
    A script's header lines turn into a form for its inputs. JavaScript scripts run in your browser.

    Given user is logged in
    And simple mode is off
    And no script named "Body mass index" is on the server
    And the Scripts gallery search is emptied at feature end
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    # caption: Scripts live in Browse, under Platform > Functions
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Functions" tree node inside browse tree
    And user clicks on "Platform > Functions > Scripts" tree node inside browse tree
    Then the "Scripts" view should be current
    # caption: New > JavaScript Script: it runs in your browser
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    And code editor should contain the text "alert('Hello World!')"
    # caption: The comment lines at the top describe the script
    When user replaces the line starting with "//name:" in code editor with "//name: Body mass index"
    And user replaces the line starting with "//description:" in code editor with "//description: BMI from height and weight"
    # caption: Inputs: type, name, default value and a caption
    And user puts "//input: double height = 175 {caption: Height, cm}" on a new line after the line starting with "//language:" in code editor
    And user puts "//input: double weight = 70 {caption: Weight, kg}" on a new line after the line starting with "//input: double height" in code editor
    # caption: The output is a variable the code sets
    And user puts "//output: double bmi" on a new line after the line starting with "//input: double weight" in code editor
    And user replaces the line starting with "alert(" in code editor with "bmi = weight / (height / 100) ** 2;"
    # caption: Run: the inputs come as a form
    And user clicks on "Run script (F5)" icon
    Then "Body mass index" dialog should be visible
    When user enters "180" into "Height, cm" input in "Body mass index" dialog
    And user clicks on OK button in "Body mass index" dialog
    Then the "Body mass index" dialog should close
    # caption: The result shows under the code
    And the script result "bmi" should read "21.60"
    # caption: Save it to run it again later
    When user saves the script
    Then 1 script named "Body mass index" should be on the server
    # caption: It is listed in Scripts with the others
    When user clicks on "Scripts" view tab
    And user types "Body mass" into gallery search
    Then "Body mass index" link in gallery should be visible
    And no errors should have been logged
