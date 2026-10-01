@tutorials @serial @realizes:tutorials.differential-equations
Feature: The Differential Equations tutorial
  Walks Scientific computing > Differential equations from its card to the end: Diff Studio run from
  Browse > Apps, the Lotka-Volterra model opened from its library, the guided tours of the view and of
  the equations, the predator equation and a new parameter typed into the editor, the model refreshed,
  and its inputs changed. Each step is claimed as ticked and as done — the app and the model open,
  the editor shown and hidden, the equations as typed, the inputs as set.
  Translated from playwright-tests/e2e/tutorials/differential-equations.test.ts, which typed into the
  editor by the index of its line. The model is solved in the browser. Opening the model adds it to
  Diff Studio's recent models, a file in the account's home that is put back.

  Fixed in the tutorial for this translation: the counter was 15 for 14 steps; the Lotka-Volterra model
  was taken as the third card of the library, whatever it was; the Prey, Delta and Finish inputs were
  taken by their position in the form, and the eta parameter the learner adds shifts it, so "Set Delta
  to 0.1" watched gamma and could never be done.
  The tutorial prepares each step before it lists it (the Apps gallery, the model view, a tour): a
  learner reads the step first, and so does this walk.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "DiffStudio" package is installed
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the file "diff-studio-recent.d42" of the user's home is put back at feature end
    And the "Differential equations" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Differential equations tutorial
    When user starts the "Differential equations" tutorial
    Then the tutorial progress should be 1 of 14
    When user clicks on Apps tree node inside browse tree
    Then the tutorial step "Open Apps" should be done
    # the tutorial prepares each step before it lists it: act once it is on the list
    Given the tutorial step "Run Diff Studio" should not be done yet
    When user double-clicks on Diff-Studio gallery card
    Then the tutorial step "Run Diff Studio" should be done
    And the "Diff Studio" view should be current
    Given the tutorial step "Run the Lotka-Volterra model" should not be done yet
    When user double-clicks on "Lotka-Volterra" model card
    Then the tutorial step "Run the Lotka-Volterra model" should be done
    And the "Lotka-Volterra" view should be current

    Given the tutorial step "Explore the interface" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore the interface" should be done

    Given the tutorial step "Open equations editor" should not be done yet
    When user clicks on Edit ribbon item
    Then the tutorial step "Open equations editor" should be done
    And code editor should be visible
    Given the tutorial step "Explore editor" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore editor" should be done

    Given the tutorial step "Complete the predator equation" should not be done yet
    When user replaces the line starting with "dy/dt" in code editor with "dy/dt = -gamma * y + delta * x * y - eta * y * y"
    Then the tutorial step "Complete the predator equation" should be done
    When user puts "eta = 0.01 {min: 0; max: 0.1; category: Parameters} [Crowding effect]" on a new line after the line starting with "#parameters" in code editor
    Then the tutorial step "Add the eta parameter" should be done
    When user clicks on Refresh ribbon item
    Then the tutorial step "Apply changes" should be done
    Given the tutorial step "Check the updates" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Check the updates" should be done

    When user clicks on Edit ribbon item
    Then the tutorial step "Close equations editor" should be done
    And code editor should be absent
    When user enters "2" into "Prey" input
    Then the tutorial step "Set \"Prey\" to 2" should be done
    And "Prey" input should have the value "2"
    When user enters "0.1" into "Delta" input
    Then the tutorial step "Set \"Delta\" to 0.1" should be done
    And "Delta" input should have the value "0.1"
    When user enters "150" into "Finish" input
    Then the tutorial step "Set \"Finish\" to 150" should be done
    And "Finish" input should have the value "150"

    And the "Differential equations" tutorial should be completed
    And the tutorial should have listed 14 steps
    And the tutorial progress should be 14 of 14
    And no hint should be shown
    And no errors should have been logged
