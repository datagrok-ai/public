@journey @realizes:public#4677
Feature: A changed flow asks before it is closed or replaced
  public#4677: closing a changed flow, or opening another one over it, lost the changes without a word. The close asks in
  the dialog every view with unsaved work closes through, whatever "Ask to save on close" says, and its "save the flow"
  link saves first; New… and Open from platform… ask in Flow's own dialog. A flow opened and left as it was asks nothing.
  What counts as a change, and the other ways to replace a flow (a template, a recent flow, a table's history), are
  Flow's package tests ("Flow: unsaved changes").

  The flow is made through the API under the run's time and removed when the feature starts and ends.

  Background:
    Given user is logged in
    And simple mode is off
    And a script "BDDFlowUnsaved{time}" is on the server:
      """
      //language: flow
      //tags: flow
      {"version":"2.0","name":"BDDFlowUnsaved{time}","description":"","author":"","created":"","modified":"","nodes":[{"id":"a","typeName":"Constants/String","label":"String","description":"","pos":{"x":100,"y":100},"properties":{"value":"x"},"inputValues":{}}],"connections":[],"annotations":[],"groups":[],"metadata":{"settings":{"scriptName":"BDDFlowUnsaved{time}","scriptDescription":"","tags":["flow"]}}}
      """
    And user opens the "Flow" app

  Scenario: Closing a changed flow asks, closing an unchanged one does not
    When user clicks on "BDDFlowUnsaved{time}" recent flow
    Then the "nodes" reading of flow editor widget should be 1
    And the "unsaved changes" reading of flow editor widget should be "false"
    When user closes "BDDFlowUnsaved{time}" view tab
    Then "BDDFlowUnsaved{time}" view tab should be absent
    Given user opens the "Flow" app
    When user clicks on "BDDFlowUnsaved{time}" recent flow
    And user clicks on "Inputs" pane header in toolbox
    And user double-clicks on flow sketcher input item
    Then the "nodes" reading of flow editor widget should be 2
    When user closes "BDDFlowUnsaved{time}" view tab
    Then "Close BDDFlowUnsaved{time}?" dialog should contain text "Edited the flow."
    When user clicks on CANCEL button in "Close BDDFlowUnsaved{time}?" dialog
    Then the "Close BDDFlowUnsaved{time}?" dialog should close
    And the "nodes" reading of flow editor widget should be 2

  Scenario: A new flow and another flow ask before they replace it
    When user picks "Flow | New…" from the top menu
    Then "Unsaved changes" dialog should contain text "BDDFlowUnsaved{time}"
    When user clicks on DISCARD button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And the "nodes" reading of flow editor widget should be 0
    When user clicks on "BDDFlowUnsaved{time}" recent flow
    Then the "nodes" reading of flow editor widget should be 1
    When user double-clicks on flow sketcher input item
    Then the "nodes" reading of flow editor widget should be 2
    When user picks "Flow | Open from platform…" from the top menu
    And user selects "BDDFlowUnsaved{time}" in Flow choice input in "Open Flow" dialog
    And user clicks on OK button in "Open Flow" dialog
    Then "Unsaved changes" dialog should contain text "BDDFlowUnsaved{time}"
    When user clicks on DISCARD button in "Unsaved changes" dialog
    Then the "nodes" reading of flow editor widget should be 1
    And the "unsaved changes" reading of flow editor widget should be "false"

  Scenario: The close dialog saves the flow first
    When user double-clicks on flow sketcher input item
    Then the "nodes" reading of flow editor widget should be 2
    When user closes "BDDFlowUnsaved{time}" view tab
    And user clicks on the "save the flow" link in "Close BDDFlowUnsaved{time}?" dialog
    Then "BDDFlowUnsaved{time}" view tab should be absent
    And the script "BDDFlowUnsaved{time}" on the server should contain "Inputs/Sketcher Input"
