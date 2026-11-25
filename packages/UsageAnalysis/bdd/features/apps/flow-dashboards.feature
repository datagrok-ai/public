@realizes:views.projects @realizes:views.space
Feature: The dashboards a flow publishes
  A flow's Save dialog publishes its result tables as a dashboard whose tables the flow computes again when it opens.
  public#4675: publishing "As a new dashboard" took the tables of the first one, which then failed to open.
  public#4676: the dashboard of a flow saved into a space called the flow by the name it had before the move, and opened
  without data until the page was reloaded; saved into the same space under the flow's name, it was taken for the flow.

  The flows are made through the API over cars.csv; they, the space and the dashboards carry the run's time and are
  removed when a scenario starts and when the feature ends (a script that is in a space when the space goes can no
  longer be deleted, so the spaces are made, and go, last).

  Background:
    Given user is logged in
    And simple mode is off

  @realizes:public#4675
  Scenario: A flow published as a new dashboard leaves the first one with its table
    Given no project named "BDDFlowFirst{time}, BDDFlowSecond{time}" is on the server
    And a script "BDDFlowPublish{time}" is on the server:
      """
      //language: flow
      //tags: flow
      //output: dataframe result
      {"version":"2.0","name":"BDDFlowPublish{time}","description":"","author":"","created":"","modified":"","nodes":[{"id":"a","typeName":"DG Functions/Uncategorized/OpenFile","label":"Open File: cars.csv","description":"","pos":{"x":40,"y":40},"properties":{},"inputValues":{"fullPath":"System:DemoFiles/cars.csv","sheetName":""}},{"id":"b","typeName":"Outputs/Table Output","label":"Table Output","description":"","pos":{"x":400,"y":40},"properties":{"paramName":"result"},"inputValues":{}}],"connections":[{"id":"c","source":"a","sourceOutput":"result","target":"b","targetInput":"table"}],"annotations":[],"groups":[],"metadata":{"settings":{"scriptName":"BDDFlowPublish{time}","scriptDescription":"","tags":["flow"]}}}
      """
    And user opens the "Flow" app
    When user clicks on "BDDFlowPublish{time}" recent flow
    Then the "nodes" reading of flow editor widget should be 2
    When user clicks on flow run button
    Then the "computed tables" reading of flow editor widget should be 1
    When user clicks on flow save button
    And user clicks on "Save & publish…" button in "Save Flow" dialog
    Then "Save project" dialog should be visible
    When user enters "BDDFlowFirst{time}" into Name text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    And user clicks on flow save button
    And user selects "As a new dashboard" in Publish choice input in "Save Flow" dialog
    And user clicks on "Save & publish…" button in "Save Flow" dialog
    Then "Save project" dialog should be visible
    When user enters "BDDFlowSecond{time}" into Name text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "cars" table of the "BDDFlowSecond{time}" project should be saved with data sync
    And the "cars" table of the "BDDFlowFirst{time}" project should be saved with data sync
    When user closes all views
    And user opens the "BDDFlowFirst{time}" project and waits for its table
    Then the table should have 30 rows

  @realizes:public#4676
  Scenario: A dashboard published while its flow is saved into a space opens with its data
    Given no project named "BDDFlowMoved{time}" is on the server
    And a script "BDDFlowMoved{time}" is on the server:
      """
      //language: flow
      //tags: flow
      //output: dataframe result
      {"version":"2.0","name":"BDDFlowMoved{time}","description":"","author":"","created":"","modified":"","nodes":[{"id":"a","typeName":"DG Functions/Uncategorized/OpenFile","label":"Open File: cars.csv","description":"","pos":{"x":40,"y":40},"properties":{},"inputValues":{"fullPath":"System:DemoFiles/cars.csv","sheetName":""}},{"id":"b","typeName":"Outputs/Table Output","label":"Table Output","description":"","pos":{"x":400,"y":40},"properties":{"paramName":"result"},"inputValues":{}}],"connections":[{"id":"c","source":"a","sourceOutput":"result","target":"b","targetInput":"table"}],"annotations":[],"groups":[],"metadata":{"settings":{"scriptName":"BDDFlowMoved{time}","scriptDescription":"","tags":["flow"]}}}
      """
    And a space named "BDDFlowHome{time}" is on the server
    And user opens the "Flow" app
    When user clicks on "BDDFlowMoved{time}" recent flow
    Then the "nodes" reading of flow editor widget should be 2
    When user clicks on flow run button
    Then the "computed tables" reading of flow editor widget should be 1
    When user clicks on flow save button
    And user clicks on the "Choose…" link in "Save Flow" dialog
    And user clicks on "BDDFlowHome{time}" tree node in "Save Flow" dialog
    And user clicks on "Save & publish…" button in "Save Flow" dialog
    Then "Save project" dialog should be visible
    And Name text input in "Save project" dialog should have value "BDDFlowMoved{time}"
    When user clicks on OK in the Save project dialog and the project uploads
    And user closes all views
    And user opens the "BDDFlowMoved{time}" project and waits for its table
    Then "Data loading error" dialog should be absent
    And the table should have 30 rows

  @realizes:public#4676
  Scenario: So does one published into the same space under the flow's name
    Given no project named "BDDFlowInSpace{time}" is on the server
    And a script "BDDFlowInSpace{time}" is on the server:
      """
      //language: flow
      //tags: flow
      //output: dataframe result
      {"version":"2.0","name":"BDDFlowInSpace{time}","description":"","author":"","created":"","modified":"","nodes":[{"id":"a","typeName":"DG Functions/Uncategorized/OpenFile","label":"Open File: cars.csv","description":"","pos":{"x":40,"y":40},"properties":{},"inputValues":{"fullPath":"System:DemoFiles/cars.csv","sheetName":""}},{"id":"b","typeName":"Outputs/Table Output","label":"Table Output","description":"","pos":{"x":400,"y":40},"properties":{"paramName":"result"},"inputValues":{}}],"connections":[{"id":"c","source":"a","sourceOutput":"result","target":"b","targetInput":"table"}],"annotations":[],"groups":[],"metadata":{"settings":{"scriptName":"BDDFlowInSpace{time}","scriptDescription":"","tags":["flow"]}}}
      """
    And a space named "BDDFlowSpace{time}" is on the server
    And user opens the "Flow" app
    When user clicks on "BDDFlowInSpace{time}" recent flow
    Then the "nodes" reading of flow editor widget should be 2
    When user clicks on flow run button
    Then the "computed tables" reading of flow editor widget should be 1
    When user clicks on flow save button
    And user clicks on the "Choose…" link in "Save Flow" dialog
    And user clicks on "BDDFlowSpace{time}" tree node in "Save Flow" dialog
    And user clicks on "Save & publish…" button in "Save Flow" dialog
    Then "Save project" dialog should be visible
    And Name text input in "Save project" dialog should have value "BDDFlowInSpace{time}"
    When user clicks on space input in "Save project" dialog
    And user clicks on "BDDFlowSpace{time}" tree node
    Then space input in "Save project" dialog should contain text "BDDFlowSpace{time}"
    When user clicks on OK in the Save project dialog and the project uploads
    And user closes all views
    And user opens the "BDDFlowInSpace{time}" project and waits for its table
    Then "Data loading error" dialog should be absent
    And the table should have 30 rows
