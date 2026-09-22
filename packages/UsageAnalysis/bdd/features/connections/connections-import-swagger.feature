@connections @full-stand @serial
Feature: Importing an OpenAPI (Swagger) file as a connection
  Opening a Swagger YAML adds a Web connection named by the file's title, with one query per
  operation, under Browse > Platform > Functions > OpenAPI. Translated from TestTrack
  Connections/import-swagger.md and import-swagger-ui.md (playwright-public
  connections/06-import-swagger.test.ts).

  The file is the package's own copy of the OpenWeatherMap spec with a title of its own
  (fixtures/bdd-swagger.yaml, "BDD-Conn-Swagger"): dev keeps other OpenWeatherMap connections and a
  shared title would name several nodes. The title is fixed, so the feature is @serial and deletes
  the connection (and the queries that go with it) before and after — checked gone. Running a query
  needs the OpenWeatherMap key in DG_OPENWEATHERMAP_API_KEY (@needs-credentials), typed from the
  environment and never printed; api.openweathermap.org is an external service (@full-stand).

  Not translated, and why: dropping the file from the desktop onto the page — a synthetic drop never
  reaches the platform's drop zone (only an OS drag does), so the file is opened with Ctrl+O, which
  runs the same import. The old spec's fallback "a grid or any balloon" and "the first child node,
  whatever it is" are not restored.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And no connection named "BDD-Conn-Swagger" is on the server
    When user opens the local file "fixtures/bdd-swagger.yaml"
    Then 1 connection named "BDD-Conn-Swagger" should be on the server
    And the "bdd-swagger.yaml" view should be current
    When user closes all views
    Given the browse panel is open
    # the tree keeps the node of a connection deleted through the API; a refresh drops it
    When user clicks on "Refresh" icon inside browse toolbar
    Given Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    And Platform---Functions---OpenAPI tree node inside browse tree is expanded

  Scenario: The imported file is a connection with a query per operation
    Then Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree should be visible
    Given Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree is expanded
    Then the following elements should be visible:
      | Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree |
      | Platform---Functions---OpenAPI---BDD-Conn-Swagger---Cities-In-Cycle tree node inside browse tree                   |
      | Platform---Functions---OpenAPI---BDD-Conn-Swagger---5-day/3-hour-Forecast-By-City-Name tree node inside browse tree |
    When user opens the context menu of Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree
    Then the open menu should list "Edit..."
    And the open menu should list "Test connection"
    When user closes the context menu
    Then no errors should have been logged

  @needs-credentials
  Scenario: With its API key the connection's query returns the weather
    When user right-clicks on Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree
    And user picks "Edit..." from the open menu
    Then "Edit Connection" dialog should be visible
    And Url input in "Edit Connection" dialog should have value "https://api.openweathermap.org/data/2.5"
    When user enters the DG_OPENWEATHERMAP_API_KEY secret into ApiKey input in "Edit Connection" dialog
    And user clicks on OK button in "Edit Connection" dialog
    Then the "Edit Connection" dialog should close
    Given Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree is expanded
    When user right-clicks on Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree
    And user picks "Run" from the open menu
    Then "Current Weather Data By City Name" dialog should be visible
    When user enters "London" into Q input in "Current Weather Data By City Name" dialog
    Given user watches the task bar
    When user clicks on OK button in "Current Weather Data By City Name" dialog
    Then the task bar should have finished "Running Current Weather Data By City Name"
    And the table should have 1 row
    And the table should have a column "weather/main"
    And no error or warning balloon should have been shown
