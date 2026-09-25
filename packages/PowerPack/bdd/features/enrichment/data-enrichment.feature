@journey @realizes:powerpack.db-explorer
Feature: Column enrichment of a database table
  DB Explorer lets a column of a query result that came from a database table be enriched: an
  enrichment joins columns of a related table onto it, is saved under a name, and applies with a
  click. Translated from TestTrack PowerPack/data-enrichment.md, against the platform's own Postgres
  database (System:Datagrok): events.session_id joined to users_sessions.id, events.event_type_id to
  event_types.id.

  The events come from two queries run from the Browse tree: a saved SQL query over the events that
  have a session, and the visual query Get Top 100 of the events table. The md writes the SQL query in
  the platform's query editor; here it is saved through the API, and running it is the claim. Every
  enrichment the feature makes starts with "bdd-enrich-" and is deleted when it ends, with the
  queries and the project.

  An enrichment is a saved join configuration (help: access/databases, Data enrichment). Applying it
  joins its columns onto the table in place, a left join; a name the table already has comes in as
  "result.<name>". Applied again after an edit, it adds the columns it now selects. Deleting an
  enrichment deletes its configuration: the pane and the server no longer have it, and the columns an
  earlier run added stay in the table. A layout keeps the look of the view, not its data: columns
  removed from the table do not come back when a layout saved before is applied. A reopened project
  runs the enrichments its table went through again, from their saved configuration; one deleted in
  the meantime is reported in an error balloon ("Could not find enrichment") and its columns do not
  come back.

  The pane of a column lists the enrichments made for that table and column, and those made for the
  column its foreign key refers to, labelled "via <table>": func_calls.session_id refers to
  users_sessions.id, so it is offered an enrichment of users_sessions.id "via users_sessions", and
  not the enrichments made for events.session_id (md 3.6).

  Background:
    Given user is logged in
    And no enrichment whose name starts with "bdd-enrich-" is on the server
    And a query "bdd-enrich-events-{time}" on "System:Datagrok" reads "select * from events where session_id is not null order by event_time desc limit 50"
    And a query "bdd-enrich-calls-{time}" on "System:Datagrok" reads "select * from func_calls where session_id is not null limit 20"
    And the context panel is open

  Scenario: Both queries over the events run from the Browse tree
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user scrolls Databases---Postgres---Datagrok---Schemas---public---events tree node inside browse tree to the middle of its list
    And user picks "Get Top 100" from the context menu of Databases---Postgres---Datagrok---Schemas---public---events tree node inside browse tree
    Then the "events" view should be current
    And the table should have 100 rows
    And the table should have a column "session_id"
    Given the browse panel is open
    When user double-clicks on Databases---Postgres---Datagrok---bdd-enrich-events-{time} tree node inside browse tree
    Then the "bdd-enrich-events-{time}" view should be current
    And the table should have 50 rows
    And the table should have the columns "id, friendly_name, session_id, event_type_id, event_time, description, error_message, error_stack_trace, exported_by"
    And "session_id" column should have no missing values
    And the context panel should show "bdd-enrich-events-{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An enrichment of users_sessions.id is made on the sessions table
    Given the browse panel is open
    When user scrolls Databases---Postgres---Datagrok---Schemas---public---users-sessions tree node inside browse tree to the middle of its list
    And user picks "Get Top 100" from the context menu of Databases---Postgres---Datagrok---Schemas---public---users-sessions tree node inside browse tree
    Then the "users_sessions" view should be current
    When user clicks on the "header id" area of grid
    Then the context panel should show "id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    When user clicks on add enrichment button
    Then "Enrich id" dialog should be visible
    When user clicks on "Add a table to join" icon in "Enrich id" dialog
    And user picks "datagrok > public > users" from the open menu
    And user picks "user_id" as the key of the main table in "Enrich id" dialog
    Then the join key in "Enrich id" dialog should read "user_id = id"
    When user opens the columns of the joined "users" table in "Enrich id" dialog
    And user toggles the "login" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user enters "bdd-enrich-users-{time}" into Name input in "Enrich id" dialog
    And user clicks on OK button in "Enrich id" dialog
    Then "bdd-enrich-users-{time}" enrichment should be visible
    And 1 enrichment named "bdd-enrich-users-{time}" should be on the server
    And no errors should have been logged
    When user closes users_sessions view
    Then the "bdd-enrich-events-{time}" view should be current

  Scenario: An enrichment of session_id is made in the editor and saved
    When user clicks on the "header session_id" area of grid
    Then the context panel should show "session_id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    Then "bdd-enrich-users-{time}" enrichment should contain text "via users_sessions"
    When user clicks on add enrichment button
    Then "Enrich session_id" dialog should be visible
    And "Enrich session_id" dialog should contain text "datagrok.public.events(1/9)"
    When user clicks on "Add a table to join" icon in "Enrich session_id" dialog
    And user picks "datagrok > public > users_sessions" from the open menu
    And user picks "session_id" as the key of the main table in "Enrich session_id" dialog
    Then the join key in "Enrich session_id" dialog should read "session_id = id"
    When user opens the columns of the joined "users_sessions" table in "Enrich session_id" dialog
    And user toggles the "ip" column in the column list of "Select columns..." dialog
    And user toggles the "started" column in the column list of "Select columns..." dialog
    And user toggles the "ended" column in the column list of "Select columns..." dialog
    And user toggles the "is_admin" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the joined "users_sessions" table in "Enrich session_id" dialog should read "datagrok.public.users_sessions(4/12)"
    When user enters "bdd-enrich-sessions-{time}" into Name input in "Enrich session_id" dialog
    And user clicks on OK button in "Enrich session_id" dialog
    Then the "Enrich session_id" dialog should close
    And "bdd-enrich-sessions-{time}" enrichment should be visible
    And 1 enrichment named "bdd-enrich-sessions-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The enrichment applied adds the four session columns, filled by the join
    Given user watches the task bar
    When user clicks on link of "bdd-enrich-sessions-{time}" enrichment
    Then the task bar should have finished "Enriching..."
    And the table should have the columns "id, friendly_name, session_id, event_type_id, event_time, description, error_message, error_stack_trace, exported_by, ended, ip, is_admin, started"
    And the table should have 50 rows
    And "ip" column should have no missing values
    And "started" column should have no missing values
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The enrichment edited to other columns adds the new column when applied again
    When user clicks on edit icon of "bdd-enrich-sessions-{time}" enrichment
    Then "Enrich session_id" dialog should be visible
    And Name input in "Enrich session_id" dialog should have value "bdd-enrich-sessions-{time}"
    And the joined "users_sessions" table in "Enrich session_id" dialog should read "datagrok.public.users_sessions(4/12)"
    When user opens the columns of the joined "users_sessions" table in "Enrich session_id" dialog
    And user toggles the "is_admin" column in the column list of "Select columns..." dialog
    And user toggles the "type" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the joined "users_sessions" table in "Enrich session_id" dialog should read "datagrok.public.users_sessions(4/12)"
    When user clicks on OK button in "Enrich session_id" dialog
    Then the "Enrich session_id" dialog should close
    And 1 enrichment named "bdd-enrich-sessions-{time}" should be on the server
    And the enrichment "bdd-enrich-sessions-{time}" on the server should select the column "type"
    Given user watches the task bar
    When user clicks on link of "bdd-enrich-sessions-{time}" enrichment
    Then the task bar should have finished "Enriching..."
    And the table should have a column "type"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A deleted enrichment is gone from the pane and the server, and its columns stay
    When user clicks on delete icon of "bdd-enrich-sessions-{time}" enrichment
    Then "bdd-enrich-sessions-{time}" enrichment should be absent
    And 0 enrichments named "bdd-enrich-sessions-{time}" should be on the server
    And the table should have a column "ip"
    And the table should have a column "type"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A second enrichment of session_id and one of event_type_id apply together
    When user closes bdd-enrich-events-{time} view
    Then the "events" view should be current
    Given the browse panel is open
    When user double-clicks on Databases---Postgres---Datagrok---bdd-enrich-events-{time} tree node inside browse tree
    Then the table should have 9 columns
    And the context panel should show "bdd-enrich-events-{time}"
    When user clicks on the "header session_id" area of grid
    Then the context panel should show "session_id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    When user clicks on add enrichment button
    And user clicks on "Add a table to join" icon in "Enrich session_id" dialog
    And user picks "datagrok > public > users_sessions" from the open menu
    And user picks "session_id" as the key of the main table in "Enrich session_id" dialog
    And user opens the columns of the joined "users_sessions" table in "Enrich session_id" dialog
    And user toggles the "token_hash" column in the column list of "Select columns..." dialog
    And user toggles the "type" column in the column list of "Select columns..." dialog
    And user toggles the "user_id" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user enters "bdd-enrich-tokens-{time}" into Name input in "Enrich session_id" dialog
    And user clicks on OK button in "Enrich session_id" dialog
    Then "bdd-enrich-tokens-{time}" enrichment should be visible
    When user clicks on the "header event_type_id" area of grid
    Then the context panel should show "event_type_id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    When user clicks on add enrichment button
    Then "Enrich event_type_id" dialog should be visible
    When user clicks on "Add a table to join" icon in "Enrich event_type_id" dialog
    And user picks "datagrok > public > event_types" from the open menu
    And user picks "event_type_id" as the key of the main table in "Enrich event_type_id" dialog
    Then the join key in "Enrich event_type_id" dialog should read "event_type_id = id"
    When user opens the columns of the joined "event_types" table in "Enrich event_type_id" dialog
    And user toggles the "friendly_name" column in the column list of "Select columns..." dialog
    And user toggles the "source" column in the column list of "Select columns..." dialog
    And user toggles the "is_error" column in the column list of "Select columns..." dialog
    And user toggles the "error_severity" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user enters "bdd-enrich-types-{time}" into Name input in "Enrich event_type_id" dialog
    And user clicks on OK button in "Enrich event_type_id" dialog
    Then "bdd-enrich-types-{time}" enrichment should be visible
    Given user watches the task bar
    When user clicks on link of "bdd-enrich-types-{time}" enrichment
    Then the task bar should have finished "Enriching..."
    And "source" column should have no missing values
    When user drags the "x scroll handle" area of grid by 1000 pixels to the left
    And user clicks on the "header session_id" area of grid
    Then the context panel should show "session_id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    And user watches the task bar
    When user clicks on link of "bdd-enrich-tokens-{time}" enrichment
    Then the task bar should have finished "Enriching..."
    And the table should have the columns "id, friendly_name, session_id, event_type_id, event_time, description, error_message, error_stack_trace, exported_by, error_severity, result.friendly_name, is_error, source, token_hash, type, user_id"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Deleting one of the applied enrichments removes its configuration only
    When user drags the "x scroll handle" area of grid by 1000 pixels to the left
    And user clicks on the "header event_type_id" area of grid
    Then the context panel should show "event_type_id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    When user clicks on delete icon of "bdd-enrich-types-{time}" enrichment
    Then "bdd-enrich-types-{time}" enrichment should be absent
    And 0 enrichments named "bdd-enrich-types-{time}" should be on the server
    And the table should have a column "source"
    And the table should have a column "token_hash"
    And no errors should have been logged

  Scenario: The enrichment of session_id is offered on the other result of the events table
    When user clicks on events view
    Then the "events" view should be current
    When user clicks on the "header session_id" area of grid
    Then the context panel should show "session_id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    Then "bdd-enrich-tokens-{time}" enrichment should be visible
    And "bdd-enrich-sessions-{time}" enrichment should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A layout saved before enriched columns are removed does not bring them back
    Given the layouts named "bdd-enrich-events-{time}" are deleted when the feature ends
    When user clicks on bdd-enrich-events-{time} view
    Then the "bdd-enrich-events-{time}" view should be current
    And the table should have a column "token_hash"
    Given the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    When user clicks on Save button in layouts pane
    Then "bdd-enrich-events-{time}" layout card should be visible
    When user drags the "x scroll handle" area of grid by 1000 pixels to the right
    And user picks "Remove" from the context menu of the "header token_hash" area of grid
    Then the table should not have a column "token_hash"
    When user clicks on "bdd-enrich-events-{time}" layout card
    Then the table should not have a column "token_hash"
    And the table should have a column "user_id"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A project with the enriched table reopens with the enrichments that still exist
    Given no project named "bdd-enrich-project-{time}" is on the server
    When user closes events view
    Then the "bdd-enrich-events-{time}" view should be current
    When user clicks on Save button in toolbar
    Then "Save project" dialog should be visible
    When user enters "bdd-enrich-project-{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then "Save project" dialog should be hidden
    And 1 project named "bdd-enrich-project-{time}" should be on the server
    When user closes all views
    Given the browse panel is open
    When user refreshes the browse tree
    And user expands "My stuff" tree node inside browse tree
    And user double-clicks on "My stuff > bdd-enrich-project-{time}" tree node inside browse tree
    Then the table should have a column "type"
    And the table should have a column "user_id"
    And the table should not have a column "token_hash"
    And the table should not have a column "source"
    And an error balloon containing "Could not find enrichment" should have been shown
    And no errors should have been logged

  Scenario: The func_calls with a session run from the Browse tree
    Given the browse panel is open
    When user double-clicks on Databases---Postgres---Datagrok---bdd-enrich-calls-{time} tree node inside browse tree
    Then the "bdd-enrich-calls-{time}" view should be current
    And the table should have 20 rows
    And "session_id" column should have no missing values

  Scenario: On func_calls.session_id the enrichment of users_sessions.id is offered via users_sessions, not those of events.session_id
    When user clicks on the "header session_id" area of grid
    Then the context panel should show "session_id"
    Given "Datagrok" accordion header in context panel is expanded
    And Enrich accordion header in context panel is expanded
    Then "bdd-enrich-users-{time}" enrichment should be visible
    And "bdd-enrich-users-{time}" enrichment should contain text "via users_sessions"
    And "bdd-enrich-tokens-{time}" enrichment should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
