@connections
Feature: A SPARQL connection
  Sparql is one of the providers the Databases tree hides until "Show more" is clicked; its New
  connection dialog asks for an endpoint and prefixes, and the connection is saved and deleted like
  any other. Translated from TestTrack Connections/sparql.md and sparql-ui.md (playwright-public
  connections/08-sparql.test.ts). The connection, BDD-Conn-Sparql-{run}, is deleted at feature end
  and checked gone.

  Not translated, and why: TEST against the md's endpoint (http://data.ontotext.com/repositories/
  data-last) — the endpoint gave no answer within a minute when probed from dev on 2026-09-22, so a
  claim on its result would test that public server, not Datagrok. The scenarios claim that TEST is
  offered (sparql-ui.md says it was missing; it is there now) and leave its result out.

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded

  Scenario: Show more reveals Sparql among the hidden providers
    Then Databases---Sparql tree node inside browse tree should be hidden
    When user clicks on "ellipsis-h" icon inside Databases---Show-more tree node inside browse tree
    Then Databases---Sparql tree node inside browse tree should be visible
    And Databases---Show-more tree node inside browse tree should be hidden
    # the revealed providers run below the fold; a click scrolls the node into view, the context-menu
    # gesture does not (library candidate)
    When user clicks on Databases---Sparql tree node inside browse tree
    And user opens the context menu of Databases---Sparql tree node inside browse tree
    Then the open menu should list "New connection..."
    And the open menu should list "Browse connections"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: A Sparql connection is saved from its dialog and deleted
    Given no connection named "BDD-Conn-Sparql-{run}" is on the server
    And the hidden providers of the Databases tree are shown
    When user clicks on Databases---Sparql tree node inside browse tree
    And user picks "New connection..." from the context menu of Databases---Sparql tree node inside browse tree
    Then "Add new connection" dialog should be visible
    And the following elements should be visible:
      | Endpoint input in "Add new connection" dialog |
      | Prefixes input in "Add new connection" dialog |
      | TEST button in "Add new connection" dialog    |
    When user enters "BDD-Conn-Sparql-{run}" into Name input in "Add new connection" dialog
    And user enters "http://data.ontotext.com/repositories/data-last" into Endpoint input in "Add new connection" dialog
    And user clicks on OK button in "Add new connection" dialog
    Then the "Add new connection" dialog should close
    And 1 connection named "BDD-Conn-Sparql-{run}" should be on the server
    Given Databases---Sparql tree node inside browse tree is expanded
    When user picks "Delete..." from the context menu of Databases---Sparql---BDD-Conn-Sparql-{run} tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 connections named "BDD-Conn-Sparql-{run}" should be on the server
    And no error or warning balloon should have been shown
