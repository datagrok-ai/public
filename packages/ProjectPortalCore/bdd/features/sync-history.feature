@portal
Feature: A portfolio's sync history
  A portfolio may declare a schedule and owners. syncDuePortfolios syncs only the portfolios whose
  schedule is due; every sync, by schedule or by hand, is recorded with what it created, updated
  and flagged and the stage spellings it could not resolve. The board's Sync history shows the
  schedule, the owners and the runs, newest first.

  Background:
    Given user opens the Project Portal app

  Scenario: The sync history of a portfolio with a schedule
    Given a sync of the "EX-SMD" portfolio has run
    When user clicks on small molecule discovery portfolio
    And user clicks on sync history button
    Then sync history dialog should be visible
    And sync schedule should contain text "Scheduled: every 6h"
    And sync owners should contain text "admin"
    And latest sync run should contain text "by hand"
