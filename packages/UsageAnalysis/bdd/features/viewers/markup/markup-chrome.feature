@journey @viewers @realizes:viewers.markup @realizes:entities.viewer.action.close-viewer
Feature: Markup chrome — Edit content, the title bar and the strip under the content
  The two ways a user changes what this viewer shows — **Edit content...** in its context menu and
  the `Content` property — and the title bar around it.
  The old spec right-clicked at `box.height - 40` to reach the menu item, a guess at where the
  content stops, and drove the dialog's editor with Control+A and typed characters. The menu is
  opened here through the viewer's own context menu, and the claim about the dialog is the one
  worth making: it opens on the content that is on screen. What it does on OK is not asserted, and
  that is a measurement, not an omission — see the note below.
  **Why the dialog is only opened and dismissed.** Adding this viewer opens a help pane for it in
  the bottom-right dock, the pane fetches `datagrok.ai/help/visualize/viewers/markup`, and on a
  stand with no route to the public internet that fetch never finishes: the pane keeps its
  `.grok-preloader` (`grok.css:1110` — `position:absolute; inset:0; z-index:100500;
  background:white; animation:fadeIn 3s`), which after three seconds is an opaque white rectangle
  over that whole panel, with `pointer-events:auto`. The Edit dialog is placed over that corner
  (measured: dialog 1211,548 534x320, its OK button at 1699,831; the preloader 1451,762 467x294),
  so every locator-driven click on the dialog's footer is intercepted and times out. The dialog is
  therefore dismissed with Escape, which does not go through the pointer. A help pane that never
  resolves should not leave a click-eating overlay behind.
  The viewer needs no special context-menu target: `constructContextMenu` is the ordinary viewer
  hook, so a right-click anywhere on it opens its own menu. An earlier draft of the status declared
  an `empty space` shortcut — the strip between the content's bottom and the viewer's — which was
  never reported, because `host` is the `.grok-help` div and it stretches to the viewer's full
  height. The area and the shortcut are gone; the last scenario is what replaced them.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a markup viewer with:
      | content | editable probe |
    Then the "text" reading of markup viewer should be "editable probe"
    And markup viewer should report no error

  Scenario: Edit content... is on the viewer's menu and opens on what is on screen
    When user opens the context menu of markup viewer
    Then the open menu should list "Edit content..."
    When user closes the context menu
    And user picks "Edit content..." from the context menu of markup viewer
    Then Edit dialog should be visible
    And input in Edit dialog should have the value "editable probe"
    When user presses Escape
    Then Edit dialog should be hidden
    And the "text" reading of markup viewer should be "editable probe"
    And "content" property of markup viewer should be "editable probe"
    And no errors should have been logged

  Scenario: The dialog opens on whatever the content is at the time, not on the default
    When user sets "content" property of markup viewer to "second probe"
    And user picks "Edit content..." from the context menu of markup viewer
    Then Edit dialog should be visible
    And input in Edit dialog should have the value "second probe"
    When user presses Escape
    Then Edit dialog should be hidden
    When user sets "content" property of markup viewer to "editable probe"
    Then the "text" reading of markup viewer should be "editable probe"
    And no errors should have been logged

  Scenario: The title the property sets is the title the bar shows
    Then title of markup viewer should not contain the text "Patient card"
    When user sets "title" property of markup viewer to "Patient card"
    Then title of markup viewer should have text "Patient card"
    When user sets "title" property of markup viewer to ""
    Then title of markup viewer should not contain the text "Patient card"
    And the "text" reading of markup viewer should be "editable probe"
    And no errors should have been logged

  Scenario: The title bar closes the markup viewer
    When user clicks on close icon of markup viewer
    Then markup viewer should be absent
    And the open tableview should have 0 markup viewers
    And no errors should have been logged

  Scenario: The viewer's own menu opens anywhere on it, and the content fills it
    Given user adds a markup viewer with:
      | content | one line |
    Then markup viewer should have a "content" area
    And markup viewer should not have an "empty space" area
    When user opens the viewer menu of markup viewer
    Then the open menu should list "Edit content..."
    When user closes the context menu
    Then no errors should have been logged
