@journey @viewers @realizes:viewers.markup
Feature: Markup content and the interpretation mode the viewer resolved
  What the markup viewer put in the DOM, and which of the four interpretation modes it actually
  used to get there.
  The mode is the point. `Mode` defaults to **Auto**, and under Auto the viewer calls `detectMode`
  — content starting with `<` is Html, everything else is Markup (`markup_viewer_core.dart:88-92`)
  — so the property says "Auto" in both cases and no property anywhere says which branch was taken.
  The old spec told the four modes apart by counting `h1` and `pre` elements in the rendered html;
  the viewer now reports the **resolved** mode as a reading, so "Auto read this as Html" is a claim
  the viewer makes about itself, and the element counts are what it produced, not how it was
  inferred.
  One old assertion is deliberately dropped: `getComputedStyle(strong).fontWeight >= 700`. That
  number comes from the user agent's stylesheet's `strong`, not from anything this viewer does; the
  claim worth making is that the markup pass produced a `strong` at all, which is `bold`.
  Fixture: demog-1000 and the default Markdown sample the viewer opens on — one `h1`, four list
  items, four links, no `pre`. `content` is scenario-scoped: every scenario below sets the content
  it reasons about, and puts `mode` and `markupEnabled` back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a markup viewer
    Then 1000 rows should pass the filter
    And the "mode" reading of markup viewer should be "Markup"
    And "mode" property of markup viewer should be "Auto"
    And the "markup enabled" reading of markup viewer should be "true"
    And markup viewer should report no error

  Scenario: The viewer opens on the Markdown sample rendered, not on its source
    Then the "heading 1" reading of markup viewer should be "What’s Markdown?"
    And the "list items" reading of markup viewer should be 4
    And the "links" reading of markup viewer should be 4
    And the "preformatted" reading of markup viewer should be 0
    And the "text" reading of markup viewer should include the text "Markdown is a lightweight markup language"
    And the "text" reading of markup viewer should not include the text "#"
    And the "text" reading of markup viewer should not include the text "https://"
    And no errors should have been logged

  Scenario: Auto resolves the mode from the content, and the property never says which
    When user sets "content" property of markup viewer to "# Heading probe"
    Then "mode" property of markup viewer should be "Auto"
    And the "mode" reading of markup viewer should be "Markup"
    And the "heading 1" reading of markup viewer should be "Heading probe"
    When user sets "content" property of markup viewer to "<b>bold probe</b> plain probe"
    Then "mode" property of markup viewer should be "Auto"
    And the "mode" reading of markup viewer should be "Html"
    And the "bold" reading of markup viewer should be 1
    And the "heading 1" reading of markup viewer should be ""
    And no errors should have been logged

  Scenario: Mode None shows the viewer's own markup verbatim, escaped, in a pre
    When user sets "content" property of markup viewer to "<b>bold probe</b> plain probe"
    And user sets "mode" property of markup viewer to "None"
    Then the "mode" reading of markup viewer should be "None"
    And the "preformatted" reading of markup viewer should be 1
    And the "bold" reading of markup viewer should be 0
    And the "text" reading of markup viewer should include the text "<b>bold probe</b> plain probe"
    When user sets "mode" property of markup viewer to "Auto"
    Then the "mode" reading of markup viewer should be "Html"
    And the "preformatted" reading of markup viewer should be 0
    And the "bold" reading of markup viewer should be 1
    And the "text" reading of markup viewer should not include the text "<b>"
    And no errors should have been logged

  Scenario: Mode Html leaves a hash a hash; Mode Markup makes it a heading whatever the content looks like
    When user sets "content" property of markup viewer to "# Heading probe"
    And user sets "mode" property of markup viewer to "Html"
    Then the "mode" reading of markup viewer should be "Html"
    And the "heading 1" reading of markup viewer should be ""
    And the "text" reading of markup viewer should include the text "# Heading probe"
    When user sets "mode" property of markup viewer to "Markup"
    Then the "mode" reading of markup viewer should be "Markup"
    And the "heading 1" reading of markup viewer should be "Heading probe"
    And the "text" reading of markup viewer should not include the text "#"
    When user sets "mode" property of markup viewer to "Auto"
    Then the "mode" reading of markup viewer should be "Markup"
    And no errors should have been logged

  Scenario: The markup pass turns the emphasis marks into a strong element and eats them
    When user sets "mode" property of markup viewer to "Markup"
    And user sets "content" property of markup viewer to "**md bold** plain tail"
    Then the "bold" reading of markup viewer should be 1
    And the "text" reading of markup viewer should be "md bold plain tail"
    And the "text" reading of markup viewer should not include the text "**"
    When user sets "mode" property of markup viewer to "Auto"
    Then the "bold" reading of markup viewer should be 1
    And no errors should have been logged

  Scenario: A Markdown list is a list of items, not four lines of source
    When user sets "content" property of markup viewer to "* one\n* two\n* three"
    Then the "mode" reading of markup viewer should be "Markup"
    And the "list items" reading of markup viewer should be 3
    And the "preformatted" reading of markup viewer should be 0
    And the "text" reading of markup viewer should not include the text "*"
    When user sets "mode" property of markup viewer to "None"
    Then the "list items" reading of markup viewer should be 0
    And the "preformatted" reading of markup viewer should be 1
    And the "text" reading of markup viewer should include the text "* one"
    When user sets "mode" property of markup viewer to "Auto"
    Then the "list items" reading of markup viewer should be 3
    And no errors should have been logged

  Scenario: Markup Enabled decides whether the table expressions are evaluated at all
    When user sets "content" property of markup viewer to "Rows: #{t.rowCount}"
    Then the "text" reading of markup viewer should be "Rows: 1000"
    When user sets "markupEnabled" property of markup viewer to "false"
    Then the "markup enabled" reading of markup viewer should be "false"
    And the "text" reading of markup viewer should be "Rows: #{t.rowCount}"
    When user sets "markupEnabled" property of markup viewer to "true"
    Then the "markup enabled" reading of markup viewer should be "true"
    And the "text" reading of markup viewer should be "Rows: 1000"
    And no errors should have been logged
