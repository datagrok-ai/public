# Home widgets features — what is missing

A gap audit (2026-10-01) of `home-widgets.feature` and `home-as-sharing-user.feature` against TestTrack
`PowerPack/Widgets/home_widgets_manual_tests.md` (15 cases) and the old Playwright specs
(`playwright-tests/e2e/Widgets`, `playwright-public/Widgets`). The features themselves were not
changed. What the library's rules leave out is not a gap: the Community links come from an outside
service, and `widgets-after-debug-delete.md` is an `apitest` case.

## 1. Gaps the existing vocabulary closes

Each block was run inside the feature it names on a local stand built from master (green twice) and
then taken out again; it can be pasted back as it is.

### Spotlight-03: "N unread" and Mark all as read (`home-as-sharing-user.feature`)

The feature says "Not translated: "N unread", the "Mark all as read" link and Mark all as read
itself". All three are claimable. After `notifications page of Spotlight home widget should contain
text "bdd-home-shared-{time}"` in "The share arrives as the one unread notification…":

```gherkin
    And notifications page of Spotlight home widget should contain text "1 unread"
    When user clicks on "Mark all as read" link in Spotlight home widget
    Then badge of Spotlight home widget should be absent
    And notifications page of Spotlight home widget should not contain text "1 unread"
    And the signed-in user should have 0 unread notifications on the server
```

The badge and the "N unread" line are removed only after `markAllAsRead()` has answered
(`spotlight-widget.ts`), so the negatives cannot read the state before it; the server count is polled.

### Controls-01: the close icon's "Remove" tooltip (`home-widgets.feature`)

Listed as not translated; it is shown on hover. In "The close icon shows on hover and removes the
widget", before the click on the icon:

```gherkin
    When user hovers over close icon of Community home widget
    Then tooltip should contain text "Remove"
```

### Spotlight-02: what the Learn tab lists (`home-widgets.feature`)

Only "Cheminformatics" is claimed. The VIDEO playlists and the WIKI links are static in PowerPack
(`learning-widget.ts`), the same on every stand. In "Spotlight has six tabs…":

```gherkin
    And the following elements should be visible:
      | "Meetings" text in Spotlight home widget         |
      | "Develop" text in Spotlight home widget          |
      | "Cheminformatics" text in Spotlight home widget  |
      | "Visualize" text in Spotlight home widget        |
      | "Explore" text in Spotlight home widget          |
    When user clicks on WIKI tab in Spotlight home widget
    Then the "WIKI" tab of Spotlight home widget should be showing
    And the following elements should be visible:
      | "Bioinformatics" text in Spotlight home widget |
      | "Overview" text in Spotlight home widget       |
      | "Access" text in Spotlight home widget         |
      | "Transform" text in Spotlight home widget      |
      | "Compute" text in Spotlight home widget        |
      | "Govern" text in Spotlight home widget         |
```

### Claims that can pass with the behaviour broken (`home-widgets.feature`)

- **A hidden widget after a reload.** In "A widget hidden in the Customize form stays hidden after a
  reload", `Community home widget should be absent` right after `user reloads the page` can read a
  Home page whose widgets are not built yet (they are appended after `getCurrentUserGroup()`
  resolves). Claim the list first:
  ```gherkin
      Then the Home page should show the widgets "Spotlight, Reports, Usage"
      And Community home widget should be absent
  ```
- **The close icon hidden before the hover** passes when the icon is not attached yet. First:
  `Then close icon of Community home widget should be present`.
- **The search hides the widgets** — claim they were shown before typing:
  `Then home widgets panel should be visible`.

## 2. What needs a new step or a signal

- **The tip of the day by weekday.** `the tip of the day of Spotlight home widget should open what it
  names` (`bindings/home.ts`) accepts whichever kind of tip the page shows; on a demo day the link text
  comes from the same list it is checked against, so a broken weekday mapping passes. The step should
  derive the expected kind from the day (a demo on Monday, Friday and the weekend, a tutorial on
  Wednesday, a plain tip on Tuesday and Thursday) and fail on another.
- **Usage Analysis finished loading.** After "Open Usage Analysis" the Overview view is current while
  its charts still load; `no errors should have been logged` reads once, and a later error is cleared
  at the next scenario's start. An `aria-busy` on the app's view, or an event when its viewers have
  data, would let the claim wait for it.
- **Pages of the account's own history.** The Favorites and My Activity pages and the Workspace hint
  ("Select a pinned item...") depend on what the account pinned and did; a part per page of `home
  widget` and a Given that pins an item for the feature (removed at its end) would make them claimable.
- **The DEMO and TUTORIALS sub-tabs** hold only where the Tutorials package is installed; a capability
  gate inside a journey skips the rest of the journey, so they need a scenario of their own behind
  `Given the "Tutorials" package is installed`.
