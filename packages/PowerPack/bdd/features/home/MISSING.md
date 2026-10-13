# Home widgets features — what is missing

A gap audit (2026-10-01) of `home-widgets.feature` and `home-as-sharing-user.feature` against TestTrack
`PowerPack/Widgets/home_widgets_manual_tests.md` (15 cases) and the old Playwright specs
(`playwright-tests/e2e/Widgets`, `playwright-public/Widgets`); triaged with the lead on 2026-10-11. What
the library's rules leave out is not a gap: the Community links come from an outside service, and
`widgets-after-debug-delete.md` is an `apitest` case.

## Resolved (2026-10-12)

- **Spotlight-03: "N unread" and Mark all as read** — `home-as-sharing-user` claims "1 unread", clicks Mark all
  as read, and claims the badge and the line gone and 0 unread notifications on the server.
- **Controls-01: the close icon's "Remove" tooltip** — `home-widgets` hovers the icon and claims the tooltip.
- **Spotlight-02: what the Learn tab lists** — the VIDEO playlists (Meetings, Develop, Cheminformatics, Visualize,
  Explore) and the WIKI links (Bioinformatics, Overview, Access, Transform, Compute, Govern), static in PowerPack
  (`learning-widget.ts`).
- **Claims that could pass with the behaviour broken** — after a reload the widget list is claimed before the
  hidden widget's absence; the close icon is claimed present before it is claimed hidden; the widgets panel is
  claimed visible before the search hides it.
- **The tip of the day by weekday** — `the tip of the day of … should open what it names` requires the kind of tip
  the weekday gives in the browser (a demo Monday, Friday and the weekend, a tutorial Wednesday, a plain tip
  Tuesday and Thursday; `spotlight-widget.ts` `TIP_TYPE_BY_DAY`), so a broken mapping fails.

## Still open

Nothing planned.

## Won't do

- **Usage Analysis finished loading.** After "Open Usage Analysis" the Overview view is current while its charts
  still load; a later error is charged to the next scenario. It needs an `aria-busy` on the app's view or an event
  when its viewers have data — a UsageAnalysis signal for one claim.
- **Pages of the account's own history.** Favorites, My Activity and the Workspace hint ("Select a pinned
  item...") depend on what the account pinned and did; each needs a part of `home widget` and a Given that pins an
  item for the feature.
- **The DEMO and TUTORIALS sub-tabs** hold only where the Tutorials package is installed; they would need a
  scenario of their own behind `Given the "Tutorials" package is installed`. The tabs themselves are claimed.
