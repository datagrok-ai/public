# UsageAnalysis behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. `features/viewers/`
holds one folder per platform viewer (every TestTrack viewer spec translated, most as `@journey`
features: the data and the viewer opened once, the scenarios in order as soft steps) plus
`viewer-chrome.feature`, the outline over the title and description every viewer shares;
`viewers/grid/grid-context-menu.feature` is not a translation but the reproduction of a bug
(2026-09-22: a right click put the current row back where it was and scrolled there), kept as the
proof that a right-clicked cell becomes current and that the Current Value actions of Chem and Helm
act on it — the stand needs those two packages;
`features/viewers/legend/` the Legend TestTrack section, translated from its manual-case md files
(seven viewers sharing one legend column, the legend under filters, its placement, molecules in
it; the section's scatter plot and line chart cases went into those viewers' legend features);
`features/spaces/` the Spaces features (the browse tree, the space view, sharing — the sharing
one shares with `DATAGROK_SHARING_LOGIN`, or with the `bddsecond` user the library's setup
creates when the variable is unset); `features/users-groups-roles/` Browse > Platform > Users,
Groups and Roles (the views, the New dialogs, memberships, disabling, favorites, global
permissions). A user can never be deleted, so the features share two fixture users made once per
stand (`bddviewed`, `bddmanaged`), and only users-create adds users: the two it makes per run,
`bdd<time>` and `bdd-svc<time>`, add up on a shared stand and go with a fresh CI database. The
group and role features, and users-manage, which makes groups and roles too, are `@serial`: their
gallery searches are fuzzy and bring up each other's fixtures, so they take turns while the rest runs
in parallel. Every TestTrack case there is
translated but Groups-19 (a group cannot be added to favorites); `features/browse/` the Browse
panel itself (its toolbar, the tree and its keyboard, browsing versus persistent views, Files, My
stuff, Platform, Databases, Apps, Dashboards, the context panel and menus, and the per-section
error matrix), translated from the manual cases, each feature naming what it left out and why.
Two of its scenarios are `@full-stand` (they name the providers and the Platform sections a full
stand carries) and one is `@compute` (the Model Hub needs the Compute package): a smaller stand
runs with `grok-bdd run --grep-invert "@full-stand|@compute"`. `features/guides/` holds the
answers to "how do I …" questions as scenarios: `grok-bdd guide features/guides/<name>.feature`
films one into `guides/<feature>/<scenario>/guide.mp4` with the numbered steps and pictures beside
it (`steps.md`), `--help-pages` re-films every `@help:`-tagged one into the help tree; `INDEX.md`
there lists the questions answered. They run with the suite, so an answer that stops being true
fails. `bindings/` keeps the steps only one
viewer can define (the bar chart's bar order and lengths, the pie chart's slices, the pivot's
aggregation against a `groupBy`, the correlation plot's coefficient against `DG.Stats`, the
Forms viewer's card rows, the tile viewer's designer, the filter panel's hierarchical card); the
rest of the vocabulary is the library's (`npx grok-bdd list-steps`).

The [known-failure audit](../../../libraries/bdd/KNOWN_FAILURES.md) records the current defects,
their observed failures and causes. The line-chart lasso scenario now passes without a tag:
checkbox menu items keep the menu open, so close it before dragging on the chart.

`features/viewers/filter-panel/` stands in for the TestTrack scenarios of
`files/TestTrack/Viewers/FilterPanel/` — `panel-core-ladder.md`, `add-remove-entry-points.md`,
`filter-type-selection-modes.md`, `hierarchical-and-combined-boolean.md`,
`compose-viewer-filtering.md`, `expression-text-filters.md`, `cloned-view-sync.md`,
`collaborative-filtering-for-linked-tables.md`, `save-and-reapply-state.md` and
`filter-summary-ui.md` (not `bio-filters.md`, which belongs to Bio) — plus Scenario 4 of
`PieChart/piechart-onclick-select-filter.md` and Scenario 1 of
`TrellisPlot/trellis-plot-click-to-filter.md` in the click-filter outline of
`compose-with-viewers.feature`. Each feature says in its description what of its md it does not
translate, and why.

From a fresh checkout of `public`, against a local stand on `http://localhost:8888` (another one:
`DATAGROK_URL=https://… npx grok-bdd run`):

```bash
cd public && grok setup                              # once per checkout: the pnpm workspace (needs `npm i -g datagrok-tools`)
cd libraries/bdd && npm run build                    # the library (a workspace dependency of this package; dist/ is not committed)
npx playwright install chromium                      # its browser, once per machine
cd ../../packages/UsageAnalysis/bdd                  # this bdd project; the workspace already links the library, no `grok-bdd link`
npx grok-bdd run --reporter=list                     # compile --check, then Playwright on 4 workers
PLAYWRIGHT_WORKERS=2 npx grok-bdd run                # a stand whose pub serve or datlas falls behind at 4 (bundle loads past 30 s, 502s)
npx grok-bdd run --workers 2 generated/viewers/box-plot   # any Playwright flag or path passes through
```

For a core dev server on port **8889**, set the URL once in the terminal session, then run from
the package directory (the localhost dev key supplies authentication):

```bash
export DATAGROK_URL=http://localhost:8889
npx grok-bdd run --workers 1 --reporter=list \
  generated/viewers/viewer-chrome.test.ts -g 'box plot shows and clears'
```

Set `DATAGROK_LOGIN` and `DATAGROK_PASSWORD` for login-form authentication when no dev key
is available. If the run reports `Requiring @playwright/test second time`, repeat
`npx grok-bdd link` here: the local library and package must use one physical Playwright
installation, even when their installed versions match.

The library maps `Control` / `Ctrl` to Command on Mac for shortcuts and selection gestures,
including typing and clearing, so existing features stay portable. `Shift+Delete` maps to
Shift+Backspace on Mac, following Datagrok's delete-command binding. Physical Control is
available as `ControlLeft`; the custom current-cell copy uses `ControlLeft+Shift+C` because
its handler specifically requires it. Rebuild `libraries/bdd` after updating those helpers.

The stand needs a platform from `core` at or after 2026-09-10, the dev key of the `localhost`
entry in `~/.grok/config.yaml` (`grok config`), and the packages used by the viewers:

```bash
grok s packages install PowerGrid GIS Charts Chem Curves PowerPack --host localhost
```

The installed packages must include the automation support in this checkout. If a viewer reports
no readings, build and publish its current package with `npx webpack && grok publish localhost`
from that package directory. Forms comes from PowerGrid and also needs its
`@datagrok-libraries/utils` dependency linked to this checkout; Map comes from GIS and Word cloud
from Charts. A registry version can be current while still predating these source changes.

Upload the 1000-row demog subset once (from the core repository root):

```bash
grok s files put public/packages/ApiTests/files/datasets/demog-1000.csv "System:DemoFiles/demog-1000.csv" --host localhost
```

Editing: change a feature, `npx grok-bdd compile`, commit the regenerated spec with it;
`npx grok-bdd compile --verbose` prints how every element phrase resolves; `npx grok-bdd run
--trace on` records a trace with DOM snapshots; `PLAYWRIGHT_JSON_OUTPUT_NAME=run.json npx grok-bdd
run --reporter=list,json` gives per-step timings. A failed step reports its feature line, the
step, the reason, and what the page shows instead — see "Reading a failure" in the library README.
