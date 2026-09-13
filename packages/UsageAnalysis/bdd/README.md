# UsageAnalysis behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. `features/viewers/`
holds one folder per platform viewer (every TestTrack viewer spec translated, most as `@journey`
features: the data and the viewer opened once, the scenarios in order as soft steps) plus
`viewer-chrome.feature`, the outline over the title and description every viewer shares;
`features/spaces/` the Spaces features (the browse tree, the space view, sharing — the sharing
one shares with `DATAGROK_SHARING_LOGIN`, or with the `bddsecond` user the library's setup
creates when the variable is unset). `bindings/` keeps the steps only one
viewer can define (the bar chart's bar order and lengths, the pie chart's slices, the pivot's
aggregation against a `groupBy`, the correlation plot's coefficient against `DG.Stats`, the
Forms viewer's card rows, the tile viewer's designer, the filter panel's hierarchical card); the
rest of the vocabulary is the library's (`npx grok-bdd list-steps`).

From a fresh checkout of `public`, against a local stand on `http://localhost:8888` (another one:
`DATAGROK_URL=https://… npx grok-bdd run`):

```bash
cd public/libraries/bdd && npm ci && npm run build   # the library (a path dependency of this package; dist/ is not committed)
npx playwright install chromium                      # its browser, once per machine (here, not in the package)
cd ../../packages/UsageAnalysis && npm ci            # the package; npm links the library in and puts grok-bdd in .bin
npx grok-bdd link                                    # ONE Playwright: the library's copy into node_modules (redo after every npm ci)
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
