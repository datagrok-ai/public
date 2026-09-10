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

The stand needs a platform from `core` at or after 2026-09-10, the `Chem` package
(`grok s packages install Chem`), the dev key of the `localhost` entry in `~/.grok/config.yaml`
(`grok config`), and the 1000-row demog subset uploaded once:

```bash
grok s files put public/packages/ApiTests/files/datasets/demog-1000.csv "System:DemoFiles/demog-1000.csv" --host localhost
```

Editing: change a feature, `npx grok-bdd compile`, commit the regenerated spec with it;
`npx grok-bdd compile --verbose` prints how every element phrase resolves; `npx grok-bdd run
--trace on` records a trace with DOM snapshots; `PLAYWRIGHT_JSON_OUTPUT_NAME=run.json npx grok-bdd
run --reporter=list,json` gives per-step timings. A failed step reports its feature line, the
step, the reason, and what the page shows instead — see "Reading a failure" in the library README.
