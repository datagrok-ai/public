# UsageAnalysis behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. Today: the viewer
features (`features/viewers/`), the TestTrack box plot spec as one `@journey`. The full guide, the
vocabulary and the stand requirements are in the library's README.

From a fresh checkout of `public`, against a local stand on `http://localhost:8888` (another one:
`DATAGROK_URL=https://… npx grok-bdd run`):

```bash
cd public/libraries/bdd && npm ci && npm run build   # the library (a path dependency of this package; dist/ is not committed)
npx playwright install chromium                      # its browser, once per machine (here, not in the package)
cd ../../packages/UsageAnalysis && npm ci            # the package; npm links the library in and puts grok-bdd in .bin
npx grok-bdd link                                    # ONE Playwright: the library's copy into node_modules (redo after every npm ci)
npx grok-bdd run --reporter=list                     # compile --check, then Playwright
```

The stand needs a platform from `core` at or after `6983855e91` (2026-09-07), the `Chem` package
(`grok s packages install Chem`), the dev key of the `localhost` entry in `~/.grok/config.yaml`
(`grok config`), and the 1000-row demog subset uploaded once:

```bash
grok s files put public/packages/ApiTests/files/datasets/demog-1000.csv "System:DemoFiles/demog-1000.csv" --host localhost
```

Editing: change a feature, `npx grok-bdd compile`, commit the regenerated spec with it;
`npx grok-bdd list-steps` prints every phrase this package can use; `npx grok-bdd compile --verbose`
prints how every element phrase resolves; `npx grok-bdd run --trace on` records a trace with DOM
snapshots; `PLAYWRIGHT_JSON_OUTPUT_NAME=run.json npx grok-bdd run --reporter=list,json` gives
per-step timings. A failed step reports its feature line, the step, the reason, and what the page
shows instead (the visible menu items, the nearest property captions) — see "Reading a failure"
in the library README.
