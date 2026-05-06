# playwright-public

End-to-end Playwright tests for the public Datagrok platform — `connections/`,
`queries/`, `scripts/`. Designed to be invoked through the standard `grok test`
runner so it shares auth, host config, and CSV reporting with the rest of the
test pipeline.

## Run on dev

From this directory:

```bash
npm install                                  # first time
npx playwright install chromium              # first time
grok test --skip-puppeteer --host dev        # full suite
grok test --skip-puppeteer --host dev --category scripts   # one suite
grok test --skip-puppeteer --host dev --category connections
grok test --skip-puppeteer --host dev --category queries
```

`--skip-puppeteer` (added to `datagrok-tools` for this directory) tells `grok
test` to skip the Puppeteer / `DG.Test` pass entirely and run only the
Playwright tests. The runner reads `~/.grok/config.yaml` to resolve the host
alias and exchanges the dev key for an auth token, then exports
`DATAGROK_URL` + `DATAGROK_AUTH_TOKEN` to the Playwright process.

`e2e/global-setup.ts` boots Chromium, drops the token into cookie +
`localStorage`, and writes the resulting storage state to the three paths the
existing specs reference (`e2e/.auth.json`, `e2e/.auth.public.json`,
`./.auth.json`). All three are gitignored.

## CSV report

Pass `--csv` to also emit `test-report.csv` and `test-report-playwright.csv`
in this directory — same format the rest of the platform pipeline ingests.

## Required env vars (test-specific)

| Var                          | Used by                                              |
|------------------------------|------------------------------------------------------|
| `DG_PG_PASSWORD`             | connections/identifiers, queries Postgres lifecycle  |
| `DG_PG_EXT_LOGIN`            | connections/external-provider                        |
| `DG_PG_EXT_PASSWORD`         | connections/external-provider                        |
| `DG_OPENWEATHERMAP_API_KEY`  | connections/import-swagger                           |

Tests that need these vars `test.skip(...)` when they are missing.
