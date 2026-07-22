# Datagrok-tools changelog

## 6.4.7 (2026-07-22)

* GROK-20452: `grok publish` — generate `dockerfiles/celery` / `dockerfiles/queue` **before** gathering files for the zip. They were generated after, so the generated Dockerfile never reached the server: `meta.queue` functions fell back to a server-side container named `<pkg>-queue-celery` while the client had built `<pkg>-queue` — image validation 404'd and the worker container never started ("Container is not started" on every call).

## 6.4.6 (2026-07-22)

* `grok test` — the whole-run Puppeteer cap (60 min) is now overridable via `GROK_TEST_INVOCATION_TIMEOUT_MS`. When the cap fires mid-pass all collected results are discarded (empty `test-report.csv`, "Passed tests: 0"), which CI reports as "no results produced"; loaded CI agents can now raise the cap instead.
* `grok test` — `--no-retry` is now honored for Playwright runs. minimist parsed `--no-retry` as `{retry:false}`, so the flag was silently dropped and failed specs were still retried once; normalized so `--retries=0` reaches Playwright.

## 6.4.5 (2026-06-23)

* `grok test` — screen recording (`--record`) is now optional: if `page.screencast()` can't start (e.g. `ffmpeg` is missing on the runner, `spawnSync ffmpeg ENOENT`), the run warns and continues instead of failing the whole Puppeteer pass with 0 tests executed.

## 6.4.4 (2026-06-22)

* `grok publish` — fixed every publish failing with a silent `exit 1` after a successful upload: a stray `fs.unlinkSync('zip')` threw `ENOENT` (the archive is streamed in-memory, no `zip` file is ever written), and the surrounding `catch` only logged under `--verbose`. The errant unlink is removed and publish errors are now always surfaced.

## 6.4.3 (2026-06-22)

* `grok api` — numeric IVP model inputs are now generated with `nullable: false`, so an emptied input field fails form validation instead of running the solver with a null.

## 6.4.2 (2026-06-18)

* `grok publish` — registry-aware Docker fallback: when a package's image isn't built locally and the target server has no compatible record, `grok publish` now checks the configured registry and Docker Hub (`docker manifest inspect`) for the expected `datagrok/<name>:<version>` (and content-hashed) tag and uses it, instead of reporting "No fallback available" and failing. Fixes dependency publishes (e.g. Bio → @datagrok/chem) on CI runners where the image exists in the registry but not locally.

## 6.4.1 (2026-06-18)

* Fixed `grok` failing with `Cannot find module './commands/build'` — the `.npmignore` `build.js` rule was unanchored and excluded the compiled `bin/commands/build.js` from the published package; anchored it to `/build.js`.
* Pinned `ignore-walk` to ^6.0.5 so the package still supports Node 18 (9.x requires Node 22+).

## 6.4.0 (2026-06-18)

* Dependencies: sanitized and updated all dependencies; `npm install` is now warning-free and `npm audit` reports 0 vulnerabilities (was 24).
* Dependencies: migrated linting to ESLint 9 flat config (`eslint.config.mjs` + `typescript-eslint` + `@stylistic`), dropping the archived `eslint-config-google`.
* Dependencies: upgraded Puppeteer to v24 and migrated to its native `page.screencast()` for `--record`, removing `puppeteer-screen-recorder` and the deprecated `fluent-ffmpeg`.
* Dependencies: replaced `archiver-promise` with `archiver` directly, and replaced `@babel/cli` with a small `@babel/core` build script (`build.js`) to drop deprecated transitive packages (glob@7, inflight).

## 6.3.3 (2026-06-16)

* Fixed Celery Docker image generation — the image wasn't built locally on publish.

## 6.3.2 (2026-06-15)

* `func-gen` webpack plugin — generated RichFunctionView model inputs now use the script-form names (argument bounds/step `_t0`/`_t1`/`_h`, loop count `_count`) instead of the deprefixed forms, so the run, fitting, and sensitivity-analysis paths share one set of input names with diff-grok's pipeline. Fixes `Inconsistent inputs: "_t0" is missing` when starting fitting/SA from a Rich Function View.
* `func-gen` webpack plugin — the generated RichFunctionView model output annotation appends the `DiffStudio Facet` viewer (last), alongside Grid and Line chart.
* `grok stresstest` — run the ApiTests Node test runner directly via tsx (`src/package-test-node.ts`) instead of a compiled `dist-node` bundle; dropped the obsolete `npm run build-node` step and the `tsconfig-paths-bootstrap.js` / `dist-node` invocation.

## 6.3.0 (2026-06-12)

* `func-gen` webpack plugin — generates RichFunctionView model wrappers from Diff Studio `#meta.role: model` `.ivp` files. Inputs/output are derived from the parsed IVP, `#meta.icon` becomes the model icon, and `#meta.inputs` lookups are emitted as real `propagateChoice: all` inputs (rendered natively by the Rich Function View, unlike `meta.inputs`). Ships a prebuilt, tree-shaken CJS bundle of diff-grok's IVP parser (`plugins/ivp-parser.bundle.cjs`); regenerate with `npm run update:ivp-parser`.

## 6.2.6 (2026-05-26)

* `grok s tables upload` — accepts `.d42` binary blobs in addition to `.csv`. Content-Type is auto-detected from the file extension (`application/octet-stream` for d42, `text/csv` otherwise); server content-negotiates and persists either form against the same `/public/v1/tables/{name}` endpoint.

## 6.2.5 (2026-05-21)

* `grok report read` — renamed `--extract-actions` to `--extract-client-log`; sidecar is now `<stem>_client_log.json`. The old flag is no longer accepted.
* `grok report read` — legacy `actions` field in pre-consolidation report zips is folded into `clientLog` at read time (stderr warning emitted), so downstream consumers see one canonical field. Companion to the platform-side merge that drops `reports_data.actions` in favor of `client_log`.

## 6.2.4 (2026-05-13)

* GROK-20097: package template — replaced deprecated `"moduleResolution": "node"` in `tsconfig.json` with `"bundler"`, suppressing the TypeScript warning for new packages.

## 6.2.3 (2026-05-06)

* `grok test` — added `--skip-puppeteer` flag (symmetric counterpart of `--skip-playwright`). Bypasses `loadPackages()` and the Puppeteer/`DG.Test` runner so Playwright-only test directories (e.g. `public/playwright-public`) can run end-to-end via `grok test --skip-puppeteer` without tripping the Dart/JS package loader.

## 6.2.2 (2026-05-05)

* `grok test` — Playwright runner now writes `test-report-playwright.csv` next to the existing merged `test-report.csv`, so CI can ship Playwright rows to a dedicated Datlas reporting bucket without disturbing the legacy `package` flow.

## 6.2.1 (2026-05-05)

* Reports: `grok report attach <ticket> <file>` — upload a file as a JIRA issue attachment via REST v2 multipart POST.

## 6.2.0 (2026-05-04)

* `grok test` — Playwright support: when a package's `package.json` declares `"playwrightTests": "<path>"`, `grok test` runs `npx playwright test` against that directory in addition to the existing Puppeteer pass and merges results into a single `test-report.csv`. Auth is unified with the Puppeteer pass (dev key from `~/.grok/config.yaml` → session token → cookie + `localStorage` injection — no login form). Optional `DATAGROK_DEV_KEY_2` env var enables a second-user identity for specs that need it (`DATAGROK_AUTH_TOKEN_2` exposed to specs). New `--skip-playwright` flag opts out of the Playwright pass for a single run.

## 6.1.14 (2026-05-01)

* Reports: `grok report comment` now converts Markdown body to JIRA wiki markup before POSTing, fixing rendered headings/list/HTML-entity mismatches in JIRA UI.

## 6.1.13 (2026-04-30)

* Reports: `grok report ticket` now uses direct JIRA REST honoring `$JIRA_PROJECT`, replacing the Datlas-mediated path that hardcoded GROK.

## 6.1.11 (2026-04-27)

* `grok report read <path | instance number>` — normalize a report zip/json into one JSON object on stdout (envelope unwrap, `_meta.json` merge, optional `--extract-screenshot` / `--extract-d42` / `--extract-actions`)
* `grok report fetch` — single round-trip via the new `/reports/by-number/<number>/zip` endpoint, with fallback to the legacy search-then-download path on 404

## 6.1.10 (2026-04-17)

* `grok s connections save` — create/update a connection from a JSON file (`--save-credentials` optional)
* `grok s connections test` — test connectivity by JSON body or by existing id/name
* `grok s users save` — create/update a user from a JSON file
* `grok s groups save` — create/update a group from a JSON file (`--save-relations` optional)
* `grok s shares add` — share an entity with one or more groups (`--access View|Edit`)
* `grok s shares list` — list who an entity is shared with
* Tools: Normalize package name to lowercase in grok create, preserve original as friendlyName

## 5.1.3 (2026-02-03)

* GROK-19407:
  * Improved developer experience for grok check, grok api, grok publish
  * Fixed warnings all over the packages when using grok check

### Bug Fixes

## 4.14.72 (2026-01-13)

### Bug Fixes

* Include ts-morph in dependencies

## 4.14.72 (2026-01-13)

### Bug Fixes

* Removed local dependency on js-api

## 4.14.71 (2026-01-10)

### Features

* Grok Check: Added validation to restrict datagrok-api imports to public API paths only (GROK-19504)
* Only allows imports from 'datagrok-api/dg', 'datagrok-api/grok', 'datagrok-api/ui'
* Rejects deep imports like 'datagrok-api/dg/events' or 'datagrok-api/src/...'
* Validates both import and export statements

### Bug Fixes

* Fixed invalid datagrok-api imports in UsageAnalysis, ChatGPT, and PowerPack packages


## 4.14.70 (2026-01-10)

### Bug Fixes

* Grok Check: Fixed `--recursive` flag to check all packages instead of stopping after the first package
* Grok Check: Fixed process exit behavior in recursive mode to continue checking remaining packages even when errors are found


## 4.14.26 (2025-07-10)

### Features

* Grok Decorators added support for ViewBase type
* Grok Decorators updated comments and line feeds
* Grok Api updated comments and line feeds
* Grok Api minor type update

## 4.14.19 (2025-06-23)

### Features

* Grok Check added header tag validation
* Grok Api added optional parameters


## 4.14.14 (2025-06-16)

### Features

* Grok Link added ability to link packages from root directory
* Grok Link added repo-only option

## 4.14.10 (2025-06-08)

### Features

* Grok Check checks update


## 4.14.9 (2025-06-06)

### Features

* Grok Check added soft mode 


## 4.14.8 (2025-06-06)

### Features

* Grok Check added collision checks 
* Grok Publish runs check before invocation 
* Grok Api creates api for the scripts in the python directory.


## 4.14.6 (2025-05-29)

### Features

* Decorators added ability to set hardcoded values
* Decorators added ability to get outputs from multiple parameters 


## 4.14.2 (2025-05-16)

### Features

* Grok Link added path mode
* Grok Link added unlink mode

## 4.14.0 (2025-05-12)

### Features

* Added annotation creation by decorators

## 4.13.70 (2025-04-14)

### Features

* Added error output

## 4.13.67 (2025-03-26)

### Features

* Updated grok link:
*   Excluded diff-grok from linking
*   Added verbose option
*   Updated help

## 4.13.66 (2025-03-21)

### Features

* Added js-flag to puppeteer options

## 4.13.65 (2025-03-18)

### Features

* Fixed func-gen-plugin for js packages 

## 4.13.64 (2025-02-27)

### Features

* Test/TestAll opens inspector in debug mode

## 4.13.63 (2025-02-27)

### Features

* Grok Link minor fixes

## 4.13.62 (2025-02-27)

### Features

* No Sandbox mode for puppeteer removed
* Added debug mode for test and testAll

## 4.13.61 (2025-01-27)

### Features

* No Sandbox mode for puppeteer added

## 4.13.60 (2025-01-09)

### Features

* Updated console output for cases "Tests not found"

## 4.13.59 (2025-01-09)

### Features

* Benchmark testing fix
* Updated csv output result(it saves benchmark and stress data now)


## 4.13.58 (2025-01-09)

### Features

* Test failed results fix


## 4.13.57 (2025-01-08)

### Features

* Test failed results fix

## 4.13.56 (2025-01-06)

### Features

* Worker to browser in test all 

## 4.13.55 (2025-01-06)

### Features

* Removed unnecessary outputs in grok check


## 4.13.54 (2025-01-06)

### Features

* Grok check regex update


## 4.13.53 (2025-01-06)

### Features

* Updated csv output for test all(added workes id)


## 4.13.52 (2024-12-30)

### Features

* Refactored workers to browsers


## 4.13.51 (2024-12-30)

### Features

* Minor fix for testToWorker order in testsAll


## 4.13.50 (2024-12-30)

### Features

* Added testToWorker order for testsAll

## 4.13.49 (2024-12-26)

### Features

* Removed load npm
* Added ability to set vue as external lib for check

## 4.13.48 (2024-12-25)

### Features

* Added load npm

## 4.13.47 (2024-12-25)

### Features

* Fixed category selection for test

## 4.13.46 (2024-12-23)

### Features

* Improved error handling:
  - Different exit codes for package errors / grok script errors
  - Graceful error handling when testing non-existing packages

## 4.13.45 (2024-12-12)

### Features

* Test all workers count variable fixes

## 4.13.44 (2024-12-06)

### Features

* Test fixed infinite testing 

## 4.13.43 (2024-12-02)

### Features

* Publish visual updates

## 4.13.42 (2024-11-29)

### Features

* Publish refresh orientates on debug versions of packages

## 4.13.41 (2024-11-29)

### Features

* Added ability to run auto tests by core variable

## 4.13.39 (2024-11-26)

### Features

* Grok publish added ability to link packages
* Grok test log for multiple workers fix

## 4.13.38 (2024-11-15)

### Features

* Minor fixes

## 4.13.37 (2024-11-15)

### Features

* Output fixes

## 4.13.36 (2024-11-15)

### Features

* Grok test reopens on each failed test

## 4.13.35 (2024-11-4)

### Features

* Grok help fixes

## 4.13.34 (2024-11-1)

### Features

* Config bug fixes

## 4.13.33 (2024-10-31)

### Features

* Added Changelog to package template
* Added ability to reload page on Execution Timeout

## 4.13.32 (2024-10-16)

### Features

* Replaced publishAll by publuish variables 
* Added global.d.ts file

## 4.13.31 (2024-10-15)

### Features

* Grok test without repository fixes


## 4.13.30 (2024-10-11)

### Features

* Grok Check fixes

## 4.13.29 (2024-10-09)

### Features

* Help fixes


## 4.13.28 (2024-10-09)

### Features

* Grok Check update(added checks on caching)

## 4.13.27 (2024-10-02)

### Features

* Publish all added

## 4.13.26 (2024-10-02)

### Features

* Create command fixes

## 4.13.25 (2024-09-30)

### Features

* Test command package publishing fix 

## 4.13.24 (2024-09-26)

### Features

* Test command fixes update

## 4.13.23 (2024-09-23)

### Features

* Test all bug fixes

## 4.13.22 (2024-09-20)

### Features

* Implemented mechanism which allows to test few packages by one function

## 4.13.21 (2024-09-17)

### Features

* Grok check updated for release candidate versions

## 4.13.20 (2024-08-28)

### Features

* Grok link command fixes 

## 4.13.19 (2024-08-27)

### Features

* Package ts data insertation fixes

## 4.13.18 (2024-08-27)

### Features

* Creation fixes

## 4.13.17 (2024-08-27)

### Features

* Added decorators for viewer functions

## 4.13.16 (2024-08-23)

### Features

* Grok link fixes

## 4.13.15 (2024-08-21)

### Features

* Grok link updated to link all dependent libs 

## 4.13.14 (2024-08-09)

### Features

* Timeout update for test runs up to 1 hour

## 4.13.13 (2024-08-09)

### Features

* Added stress test flag

## 4.13.12 (2024-08-05)

### Features

* Puppeteer-screen-recorder version fixed to 3.0.3

## 4.13.11 (2024-07-30)

### Features

* Added commit label to context panel for 'manually' published versions

## 4.13.10 (2024-07-29)

### Features

* Added ability to select test by args variable 

## 4.13.9 (2024-07-29)

### Features

* Added ability to select category by args variable 

## 4.13.8 (2024-07-18)

### Bug Fixes

* Fixed path chacks for npmignore

## 4.13.6 (2024-07-15)

### Features

* Splitted warnings and errors in grok check command

## 4.13.4 (2024-06-25)

### Bug Fixes

* Check source maps fixes 

## 4.13.3 (2024-06-21)

### Features

* Added check on source maps in grock check command

## 4.13.1 (2024-06-04)

### Features

* Added argument `package`, which allows you to choose a package for the `test` command


## 4.12.13 (2023-08-17)

### Bug Fixes

* Datagrok API version check fix
* Check Changelog fix
* Sync --csv and --verbose flags

## 4.12.12 (2023-08-07)

### Features

* Check for datagrok-api dependency

### Bug Fixes

* Latest package version in CHANGELOG check fix

## 4.12.11 (2023-08-04)

### Features

* GROK-13643 Check improvements:
  * There is no beta property in package.json
  * No datagrok-tools in dependencies (or latest version)
  * Latest version from package.json is in CHANGELOG (warning)
  * Ignore CHANGELOG checks for service packages ("servicePackage" property in package.json)
  * Change supported h2 formats (1.7.9 (2023-07-24) and 1.7.9 (WIP))
  * For packages < 1.0.0 exit with exit code 0, and only show warnings. And for packages >= 1.0.0, exit with a non-zero code (only for check command)
  * If an invalid flag/command is specified, output the help and exit with exit code 1

## 4.12.10 (2023-08-01)

### Features

* Video recording enhancements

### Bug Fixes

* FuncSignatures check fix

## 4.12.7 (2023-07-28)

### Features

* Tools: Changelog h2 new format

## 4.12.4 (2023-07-24)

### Features

* GROK-13573 Tools: simplify output (add --verbose flag)
* GROK-13573 Tools: checks for changelog
