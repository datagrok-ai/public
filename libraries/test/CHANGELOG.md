# test changelog

## 1.5.0 (WIP)

* A benchmark run does not enter a category that has no benchmark test, so its `before()` and `after()` no longer run for nothing (in ApiTests' benchmark pass every category's hooks ran, for 13 tests)
* The garbage collection `grok test` forces in the browser runs after a test that was executed, at most once in 5 seconds, instead of after every test, executed or not: one collection of the platform's heap takes 0.1-0.3 s, which for a suite of small tests was more than the tests themselves (ApiTests: 1717 tests, ~10 of 19 minutes)
* Test results are sent to the server behind the tests, one at a time, and are awaited at the end of the run; a result the server did not take no longer fails the run
* When a category's `before()` fails, each of its tests reports the reason (`before() failed: <reason>`) instead of the bare `before() failed`
* Added `skipCategoryWithoutSpawner`, which skips a category, with the reason, when `grok.dapi.docker.getStatus()` reports Grok Spawner switched off or unreachable (a test's own skip reason is kept), and fails it when the status cannot be read or the platform has no `getStatus`
* [core#1050](https://github.com/datagrok-ai/core/issues/1050): Added `waitUntil`, `waitUntilAsync` and `settled(viewer)`; `testViewer` and its helpers and the demo auto-tests wait for viewers to settle instead of fixed delays
* The `BinaryImage` column of the detectors' test data holds a PNG instead of zero bytes, which PowerGrid's detector no longer treats as an image
* Added the TestTrack shared-page lanes (`shared-page`: worker-shared `test`/`localTest` with `resetShell`), the `known-open-bug`, `filter-panel`, `forms` and `balloons` helpers, and the viewer, property-grid, canvas and event-wait helpers in `viewers`
* Added `openDatagrok`, local/server lanes, the timing ledger, the `__readCsv` bridge and per-context token minting to `spec-login`, and `saveProjectViaApi`/`saveProjectViaUI`/`drainPendingDeletes` to `projects`, keeping the CI fixes and existing signatures
* `uiSettled` resolves at once under Node.js, so the Node runner no longer reports an `Unhandled exceptions` failure after every category

## 1.4.1 (2026-09-17)

* `grok test` runs a garbage collection after every test in the browser

## 1.4.0 (2026-07-23)

* Added the `node` option to `TestOptions`/`CategoryOptions` (test runs headless under the js-api Node runtime, no browser) and `nodeOnly`/`excludeNodeTests` to `TestExecutionOptions`; categories with no tests in the requested node/browser target are skipped entirely, including `before()`/`after()`
* Auto tests generated from query `-- test:` annotations are marked `node: true` — they evaluate server-side and, with the shared `OpenFile` fallback in the Node runtime, run headless (moves ~200 DBTests query tests off the browser)

## 1.3.1 (WIP)

* Compiled `src/playwright/` to CommonJS `.js` + `.d.ts` (via `src/playwright/tsconfig.build.json`) so consumers can `require()` the helpers from `node_modules` — Playwright does not transpile `.ts` inside `node_modules`, so 1.3.0 shipped without usable helper output

## 1.3.0 (2026-07-07)

* Added `src/playwright/` — shared Playwright-driver E2E helpers (spec-login, viewers, projects, openers, bio, chem, models-helpers, session), `global-setup`, and `base-config` for playwright-public and package-owned `playwright/` suites; `@playwright/test` declared as an optional peer

## 1.0.0 (2021-03-18)