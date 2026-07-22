# test changelog

## 1.4.0 (WIP)

* Added the `node` option to `TestOptions`/`CategoryOptions` (test runs headless under the js-api Node runtime, no browser) and `nodeOnly`/`excludeNodeTests` to `TestExecutionOptions`; categories with no tests in the requested node/browser target are skipped entirely, including `before()`/`after()`
* Auto tests generated from query `-- test:` annotations are marked `node: true` — they evaluate server-side and, with the shared `OpenFile` fallback in the Node runtime, run headless (moves ~200 DBTests query tests off the browser)

## 1.3.1 (WIP)

* Compiled `src/playwright/` to CommonJS `.js` + `.d.ts` (via `src/playwright/tsconfig.build.json`) so consumers can `require()` the helpers from `node_modules` — Playwright does not transpile `.ts` inside `node_modules`, so 1.3.0 shipped without usable helper output

## 1.3.0 (2026-07-07)

* Added `src/playwright/` — shared Playwright-driver E2E helpers (spec-login, viewers, projects, openers, bio, chem, models-helpers, session), `global-setup`, and `base-config` for playwright-public and package-owned `playwright/` suites; `@playwright/test` declared as an optional peer

## 1.0.0 (2021-03-18)