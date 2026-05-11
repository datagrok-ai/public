---
name: test-package
description: Run a Datagrok package's tests via `grok test`, Test Manager, the platform console, or GitHub Actions — and write robust async + viewer assertions
harness-authored: true
---

# test-package

## When to use

Your package's test harness is already wired (see `add-package-tests`)
and you need to RUN the tests — locally via `grok test`, in the
platform via Test Manager / console, or in CI via the Packages
GitHub Action. Also covers the runtime utilities (`awaitCheck`,
`delay`, `expectExceptionAsync`, `testViewer`) and how to triage
failures from the CI artifact bundle.

## Prerequisites

- `datagrok-tools` installed (`npm i -g datagrok-tools`).
- Server alias in `~/.grok/config.yaml` (`grok config add`) for
  `grok test --host <alias>`.
- `src/package-test.ts` exports `tests`, declares the `//name: test`
  function, and side-effect-imports at least one category file. If
  not, run `add-package-tests` first.
- Imports come from `@datagrok-libraries/test/src/test`, NOT
  `…/utils/src/test` (`DG-FACT-DRIFT-APT-001`, `DG-FACT-DRIFT-TP-001`).

## Steps

1. **Smoke-run the whole suite locally** against a configured host
   alias (`DG-FACT-356`).
   ```bash
   cd <package-dir>
   grok test --host <alias>
   ```
   Expected: `grok` builds the package, publishes it in debug, runs
   every category, prints a summary like `Passed N, Failed 0, Skipped K`.
   Exit code is `0` iff no test failed.

2. **Iterate fast.** Once the deployed package is current, skip
   build+publish and re-run a single category (`DG-FACT-356`).
   ```bash
   grok test --host <alias> --skip-build --skip-publish \
     --category=Examples --csv
   ```
   Expected: `--csv` writes `test-report.csv` next to `package.json`.

3. **Debug visually.** `--gui` disables headless mode; `--debug`
   pauses at a breakpoint before each test (GUI-only, `DG-FACT-356`).
   ```bash
   grok test --host <alias> --skip-build --skip-publish \
     --category=Examples --gui --debug
   ```
   Expected: a Chromium window opens and pauses at each test's entry
   for DevTools inspection.

4. **Run via Test Manager.** Open `Tools | Dev | Test Manager`
   (from `DevTools`, `DG-FACT-357`). Right-click → Run, Enter, or
   double-click a leaf. Deep-link (auto-runs with `?run=true`):
   ```
   <host>/apps/DevTools/TestManager/<PackageName>/<CategoryName>/<TestName>?run=true
   ```
   Expected: progress icons flip to result icons; context panel
   shows per-test result; bottom progress bar tracks percent.

5. **Run from the platform console.** Press `~` (or `Windows |
   Console`) and call the package's `test` function (`DG-FACT-358`):
   ```text
   <PackageName>:test(category="<CategoryName>")
   <PackageName>:test(category="<CategoryName>", test="<TestName>")
   ```
   Expected: a DataFrame of results inline. Example —
   `ApiTests:test(category="Layouts", test="ViewLayout.toJson()")`.

6. **Write async-safe waits.** `awaitCheck` polls a predicate every
   `interval` ms (default `50`) until it returns true OR `wait` ms
   (default `500`) elapses (`DG-FACT-352`) — the 500 ms default is
   rarely enough; pass an explicit timeout.
   ```typescript
   import {awaitCheck, delay} from '@datagrok-libraries/test/src/test';

   const smiles = grok.data.demo.molecules(20);
   const v = grok.shell.addTableView(smiles);
   await awaitCheck(() => document.querySelector('canvas') !== null,
                    'cannot load table', 3000);
   await delay(1000);                  // fixed pause (DG-FACT-353)
   v.close();
   ```
   Expected: the test only proceeds once the canvas mounts; otherwise
   it fails with `Error('cannot load table')` after 3 s. Reach for
   `delay` ONLY when the wait has no checkable end state.

7. **Assert async exceptions.** `expectExceptionAsync` passes on any
   throw, or on a throw satisfying the optional `check` predicate
   (`DG-FACT-354`). No sync counterpart — wrap sync code in an async
   closure.
   ```typescript
   import {expectExceptionAsync} from '@datagrok-libraries/test/src/test';

   await expectExceptionAsync(
     () => grok.functions.call('nonExistingFunction'));
   ```
   Expected: passes when the call rejects; otherwise fails with
   `'An exception is expected but not thrown'`.

8. **Cover viewers with `testViewer`** — runs four phases (open/close;
   selection+filter+current-row; options+save-layout; load-layout)
   on your dataframe AND on a built-in categorical fallback when
   `arbitraryDfTest !== false` (default `true`); the viewer must
   survive both (`DG-FACT-355`). Input needs ≥3 rows.
   ```typescript
   import {testViewer, awaitCheck} from '@datagrok-libraries/test/src/test';

   const smiles = grok.data.demo.molecules(100);
   await testViewer('Chem Similarity Search', smiles,
                    {detectSemanticTypes: true});

   // For viewers with async rendering, gate on a DOM marker:
   await testViewer('Chem Diversity Search', smiles, {
     detectSemanticTypes: true,
     awaitViewer: async (v) => {
       await awaitCheck(() => !!v.root.querySelector('.chem-diversity-search'),
                        'Viewer not rendered', 1000);
     },
   });
   ```
   Expected: every phase passes silently; any phase that throws fails
   the test with the offending phase's error.

9. **Trigger the Packages workflow manually + triage CI** when a
   commit-driven run fails or a flaky run needs a retry
   (`DG-FACT-359`). In the GitHub UI: `Actions → Packages → Run
   workflow` → space-separated package list (e.g., `Demo Tutorials`)
   + target branch. NPM publish runs only on `master`. Each run's
   Summary → Artifacts panel ships a zip with three files:
   - `test-console-output.log` — stack traces and stdout.
   - `test-record.mp4` — headless-browser screen capture.
   - `test-report.csv` — one row per test (status, duration, error).

## Common failure modes

- **`awaitCheck` times out after 500 ms.** You omitted the third arg
  and got the default (`DG-FACT-352`). Pass an explicit timeout
  (1000–30000 typical) and a descriptive `error` string.
- **`testViewer` fails on the categorical fallback frame.** Phases
  repeat against a built-in categorical df (`DG-FACT-355`). Either
  widen the viewer's input handling or set `{arbitraryDfTest: false}`
  with a justifying comment.
- **`expectExceptionAsync` reports "exception is expected but not
  thrown".** The action resolved instead of rejecting (`DG-FACT-354`)
  — confirm the failing path is actually `await`ed and not caught.
- **`grok test --host <alias>` errors `unknown host`.** Alias missing
  from `~/.grok/config.yaml`. Run `grok config add` once or omit
  `--host`.
- **`--debug` does nothing.** Requires `--gui` (`DG-FACT-356`) —
  headless Chromium can't honor breakpoints.
- **Console call returns 0 rows.** Category file not side-effect-
  imported from `src/package-test.ts` (`DG-FACT-289` from
  `add-package-tests`). Add `import './tests/<file>';`.

## Verification

- `grok test --host <alias>` exits 0 and the printed summary lists
  every category the package declares.
- The same suite runs cleanly from Test Manager (`Tools | Dev |
  Test Manager`) — drilling into a category shows green icons for
  every leaf.
- `<PackageName>:test(category="<CategoryName>")` from the platform
  console returns the same pass/fail counts.
- For CI: the GitHub Actions Packages workflow run for the latest
  push is green, and the artifact zip's `test-report.csv` shows no
  rows with `Failed`.

## See also

- Source articles:
  - `help/develop/how-to/tests/test-packages.md` (mirrored at
    `docs/_internal/articles-mirror/how-to/tests/test-packages.md`)
  - `help/develop/how-to/tests/add-package-tests.md` — prerequisite
    harness wiring (mirrored alongside).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-352` (`awaitCheck`), `DG-FACT-353` (`delay`),
  `DG-FACT-354` (`expectExceptionAsync`), `DG-FACT-355` (`testViewer`
  four-phase lifecycle), `DG-FACT-356` (`grok test` flag matrix),
  `DG-FACT-357` (Test Manager URL form), `DG-FACT-358` (console
  syntax), `DG-FACT-359` (CI artifact triad), `DG-FACT-360` (ApiTests
  reference); drift `DG-FACT-DRIFT-TP-001` (stale `utils/src/test.ts`
  deep-links — use `@datagrok-libraries/test/src/test`).
- Reference packages:
  - `packages/ApiTests/src/package-test.ts:1-128` — 49 side-effect
    category imports + `testPlatform()` / `testPackages()` discovery.
  - `packages/Chem/src/tests/chemprop-tests.ts:4-40` — `category` +
    `test(..., {timeout, skipReason: 'GROK-19084'})`.
  - `packages/DevTools/src/package-testing.ts:66-150` — Test Manager
    URL parser + `?run=true` auto-run hook.
- Related skills: `add-package-tests` (prerequisite).
