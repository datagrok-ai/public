---
name: test-packages
version: 0.1.0
description: |
  Run an existing Datagrok package test suite — locally before pushing,
  in the platform UI, from the console, or via CI workflow dispatch —
  and diagnose failures using the runner's CSV/log/video artifacts.
  For plugin authors validating a change (regression check, flaky CI
  triage, interactive single-test debug). Produces a pass/fail summary,
  an optional `test-report.csv`, and a clear next step when a test
  fails.
  Use when asked to "regression-check my plugin before pushing",
  "diagnose a flaky CI failure on an async viewer", or "get a CSV
  pass/fail report I can diff against the previous run".
triggers:
  - regression check before pushing
  - diagnose a flaky ci failure
  - pass fail csv report for plugin
  - run my plugin suite locally
  - step through one test interactively
  - inspect github actions test artifacts
allowed-tools:
  - Read
  - Bash
harness-authored: true
---

# test-packages

## When to use

Your Datagrok package already has tests (scaffolded via `add-package-tests`)
and you need to *run* them — before a push, inside the platform UI, or
to triage a CI failure. Use this skill to pick the right runner for
the situation and to read the artifacts when a test goes red.

## Prerequisites

- A package with a working test harness — `src/package-test.ts` plus
  one or more `src/tests/*-tests.ts` files. If absent, run the
  `add-package-tests` skill first.
- `datagrok-tools` installed (`npm i -g datagrok-tools`) for the
  `grok test` CLI.
- A server alias in `~/.grok/config.yaml` for the host you want to
  run against (created once with `grok config add`).
- For the Test Manager and console runners: a Datagrok instance with
  the `DevTools` package installed (`DG-FACT-357`).

## Steps

1. **Run the suite locally against your dev server.**
   ```bash
   cd <package-dir>
   grok test --host <alias> --csv
   ```
   The runner builds the package, publishes it to `<alias>` in debug
   mode, then executes every `category(...)` registered through
   `src/package-test.ts` (`DG-FACT-356`).
   Expected: stdout summary `Passed N, Failed M, Skipped K`; a
   `test-report.csv` next to `package.json` with one row per test
   (name, category, pass/fail/skip, duration, error).

2. **Iterate fast — re-run against the already-deployed build.**
   When you're debugging logic and the package is already published,
   skip the slow steps:
   ```bash
   grok test --host <alias> --skip-build --skip-publish \
     --category "<CategoryName>"
   ```
   `--category` narrows the run to one `category(...)` block; combine
   with `--skip-build --skip-publish` to re-run in seconds (`DG-FACT-356`).
   Expected: only the named category is exercised; build + publish
   logs do not appear.

3. **Step through a single failing test interactively.**
   `--debug` is only honoured in non-headless mode — pair it with
   `--gui`, otherwise the breakpoint is silently dropped (`DG-FACT-356`).
   ```bash
   grok test --host <alias> --gui --debug \
     --category "<CategoryName>"
   ```
   Expected: a Chromium window opens; execution pauses at a
   breakpoint before each test; DevTools is attachable for live
   stepping.

4. **Record the run for later inspection (flaky test triage).**
   ```bash
   grok test --host <alias> --record --csv --verbose
   ```
   `--record` produces an mp4 of the headless browser; `--verbose`
   adds per-test passed/skipped detail to stdout (`DG-FACT-356`).
   Useful locally when you suspect the same visual/timing issue
   you saw in CI.

5. **Run a single test from the Datagrok UI (Test Manager).**
   In the platform: **Top menu → Tools → Dev → Test Manager**
   (`DG-FACT-357`). The tree lists every package whose
   `src/package-test.ts` is loaded; expand → category → test.
   - Right-click → `Run`, or select + `Enter`, or double-click a leaf.
   - `Run all` on the ribbon executes every package.

   Deep-link a test directly from the URL bar:
   ```
   <host>/apps/DevTools/TestManager/<PackageName>/<CategoryName>/<TestName>
   ```
   Append `?run=true` to auto-execute on load (`DG-FACT-357`).
   Expected: a pass/fail icon appears next to the selected node;
   results land in the context panel and a Workspace-addable grid.

6. **Run from the platform console.**
   Press <kbd>~</kbd> (or **Windows → Console**), then:
   ```
   <PackageName>:test(category="<CategoryName>")
   <PackageName>:test(category="<CategoryName>", test="<TestName>")
   ```
   This invokes the same `//name: test` function the runner uses
   (`DG-FACT-358`). The console prints a result DataFrame.

7. **Re-run the GitHub Actions suite for a package on demand.**
   Open `https://github.com/datagrok-ai/public/actions/workflows/packages.yml`
   → **Run workflow** → enter packages space-separated
   (e.g. `Chem ApiTests`) and a target branch. Publish-to-NPM only
   fires on `master` (`DG-FACT-359`).

8. **Diagnose a CI failure from its artifact zip.**
   On the failed Action's **Summary** page, download the artifact
   from the **Artifacts** pane. The zip contains exactly three
   files (`DG-FACT-359`):
   ```text
   test-console-output.log   # raw console + stack traces
   test-record.mp4           # full video of the headless run
   test-report.csv           # per-test pass/fail/duration/error
   ```
   Triage order: open `test-report.csv` → find the failing row →
   search its name in `test-console-output.log` for the stack trace
   → if the failure is timing-dependent or DOM-related, scrub
   `test-record.mp4` to the same elapsed time.

## Common failure modes

- **`--debug` does nothing — execution runs through.** `--debug`
  requires `--gui`; without it the breakpoint is silently dropped
  (`DG-FACT-356`, tools/README.md:141). Add `--gui`.
- **Async viewer test times out only in CI.** Default `awaitCheck`
  wait is `500` ms — far too short for headless rendering
  (`DG-FACT-352`). Pass an explicit wait
  (`await awaitCheck(() => …, 'msg', 30000)`); confirm by scrubbing
  `test-record.mp4` to the failure point.
- **`grok test --host <alias>` errors `unknown host`.** `<alias>`
  isn't in `~/.grok/config.yaml`. Run `grok config add`, or omit
  `--host` to use the configured default.
- **Test Manager shows zero packages or none of mine.** Either
  `DevTools` isn't installed on that server (`DG-FACT-357`), or your
  package's `src/package-test.ts` doesn't declare `//name: test` /
  isn't deployed in debug mode — `grok test` won't find it either.
- **Console call `<Package>:test(category="X")` returns nothing.** The
  category file isn't side-effect-imported from `src/package-test.ts`
  (see `add-package-tests` for the import contract). Confirm with
  `grep "import './tests/" src/package-test.ts`.
- **CI-only `expectExceptionAsync` failure.** Reads `An exception is
  expected but not thrown` — the action under test resolved instead of
  rejected (`DG-FACT-354`). Likely a race: the rejection path didn't
  fire before the next tick. Wrap with `awaitCheck` before the call.

## Verification

- `grok test --host <alias> --csv` exits `0` and writes a
  `test-report.csv` whose row count matches the stdout summary.
- The same test reached three ways gives the same verdict: the
  CLI, **Tools → Dev → Test Manager**, and the console
  `<Package>:test(category="…", test="…")` form.
- For a CI failure, the failing test in `test-report.csv` is also
  the failing node in `test-console-output.log` (matched by name)
  and its failure timestamp aligns with the visible problem in
  `test-record.mp4`.

## See also

- Source articles:
  - `help/develop/how-to/tests/test-packages.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-352` (`awaitCheck` defaults), `DG-FACT-353` (`delay`
  semantics), `DG-FACT-354` (`expectExceptionAsync` failure
  messages), `DG-FACT-355` (`testViewer` four-phase lifecycle),
  `DG-FACT-356` (`grok test` CLI flag matrix), `DG-FACT-357`
  (Test Manager UI + deep-link URL format), `DG-FACT-358` (platform
  console syntax), `DG-FACT-359` (GitHub Actions artifact contents).
- Related skills: `add-package-tests` (prerequisite — scaffolds
  the harness this skill runs); `publish-packages` (`grok test`
  invokes `grok publish` internally unless `--skip-publish`).
