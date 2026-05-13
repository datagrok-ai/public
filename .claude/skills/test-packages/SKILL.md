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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

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
   The runner builds, publishes, and executes every registered `category(...)` (see DG-FACT-356). Output: stdout summary + `test-report.csv` next to `package.json`.

2. **Iterate fast — re-run against the already-deployed build.** Combine `--skip-build --skip-publish` with `--category` to re-run in seconds (see DG-FACT-356).
   ```bash
   grok test --host <alias> --skip-build --skip-publish \
     --category "<CategoryName>"
   ```

3. **Step through a single failing test interactively.** `--debug` requires `--gui`; otherwise the breakpoint is silently dropped (see DG-FACT-356).
   ```bash
   grok test --host <alias> --gui --debug \
     --category "<CategoryName>"
   ```

4. **Record the run for later inspection (flaky test triage).**
   ```bash
   grok test --host <alias> --record --csv --verbose
   ```

5. **Run a single test from the Datagrok UI (Test Manager).** **Top menu → Tools → Dev → Test Manager** (see DG-FACT-357). Expand package → category → test; right-click → `Run`, or double-click. Deep-link:
   ```
   <host>/apps/DevTools/TestManager/<PackageName>/<CategoryName>/<TestName>
   ```
   Append `?run=true` to auto-execute on load.

6. **Run from the platform console.** Press <kbd>~</kbd> (or **Windows → Console**), then (see DG-FACT-358):
   ```
   <PackageName>:test(category="<CategoryName>")
   <PackageName>:test(category="<CategoryName>", test="<TestName>")
   ```

7. **Re-run the GitHub Actions suite for a package on demand.** Open `https://github.com/datagrok-ai/public/actions/workflows/packages.yml` → **Run workflow** → enter packages (space-separated) and a target branch (see DG-FACT-359).

8. **Diagnose a CI failure from its artifact zip.** The zip contains exactly three files (see DG-FACT-359):
   ```text
   test-console-output.log   # raw console + stack traces
   test-record.mp4           # full video of the headless run
   test-report.csv           # per-test pass/fail/duration/error
   ```
   Triage order: `test-report.csv` (find failing row) → search name in `test-console-output.log` (stack trace) → scrub `test-record.mp4` to that elapsed time for visual/timing failures.

## Common failure modes

- **`--debug` does nothing.** Add `--gui` (see DG-FACT-356).
- **Async viewer test times out only in CI.** Default `awaitCheck` wait is 500ms — pass an explicit longer `wait` (see DG-FACT-352).
- **`grok test --host <alias>` errors `unknown host`.** Run `grok config add`, or omit `--host`.
- **Test Manager shows zero packages.** Either `DevTools` not installed (see DG-FACT-357), or `src/package-test.ts` missing `//name: test` / not deployed in debug mode.
- **Console call returns nothing.** Category file isn't side-effect-imported from `src/package-test.ts` — `grep "import './tests/" src/package-test.ts`.
- **CI-only `expectExceptionAsync` failure** (`An exception is expected but not thrown`). Likely a race — wrap with `awaitCheck` before the call (see DG-FACT-354).

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
