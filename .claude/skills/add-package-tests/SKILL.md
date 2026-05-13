---
name: add-package-tests
version: 0.1.0
description: |
  Wire the Datagrok test harness into a package so unit specs and
  function-level `//test:` assertions are discoverable by `grok test`,
  the platform's Test Manager, and the console `<Package>:test(...)`
  call. For plugin authors who need an automated check suite CI can
  run on every PR — basic smoke plus a few integration scenarios.
  Produces `src/package-test.ts`, one `src/tests/<topic>-tests.ts`
  per category, and optional `//test:` lines on exported functions.
  Use when asked to "set up an automated check suite for my plugin",
  "wire CI checks for my package", "where do package tests live so
  the platform finds them", or "auto-verify a function with an inline
  assertion".
triggers:
  - automated check suite for plugin
  - wire ci checks for package
  - where do package tests live
  - inline assertion on a function
  - run my plugin through test manager
  - smoke test a datagrok package
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# add-package-tests

## When to use

Your Datagrok package has no automated checks yet, or has only ad-hoc
scripts. You want a category/test scaffold the platform discovers,
a `grok test` command that builds + deploys + runs the suite, and
optionally per-function `//test:` assertions wired into the same run.

## Prerequisites

- A package directory (fresh from `grok create`, or existing without tests).
- A server alias in `~/.grok/config.yaml` if you plan to use
  `grok test --host <alias>` (set once via `grok config add`).
- Familiarity with the canonical test-library exports — `category`,
  `test`, `expect`, `expectFloat`, `expectTable`, `before`, `after`,
  `delay`, `awaitCheck` (`DG-FACT-284`, `DG-FACT-286`).

## Steps

1. **Scaffold the test harness.** Pick ONE path (`DG-FACT-285`).
   ```bash
   # (a) NEW package — harness included from day one
   grok create MyPackage --test

   # (b) EXISTING package — add harness in place
   cd MyPackage
   grok add tests
   ```
   Expected: `src/package-test.ts` exists, `src/tests/test-examples.ts`
   exists, `package.json` has `"test": "grok test"`, and
   `webpack.config.js` has a `test:` entry. Shape matches
   `tools/package-template/src/package-test.ts:1-21`.

2. **Confirm the canonical import path** is `@datagrok-libraries/test/src/test`, not the legacy `utils/src/test` duplicate (see `DG-FACT-284`). If legacy, rewrite imports and ensure the dep is present:
   ```bash
   npm i -S @datagrok-libraries/test     # only if missing
   npm install
   ```
   Expected: `grep -r "@datagrok-libraries/utils/src/test" src/`
   returns nothing.

3. **Add a test file under `src/tests/`** — one file per category (`DG-FACT-289`). `expect` is STRICT `!==` with default `expected = true` (`DG-FACT-287`).
   ```typescript
   // src/tests/examples-tests.ts
   import {category, expect, test} from '@datagrok-libraries/test/src/test';

   category('Examples', () => {
     test('Success', async () => { expect(1, 1); });
     test('Float window', async () => { expect(10 < 12.5 && 12.5 < 20, true); });
     test('Skipped — tracked', async () => {
       expect(1, 11);
     }, {skipReason: 'GROK-99999'});       // DG-FACT-286, DG-FACT-288
   });
   ```

4. **Side-effect-import the file from `src/package-test.ts`** (`DG-FACT-289`, `DG-FACT-290`) — without this, `runTests` finds zero tests for the category.
   ```typescript
   // src/package-test.ts — ADD near the other imports
   import './tests/examples-tests';
   ```
   Keep the template shape: `_package`, `tests`, plus both `//name: test` and `//name: initAutoTests` entries (`DG-FACT-290`).

5. **(Optional) Annotate functions with `//test:` for auto-tests** — boolean Grok-script expression per line, inline sub-options `skip:`, `wait:<ms>`, `cat:<category>`, `timeout:<ms>` (`DG-FACT-291`, `DG-FACT-292`). Parameter-type matrix in `DG-FACT-293`.
   ```typescript
   //name: square
   //input: int x
   //output: int y
   //test: square(1) == 1
   //test: square(2) == 4
   //test: square(3) == 9, cat: Math
   export function square(x: number): number { return x ** 2; }
   ```

6. **Build, publish, and run.**
   ```bash
   grok publish <host> --release       # or omit --release for dev
   grok test --host <host>             # builds, publishes in debug, runs
   ```
   Useful flags: `--skip-build`, `--skip-publish`, `--csv`
   (`test-report.csv` in the package folder), `--gui` (non-headless
   browser for debugging).
   Expected: a summary like `Passed N, Failed 0, Skipped 1` on stdout;
   with `--csv`, a `test-report.csv` next to `package.json`.

## Common failure modes

- **`Cannot find module '@datagrok-libraries/test/src/test'`** — legacy `utils/src/test` path or missing dep. Fix per `DG-FACT-284`.
- **Tests don't appear in Test Manager / `runTests` returns 0** — missing side-effect import in `src/package-test.ts` (`DG-FACT-289`).
- **`expect(someBool)` fails for a truthy value** — `expect` is strict `!==` with default `true` (`DG-FACT-287`). Use `expect(!!val, true)`.
- **`//test:` line never runs** — `src/package-test.ts` is missing `//name: initAutoTests` (`DG-FACT-290`).
- **Test killed at 30 s** — default timeout (`DG-FACT-286`). Bump via `{timeout: 90000}` per-test or per-category.
- **`grok test --host <alias>` errors `unknown host`** — add via `grok config add` or omit `--host`.

## See also

- Source: `help/develop/how-to/tests/add-package-tests.md`;
  `help/develop/how-to/tests/test-packages.md` (Test Manager, `grok
  test` flags, console runs, GitHub Actions retry).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-284` (canonical import path), `DG-FACT-285` (CLI
  scaffolds), `DG-FACT-286` (`TestOptions`, `STANDART_TIMEOUT`),
  `DG-FACT-287` (`expect` strict semantics), `DG-FACT-288`
  (`skipReason`), `DG-FACT-289` (per-file category + side-effect
  import), `DG-FACT-290` (`package-test.ts` entry contract),
  `DG-FACT-291`–`DG-FACT-293` (`//test:` expression, inline
  sub-options, parameter-type matrix).
- Reference: `packages/Chem/src/package-test.ts:1-50` (30+
  side-effect imports); `packages/Chem/src/tests/chemprop-tests.ts:11-40`
  (`before`, `test(..., {timeout, skipReason})`);
  `packages/ApiTests/src/package-test.ts:1-40` (large-package layout).
- Related skills: `publish-packages` (prerequisite — `grok publish`
  is invoked by `grok test`).
