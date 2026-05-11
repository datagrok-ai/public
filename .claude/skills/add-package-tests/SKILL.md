---
name: add-package-tests
description: Add the Datagrok test harness to a package so unit tests and function-level `//test:` assertions run via `grok test` and Test Manager
harness-authored: true
---

# add-package-tests

## When to use

Your Datagrok package needs unit tests — either a fresh scaffold with the
test harness wired in, or an existing package that has none. Triggers:
"set up tests for my package", "wire up `grok test`", "add `category/test`
specs", "run my package in Test Manager", "annotate a function with a
`//test:` assertion so it's auto-tested."

## Prerequisites

- `datagrok-tools` installed globally (`npm i -g datagrok-tools`) — provides
  the `grok` CLI used below.
- A package scaffold (either fresh from `grok create`, or an existing
  package directory).
- A Datagrok server alias in `~/.grok/config.yaml` (`grok config add`); used
  by `grok test --host <alias>`.
- Familiarity with the canonical test-library exports (`category`, `test`,
  `expect`, `expectFloat`, `expectTable`, `before`, `after`, `delay`,
  `awaitCheck` — see `DG-FACT-284`, `DG-FACT-286`).

## Steps

1. **Scaffold the test harness.** Pick ONE path (`DG-FACT-285`):
   ```bash
   # (a) NEW package — harness included from day one
   grok create MyPackage --test

   # (b) EXISTING package — add harness in place
   cd MyPackage
   grok add tests
   ```
   Expected: `src/package-test.ts` exists; `src/tests/test-examples.ts`
   exists; `package.json` has `"test": "grok test"`; `webpack.config.js`
   has a `test:` entry. Compare `tools/package-template/src/package-test.ts`.

2. **Fix the DRIFT introduced by `grok add tests`** (skip on path (a)).
   `grok add tests` writes a `package-test.ts` that imports from
   `@datagrok-libraries/test/src/test` but only adds
   `@datagrok-libraries/utils` to `dependencies`
   (`DG-FACT-DRIFT-APT-002`). Add the missing dep and normalize the
   example file to the canonical path (`DG-FACT-DRIFT-APT-001`,
   `DG-FACT-DRIFT-APT-003`):
   ```bash
   npm i -S @datagrok-libraries/test
   ```
   Then edit `src/tests/test-examples.ts` so its first line reads
   `import {category, expect, test} from '@datagrok-libraries/test/src/test';`
   (NOT `…/utils/src/test` as the upstream article and entity template
   still teach).

3. **Install deps.**
   ```bash
   npm install
   ```
   Expected: clean install; `node_modules/@datagrok-libraries/test/src/test.ts`
   resolves (the file exporting `category`, `test`, `expect`, … —
   `DG-FACT-284`).

4. **Add a test file** under `src/tests/`, one file per category
   (`DG-FACT-289`).
   ```typescript
   // src/tests/examples-tests.ts
   import {category, expect, test} from '@datagrok-libraries/test/src/test';

   category('Examples', () => {
     test('Success', async () => {
       expect(1, 1);                            // strict !== compare (DG-FACT-287)
     });

     test('Float window', async () => {
       expect(10 < 12.5 && 12.5 < 20, true);    // expect default is `true`
     });

     test('Skipped — tracked', async () => {
       expect(1, 11);
     }, {skipReason: 'GROK-99999'});            // DG-FACT-286, DG-FACT-288
   });
   ```
   Expected: file compiles; `category(...)` runs at module load and
   registers tests in the library's shared `tests` registry
   (`DG-FACT-289`).

5. **Side-effect import the file from `src/package-test.ts`**
   (`DG-FACT-289`, `DG-FACT-290`). Without this line, `runTests` finds
   zero tests for the category.
   ```typescript
   // src/package-test.ts — ADD this line near the other imports
   import './tests/examples-tests';
   ```
   Expected: `src/package-test.ts` retains its template shape — exports
   `_package = new DG.Package()` and `tests`; declares both
   `//name: test` (calls `runTests(...)`) and `//name: initAutoTests`
   (calls `initAutoTests(_package, _package.getModule('package-test.js'))`).
   Compare `packages/Chem/src/package-test.ts:1-40`.

6. **(Optional) Annotate functions with `//test:` for auto-tests.**
   One line per assertion; the expression is a Grok-script boolean
   (`DG-FACT-291`). Multiple lines create numbered tests (`square 1`,
   `square 2`, …). Inline sub-options after the expression
   (`DG-FACT-292`): `skip:`, `wait:<ms>`, `cat:<category>`,
   `timeout:<ms>`.
   ```typescript
   //name: square
   //input: int x
   //output: int y
   //test: square(1) == 1
   //test: square(2) == 4
   //test: square(3) == 9, cat: Math
   export function square(x: number): number { return x ** 2; }
   ```
   Parameter-type support (`DG-FACT-293`): `int`, `double`, `bool`,
   `string` (use `""`), `datetime` (`Date(YYYY, MM, DD)` /
   `DateDiff(DateTime(...), "<iso>")`), `map`. For `dataframe`,
   `column_list`, `column`, `file`, `blob` — wrap via a helper call
   like `f(getDataframe())`.

7. **Build, publish, and run.** From the package root, with a configured
   host alias:
   ```bash
   grok publish <host> --release       # or omit --release for dev
   grok test --host <host>             # builds, publishes in debug, runs
   ```
   Useful flags: `--skip-build`, `--skip-publish`, `--csv` (report file
   in package folder), `--gui` (non-headless browser for debugging).
   Expected: a summary like `Passed N, Failed 0, Skipped 1` on stdout
   and (with `--csv`) a `test-report.csv` next to `package.json`.

## Common failure modes

- **`Cannot find module '@datagrok-libraries/test/src/test'`.** Hit the
  `grok add tests` drift — the command wrote the import but skipped the
  dep (`DG-FACT-DRIFT-APT-002`). Fix: `npm i -S @datagrok-libraries/test`.
- **Tests don't appear in Test Manager / `runTests` returns 0 for a
  category.** The test file isn't imported from `src/package-test.ts`,
  so its `category(...)` never executed (`DG-FACT-289`). Add the
  side-effect `import './tests/<file>';` line.
- **`expect(someBool)` fails for a truthy value.** `expect` is STRICT
  (`!==`), not JS-truthy, and defaults `expected = true`
  (`DG-FACT-287`). `expect(1)` fails because `1 !== true`. Use
  `expect(!!val, true)` or `expect(val, true)` when `val` is already a
  boolean.
- **`//test: f(123) == 246` line never runs.** No
  `initAutoTests` codegen — `src/package-test.ts` is missing the
  `//name: initAutoTests` annotation (`DG-FACT-290`). Restore it from
  `tools/package-template/src/package-test.ts`.
- **Test hangs and is killed at 30 s.** `STANDART_TIMEOUT` (sic) is
  30000 ms (`DG-FACT-286`). Bump per-test via the third arg:
  `test('long', async () => {...}, {timeout: 90000});` or set
  `category(..., () => {...}, {timeout: 90000})`.
- **`grok test --host <alias>` errors `unknown host`.** The alias
  isn't in `~/.grok/config.yaml`. Run `grok config add` once, or omit
  `--host` to use the default.

## Verification

- `grok test --host <host>` exits 0 and prints a summary with all the
  tests in your new category counted.
- In the platform: **Tools → Dev → Test manager** lists your package;
  expanding it shows your category with its tests. Running them inline
  shows the same pass/fail/skip counts as the CLI.
- Open the platform Console (`~`) and run
  `<PackageName>:test(category="Examples")` — the same results come
  back as a DataFrame.

## See also

- Source articles:
  - `help/develop/how-to/tests/add-package-tests.md` (mirror in
    `docs/_internal/articles-mirror/how-to/tests/add-package-tests.md`)
  - `help/develop/how-to/tests/test-packages.md` (mirror in
    `docs/_internal/articles-mirror/how-to/tests/test-packages.md`) —
    Test Manager, `grok test` flags, console-based runs, GitHub
    Actions retry.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-284` (canonical import path) … `DG-FACT-293` (`//test:`
  parameter-type matrix); drifts
  `DG-FACT-DRIFT-APT-001` (article uses `utils/`, code uses `test/`),
  `DG-FACT-DRIFT-APT-002` (`grok add tests` injects the wrong dep),
  `DG-FACT-DRIFT-APT-003` (`entity-template/test.ts` mismatched with
  orchestrator).
- Reference packages:
  - `packages/Chem/src/package-test.ts:1-40` — 30+ side-effect imports
    binding category files to the entry point.
  - `packages/Chem/src/tests/chemprop-tests.ts:4-40` — `category` with
    `before`, `test(..., {timeout, skipReason: 'GROK-19084'})`.
  - `packages/Bio/src/tests/msa-tests.ts:85-95` — free-form
    `skipReason: 'Fails in docker'`.
  - `packages/ApiTests/src/package-test.ts:1-30` — large-package
    side-effect-import layout.
- Related skills:
  - `publish-packages` (prerequisite — `grok publish` is invoked by
    `grok test` under the hood).
