# Testing

**Analysis Date:** 2026-05-11

## Framework

The package uses `@datagrok-libraries/test` — Datagrok's in-house TypeScript test framework. Tests are bundled with the package and run inside a Datagrok session via the `grok test` CLI.

```ts
import {category, test, expect} from '@datagrok-libraries/test/src/test';
```

Three core primitives:
- `category(name, body)` — groups tests
- `test(name, async body)` — declares a single test
- `expect(actual, expected)` — positional assertion (not `.toBe()` chain style)

Example from `src/tests/parsers.ts`:

```ts
category('Parsers: MaxQuant', () => {
  test('parses proteinGroups.txt', async () => {
    const text = makeTsv(...);
    const df = parseMaxQuantText(text);
    expect(df.rowCount, 27);
    expect(df.getTag('proteomics.source'), 'maxquant');
  });
});
```

## Entry Point and Discovery

`src/package-test.ts` is the test entry — it imports each test file via side-effect import, which triggers the `category()/test()` registrations:

```ts
import './tests/parsers';
import './tests/spectronaut-parser';
import './tests/generic-parser';
import './tests/analysis';
import './tests/qc-dashboard';
import './tests/enrichment';
import './tests/enrichment-visualization';

export const tests = (DG as any).Test?.tests ?? [];
```

The exported `tests` array is what `grok test` discovers and runs.

## Test File Layout

Test files live under `src/tests/`, one per domain area:

| File | Covers |
|------|--------|
| `src/tests/parsers.ts` | MaxQuant and Spectronaut parser behavior |
| `src/tests/spectronaut-parser.ts` | Spectronaut-specific edge cases (long-format pivot, pre-normalized detection) |
| `src/tests/generic-parser.ts` | Generic matrix parser; delimiter detection |
| `src/tests/analysis.ts` | Normalization, imputation, DE (client-side paths) |
| `src/tests/qc-dashboard.ts` | QC computations — MA, CV, missingness |
| `src/tests/enrichment.ts` | g:Profiler integration (mocked at the function level) |
| `src/tests/enrichment-visualization.ts` | Enrichment viewer factories and dot-plot construction |

## Test Fixtures

Tests build DataFrames inline using small helper factories defined at file scope (not `beforeEach`):

```ts
function makeTsv(headers: string[], rows: string[][]): string {
  return [headers.join('\t'), ...rows.map(r => r.join('\t'))].join('\n');
}

function makeTestDf(): DG.DataFrame {
  return DG.DataFrame.fromCsv(`Protein IDs,Sample1,Sample2
P00001,100,200
P00002,150,250`);
}
```

For larger / realistic fixtures, the package's demo files double as test data:

| File | Use |
|------|-----|
| `files/demo/proteinGroups.txt` | Small (27 proteins, 6 samples) — fast parser tests |
| `files/demo/cptac-spike-in.txt` | Large (1,569 proteins) — full-pipeline integration |
| `files/demo/spectronaut-hye-mix.tsv` | Spectronaut long-format (93 proteins, 8 runs) — DIA parser tests |

Tests load these via `_package.files.readAsText('demo/...')` when they need realistic data.

## Async Convention

**All tests are `async`** even when the body is synchronous. This is a framework requirement — `test(name, fn)` expects `fn: () => Promise<void>`.

```ts
test('extracts log2 columns', async () => {
  const cols = log2TransformColumns(df, ['Sample1', 'Sample2']);
  expect(cols.length, 2);
});
```

## Assertion Style

`expect()` takes two positional arguments: `expect(actual, expected)`. There is no `.toBe()` / `.toEqual()` chain style.

```ts
expect(df.rowCount, 27);                      // value equality
expect(df.getTag('proteomics.source'), 'maxquant');  // tag check
expect(df.columns.contains('log2FC'), true);  // boolean
```

Tag assertions are the **primary pattern** for verifying workflow state — the package's design (tag-based state) makes "did this step run?" trivially testable.

## Mocking

**No mocking framework.** Tests call real implementation functions directly. The R-dependent paths (`runLimmaDE`, `runDeqmsDE`, `vsnNormalize`) are tested only by triggering their JS-fallback branches — when R is not available in the test environment, the catch block fires and the client-side equivalent runs.

Consequence: the R output parsing logic (column renaming, type conversion in `src/analysis/differential-expression.ts:148-151`) is **not covered**. See CONCERNS.md "R-Dependent Paths Not Tested" for details and the risk argument.

## Running Tests

```bash
grok test                          # default server
grok test --host localhost         # local Datagrok instance
grok test --host dev               # dev server
grok test --gui                    # visible browser (debugging)
grok test --gui --debug            # with breakpoints
grok test --verbose                # detailed output
grok test --test "TestName"        # single test
grok test --category "CategoryName"  # one category
grok test --skip-build             # use existing dist/
grok test --skip-publish           # use already-published version
```

Common workflow during development:

```bash
npm install
grok publish local            # publish current build
grok test --host localhost    # run against it
```

The test runner spins up a headless Puppeteer browser connected to the target Datagrok server, loads the package's test bundle (`dist/package-test.js`), and executes each registered `test(...)`.

Convenience npm scripts:

```bash
npm test          # = grok test
npm run test-local  # = grok test --host localhost
npm run test-dev    # = grok test --host dev
```

## Coverage Gaps

See CONCERNS.md "Test Coverage Gaps" for the detailed audit. Highest-priority gaps:

- **PCA computation** (`src/analysis/pca.ts` — hand-rolled Jacobi eigendecomposition, zero tests)
- **R output parsing** (column renaming logic in DE — runs only with R available, untested without it)
- **Viewer rendering** (no test calls `createExpressionHeatmap`, `createVolcanoPlot`, or `createPcaPlot`)
- **External API JSON parsing** (UniProt and g:Profiler response → DataFrame logic untested)

## Test Output

`grok test` reports per-category pass/fail counts and detailed failure messages. CSV export available via `grok test --csv`. Test runs are recorded with `grok test --record` for replay debugging.

---

*Testing analysis: 2026-05-11*
