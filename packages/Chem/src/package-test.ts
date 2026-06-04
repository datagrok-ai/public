import * as DG from 'datagrok-api/dg';
import {runTests, TestContext, tests, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';

import './tests/calculate';
import './tests/menu-tests-chem-space';
import './tests/menu-tests-cliffs';
import './tests/menu-tests-similarity-diversity';
import './tests/menu-tests-rgroups';
import './tests/menu-tests-script-based';

import './tests/col-panel-tests';
import './tests/cell-panel-tests';

import './tests/substructure-search-tests';
import './tests/rendering-tests';
import './tests/rendering-scatter-plot-tooltip-tests';
import './tests/sketcher-tests';

import './tests/detector-tests';
import './tests/api-based-tests';
import './tests/notation-converter-tests';
import './tests/screening-tools';
import './tests/pharmacophore-features-tests';

import './tests/save-as-sdf-tests';
import './tests/substructure-filter-tests';

import './tests/mol2-importer-tests';
import './tests/chemical-table-parsing';
import './tests/is-smarts-tests';
import './tests/fingerprints';
import './tests/scaffold-tree-tests';
import './tests/projects-tests';
import './tests/clone-layout-tests';
import './tests/mmpa-tests';
import './tests/chemprop-tests';
import './tests/vector-funcs-tests';
import './tests/synthon-search-tests';

import './tests/atom-index-mapper-tests';
import './tests/atom-picker-tests';
import './tests/reaction-enumeration-tests';
// import './tests/atom-picker-3d-hover-tests';
import './tests/atom-picker-escape-tests';
import './tests/renderer-hover-tests';

import './tests/viewers';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool stressTest {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext, stressTest?: boolean): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, stressTest});
  printTimingSummary(data);
  return DG.DataFrame.fromObjects(data)!;
}

// Temporary diagnostics: prints per-test durations and the grand total so a run's log ends
// with a sortable timing table. The suite hang floats between viewer tests run-to-run, so this
// shows which test ate the time regardless of where the wedge landed. Remove once diagnosed.
function printTimingSummary(data: {category?: string, name?: string, ms?: number, skipped?: boolean}[]): void {
  const log = (msg: string): void => console.log(`[TIMING-SUMMARY] ${msg}`);
  if (!Array.isArray(data) || data.length === 0) {
    log('no test results to summarize');
    return;
  }
  const run = data.filter((d) => !d.skipped && typeof d.ms === 'number');
  const totalMs = run.reduce((s, d) => s + (d.ms ?? 0), 0);
  const skipped = data.length - run.length;

  const byCategory = new Map<string, number>();
  for (const d of run)
    byCategory.set(d.category ?? '?', (byCategory.get(d.category ?? '?') ?? 0) + (d.ms ?? 0));

  log(`==== TEST TIMING SUMMARY: ${run.length} tests run, ${skipped} skipped, ` +
    `total ${(totalMs / 1000).toFixed(1)}s (${totalMs} ms) ====`);

  log('---- per category (slowest first) ----');
  for (const [cat, ms] of [...byCategory.entries()].sort((a, b) => b[1] - a[1]))
    log(`  ${(ms / 1000).toFixed(1)}s  ${cat}`);

  log('---- per test (slowest first) ----');
  for (const d of [...run].sort((a, b) => (b.ms ?? 0) - (a.ms ?? 0)))
    log(`  ${((d.ms ?? 0) / 1000).toFixed(1)}s  ${d.category} / ${d.name}`);

  log(`==== END TEST TIMING SUMMARY (total ${(totalMs / 1000).toFixed(1)}s) ====`);
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
