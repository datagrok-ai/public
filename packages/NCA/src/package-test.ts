import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {registerNcaRepositoryFactory, EmsNcaRepository} from './adapter';

// The live categories use the real store, as the app does (package.ts registers the same factory).
// On a shared server it refuses writes: a store-writing test missing from the skip list fails loudly.
const sharedStand = !isLoopbackServer(grok.dapi.root);
registerNcaRepositoryFactory(() =>
  sharedStand ? refuseStoreWrites(new EmsNcaRepository()) : new EmsNcaRepository());

// Each tests file registers its own `category(...)` blocks at module load.
import './tests/adapter-column-validators-tests';
import './tests/adapter-profile-index-tests';
import './tests/adapter-profile-inputs-tests';
import './tests/adapter-results-writer-tests';
import './tests/adapter-results-aggregator-tests';
import './tests/adapter-grouping-columns-tests';
import './tests/adapter-autodetect-tests';
import './tests/adapter-route-normalize-tests';
import './tests/adapter-event-driven-tests';
import './tests/adapter-mean-sd-profile-tests';
import './tests/dg-null-sentinel-tests';
import './tests/dg-null-detector-tests';
import './tests/bioavailability-tests';
import './tests/toxicokinetics-tests';
import './tests/dose-proportionality-tests';
import './tests/sparse-nca-tests';
import './tests/sparse-strata-tests';
import './tests/sparse-iv-anchor-tests';
import './tests/adapter-na-normalize-tests';
import './tests/adapter-csv-float64-tests';
import './tests/adapter-unit-conversion-tests';
import './tests/controller-pipeline-tests';
import './tests/cl-plausibility-tests';
import './tests/conc-misscale-guard-tests';
import './tests/sex-mapping-tests';
import './tests/sex-stratifier-tests';
import './tests/reference-suite-tests';
import './tests/reference-nonmem-tests';
import './tests/results-view-presentation-tests';
import './tests/results-grid-platform-tests';
import './tests/dock-layout-tests';
import './tests/results-grid-columns-tests';
import './tests/compute-notices-tests';
import './tests/study-schema-tests';
import './tests/studies-router-tests';
import './tests/recipe-round-trip-tests';
import './tests/summary-cache-tests';
import './tests/study-grid-tests';
import './tests/multi-column-group-by-tests';
import './tests/passthrough-columns-tests';
import './tests/context-panel-forms-toggle-tests';
import './tests/hub-onboarding-tests';
import './tests/lambdaz-manual-points-tests';
import './tests/nca-controller-override-tests';
import './tests/lambdaz-service-columns-tests';
import './tests/mapping-unit-dual-mode-tests';
import './tests/workflow-stages-adapter-tests';
import './tests/workflow-stages-exclusion-tests';
import './tests/workflow-stages-controller-tests';
import './tests/workflow-stages-storage-tests';
import './tests/workflow-stages-dock-probe-tests';
import './tests/workflow-stages-ui-tests';
import './tests/workflow-stages-live-tests';
import './tests/pipeline-functions-tests';
import './tests/pipeline-save-tests';
import './tests/ems-goldens-tests';
import './tests/ems-platform-tests';
import './tests/ems-nca-rows-tests';
import './tests/nca-repository-contract-tests';
import './tests/shared-stand-guard-tests';
import './tests/demo-studies-tests';
import './tests/occasion-tests';
import {skipStoreWritingTestsOnSharedStand, isLoopbackServer, refuseStoreWrites} from './tests/shared-stand-guard';

// After every category registered: the store-writing tests run only against a local / CI server.
skipStoreWritingTestsOnSharedStand(tests, grok.dapi.root);

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, verbose: true});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
