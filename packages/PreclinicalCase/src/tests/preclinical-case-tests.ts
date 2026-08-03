import {before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {createStudyWithConfig, readClinicalData, studies} from '../utils/app-utils';
import {CONFIGURATION_VIEW_NAME, MATRIX_TABLE_VIEW_NAME, MEASUREMENT_PROFILE_TABLE_VIEW_NAME,
  MICROSCOPIC_FINDINGS_TABLE_VIEW_NAME, OBSERVATION_TIMELINES_VIEW_NAME,
  SUMMARY_VIEW_NAME, VALIDATION_VIEW_NAME} from '../constants/view-names-constants';
import {calculateLBBaselineColumns, createAllMeasurementsDf} from '../data-preparation/data-preparation';
import { _package } from '../package-test';
import { funcs } from '../package-api';

const TEST_STUDY_ID = 'PC201708';
const VALIDATION_TIMEOUT_MS = 30000;

async function ensureStudyLoaded(): Promise<void> {
  if (studies[TEST_STUDY_ID]?.initCompleted)
    return;

  const files = await _package.files.list(`SEND/${TEST_STUDY_ID}`);
  const stubNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('PreclinicalCase tests');
  await createStudyWithConfig(files, stubNode, true);
  await readClinicalData(studies[TEST_STUDY_ID], files);
  studies[TEST_STUDY_ID].init();

  // process() is sync; validate() runs async — wait for it (or time out gracefully).
  await new Promise<void>((resolve) => {
    if (studies[TEST_STUDY_ID].validated) {
      resolve();
      return;
    }
    const sub = studies[TEST_STUDY_ID].validationCompleted.subscribe(() => {
      sub.unsubscribe();
      resolve();
    });
    setTimeout(() => {
      sub.unsubscribe();
      resolve();
    }, VALIDATION_TIMEOUT_MS);
  });
}

category('preclinicalCase app', () => {

  test('appOpensWithoutExceptions', async () => {
    const view = await funcs.preclinicalCaseApp();
    expect(view != null, true);
  });
 
});

category('preclinicalCase views and functions', () => {
  before(async () => {
    // Loaded into the test bundle so the data-preparation tests below can read
    // studies[TEST_STUDY_ID].domains directly. Views go through the registered
    // openPreclinicalCaseView function so they construct in the bound production bundle.
    await ensureStudyLoaded();
  });

  test('studyLoaded', async () => {
    const study = studies[TEST_STUDY_ID];
    expect(study != null, true, 'PC201708 study should be registered after load');
    expect(study.initCompleted, true);
    expect(study.domains.dm != null, true, 'dm domain should be present');
    expect(study.subjectsCount > 0, true, 'study should have at least one subject');
  });

  // Table views kick off async work in onTableViewAdded (dock, setOrder, hideValidationColumns).
  // The next test's addPreview tears the previous view down before that work finishes, which
  // triggers post-test errors against a disposed grid/dockManager. A short settle here lets
  // each view's deferred work drain before the framework moves on.
  const VIEW_SETTLE_MS = 2000;

  test('view: Summary', async () => {
    const view = await funcs.openPreclinicalCaseView(TEST_STUDY_ID, SUMMARY_VIEW_NAME);
    expect(view != null, true);
    await delay(VIEW_SETTLE_MS);
  });

  test('view: Validation', async () => {
    const view = await funcs.openPreclinicalCaseView(TEST_STUDY_ID, VALIDATION_VIEW_NAME);
    expect(view != null, true);
    await delay(VIEW_SETTLE_MS);
  });

  test('view: Matrix', async () => {
    const view = await funcs.openPreclinicalCaseView(TEST_STUDY_ID, MATRIX_TABLE_VIEW_NAME);
    expect(view != null, true);
    await delay(VIEW_SETTLE_MS);
  });

  test('view: Measurements', async () => {
    const view = await funcs.openPreclinicalCaseView(TEST_STUDY_ID, MEASUREMENT_PROFILE_TABLE_VIEW_NAME);
    expect(view != null, true);
    await delay(VIEW_SETTLE_MS);
  });

  test('view: Microscopic Findings', async () => {
    const view = await funcs.openPreclinicalCaseView(TEST_STUDY_ID, MICROSCOPIC_FINDINGS_TABLE_VIEW_NAME);
    expect(view != null, true);
    await delay(VIEW_SETTLE_MS);
  });

  test('view: Observation timelines', async () => {
    const view = await funcs.openPreclinicalCaseView(TEST_STUDY_ID, OBSERVATION_TIMELINES_VIEW_NAME);
    expect(view != null, true);
    await delay(VIEW_SETTLE_MS);
  });

  test('view: Configuration', async () => {
    const view = await funcs.openPreclinicalCaseView(TEST_STUDY_ID, CONFIGURATION_VIEW_NAME);
    expect(view != null, true);
    await delay(VIEW_SETTLE_MS);
  });

  test('calculateLBBaselineColumns adds baseline / change / pct columns', async () => {
    const lb = studies[TEST_STUDY_ID].domains.lb;
    expect(lb != null, true, 'lb domain should be loaded');
    calculateLBBaselineColumns(lb!);
    for (const colName of ['LB_BASELINE', 'LB_CHG', 'LB_PCT_CHG', 'MAX_POST_VALUE', 'MIN_PCT_CHG', 'MAX_PCT_CHG']) {
      expect(lb!.col(colName) != null, true, `expected column ${colName} after calculateLBBaselineColumns`);
    }
  });

  test('createAllMeasurementsDf builds a combined frame and excludes BG', async () => {
    const measDf = createAllMeasurementsDf(TEST_STUDY_ID);
    expect(measDf != null, true, 'createAllMeasurementsDf should return a DataFrame');
    expect(measDf!.rowCount > 0, true, 'combined measurements should have rows');

    const domainCol = measDf!.col('DOMAIN');
    expect(domainCol != null, true, 'combined measurements should have a DOMAIN column');
    if (domainCol) {
      const domainValues = new Set<string>(domainCol.toList().map((v: any) => String(v).toLowerCase()));
      expect(domainValues.has('bg'), false, 'BG (body weight gain) must be excluded from combined measurements');
    }

    for (const colName of ['USUBJID', 'test', 'result', 'visit_day']) {
      expect(measDf!.col(colName) != null, true, `expected combined-measurements column ${colName}`);
    }
  });
});
