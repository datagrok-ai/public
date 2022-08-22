import { after, before, category, expect, expectFloat, test } from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { _testDimensionalityReducer } from './dimensionality-reduce-utils';
import { createTableView, readDataframe } from './utils';
import { _package } from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import { _testActivityCliffsOpen } from './activity-cliffs-utils';


category('activityCliffs', async () => {


  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('activityCliffsOpen', async () => {
    await _testActivityCliffsOpen('activity_cliffs.csv', 92);
  });

  test('activityCliffsWithEmptyRows', async () => {
    await _testActivityCliffsOpen('tests/activity_cliffs_empty_rows.csv', 91);
  });

  after(async () => {
  });

});