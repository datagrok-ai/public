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

  let actCliffsTableView: DG.TableView;
  let actCliffsDf: DG.DataFrame;
  let actCliffsEmptyRowsTableView: DG.TableView;
  let actCliffsEmptyRowsDf: DG.DataFrame

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    actCliffsTableView = await createTableView('activity_cliffs.csv');
    actCliffsDf = await readDataframe('activity_cliffs.csv');
    actCliffsEmptyRowsTableView = await createTableView('tests/activity_cliffs_empty_rows.csv');
    actCliffsEmptyRowsDf = await readDataframe('tests/activity_cliffs_empty_rows.csv');
  });

  test('activityCliffsOpen', async () => {
    await _testActivityCliffsOpen(actCliffsDf, 92);
  });

  test('activityCliffsWithEmptyRows', async () => {
    await _testActivityCliffsOpen(actCliffsEmptyRowsDf, 91);
  });

  after(async () => {
    actCliffsTableView.close();
    actCliffsEmptyRowsTableView.close();
  });

});