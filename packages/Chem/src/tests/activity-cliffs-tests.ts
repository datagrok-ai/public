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

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    actCliffsTableView = await createTableView('activity_cliffs.csv');
    actCliffsDf = await readDataframe('activity_cliffs.csv');
  });

  test('activityCliffsOpen', async () => {
    await _testActivityCliffsOpen(actCliffsDf);
  });

  after(async () => {
    actCliffsTableView.close();
  });

});