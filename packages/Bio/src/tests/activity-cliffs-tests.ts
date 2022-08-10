import {after, before, category, test} from '@datagrok-libraries/utils/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {readDataframe} from './utils';
import { _testActivityCliffsOpen } from './activity-cliffs-utils';


category('activityCliffs', async () => {
  let actCliffsTableView: DG.TableView;
  let actCliffsDf: DG.DataFrame;
  let actCliffsTableViewWithEmptyRows: DG.TableView;
  let actCliffsDfWithEmptyRows: DG.DataFrame;
  

  before(async () => {
    actCliffsDf = await readDataframe('samples/sample_MSA.csv');
    actCliffsTableView = grok.shell.addTableView(actCliffsDf);
    actCliffsDfWithEmptyRows = await readDataframe('samples/sample_HELM_empty_vals.csv');
    actCliffsTableViewWithEmptyRows = grok.shell.addTableView(actCliffsDfWithEmptyRows);
  });

  after(async () => {
    grok.shell.closeTable(actCliffsDf);
    actCliffsTableView.close();
    grok.shell.closeTable(actCliffsDfWithEmptyRows);
    actCliffsTableViewWithEmptyRows.close();
  });

  test('activityCliffsOpen', async () => {
    await _testActivityCliffsOpen(actCliffsDf, 53, 'UMAP', 'MSA');
  });

  test('activityCliffsOpenWithEmptyRows', async () => {
    await _testActivityCliffsOpen(actCliffsDfWithEmptyRows, 53, 'UMAP', 'HELM');
  });
  
});
