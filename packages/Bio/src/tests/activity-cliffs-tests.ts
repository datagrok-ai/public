import {after, before, category, test} from '@datagrok-libraries/utils/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {readDataframe} from './utils';
import {_testActivityCliffsOpen} from './activity-cliffs-utils';


category('activityCliffs', async () => {
  let actCliffsTableView: DG.TableView;
  let actCliffsDf: DG.DataFrame;
  let actCliffsTableViewWithEmptyRows: DG.TableView;
  let actCliffsDfWithEmptyRows: DG.DataFrame;

  let viewList: DG.ViewBase[] = [];
  let dfList: DG.DataFrame[] = [];

  before(async () => {
    viewList = [];
    dfList = [];
  });

  after(async () => {
    for (const view of viewList) view.close();
    for (const df of dfList) grok.shell.closeTable(df);
  });

  test('activityCliffsOpens', async () => {
    actCliffsDf = await readDataframe('tests/sample_MSA_data.csv');
    dfList.push(actCliffsDf);
    actCliffsTableView = grok.shell.addTableView(actCliffsDf);
    viewList.push(actCliffsTableView);

    await _testActivityCliffsOpen(actCliffsDf, 57, 'UMAP', 'MSA');
  });

  test('activityCliffsWithEmptyRows', async () => {
    actCliffsDfWithEmptyRows = await readDataframe('tests/sample_MSA_data_empty_vals.csv');
    dfList.push(actCliffsDfWithEmptyRows);
    actCliffsTableViewWithEmptyRows = grok.shell.addTableView(actCliffsDfWithEmptyRows);
    viewList.push(actCliffsTableViewWithEmptyRows);

    await _testActivityCliffsOpen(actCliffsDfWithEmptyRows, 57, 'UMAP', 'MSA');
  });
});
