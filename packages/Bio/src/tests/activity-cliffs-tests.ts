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

  test('activityCliffsOpens', async () => {
    actCliffsDf = await readDataframe('tests/sample_MSA_data.csv');
    actCliffsTableView = grok.shell.addTableView(actCliffsDf);
    await _testActivityCliffsOpen(actCliffsDf, 57, 'UMAP', 'MSA');
    grok.shell.closeTable(actCliffsDf);
    actCliffsTableView.close();
  });

  test('activityCliffsWithEmptyRows', async () => {
    actCliffsDfWithEmptyRows = await readDataframe('tests/sample_MSA_data_empty_vals.csv');
    actCliffsTableViewWithEmptyRows = grok.shell.addTableView(actCliffsDfWithEmptyRows);
    await _testActivityCliffsOpen(actCliffsDfWithEmptyRows, 57, 'UMAP', 'MSA');
    grok.shell.closeTable(actCliffsDfWithEmptyRows);
    actCliffsTableViewWithEmptyRows.close();
  });
});
