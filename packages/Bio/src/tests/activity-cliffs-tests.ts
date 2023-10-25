import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test} from '@datagrok-libraries/utils/src/test';

import {readDataframe} from './utils';
import {_testActivityCliffsOpen} from './activity-cliffs-utils';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';

import {_package} from '../package-test';


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
    // for (const df of dfList) grok.shell.closeTable(df);
    // for (const view of viewList) view.close();
  });

  test('activityCliffsOpens', async () => {
    actCliffsDf = await readDataframe(
      DG.Test.isInBenchmark ? 'test/peptides_motif-with-random_10000.csv' : 'tests/100_3_clustests.csv',
    );
    dfList.push(actCliffsDf);
    actCliffsTableView = grok.shell.addTableView(actCliffsDf);
    viewList.push(actCliffsTableView);
    const cliffsNum = DG.Test.isInBenchmark ? 6 : 3;

    await _testActivityCliffsOpen(actCliffsDf, DimReductionMethods.UMAP,
      'sequence', 'Activity', 90, cliffsNum);
  });

  test('activityCliffsWithEmptyRows', async () => {
    actCliffsDfWithEmptyRows = await readDataframe('tests/100_3_clustests_empty_vals.csv');
    dfList.push(actCliffsDfWithEmptyRows);
    actCliffsTableViewWithEmptyRows = grok.shell.addTableView(actCliffsDfWithEmptyRows);
    viewList.push(actCliffsTableViewWithEmptyRows);

    await _testActivityCliffsOpen(actCliffsDfWithEmptyRows, DimReductionMethods.UMAP,
      'sequence', 'Activity', 90, 3);
  });

  test('Helm', async () => {
    const df = await _package.files.readCsv('samples/sample_HELM.csv');
    const view = grok.shell.addTableView(df);

    await _testActivityCliffsOpen(df, DimReductionMethods.UMAP,
      'HELM', 'Activity', 90, 53);
  });
});
