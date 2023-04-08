import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, delay} from '@datagrok-libraries/utils/src/test';
import {readDataframe} from './utils';
import {_testSequenceSpaceReturnsResult} from './sequence-space-utils';

category('sequenceSpace', async () => {
  let testFastaDf: DG.DataFrame;
  let testFastaTableView: DG.TableView;
  let testHelmWithEmptyRows: DG.DataFrame;
  let testHelmWithEmptyRowsTableView: DG.TableView;

  test('sequenceSpaceOpens', async () => {
    testFastaDf = await readDataframe('tests/sample_MSA_data.csv');
    testFastaTableView = grok.shell.addTableView(testFastaDf);
    await _testSequenceSpaceReturnsResult(testFastaDf, 'UMAP', 'MSA');
    grok.shell.closeTable(testFastaDf);
    testFastaTableView.close();
  }, {skipReason: 'GROK-12136'});

  test('sequenceSpaceWithEmptyRows', async () => {
    testHelmWithEmptyRows = await readDataframe('tests/sample_MSA_data_empty_vals.csv');
    testHelmWithEmptyRowsTableView = grok.shell.addTableView(testHelmWithEmptyRows);
    await _testSequenceSpaceReturnsResult(testHelmWithEmptyRows, 'UMAP', 'MSA');
    grok.shell.closeTable(testHelmWithEmptyRows);
    testHelmWithEmptyRowsTableView.close();
  }, {skipReason: 'GROK-12136'});
});
