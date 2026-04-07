import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {category/*, expect*/, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {PdbResDataFrame} from '../utils/pdb-helper';

import {_package} from '../package-test';

category('pdbHelper', () => {
  test('pdbToDf', async () => {
    await _testPdbToDf();
  });
});

async function _testPdbToDf(): Promise<void> {
  const ph: IPdbHelper = await grok.functions.call('BiostructureViewer:getPdbHelper');
  const pdbCnt: string = await grok.dapi.files.readAsText(`System:AppData/${_package.name}/samples/1bdq.pdb`);
  const df: DG.DataFrame = await ph.pdbToDf(pdbCnt, '1bdq');
  expect(df.rowCount, 198);
  expectArray(df.columns.names(), Object.values(PdbResDataFrame.ColNames));
}
