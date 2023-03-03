import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {_packageName} from './utils';
import {category/*, expect*/, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {getPdbHelper} from '../package';
import {PdbResDataFrame} from '../utils/pdb-helper';


category('pdbHelper', () => {
  test('pdbToDf', async () => {
    await _testPdbToDf();
  });
});

async function _testPdbToDf(): Promise<void> {
  const ph: IPdbHelper = await getPdbHelper();
  const pdbCnt: string = await grok.dapi.files.readAsText(`System:AppData/${_packageName}/samples/1bdq.pdb`);
  const df: DG.DataFrame = await ph.pdbToDf(pdbCnt, '1bdq');
  expect(df.rowCount, 198);
  expectArray(df.columns.names(), Object.values(PdbResDataFrame.ColNames));
}
