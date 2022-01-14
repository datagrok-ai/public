import {category/*, expect*/, test} from '@datagrok-libraries/utils/src/test';
import {
  _testRCSBAlive,
  /*requireText,*/
} from './utils';
/*
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
*/

category('biosv', () => {
  test('biosv.pdb_entry.RSCB_is_alive', async () => {
    await _testRCSBAlive();
  });
});
