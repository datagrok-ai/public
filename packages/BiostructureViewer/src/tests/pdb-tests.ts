import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {category/*, expect*/, test} from '@datagrok-libraries/test/src/test';
import {
  _testRCSBAlive,
  _testParseExamplePDBFile,
  _testMolstarViewerIsOpening,
  /*requireText,*/
} from './utils';

category('biosv', () => {
  test('pdb_entry.RCSB_is_alive', async () => {
    await _testRCSBAlive();
  });

  test('pdb_entry.Parse_sample_2V0A', async () => {
    await _testParseExamplePDBFile();
  });

  test('molstar_viewer.is_opening', async () => {
    await _testMolstarViewerIsOpening();
  });
});
