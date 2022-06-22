import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_testPaletteN, _testPaletteAA} from '@datagrok-libraries/bio/src/tests/palettes.test';

category('Palettes', () => {
  test('testPaletteN', async () => { _testPaletteN(); });
  test('testPaletteAA', async () => { _testPaletteAA(); });
});
