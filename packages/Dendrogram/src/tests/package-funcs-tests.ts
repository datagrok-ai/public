import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package-test';
import {
  after,
  before,
  category,
  test,
  expect,
  expectArray,
  expectObject,
  awaitCheck
} from '@datagrok-libraries/utils/src/test';

/** Tests for package functions, test apps, file previews, file handlers, ... */
category('packageFuncs', () => {
  test('previewNewick', async () => {
    const nwkFi: DG.FileInfo = (await grok.dapi.files.list(
      `System:AppData/${_package.name}/data`, false, 'tree95.nwk'))[0];
    await grok.functions.call(`${_package.name}:previewNewick`, {file: nwkFi});
    await awaitCheck(() => true, 'Error', 200);
  });
});
