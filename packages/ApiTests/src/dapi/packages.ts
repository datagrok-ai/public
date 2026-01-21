import type * as _grok from 'datagrok-api/grok';
declare let grok: typeof _grok;

import {category, expect, test, expectExceptionAsync} from '@datagrok-libraries/test/src/test';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

category('Dapi: packages', () => {
  test('list', async () => {
    const list = (await grok.dapi.packages.list());
    expect(list.some((pack) => pack.name === 'ApiTests'), true);
  }, {stressTest: true});

  test('find', async () => {
    const apiTestsPackage = await grok.dapi.packages.find(_package.id);
    expect(apiTestsPackage.version, _package.version);
    const res = await grok.dapi.packages.find('00000000-0000-0000-0000-000000000000');
    expect(typeof res === 'undefined');
  }, {stressTest: true});

  test('webRoot content', async () => {
    const apiTestsPackage = await grok.dapi.packages.find(_package.id);
    expect(apiTestsPackage.webRoot + '/', _package.webRoot);
  });

  test('readCsv error', async () => {
    await expectExceptionAsync(() => _package.files.readCsv('datasets/noFile.csv').then());
  });
});
