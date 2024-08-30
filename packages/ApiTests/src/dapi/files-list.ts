import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, expectArray, expectTable, test} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

category('Dapi: files: list', () => {
  const filePath = 'System:AppData/ApiTests/list';
  const fileList = ['file1.txt', 'file2.csv', 'file3.json'];

  test('all', async () => {
    const allFileList = await grok.dapi.files.list(filePath);
    expectArray(allFileList.map((fi) => fi.name), fileList);
    expectArray(allFileList.map((fi) => fi.fullPath), fileList.map((fn) => `${filePath}/${fn}`));
  });

  test('with-name', async () => {
    const logPrefix = `ApiTests: Tests.Dapi.files.list.with-name`;
    for (const fn of fileList) {
      const fiList = await grok.dapi.files.list(filePath, false, fn);
      _package.logger.debug(`${logPrefix}, fn: '${fn}', fiList (#${fiList.length}):\n${fiList.map((fi) => fi.fullPath).join('\n')}`);
      expect(fiList.length, 1);
      const fi = fiList[0];
      expect(fi.name, fn);
      expect(fi.fullPath, `${filePath}/${fn}`);
    }
  });

  test('empty', async () => {
    const logPrefix = `ApiTests: Tests.Dapi.files.list.empty`;
    const fiList = await grok.dapi.files.list(filePath, false, '');
    expect(fiList.length, 3);
  });

  test('with-name-not-exist', async () => {
    const logPrefix = `ApiTests: Tests.Dapi.files.list.with-name-not-exist`;
    const fn = `file-name-not-exist.ext`;
    const fiList = await grok.dapi.files.list(filePath, false, fn);
    expect(fiList.length, 0);
  });
});

