import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, expectTable, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';

category('Dapi: files', () => {
  const filePrefix = 'System:AppData/ApiTests/';

  test('exists', async () => {
    const path = filePrefix + DG.Utils.randomString(6);
    try {
      await grok.dapi.files.writeAsText(path, 'testString');
      if (!await grok.dapi.files.exists(path))
        throw new Error('File doesn\'t exist');
    } finally {
      await safeDeleteFile(path);
    }
  }, {stressTest: true});

  test('write/read text', async () => {
    const filePath = filePrefix + `${DG.Utils.randomString(6)}.txt`;
    const fileText = 'testString';

    try {
      await grok.dapi.files.writeAsText(filePath, fileText);
      if (fileText !== await grok.dapi.files.readAsText(filePath))
        throw new Error('Сontent is wrong');
    } finally {
      await safeDeleteFile(filePath);
    }
  }, {stressTest: true});

  test('write/read blob', async () => {
    const filePath = filePrefix + `${DG.Utils.randomString(6)}.txt`;
    const content = [0, 1, 2, 3];

    try {
      await grok.dapi.files.write(filePath, content);
      if (content.toString() !== (await grok.dapi.files.readAsBytes(filePath)).toString())
        throw new Error('Сontent is wrong');
    } finally {
      await safeDeleteFile(filePath);
    }
  }, {stressTest: true});

  test('search', async () => {
    const fileName = DG.Utils.randomString(6);
    const path = filePrefix + fileName;
    try {
      await grok.dapi.files.writeAsText(path, 'testString');
      if ((await grok.dapi.files.list(filePrefix, false, fileName)).length !== 1)
        throw new Error('Can\'t find the file');
    } finally {
      await safeDeleteFile(path);
    }

  }, {stressTest: true});

  test('package files', async () => {
    const files = await _package.files.list('datasets', false, 'csv');
    expect(files.length > 0, true);
    files.every((f) => expect(f.extension, 'csv'));
  }, {stressTest: true});

  // test('move', async () => {
  //     let fileName = 'move.txt';
  //     let filePath = filePrefix + fileName;
  //     let newFilePrefix = 'texts';
  //
  //     await grok.dapi.files.writeAsText(filePath, 'testString');
  //     await grok.dapi.files.move([filePath], newFilePrefix);
  //
  //     if(await grok.dapi.files.exists(filePath))
  //         throw "Old file exists";
  //
  //     if(!await grok.dapi.files.exists(filePrefix + newFilePrefix + '/' + fileName))
  //         throw "New file doesn't exist";
  //
  //     await grok.dapi.files.delete(filePath);
  // });

  test('delete', async () => {
    const filePath = filePrefix + `${DG.Utils.randomString(6)}.txt`;
    try {
      await grok.dapi.files.writeAsText(filePath, 'testString');
      await grok.dapi.files.delete(filePath);

      if (await grok.dapi.files.exists(filePath))
        throw new Error('File exists');
    } finally {
      await safeDeleteFile(filePath);
    }
  }, {stressTest: true});

  test('readBinaryDataFrames', async () => {
    const dfList = await _package.files.readBinaryDataFrames('datasets/country-languages.d42');
    expect(dfList.length, 1);
    expect(dfList[0] instanceof DG.DataFrame, true);
  }, {stressTest: true});

  test('writeBinaryDataFrames', async () => {
    const df = grok.data.demo.demog(10);
    const filePath = `${DG.Utils.randomString(6)}.d42`;
    try {
      //@ts-ignore
      await _package.files.writeBinaryDataFrames(filePath, [df]);
      const dfList = await _package.files.readBinaryDataFrames(filePath);
      expect(dfList.length, 1, `Saved ${dfList.length} dataframes instead of 1`);
      expectTable(dfList[0], df, 'Saved dataframe has wrong data');
      await _package.files.delete(filePath);
    } finally {
      await safeDeleteFile(filePrefix + filePath);
    }
  });

  test('readAsText', async () => {
    const files = await _package.files.list('datasets', true, 'demog.csv');
    const res = await _package.files.readAsText(files[0]);
    expect(!!res);
  }, {stressTest: true});

}, {owner: 'aparamonov@datagrok.ai'});

category('Dapi: files: formats', () => {
  const extensions = ['csv', 'd42', 'json', 'tar', 'tar.gz', 'tsv', 'txt', 'xml', 'zip', 'kmz', 'kml'];

  for (const ext of extensions) {
    test(ext, async () => {
      const df = await grok.data.files.openTable('System:AppData/ApiTests/datasets/formats/cars.' + ext);
      expect(df.rowCount, 10, 'wrong rows number');
      expect(df.columns.length, 10, 'wrong columns number');
    }, {skipReason: typeof process !== 'undefined' ? 'NodeJS environment'
          : ['kmz', 'kml'].includes(ext) ? 'GROK-13263'
              : undefined});
  }
}, {owner: 'aparamonov@datagrok.ai'});

async function safeDeleteFile(path: string): Promise<void> {
  try {
    await grok.dapi.files.delete(path);
  } catch (_) {}
}