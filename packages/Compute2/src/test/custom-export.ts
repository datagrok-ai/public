import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectExceptionAsync} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {replaceForWindowsPath} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {exportedFile} from '../call-export';
import {DEFAULT_EXPORT, exportCall} from '../run-exports';

// Headless check for the custom-export "report handler" contract: the function named in
// `meta.customExports` (NOT the model itself) must run, receive the step's funcCall and return a file.
// Fixtures: Compute2:TestCustomExportModel (declares the exports), Compute2:TestCustomExportRecorder
// (returns a file echoing the funcCall it is handed) and Compute2:TestCustomExportNoFile — see package.ts.
category('Custom export: report handler', () => {
  test('passes the step funcCall to the declared export function and returns its file', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 777});

    const file = await exportCall(fc, 'rec');

    expect(file instanceof DG.FileInfo, true, 'the export must return a file');
    expect(file.fileName, 'recorded.txt');
    expect(new TextDecoder().decode(file.data), 'Compute2:TestCustomExportModel|777',
      'the declared export must run and receive the model funcCall with its inputs');
  });

  test('throws for an unknown export name', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 1});
    await expectExceptionAsync(async () => {await exportCall(fc, 'missing');});
  });

  test('makes a file of the Blob a legacy export returns, telling it not to download', async () => {
    const fc = DG.Func.byName('Compute2:TestBlobExportModel').prepare({a: 1});
    const file = await exportCall(fc, 'blob');
    expect(file.fileName, 'Compute2_TestBlobExportModel - Not completed');
    const ran = await fc.call();
    const completed = await exportCall(ran, 'blob', {completed: true});
    expect(completed.fileName, replaceForWindowsPath(`Compute2:TestBlobExportModel - ${ran.started}`));
    expect(new TextDecoder().decode(file.data), 'Compute2:TestBlobExportModel|false');
  });

  test('names a legacy Blob .xlsx only when it is Excel', async () => {
    const excel = new Blob(['x'], {type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'});
    expect((await exportedFile(excel, 'e', 'run')).fileName, 'run.xlsx');
    expect((await exportedFile(new Blob(['x'], {type: 'application/pdf'}), 'e', 'run')).fileName, 'run');
  });

  test('throws for an export that returns no file', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 1});
    await expectExceptionAsync(async () => {await exportCall(fc, 'none');});
  });

  test('CallExports lists the Export menu, without Default Excel when the function hides it', async () => {
    const ids = async (nqName: string) => (await DG.Func.byName('Compute2:CallExports')
      .apply({func: DG.Func.byName(nqName)})).map(({id}: {id: string}) => id);
    expectDeepEqual(await ids('Compute2:TestPreferredExportModel'), ['none', 'hook', DEFAULT_EXPORT, 'rec']);
    expectDeepEqual(await ids('Compute2:TestNoDefaultExportModel'), ['rec']);
  });

  test('ExportCall runs an export as the menu does and passes the args on', async () => {
    const run = (funcCall: DG.FuncCall, exportId: string, args?: object) =>
      DG.Func.byName('Compute2:ExportCall').apply({funcCall, exportId, args});
    const fc = DG.Func.byName('Compute2:TestCustomExportModel').prepare({a: 1});
    const hooked = await run(fc, 'hook', {metaState: {a: {value: 'm'}}});
    expect(new TextDecoder().decode(hooked.data), 'undefined|m');
    const excel = await run(await fc.call(), DEFAULT_EXPORT, {isOutputOutdated: true});
    expect(excel.fileName, 'Compute2_TestCustomExportModel - Not completed.xlsx');
    await expectExceptionAsync(async () => {await run(fc, 'missing');});
  });
});
