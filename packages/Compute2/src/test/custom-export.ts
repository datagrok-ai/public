import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectExceptionAsync} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {replaceForWindowsPath} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {exportedFile} from '../components/TreeWizard/tree-export';
import {exportCall, viewExports} from '../run-exports';

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

  test('a custom view lists Default Excel and gives its exports the hook of the run\'s function', async () => {
    const exports = viewExports(DG.Func.byName('Compute2:TestCustomView'));
    expectDeepEqual(exports.map(({name}) => name), ['Default Excel', 'rec', 'hook']);
    const call = DG.Func.byName('Compute2:ViewersHookModel').prepare({n: 2});
    const file = await exports.find(({name}) => name === 'hook')!.export(call);
    expect(new TextDecoder().decode(file.data), 'function|undefined');
  });
});
