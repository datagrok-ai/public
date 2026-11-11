import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectExceptionAsync} from '@datagrok-libraries/test/src/test';
import {applyCustomExport} from '../components/TreeWizard/tree-export';

// Headless check for the custom-export "report handler" contract: the function named in
// `meta.customExports` (NOT the model itself) must run, receive the step's funcCall and return a file.
// Fixtures: Compute2:TestCustomExportModel (declares the exports), Compute2:TestCustomExportRecorder
// (returns a file echoing the funcCall it is handed) and Compute2:TestCustomExportNoFile — see package.ts.
category('Custom export: report handler', () => {
  test('passes the step funcCall to the declared export function and returns its file', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 777});

    const file = await applyCustomExport(fc, 'rec');

    expect(file instanceof DG.FileInfo, true, 'the export must return a file');
    expect(file.fileName, 'recorded.txt');
    expect(new TextDecoder().decode(file.data), 'Compute2:TestCustomExportModel|777',
      'the declared export must run and receive the model funcCall with its inputs');
  });

  test('throws for an unknown export name', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 1});
    await expectExceptionAsync(async () => {await applyCustomExport(fc, 'missing', {});});
  });

  test('throws for an export that returns no file', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 1});
    await expectExceptionAsync(async () => {await applyCustomExport(fc, 'none');});
  });
});
