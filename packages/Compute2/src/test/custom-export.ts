import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectExceptionAsync} from '@datagrok-libraries/test/src/test';
import {applyCustomExport} from '../utils';

// Headless check for the custom-export "report handler" contract: the function named in
// `meta.customExports` (NOT the model itself) must run, and it must receive the step's funcCall.
// Fixtures: Compute2:TestCustomExportModel (declares the export) and
// Compute2:TestCustomExportRecorder (echoes back the funcCall it is handed) — see package.ts.
category('Custom export: report handler', () => {
  test('passes the step funcCall to the declared export function', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 777});

    // 'rec' resolves to TestCustomExportRecorder (a different function than the model),
    // which returns `<nqName>|<input a>|<startDownload>` of whatever funcCall it received.
    const res = await applyCustomExport(fc, 'rec', {startDownload: false});

    expect(res, 'Compute2:TestCustomExportModel|777|false',
      'the declared export must run and receive the model funcCall with its inputs');
  });

  test('throws for an unknown export name', async () => {
    const model = DG.Func.byName('Compute2:TestCustomExportModel');
    const fc = model.prepare({a: 1});
    await expectExceptionAsync(async () => { await applyCustomExport(fc, 'missing', {}); });
  });
});
