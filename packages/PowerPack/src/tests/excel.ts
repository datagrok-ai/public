import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, delay, test} from '@datagrok-libraries/test/src/test';

// Test xlsx files live in the public data.datagrok.ai bucket (s3://datagrok-data/tests/excel),
// read through an anonymous S3 connection created in before().
let dc: DG.DataConnection;

async function testExcelImport(fileName: string, isBenchmarkTest: boolean) {
  if (isBenchmarkTest && !DG.Test.isInBenchmark)
    return;
  let error = '';
  let views: DG.TableView[] = [];
  try {
    const tables = await grok.data.files.openTables(`${dc.nqName}/${fileName}`);
    views = tables.map((table) => grok.shell.addTableView(table));
    await delay(10);
  } catch (e) {
    error = (e as Error).message;
  } finally {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
    views.forEach((view) => view.close());
    expect(error, '', `'${error}' is shown for the correct input`);
  }
}

category('Excel', () => {
  category('Excel: Import', () => {
    before(async () => {
      // Random suffix: parallel CI runs with a same-named connection cross-delete each other
      dc = DG.DataConnection.create(`PowerPack Excel Tests ${DG.Utils.randomString(8)}`, {
        dataSource: 'S3',
        region: 'us-east-2',
        bucket: 'datagrok-data',
        dir: 'tests/excel',
        anonymous: true,
      } as any);
      dc = await grok.dapi.connections.save(dc);
    });
    after(async () => {
      await grok.dapi.connections.delete(dc);
    });
    test('rich text test', async () => await testExcelImport('excel-rich-text-test.xlsx', false));
    test('1MB', async () => await testExcelImport('excel-1mb.xlsx', false));
    test('5MB', async () => await testExcelImport('excel-5mb.xlsx', true), {benchmark: true});
    test('10MB', async () => await testExcelImport('excel-10mb.xlsx', true), {benchmark: true});
    test('40MB 2 spreadsheets', async () => await testExcelImport('excel-40mb-2-spreadsheets.xlsx', true), {benchmark: true, benchmarkTimeout: 70000, timeout: 70000});
    test('50.2MB', async () => await testExcelImport('excel-50.2mb.xlsx', true), {benchmark: true, benchmarkTimeout: 90000, timeout: 90000});
    test('80MB MAX', async () => await testExcelImport('excel-80mb-max.xlsx', true), {benchmark: true, benchmarkTimeout: 120000, timeout: 120000});
  });
}, {owner: 'dkovalyov@datagrok.ai', clear: false});
