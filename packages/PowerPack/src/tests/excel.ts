import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, delay, test} from '@datagrok-libraries/utils/src/test';


async function testExcelImport(path: string, isBenchmarkTest: boolean) {
  if (isBenchmarkTest && !DG.Test.isInBenchmark)
    return;
  let error = '';
  let views: DG.TableView[] = [];
  try {
    const tables = await grok.data.files.openTables(path);
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
    test('rich text test', async () => await testExcelImport('System:DemoFiles/test/excel/excel-rich-text-test.xlsx', false));
    test('1MB', async () => await testExcelImport('System:DemoFiles/test/excel/excel-1mb.xlsx', false));
    test('5MB', async () => await testExcelImport('System:DemoFiles/test/excel/excel-5mb.xlsx', true), {benchmark: true});
    test('10MB', async () => await testExcelImport('System:DemoFiles/test/excel/excel-10mb.xlsx', true), {benchmark: true});
    test('40MB 2 spreadsheets', async () => await testExcelImport('System:DemoFiles/test/excel/excel-40mb-2-spreadsheets.xlsx', true), {benchmark: true, benchmarkTimeout: 60000, timeout: 60000});
    test('50.2MB', async () => await testExcelImport('System:DemoFiles/test/excel/excel-50.2mb.xlsx', true), {benchmark: true, benchmarkTimeout: 75000, timeout: 75000});
    test('80MB MAX', async () => await testExcelImport('System:DemoFiles/test/excel/excel-80mb-max.xlsx', true), {benchmark: true, benchmarkTimeout: 110000, timeout: 110000});
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
}, {owner: 'dkovalyov@datagrok.ai', clear: false});
