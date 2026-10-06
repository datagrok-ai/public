import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, delay, test} from '@datagrok-libraries/test/src/test';

const FIXTURES = 'System:AppData/PowerPack/tests';

async function testExcelImport(fileName: string, expectedTables: number) {
  let error = '';
  let views: DG.TableView[] = [];
  let tables: DG.DataFrame[] = [];
  try {
    tables = await grok.data.files.openTables(`${FIXTURES}/${fileName}`);
    views = tables.map((table) => grok.shell.addTableView(table));
    await delay(10);
  } catch (e) {
    error = (e as Error).message;
  } finally {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
    views.forEach((view) => view.close());
  }
  expect(error, '', `'${error}' is shown for the correct input`);
  expect(tables.length, expectedTables, `${fileName} should open as ${expectedTables} table(s)`);
}

category('Excel', () => {
  category('Excel: Import', () => {
    test('rich text test', async () => await testExcelImport('excel-rich-text-test.xlsx', 1));
    // The only multi-sheet coverage: one dataframe per worksheet.
    test('two sheets', async () => await testExcelImport('excel-two-sheets.xlsx', 2));
  });
}, {owner: 'dkovalyov@datagrok.ai', clear: false});
