import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectFloat, delay, test} from '@datagrok-libraries/test/src/test';

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

/** The table's columns (name: type) and its rows, a date as its UTC day. */
function expectSheet(df: DG.DataFrame, columns: {[name: string]: string}, rows: (string | number)[][]): void {
  expectArray(df.columns.names(), Object.keys(columns));
  expectArray(df.columns.toList().map((c) => c.type), Object.values(columns));
  expect(df.rowCount, rows.length, `rows of ${df.name}`);
  for (let i = 0; i < rows.length; i++) {
    for (let j = 0; j < rows[i].length; j++) {
      const col = df.columns.byIndex(j);
      const where = `${df.name}: "${col.name}" of row ${i + 1}`;
      const value = col.type === DG.COLUMN_TYPE.DATE_TIME ? col.get(i).toISOString().slice(0, 10) : col.get(i);
      if (col.type === DG.COLUMN_TYPE.FLOAT)
        expectFloat(value, rows[i][j] as number, 0.0001, where);
      else
        expect(value, rows[i][j], where);
    }
  }
}

category('Excel', () => {
  category('Excel: Import', () => {
    test('rich text test', async () => await testExcelImport('excel-rich-text-test.xlsx', 1));
    // The only multi-sheet coverage: one dataframe per worksheet.
    test('two sheets', async () => await testExcelImport('excel-two-sheets.xlsx', 2));
    // The workbook the bdd features open (bdd/fixtures/xlsx-open-test.xlsx): a table per sheet, in their order.
    test('three sheets', async () => {
      const tables = await grok.data.files.openTables(`${FIXTURES}/xlsx-open-test.xlsx`);
      expectArray(tables.map((t) => t.name), ['Customers', 'Orders', 'Products']);
      expectSheet(tables[0], {CustomerID: 'int', Name: 'string', Country: 'string', Since: 'datetime'}, [
        [1, 'Alfreds', 'Germany', '2020-03-15'],
        [2, 'Bottom Dollar', 'Canada', '2019-11-02'],
        [3, 'Chop-suey', 'Switzerland', '2021-06-30'],
        [4, 'Du monde entier', 'France', '2018-01-20'],
        [5, 'Ernst Handel', 'Austria', '2022-09-09'],
      ]);
      expectSheet(tables[1], {OrderID: 'int', CustomerID: 'int', Amount: 'double', Shipped: 'string'}, [
        [101, 1, 250.5, 'yes'],
        [102, 3, 99.99, 'no'],
        [103, 2, 1200, 'yes'],
        [104, 5, 47.25, 'yes'],
        [105, 4, 310, 'no'],
        [106, 1, 18.75, 'yes'],
      ]);
      expectSheet(tables[2], {ProductID: 'int', Product: 'string', Price: 'double', InStock: 'int'}, [
        [11, 'Chai', 18, 39],
        [12, 'Chang', 19, 17],
        [13, 'Aniseed Syrup', 10, 13],
        [14, 'Tofu', 23.25, 35],
      ]);
    });
  });
}, {owner: 'dkovalyov@datagrok.ai', clear: false});
