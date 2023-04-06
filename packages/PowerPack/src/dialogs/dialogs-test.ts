import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, expect, expectArray, expectFloat,
  expectObject, isDialogPresent, test} from '@datagrok-libraries/utils/src/test';
import {AddNewColumnDialog} from './add-new-column';
import {FormulaLinesDialog} from './formula-lines';
import {FUNC_TESTS, EXCLUDED} from './func-tests';


category('Dialogs', () => {
  let tv: DG.TableView;
  let df: DG.DataFrame;
  const dialogs: DG.Dialog[] = [];

  before(async () => {
    df = grok.data.demo.demog(10);
    tv = grok.shell.addTableView(df);
  });

  test('AddNewColumn', async () => {
    const dlg = new AddNewColumnDialog();
    dialogs.push(dlg.uiDialog!);
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    const funcs = DG.Func.find();
    function getUniqueColumnName() {
      let counter = 0;
      function increment() {
        counter++;
        return `Column ${counter}`;
      }
      return increment;
   }
   const getColumnName = getUniqueColumnName();

    for (const f of funcs) {
      if (!(f.name in FUNC_TESTS) || EXCLUDED.includes(f.name))
        continue;
      for (const [expression, result] of Object.entries(FUNC_TESTS[f.name])) {
        const columnName = getColumnName();
        dlg.inputName!.value = columnName;
        dlg.inputExpression!.value = expression;
        dlg.uiDialog!.getButton('OK').click();
        try {
          await awaitCheck(() => df.columns.contains(columnName));
          const column = df.col(columnName)!;
          expectTyped(column.get(0), result == null && column.type === DG.TYPE.STRING ? '' : result);
        } catch(e) {
          throw new Error(`Expression: ${expression}: ${(e as Error).message}`);
        }
      }
    }
  });

  test('FormulaLines', async () => {
    const dlg = new FormulaLinesDialog(df);
    dialogs.push(dlg.dialog);
    await awaitCheck(() => isDialogPresent(dlg.dialog.title));
  });

  after(async () => {
    dialogs.forEach((d) => d.close());
    tv.close();
    grok.shell.closeTable(df);
  });
});

function expectTyped(actual: any, expected: any) {
  if (expected == null)
    expect(actual == expected, true);
  if (Array.isArray(expected))
    expectArray(actual, expected);
  else if (typeof expected === 'object' && expected instanceof dayjs)
    expect((<dayjs.Dayjs>expected).isSame(actual), true);
  else if (typeof expected === 'object')
    expectObject(actual, expected);
  else if (typeof expected === 'number' && Number.isFinite(expected) && !Number.isInteger(expected))
    expectFloat(actual, expected);
  else
    expect(actual == expected, true);
}
