import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, expect, expectArray, expectFloat,
  expectObject, isDialogPresent, test} from '@datagrok-libraries/utils/src/test';
import {AddNewColumnDialog} from './add-new-column';
import {FormulaLinesDialog} from './formula-lines';
import {FUNC_TESTS} from './func-tests';


category('Dialogs', () => {
  let tv: DG.TableView;
  let df: DG.DataFrame;
  const dialogs: DG.Dialog[] = [];

  before(async () => {
    df = grok.data.demo.demog(1000);
    tv = grok.shell.addTableView(df);
  });

  test('AddNewColumn', async () => {
    const dlg = new AddNewColumnDialog();
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    dialogs.push(dlg.uiDialog!);
    const funcs = DG.Func.find();

    for (const f of funcs) {
      if (!(f.name in FUNC_TESTS))
        continue;
      for (const [expression, result] of Object.entries(FUNC_TESTS[f.name])) {
        dlg.inputName!.value = expression;
        dlg.inputExpression!.value = expression;
        dlg.uiDialog!.getButton('OK').click();
        await awaitCheck(() => df.col(expression) != null);
        expectTyped(df.col(expression)!.get(0), result);
      }
    }
  });

  test('FormulaLines', async () => {
    const dlg = new FormulaLinesDialog(df);
    await awaitCheck(() => isDialogPresent(dlg.dialog.title));
  });

  after(async () => {
    dialogs.forEach((d) => d.close());
    tv.close();
    grok.shell.closeTable(df);
  });
});

function expectTyped(actual: any, expected: any) {
  if (Array.isArray(expected))
    expectArray(actual, expected);
  else if (typeof expected === 'object' && expected instanceof dayjs)
    expect((<dayjs.Dayjs>expected).isSame(actual), true);
  else if (typeof expected === 'object' && expected !== null)
    expectObject(actual, expected);
  else if (typeof expected === 'number' && Number.isFinite(expected) && !Number.isInteger(expected))
    expectFloat(actual, expected);
  else
    expect(actual == expected, true);
}
