import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, expect, expectArray, expectFloat,
  expectObject, isDialogPresent, test} from '@datagrok-libraries/utils/src/test';
import {AddNewColumnDialog} from '../dialogs/add-new-column';
import {FormulaLinesDialog} from '../dialogs/formula-lines';
import {FUNC_TESTS} from './utils';


category('Dialogs', () => {
  let df: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(2);
    grok.shell.addTableView(df);
  });

  test('AddNewColumn', async () => {
    const dlg = new AddNewColumnDialog();
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    const funcs = Object.keys(FUNC_TESTS).map((name) => DG.Func.find({ name: name })[0]);

    for (const f of funcs) {
      for (const [expression, result] of Object.entries(FUNC_TESTS[f.name])) {
        const columnName = df.columns.getUnusedName(expression);
        dlg.inputName!.value = columnName;
        dlg.codeMirror!.dispatch({
          changes: {
            from: 0,
            to: dlg.codeMirror!.state.doc.length,
            insert: expression
          }
        });
        await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === expression);
        dlg.uiDialog!.getButton('OK').click();
        try {
          await awaitCheck(() => df.columns.contains(columnName));
          const column = df.col(columnName)!;
          expectTyped(column.get(0), result == null && column.type === DG.TYPE.STRING ? '' : result);
        } catch (e) {
          throw new Error(`Expression: ${expression}: ${(e as Error).message}`);
        }
      }
    }
  });

  test('FormulaLines', async () => {
    const dlg = new FormulaLinesDialog(df);
    await awaitCheck(() => isDialogPresent(dlg.dialog.title));
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
}, {clear: false});

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
