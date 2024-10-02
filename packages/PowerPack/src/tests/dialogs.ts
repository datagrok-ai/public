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
  
  test('FormulaLines', async () => {
    const dlg = new FormulaLinesDialog(df);
    await awaitCheck(() => isDialogPresent(dlg.dialog.title));
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
}, {clear: false});

export function expectTyped(actual: any, expected: any) {
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
