import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, waitUntil, before, category, expect, expectArray, expectFloat,
  expectObject, isDialogPresent, test} from '@datagrok-libraries/test/src/test';
import {AddNewColumnDialog} from '../dialogs/add-new-column';
import {FormulaLinesDialog} from '../dialogs/formula-lines';
import {FUNC_TESTS} from './utils';
import {viewersDialog} from '../viewers-gallery';


category('Dialogs', () => {
  let df: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(2);
    grok.shell.addTableView(df);
  });
  
  test('FormulaLines', async () => {
    const dlg = new FormulaLinesDialog(df);
    await waitUntil(() => isDialogPresent(dlg.dialog.title), 'timed out waiting for isDialogPresent(dlg.dialog.title)');
  });

  test('Viewer gallery with a package-less viewer', async () => {
    grok.functions.register({signature: 'viewer PowerPackTestNoPackageViewer()',
      run: () => DG.Viewer.grid(df), options: {role: DG.FUNC_TYPES.VIEWER}});
    viewersDialog(grok.shell.tv, df);
    await waitUntil(() => isDialogPresent('Add Viewer'), 'timed out waiting for isDialogPresent(\'Add Viewer\')');
    DG.Dialog.getOpenDialogs().find((d) => d.title === 'Add Viewer')?.close();
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
