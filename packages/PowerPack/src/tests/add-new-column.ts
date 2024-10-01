import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, expect,
  isDialogPresent, test} from '@datagrok-libraries/utils/src/test';
import {AddNewColumnDialog} from '../dialogs/add-new-column';
import {FUNC_HINTS, FUNC_TESTS, FUNC_VALIDATION} from './utils';
import { expectTyped } from './dialogs';


category('Add new column', () => {
  let df: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(10);
    grok.shell.addTableView(df);
  });

  test('functions without errors', async () => {
    const funcs = Object.keys(FUNC_TESTS).map((name) => DG.Func.find({ name: name })[0]);

    for (const f of funcs) {
      const dlg = new AddNewColumnDialog();
      await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
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
        expect(dlg.error, '', `Error '${dlg.error}' is shown for the correct input`);
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

  test('validation', async () => {
    const dlg = new AddNewColumnDialog();
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    for (const f of Object.keys(FUNC_VALIDATION)) {
        dlg.codeMirror!.dispatch({
          changes: {
            from: 0,
            to: dlg.codeMirror!.state.doc.length,
            insert: f
          }
        });
        await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === f, 'expression has\'t been set');
        await awaitCheck(() => dlg.error === FUNC_VALIDATION[f], 'incorrect validation error');
    }
  });

  test('hints', async () => {
    const dlg = new AddNewColumnDialog();
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    for (const f of Object.keys(FUNC_HINTS)) {
      dlg.codeMirror!.dispatch({
        changes: {
          from: 0,
          to: dlg.codeMirror!.state.doc.length,
          insert: f
        }
      });
      await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === f, 'expression has\'t been set');
      dlg.codeMirror!.dispatch({selection: {anchor: 5, head: 5}})
      await awaitCheck(() => dlg.hintDiv.children[0].innerHTML === FUNC_HINTS[f], 'incorrect hint');
    }
  });

  test('insert function on click', async () => {
    const clear = async () => {
      dlg.codeMirror!.dispatch({
        changes: {
          from: 0,
          to: dlg.codeMirror!.state.doc.length,
          insert: ''
        }
      });
      await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === '', 'code mirror has\'t been cleared');
    }      
    const dlg = new AddNewColumnDialog();
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    const absFuncLink = dlg.widgetFunctions?.root.querySelector('[name="span-Abs"]') as HTMLElement;
    //check function is added on click
    absFuncLink.click();
    await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === 'Abs(num)', 'expression has\'t been set');
    await clear();
    //check function is added on click and selected column with matching type is added automatically
    dlg.columnsDf!.currentRowIdx = 3;
    await awaitCheck(() => dlg.selectedColumn!.name === 'age', 'column has\'t been set');
    absFuncLink.click();
    await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === 'Abs(${age})', 'expression has\'t been set', 3000);
    await awaitCheck(() => dlg.gridPreview!.dataFrame.col('Abs(${age})') ? dlg.gridPreview!.dataFrame.get('Abs(${age})', 0) === 61 : false, 'incorrect preview data', 3000);
    await awaitCheck(() => dlg.gridPreview!.dataFrame.get('Abs(${age})', 9) === 26, 'incorrect preview data', 1000);
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
}, {clear: false});