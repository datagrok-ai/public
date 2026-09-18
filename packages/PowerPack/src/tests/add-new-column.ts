import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, expect,
  isDialogPresent, test} from '@datagrok-libraries/test/src/test';
import {AddNewColumnDialog} from '../dialogs/add-new-column';
import {FUNC_HINTS, FUNC_TESTS, FUNC_VALIDATION} from './utils';
import { expectTyped } from './dialogs';


category('Add new column', () => {
  let df: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(10);
    df.name = 'demog';
    grok.shell.addTableView(df);
  });

  test('functions without errors', async () => {
    const funcs = Object.keys(FUNC_TESTS).map((name) => DG.Func.find({ name: name })[0]);

    for (const f of funcs) {
      const call = DG.Func.find({ name: 'AddNewColumn' })[0].prepare({'table': df});
      const dlg = new AddNewColumnDialog(call);
      await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
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
    const call = DG.Func.find({ name: 'AddNewColumn' })[0].prepare({'table': df});
    const dlg = new AddNewColumnDialog(call);
    await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
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

  // The editor also serves hosts that store a FORMULA as a value (a pipeline
  // node's filter condition, a viewer's `filter` property) — `expressionEditorOnly`.
  // Deciding whether a formula is a condition means evaluating it, so the check
  // is asynchronous and the verdict has to reach a host in another package.
  test('expression mode validates and never touches the table', async () => {
    const call = DG.Func.find({name: 'AddNewColumn'})[0]
      .prepare({table: df, expression: '', name: 'Condition', type: DG.COLUMN_TYPE.BOOL});
    call.setAuxValue('expressionEditorOnly', true);
    call.setAuxValue('filterFormulaEditor', true);
    const widget = new DG.Widget(ui.div());
    const dlg = new AddNewColumnDialog(call, widget);
    await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
    const columnsBefore = df.columns.length;

    const type = (expression: string) => dlg.codeMirror!.dispatch({changes: {
      from: 0, to: dlg.codeMirror!.state.doc.length, insert: expression}});
    // Whatever this table actually calls its numbers — a hard-coded name that
    // is merely missing would fail the same assertions for the wrong reason.
    const numeric = df.columns.toList().find((c) => c.matches('numerical'))!.name;

    // Wait for the verdict's CONTENT, not merely for an attribute — checking
    // presence alone is satisfied by whatever the previous formula left behind.
    type(`\${${numeric}} + 1`);
    await awaitCheck(() => (widget.root.getAttribute('data-expression-error') ?? '').includes('true/false'),
      'an arithmetic formula was accepted as a condition', 20000);
    expect(dlg.error.includes('true/false'), true, `unhelpful message: "${dlg.error}"`);

    type(`\${${numeric}} > 30`);
    await awaitCheck(() => !widget.root.hasAttribute('data-expression-error'),
      'the error outlived the formula that caused it', 20000);

    // Accepting publishes the text back onto the call — a host is storing a
    // value, so appending a column to its table would be a side effect nobody
    // asked for.
    await dlg.addNewColumnAction();
    expect(call.getParamValue('expression'), `\${${numeric}} > 30`);
    expect(df.columns.length, columnsBefore, 'the table gained a column');
  });

  test('hints', async () => {
    const call = DG.Func.find({ name: 'AddNewColumn' })[0].prepare({'table': df});
    const dlg = new AddNewColumnDialog(call);
    await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
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

  test('rename during preview refresh', async () => {
    const call = DG.Func.find({ name: 'AddNewColumn' })[0]
      .prepare({'table': df, 'name': 'My new column', 'expression': '${age}*2', 'type': 'auto'});
    const dlg = new AddNewColumnDialog(call);
    await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    // pre-filled expression triggers the initial preview via a 1s debounce
    await awaitCheck(() => dlg.gridPreview?.dataFrame.col('My new column') != null,
      'initial preview has not been computed', 10000);
    // start a preview refresh and rename while it is still computing
    const refresh = dlg.updatePreview('${age}*2', false);
    dlg.inputName!.value = 'first rename';
    await dlg.updatePreview('${age}*2', true);
    await refresh;
    expect(dlg.gridPreview!.dataFrame.col('first rename') != null, true,
      `preview grid shows '${dlg.gridPreview!.dataFrame.columns.names().join(', ')}' instead of renamed column`);
    // subsequent renames must keep working
    dlg.inputName!.value = 'second rename';
    await dlg.updatePreview('${age}*2', true);
    expect(dlg.gridPreview!.dataFrame.col('second rename') != null, true,
      'rename after a preview refresh has no effect');
    dlg.uiDialog!.close();
  });

  test('function at caret becomes current object', async () => {
    const call = DG.Func.find({name: 'AddNewColumn'})[0].prepare({'table': df});
    const dlg = new AddNewColumnDialog(call);
    await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    await awaitCheck(() => dlg.findFunc('Abs') != null, 'functions are not registered', 5000);
    dlg.codeMirror!.dispatch({changes: {from: 0, to: 0, insert: 'Abs(${age})'}});
    await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === 'Abs(${age})');
    dlg.codeMirror!.dispatch({selection: {anchor: 2}, userEvent: 'select.pointer'});
    await awaitCheck(() => grok.shell.o instanceof DG.Func && (grok.shell.o as DG.Func).name === 'Abs',
      'Abs is not the current object', 3000);
    // inside the argument list the hint still describes the enclosing call
    dlg.codeMirror!.dispatch({selection: {anchor: 6}, userEvent: 'select'});
    await awaitCheck(() => dlg.hintDiv.children[0]?.textContent === 'Abs(x:num): num', 'no signature inside the call', 3000);
    expect(dlg.hintDiv.children[1]?.textContent ?? '', dlg.findFunc('Abs')!.description ?? '', 'wrong description in the hint');
    expect(dlg.functionInfo('Abs')[0]?.textContent, 'Abs(x:num): num', 'no completion info');
    const table = dlg.findFunc('Table');
    if (table)
      expect(dlg.functionInfo('Table')[1]?.textContent, table.description, 'no description in completion info');
    dlg.uiDialog!.close();
  });

  test('table and column argument selectors', async () => {
    const {CompletionContext} = await import('@codemirror/autocomplete');
    const {EditorState} = await import('@codemirror/state');
    const call = DG.Func.find({name: 'AddNewColumn'})[0].prepare({'table': df});
    const dlg = new AddNewColumnDialog(call);
    await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    dlg.coreFunctionsParams['TestTableFunc'] = {isVectorFunc: false, params: [
      {propName: 'table', propertyType: DG.TYPE.DATA_FRAME},
      {propName: 'column', propertyType: DG.TYPE.COLUMN},
      {propName: 'output', propertyType: DG.TYPE.INT},
    ]};
    // String parameters annotated the way Column(columnName, [tableName]) is: the
    // semtype names the selector, `options.table` names the table parameter.
    dlg.coreFunctionsParams['TestNamesFunc'] = {isVectorFunc: false, params: [
      {propName: 'columnName', propertyType: DG.TYPE.STRING,
        semType: DG.SEMTYPE.COLUMN_NAME, tableParam: 'tableName', columnTypeFilter: DG.COLUMN_TYPE_FILTER.NUMERICAL},
      {propName: 'row', propertyType: DG.TYPE.INT},
      {propName: 'tableName', propertyType: DG.TYPE.STRING, semType: DG.SEMTYPE.TABLE_NAME},
      {propName: 'output', propertyType: DG.TYPE.INT},
    ]};
    const other = DG.DataFrame.fromCsv('x,y\n1,2');
    other.name = 'other';
    const otherView = grok.shell.addTableView(other);
    const complete = (text: string, pos?: number) =>
      dlg.argumentCompletions(new CompletionContext(EditorState.create({doc: text}), pos ?? text.length, true));
    const labels = (text: string, pos?: number) => complete(text, pos)?.options.map((o) => o.label) ?? [];
    const applied = (text: string, label: string) =>
      complete(text)?.options.find((o) => o.label === label)?.apply;
    try {
      const ctx = dlg.getCallContext('Foo("a(b", Bar(1, ', 18);
      expect(ctx?.funcName, 'Bar');
      expect(ctx?.argIndex, 1);
      expect(dlg.getCallContext('Foo("a(b", ', 11)?.argIndex, 1);

      expect(labels('TestTableFunc(').includes(df.name), true, 'table names are not offered');
      expect(labels('TestTableFunc(').includes('other'), true, 'the second table is not offered');
      expect(applied('TestTableFunc(', 'other'), '"other"');
      expect(applied('TestTableFunc("', 'other'), 'other"');
      expect(complete('TestTableFunc("oth')?.from, 15);

      expect(labels(`TestTableFunc("${df.name}", `).includes('age'), true, 'source columns are not offered');
      expect(applied(`TestTableFunc("${df.name}", `, 'age'), '${age}');
      expect(labels('TestTableFunc("other", ').join(','), 'x,y');
      expect(applied('TestTableFunc("other", ', 'x'), '"x"');

      const placeholder = complete('TestTableFunc(table, column)', 19)!;
      expect(placeholder.from, 14);
      expect(placeholder.filter, false);

      expect(complete('Abs('), null);
      expect(complete('TestTableFunc("other", ${'), null);

      // semtype-driven selectors on string parameters
      // earlier tests in this category add calculated columns to `df`, so derive the expectation
      const numerical = df.columns.toList().filter((c) => c.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL)).map((c) => c.name);
      expect(labels('TestNamesFunc(').join(','), numerical.join(','), 'numerical filter is not applied');
      expect(applied('TestNamesFunc(', 'age'), '"age"');
      expect(labels('TestNamesFunc("age", 0, ').includes('other'), true, 'table names are not offered for TableName');
      expect(labels('TestNamesFunc("x", 0, "other")', 15).join(','), 'x,y', 'options.table is not honoured');
      expect(complete('TestNamesFunc("age", '), null);
    } finally {
      otherView.close();
      grok.shell.closeTable(other);
      dlg.uiDialog!.close();
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
    const call = DG.Func.find({ name: 'AddNewColumn' })[0].prepare({'table': df});
    const dlg = new AddNewColumnDialog(call);
    await awaitCheck(() => dlg.codeMirror != null, 'cannot load CodeMirror', 5000);
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    const absFuncLink = dlg.widgetFunctions?.root.querySelector('[name="span-Abs"]') as HTMLElement;
    //check function is added on click
    absFuncLink.click();
    await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === 'Abs(x)', 'expression has\'t been set');
    await clear();
    //check function is added on click and selected column with matching type is added automatically
    dlg.columnsDf!.currentRowIdx = 3;
    await awaitCheck(() => dlg.selectedColumn?.name === 'age', 'column has\'t been set');
    absFuncLink.click();
    await awaitCheck(() => dlg.codeMirror!.state.doc.toString() === 'Abs(${age})', 'expression has\'t been set', 3000);
    await awaitCheck(() => dlg.gridPreview!.dataFrame.col('Abs(${age})') ? dlg.gridPreview!.dataFrame.get('Abs(${age})', 0) === 61 : false, 'incorrect preview data', 3000);
    await awaitCheck(() => dlg.gridPreview!.dataFrame.get('Abs(${age})', 9) === 26, 'incorrect preview data', 1000);
  });
}, {clear: false});
