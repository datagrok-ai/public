/** Tests for the custom function-editor integration: routing
 *  (`shouldUseFunctionEditor`), the FuncCall→panel value conversion, the
 *  write-back rules (connected inputs win, columns come back as names), and
 *  the Input Parameters header icon. The dialog round-trip itself
 *  (`createFuncCallEditor`) is interactive and exercised manually. */
import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {PropertyPanel} from '../panel/property-panel';
import {
  FuncEditorLauncher, applyEditorResult, editorValueToPanelValue, tableParamForColumn,
} from '../panel/func-editor-launcher';
import {shouldUseFunctionEditor, pollDialogCreation} from '../utils/func-editor-utils';
import {ExecutionController} from '../execution/execution-controller';
import {NodeExecStatus} from '../execution/execution-state';
import {FuncNode} from '../rete/nodes/func-node';
import {FlowNode} from '../rete/scheme';
import {safeGet} from '../utils/dart-proxy-utils';
import {makeEditor, destroyEditor, addNode} from './test-utils';

/** A registered func by exact name, or null. */
function registeredByName(name: string): {func: DG.Func; nodeTypeName: string} | null {
  const info = getRegisteredFuncs().find((f) => f.func.name === name);
  return info ? {func: info.func, nodeTypeName: info.nodeTypeName} : null;
}

category('Flow: func editor', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('shouldUseFunctionEditor: allowlist, editor meta, and the plain majority', async () => {
    // AddNewColumn is explicitly allowlisted (core:AddNewColumn).
    const anc = DG.Func.find({name: 'AddNewColumn'})[0];
    if (anc) expect(shouldUseFunctionEditor(anc), true, 'AddNewColumn is allowlisted');

    // Any function declaring `editor:` meta qualifies.
    const withEditor = DG.Func.find({package: 'Chem'}).find((f) => !!safeGet(f.options, 'editor'));
    if (withEditor)
      expect(shouldUseFunctionEditor(withEditor), true, `${withEditor.name} declares an editor`);

    // A function without editor meta and off the allowlist does not.
    const plain = getRegisteredFuncs().find((f) => {
      try {
        return !safeGet(f.func.options, 'editor') && f.func.nqName !== 'core:AddNewColumn';
      } catch {return false;}
    });
    expect(!!plain, true, 'a plain function exists');
    expect(shouldUseFunctionEditor(plain!.func), false, `${plain!.func.name} has no custom editor`);
  });

  test('editorValueToPanelValue converts editor values to panel representations', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.INT, 'a', [1, 2]),
      DG.Column.fromList(DG.TYPE.INT, 'b', [3, 4]),
    ]);
    // A live column → its name (the panel edits name strings).
    expect(editorValueToPanelValue(df.col('a'), 'column'), 'a');
    // A column list → comma-separated names.
    expect(editorValueToPanelValue([df.col('a'), df.col('b')], 'column_list'), 'a, b');
    // string_list → comma-separated (the panel's comma-separated field).
    expect(editorValueToPanelValue(['x', 'y'], 'string_list'), 'x, y');
    // A plain `list` stays an array (DG's List input holds a JS array).
    const arr = editorValueToPanelValue(['x', 'y'], 'list');
    expect(Array.isArray(arr), true, 'list stays an array');
    // A dataframe can't live in the panel.
    expect(editorValueToPanelValue(df, 'dataframe') === undefined, true);
    // Primitives pass through.
    expect(editorValueToPanelValue(42, 'int'), 42);
    expect(editorValueToPanelValue('text', 'string'), 'text');
    expect(editorValueToPanelValue(false, 'bool'), false);
  });

  test('tableParamForColumn: explicit association → suffix pairing → first table', async () => {
    const fakeNode = (columnTables?: Record<string, string>): FlowNode =>
      ({properties: columnTables ? {columnTables} : {}} as unknown as FlowNode);
    // Explicit association wins.
    expect(tableParamForColumn(fakeNode({keys2: 'table1'}), 'keys2', ['table1', 'table2']), 'table1');
    // No association → shared numeric suffix.
    expect(tableParamForColumn(fakeNode(), 'keys2', ['table1', 'table2']), 'table2');
    // No suffix match → the first dataframe input.
    expect(tableParamForColumn(fakeNode(), 'column', ['table1', 'table2']), 'table1');
    // A stale association pointing at a non-existent table falls through.
    expect(tableParamForColumn(fakeNode({keys2: 'gone'}), 'keys2', ['table1', 'table2']), 'table2');
    // No dataframe inputs at all → undefined.
    expect(tableParamForColumn(fakeNode(), 'column', []) === undefined, true);
  });

  test('applyEditorResult writes back editable inputs; connected inputs win', async () => {
    const anc = DG.Func.find({name: 'AddNewColumn'})[0];
    if (!anc) return; // not on this stand — skip
    const node = new FuncNode(anc);
    // Pick two primitive editable params to exercise (name/expression on ANC).
    const prims = anc.inputs.filter((p) =>
      ['string', 'int', 'double', 'bool'].includes(String(p.propertyType)) && p.name in node.inputValues);
    if (prims.length < 2) return;
    const [p1, p2] = [prims[0].name, prims[1].name];

    const fc = anc.prepare({});
    fc.setParamValue(p1, 'fromEditor1');
    fc.setParamValue(p2, 'fromEditor2');

    // Nothing connected → both written back.
    const applied = applyEditorResult(node, anc, fc, () => false);
    expect(applied.includes(p1) && applied.includes(p2), true, 'both params applied');
    expect(String(node.inputValues[p1]), 'fromEditor1');
    expect(String(node.inputValues[p2]), 'fromEditor2');

    // p1 connected → the wire wins, only p2 is updated.
    fc.setParamValue(p1, 'shouldNotLand');
    fc.setParamValue(p2, 'secondPass');
    const applied2 = applyEditorResult(node, anc, fc, (n) => n === p1);
    expect(applied2.includes(p1), false, 'connected param skipped');
    expect(String(node.inputValues[p1]), 'fromEditor1', 'connected value untouched');
    expect(String(node.inputValues[p2]), 'secondPass');
  });

  test('applyEditorResult converts a column value back to a name string', async () => {
    // Any registered func with an editable `column` input.
    const info = getRegisteredFuncs().find((f) =>
      f.func.inputs.some((p) => String(p.propertyType) === 'column'));
    if (!info) return; // none on this stand — skip
    const func = info.func;
    const colParam = func.inputs.find((p) => String(p.propertyType) === 'column')!.name;
    const node = new FuncNode(func);
    if (!(colParam in node.inputValues)) return; // connection-only on this func — skip

    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.TYPE.INT, 'picked', [1, 2, 3])]);
    const fc = func.prepare({});
    fc.setParamValue(colParam, df.col('picked'));
    applyEditorResult(node, func, fc, () => false);
    expect(String(node.inputValues[colParam]), 'picked', `${func.name}.${colParam} came back as a name`);
  });

  test('launcher: gate refuses without a table; a captured table opens the editor', async () => {
    // End-to-end against the live platform: unconnected → refused (balloon);
    // connected with a captured upstream table → the function's OWN dialog
    // ("Add New Column", not a generic form) opens seeded with it; closing it
    // resolves the round-trip and writes the values back into the node.
    const anc = registeredByName('AddNewColumn');
    if (!anc) return;
    const e = makeEditor();
    try {
      const ctrl = new ExecutionController(e.flow);
      const launcher = new FuncEditorLauncher(e.flow, ctrl, () => ({name: 'T', description: '', tags: []}));
      const node = await addNode(e.flow, anc.nodeTypeName);

      // 1. Table not connected → gate refuses without opening any dialog.
      expect(await launcher.open(node), false, 'unconnected table input → refused');

      // 2. Connect a table input and pre-populate its captured result, so the
      //    launcher reuses it (no confirm / no slice run — the stash path).
      const src = await addNode(e.flow, 'Inputs/Table Input', 0, 200);
      await e.flow.addConnectionByKeys(src.id, 'table', node.id, 'table');
      const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.TYPE.INT, 'a', [1, 2, 3])]);
      ctrl.state.setNodeStatus(src.id, NodeExecStatus.completed,
        {outputs: {table: {type: 'dataframe', clone: df}}});
      node.inputValues['name'] = 'myNewCol';

      // 3. Launch; the function's own editor dialog opens; close it → the
      //    round-trip resolves and the seeded value survives the write-back.
      const openPromise = launcher.open(node);
      const dlg = await pollDialogCreation(10_000);
      expect(!!dlg, true, 'the editor dialog opened');
      expect(dlg!.title, 'Add New Column', 'the function’s own editor, not a generic form');
      dlg!.close();
      expect(await openPromise, true, 'the round-trip resolved');
      expect(String(node.inputValues['name']), 'myNewCol', 'seeded value survives the round-trip');
    } finally {
      destroyEditor(e);
    }
  });

  test('the Input Parameters header shows the editor icon only when applicable', async () => {
    const anc = registeredByName('AddNewColumn');
    if (!anc) return; // not in the catalog on this stand — skip
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, anc.nodeTypeName);

      // Callback not wired → no icon.
      panel.showNode(node);
      expect(!!panel.root.querySelector('[data-testid="ff-prop-func-editor"]'), false,
        'no icon without onEditFuncParams');

      // Callback wired + editor-capable func → icon in the pane header.
      let clicked: FlowNode | null = null;
      panel.onEditFuncParams = (n): void => {clicked = n;};
      panel.showNode(node);
      const icon = panel.root.querySelector('[data-testid="ff-prop-func-editor"]') as HTMLElement | null;
      expect(!!icon, true, 'icon rendered for an editor-capable function');
      icon!.click();
      expect(clicked === node, true, 'clicking the icon passes the node to the callback');

      // A function without a custom editor → no icon even with the callback.
      const plain = getRegisteredFuncs().find((f) => {
        try {
          return !shouldUseFunctionEditor(f.func) && f.func.inputs.length > 0;
        } catch {return false;}
      });
      if (plain) {
        const plainNode = await addNode(e.flow, plain.nodeTypeName, 100, 100);
        panel.showNode(plainNode);
        expect(!!panel.root.querySelector('[data-testid="ff-prop-func-editor"]'), false,
          `no icon for ${plain.func.name} (no custom editor)`);
      }
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });
});
