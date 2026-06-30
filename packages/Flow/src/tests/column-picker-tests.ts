/** Column picker on func-node column inputs: the property panel renders a
 *  picker icon next to each `column` / `column_list` input, and the request it
 *  fires resolves to the correct dataframe input (per-column for multi-table
 *  funcs like JoinTables: keys1→table1, keys2→table2). */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, createNode} from '../rete/node-factory';
import {PropertyPanel} from '../panel/property-panel';
import {ColumnPickRequest} from '../panel/column-picker';
import {ExecutionController} from '../execution/execution-controller';
import {NodeExecStatus} from '../execution/execution-state';
import {makeEditor, destroyEditor, addNode} from './test-utils';

/** Registered factory name for JoinTables, or null if absent on this server. */
function joinTablesTypeName(): string | null {
  const info = getRegisteredFuncs().find((f) => f.func.name === 'JoinTables');
  return info?.nodeTypeName ?? null;
}

category('Flow: column picker', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('column inputs render a picker that resolves the right table per column', async () => {
    const typeName = joinTablesTypeName();
    if (!typeName) return; // JoinTables not registered on this server — skip
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const join = await addNode(e.flow, typeName);
      const reqs: ColumnPickRequest[] = [];
      panel.onPickColumns = (r) => reqs.push(r);
      panel.showNode(join);

      const btn1 = panel.root.querySelector('[data-testid="ff-prop-pick-columns-keys1"]') as HTMLElement | null;
      const btn2 = panel.root.querySelector('[data-testid="ff-prop-pick-columns-keys2"]') as HTMLElement | null;
      expect(!!btn1, true, 'keys1 has a picker button');
      expect(!!btn2, true, 'keys2 has a picker button');

      btn1!.click();
      expect(reqs.length, 1, 'click fires exactly one request');
      expect(reqs[0].paramName, 'keys1');
      expect(reqs[0].isList, true, 'keys1 is a column_list (multi-select)');
      expect(reqs[0].tableParam, 'table1', 'keys1 resolves against table1');

      btn2!.click();
      expect(reqs.length, 2);
      expect(reqs[1].tableParam, 'table2', 'keys2 resolves against table2');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });

  test('cloneForNode reuses a completed result but never a stale one', async () => {
    const e = makeEditor();
    try {
      const ctrl = new ExecutionController(e.flow);
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('c', ['a', 'b'])]);
      // A completed slice run leaves a captured clone — the picker reuses it
      // (so a second pick against the same table doesn't re-run it).
      ctrl.state.setNodeStatus('n1', NodeExecStatus.completed, {outputs: {res: {type: 'dataframe', clone: df}}});
      expect(ctrl.cloneForNode('n1') === df, true, 'completed result is reused');
      // A graph edit marks results stale → the picker must recompute, not reuse
      // a possibly-outdated table.
      ctrl.state.setNodeStatus('n1', NodeExecStatus.stale, {});
      expect(ctrl.cloneForNode('n1'), null, 'stale result is not reused');
      expect(ctrl.cloneForNode('missing'), null, 'unknown node → null');
    } finally {
      destroyEditor(e);
    }
  });

  test('no picker is rendered when the panel has no onPickColumns handler', async () => {
    const typeName = joinTablesTypeName();
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    document.body.appendChild(panel.root);
    try {
      const join = await addNode(e.flow, typeName);
      panel.showNode(join);   // onPickColumns left unset
      const btn = panel.root.querySelector('[data-testid="ff-prop-pick-columns-keys1"]');
      expect(!!btn, false, 'no picker icon without a handler');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });
});
