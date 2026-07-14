/** Suggestion engine + toolbox pane: context-driven ranking (semType columns
 *  → domain ops, clicked Molecule cell → prefilled searches, two selected
 *  tables → auto-wired Join), suppression of already-wired steps, the ≤10 cap,
 *  and the end-to-end pane flow (render → accept → node created and wired). */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, FuncInfo} from '../rete/node-factory';
import {
  computeSuggestions, collectSuggestContext, SuggestContext, TableSignal, MAX_SUGGESTIONS,
} from '../suggest/suggestion-engine';
import {FuncFlowView} from '../funcflow-view';
import {FF_SUGGEST_MIME} from '../panel/suggestion-pane';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

function ctxOf(over: Partial<SuggestContext>): SuggestContext {
  return {
    nodeCount: 1, selectedCount: 0, tables: [], scalars: [], cell: null,
    canvasFuncNames: new Set(), canvasDomains: new Set(), wiredTargets: new Set(),
    ...over,
  };
}

function tableSignal(over: Partial<TableSignal>): TableSignal {
  return {nodeId: 'n1', nodeLabel: 'Open File', outputKey: 'result', passthrough: false,
    selected: true, columns: [], ...over};
}

function findFunc(pred: (f: FuncInfo) => boolean): FuncInfo | undefined {
  return getRegisteredFuncs().find((f) => {
    try {
      return pred(f);
    } catch {
      return false;
    }
  });
}

const byName = (name: string): FuncInfo | undefined =>
  findFunc((f) => f.func.name.toLowerCase() === name.toLowerCase());

category('Flow: suggestions', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('empty canvas suggests starting points', async () => {
    const items = computeSuggestions(ctxOf({nodeCount: 0}));
    expect(items.some((s) => s.typeName === 'Inputs/Table Input'), true, 'Table Input offered');
    if (byName('OpenFile'))
      expect(items.some((s) => s.label.toLowerCase().includes('file')), true, 'Open File offered');
  });

  test('a Molecule column drives chem suggestions, wired and prefilled', async () => {
    // Data-driven guard: the rule only exists if the catalog has functions
    // declaring a Molecule-qualified column input (Chem installed).
    const molFunc = findFunc((f) =>
      f.func.inputs.some((p) => String(p.propertyType) === 'column' && String(p.semType ?? '') === 'Molecule') &&
      f.func.inputs.some((p) => String(p.propertyType) === 'dataframe'));
    if (!molFunc) return;

    const items = computeSuggestions(ctxOf({
      tables: [tableSignal({columns: [
        {name: 'smiles', type: 'string', semType: 'Molecule'},
        {name: 'activity', type: 'double', semType: null},
      ]})],
    }));
    const chem = items.filter((s) => s.prefill && Object.values(s.prefill).includes('smiles'));
    expect(chem.length > 0, true, 'semType-matched functions suggested');
    expect(chem[0].wire.length, 1, 'wired to the source table');
    expect(chem[0].wire[0].fromNodeId, 'n1');
    expect(items[0].reason.includes('Molecule'), true, 'semType suggestions rank first');
    if (byName('descriptors'))
      expect(items.some((s) => s.label.toLowerCase().includes('descriptor')), true, 'featured op present');
  });

  test('a Macromolecule column drives bio suggestions', async () => {
    const bioFunc = findFunc((f) =>
      f.func.inputs.some((p) => String(p.propertyType) === 'column' && String(p.semType ?? '') === 'Macromolecule'));
    if (!bioFunc) return;
    const items = computeSuggestions(ctxOf({
      tables: [tableSignal({columns: [{name: 'sequence', type: 'string', semType: 'Macromolecule'}]})],
    }));
    expect(items.some((s) => s.reason.includes('Macromolecule')), true, 'bio ops suggested');
  });

  test('two selected tables suggest Join Tables wired to both', async () => {
    const join = byName('JoinTables');
    if (!join) return;
    const items = computeSuggestions(ctxOf({
      nodeCount: 2, selectedCount: 2,
      tables: [
        tableSignal({nodeId: 'a', nodeLabel: 'A', outputKey: 'table'}),
        tableSignal({nodeId: 'b', nodeLabel: 'B', outputKey: 'table'}),
      ],
    }));
    const s = items.find((x) => x.typeName === join.nodeTypeName);
    expect(s != null, true, 'JoinTables suggested');
    expect(s!.wire.length, 2, 'both tables wired');
    expect(s!.wire[0].fromNodeId, 'a');
    expect(s!.wire[1].fromNodeId, 'b');
    expect(s!.wire[0].toInput !== s!.wire[1].toInput, true, 'distinct dataframe inputs');
  });

  test('a single table suggests common next steps and matching viewers', async () => {
    const items = computeSuggestions(ctxOf({
      tables: [tableSignal({columns: [
        {name: 'x', type: 'int', semType: null},
        {name: 'y', type: 'double', semType: null},
        {name: 'group', type: 'string', semType: null},
      ]})],
    }));
    expect(items.some((s) => s.typeName === 'Viewers/Scatter Plot'), true, '2 numeric → scatter plot');
    expect(items.some((s) => s.typeName === 'Viewers/Bar Chart' || s.typeName === 'Viewers/Histogram'),
      true, 'category/numeric viewer offered');
    for (const s of items)
      if (s.wire.length > 0) expect(s.wire[0].fromNodeId, 'n1', 'everything wires to the context table');
  });

  test('clicked Molecule cell suggests single-molecule searches with the value prefilled', async () => {
    const cellFunc = findFunc((f) =>
      f.func.inputs.some((p) => String(p.propertyType) !== 'column' &&
        String(p.propertyType) !== 'dataframe' && String(p.semType ?? '') === 'Molecule'));
    if (!cellFunc) return;
    const items = computeSuggestions(ctxOf({
      cell: {semType: 'Molecule', column: 'smiles', value: 'CCO'},
    }));
    const s = items.find((x) => x.prefill && Object.values(x.prefill).includes('CCO'));
    expect(s != null, true, 'clicked value prefills a search function');
    expect(s!.wire.length, 0, 'value functions need no wire');
  });

  test('already-wired steps are suppressed and the cap holds', async () => {
    const withOut = computeSuggestions(ctxOf({tables: [tableSignal({})]}));
    expect(withOut.some((s) => s.typeName === 'Outputs/Table Output'), true);
    const suppressed = computeSuggestions(ctxOf({
      tables: [tableSignal({})],
      wiredTargets: new Set(['n1|Outputs/Table Output']),
    }));
    expect(suppressed.some((s) => s.typeName === 'Outputs/Table Output'), false, 'wired → gone');

    const busy = computeSuggestions(ctxOf({
      tables: [tableSignal({columns: [
        {name: 'smiles', type: 'string', semType: 'Molecule'},
        {name: 'seq', type: 'string', semType: 'Macromolecule'},
        {name: 'x', type: 'int', semType: null}, {name: 'y', type: 'double', semType: null},
      ]})],
      cell: {semType: 'Molecule', column: 'smiles', value: 'CCO'},
    }));
    expect(busy.length <= MAX_SUGGESTIONS, true, 'never more than 10');
  });

  test('collectSuggestContext reads selection, focus order, and wiring', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const b = await addNode(e.flow, 'Inputs/Table Input', 300, 0);
      const out = await addNode(e.flow, 'Outputs/Table Output', 600, 0);
      await e.flow.addConnectionByKeys(a.id, 'table', out.id, 'table');

      // Nothing selected → canvas-wide fallback: both tables, none marked selected.
      let ctx = await collectSuggestContext(e.flow, null, null, null);
      expect(ctx.nodeCount, 3);
      expect(ctx.tables.length, 2, 'both table-bearing nodes scanned');
      expect(ctx.tables.every((t) => !t.selected), true);
      expect(ctx.wiredTargets.has(`${a.id}|Outputs/Table Output`), true, 'existing wiring indexed');

      // Selection: only selected nodes enter, focus first.
      await e.flow.selectNode(b.id);
      await e.flow.selectNode(a.id, true);
      ctx = await collectSuggestContext(e.flow, null, a.id, null);
      expect(ctx.tables.length, 2);
      expect(ctx.tables[0].nodeId, a.id, 'focus node first');
      expect(ctx.tables.every((t) => t.selected), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('pane end-to-end: renders suggestions and accepting Join wires both tables', async () => {
    const join = byName('JoinTables');
    if (!join) return;
    try {
      localStorage.removeItem('funcflow-suggestions-collapsed');
    } catch {/* blocked storage */}
    const view = new FuncFlowView();
    try {
      await until(() => (view as never as {flow?: unknown}).flow != null);
      const flow = (view as never as {flow: import('../rete/flow-editor').FlowEditor}).flow;
      const a = await addNode(flow, 'Inputs/Table Input', 0, 0);
      const b = await addNode(flow, 'Inputs/Table Input', 0, 200);
      await flow.selectNode(a.id);
      await flow.selectNode(b.id, true);

      await view.suggestionPane.refreshNow();
      expect(view.suggestionPane.suggestions.length > 0, true, 'pane has suggestions');
      const list = view.suggestionPane.root.querySelector('[data-testid="ff-suggest-pane-list"]')!;
      expect(list.children.length > 0, true, 'items rendered');

      // Double-click the Join item — same gesture as the toolbox catalog.
      const idx = view.suggestionPane.suggestions.findIndex((x) => x.typeName === join.nodeTypeName);
      expect(idx >= 0, true, 'Join Tables suggested for the two selected tables');
      const joinItem = list.children[idx] as HTMLElement;
      expect(joinItem.draggable, true, 'items are draggable');
      joinItem.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
      expect(await until(() => flow.getNodeCount() === 3), true, 'double-click created the join node');
      expect(await until(() => flow.getConnectionCount() === 2), true, 'both tables auto-connected');

      // Drag-drop: a suggestion travels as JSON under FF_SUGGEST_MIME; the drop
      // applies it (wiring included) at the drop point.
      const canvas = (view as never as {canvasContainer: HTMLElement}).canvasContainer;
      const r = canvas.getBoundingClientRect();
      const dt = new DataTransfer();
      dt.setData(FF_SUGGEST_MIME, JSON.stringify({
        typeName: 'Outputs/Table Output', label: 'Table Output', reason: 'test', score: 1,
        wire: [{fromNodeId: a.id, fromOutputKey: 'table', toInput: 'table'}],
      }));
      canvas.dispatchEvent(new DragEvent('drop', {
        dataTransfer: dt, bubbles: true, cancelable: true,
        clientX: r.left + 400, clientY: r.top + 300,
      }));
      expect(await until(() => flow.getNodeCount() === 4), true, 'drop created the node');
      expect(await until(() => flow.getConnectionCount() === 3), true, 'drop applied the suggested wiring');
      const dropped = flow.getNodes().find((n) => n.dgTypeName === 'Outputs/Table Output')!;
      const want = flow.screenToCanvas(r.left + 400, r.top + 300);
      expect(Math.abs(dropped.pos.x - want.x) < 1 && Math.abs(dropped.pos.y - want.y) < 1, true,
        'dropped node lands at the drop point');

      // Collapse via caret: list hidden, state persisted.
      const caret = view.suggestionPane.root.querySelector<HTMLElement>('[data-testid="ff-suggest-pane-caret"]')!;
      caret.click();
      expect(view.suggestionPane.root.getAttribute('data-collapsed'), 'true');
      caret.click();
      expect(view.suggestionPane.root.getAttribute('data-collapsed'), 'false');
    } finally {
      ((view as never as {flow?: {destroy?: () => void}}).flow)?.destroy?.();
      view.root.remove();
    }
  }, {timeout: 60000});
});
