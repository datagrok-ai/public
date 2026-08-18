/** Suggestion engine + toolbox pane: context-driven ranking, suppression of
 *  already-wired steps, the ≤10 cap, and the end-to-end pane flow. */
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
    // The rule only exists if the catalog declares Molecule-qualified column inputs (Chem installed).
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

    const append = byName('AppendTables');
    if (append) {
      const ap = items.find((x) => x.typeName === append.nodeTypeName);
      expect(ap != null, true, 'wrapped AppendTables suggested for two tables');
      expect(ap!.wire.length, 2, 'both tables wired to the exposed inputs');
      expect(ap!.wire.map((w) => w.toInput).sort().join(','), 'table1,table2',
        'wired to the wrapper-exposed sockets');
    }

    const withMol = computeSuggestions(ctxOf({
      nodeCount: 2, selectedCount: 2,
      tables: [
        tableSignal({nodeId: 'a', nodeLabel: 'A', outputKey: 'table',
          columns: [{name: 'smiles', type: 'string', semType: 'Molecule'}]}),
        tableSignal({nodeId: 'b', nodeLabel: 'B', outputKey: 'table'}),
      ],
    }));
    const joinIdx = withMol.findIndex((x) => x.typeName === join.nodeTypeName);
    expect(joinIdx > -1, true, 'JoinTables still suggested next to chem matches');
    const firstSingle = withMol.findIndex((x) => x.wire.length === 1);
    expect(firstSingle === -1 || joinIdx < firstSingle, true,
      'two-table suggestions rank above single-table (chem) ones');
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

      let ctx = await collectSuggestContext(e.flow, null, null, null);
      expect(ctx.nodeCount, 3);
      expect(ctx.tables.length, 2, 'both table-bearing nodes scanned');
      expect(ctx.tables.every((t) => !t.selected), true);
      expect(ctx.wiredTargets.has(`${a.id}|Outputs/Table Output`), true, 'existing wiring indexed');

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

      const idx = view.suggestionPane.suggestions.findIndex((x) => x.typeName === join.nodeTypeName);
      expect(idx >= 0, true, 'Join Tables suggested for the two selected tables');
      const joinItem = list.children[idx] as HTMLElement;
      expect(joinItem.draggable, true, 'items are draggable');
      joinItem.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
      expect(await until(() => flow.getNodeCount() === 3), true, 'double-click created the join node');
      expect(await until(() => flow.getConnectionCount() === 2), true, 'both tables auto-connected');

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
      // This view is detached (rects are all zero), so only the strip-row presence is asserted.
      const dropped = flow.getNodes().find((n) => n.dgTypeName === 'Outputs/Table Output')!;
      const asRow = await until(() =>
        canvas.querySelector(`.ff-output-row[data-node-id="${dropped.id}"]`) != null);
      expect(asRow, true, 'dropped output node renders as an Outputs-strip row');

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

  test('drag-out menu leads with the pane\'s picks and searches descriptions', async () => {
    const view = new FuncFlowView();
    try {
      await until(() => (view as never as {flow?: unknown}).flow != null);
      const flow = (view as never as {flow: import('../rete/flow-editor').FlowEditor}).flow;
      const a = await addNode(flow, 'Inputs/Table Input', 0, 0);

      const expected = computeSuggestions(await collectSuggestContext(flow, null, a.id, null))
        .filter((s) => s.wire.some((w) => w.fromNodeId === a.id));
      expect(expected.length > 0, true, 'the engine has picks for a table node');

      const openMenu = (flow as never as {openSuggestionMenu: (
        x: number, y: number, src: {nodeId: string; outputKey: string; dgType: string}) => Promise<void>})
        .openSuggestionMenu.bind(flow);
      const menuDone = openMenu(200, 200, {nodeId: a.id, outputKey: 'table', dgType: 'dataframe'});

      expect(await until(() => document.querySelector('.ff-suggest-popup') != null), true, 'popup opened');
      const popup = document.querySelector('.ff-suggest-popup')!;
      const rows = (): HTMLElement[] => Array.from(popup.querySelectorAll('.ff-suggest-item'));

      const first = rows()[0];
      expect(first.classList.contains('ff-suggest-item-suggested'), true, 'first row is an engine pick');
      expect(first.dataset.nodeTypeName, expected[0].typeName, 'same top pick as the Suggestions pane');
      const reason = first.querySelector('.ff-suggest-item-reason');
      expect(reason?.textContent, expected[0].reason, 'the pane\'s reason shows inline');
      for (let i = 0; i < expected.length && i < rows().length; i++) {
        expect(rows()[i].dataset.nodeTypeName, expected[i].typeName,
          `engine pick #${i} leads the menu in engine order`);
      }

      // "marks a dataframe" is Table Output's description — its name says nothing of the sort.
      const search = popup.querySelector<HTMLInputElement>('.ff-suggest-search')!;
      search.value = 'marks a dataframe';
      search.dispatchEvent(new Event('input', {bubbles: true}));
      expect(await until(() =>
        rows().length > 0 && rows().every((r) => r.dataset.nodeTypeName === 'Outputs/Table Output')),
      true, 'description-only query narrows to Table Output');

      rows()[0].dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true, button: 0}));
      await menuDone;
      expect(flow.getNodeCount(), 2, 'chosen node created');
      const conn = flow.getConnections().find((c) => c.source === a.id && String(c.sourceOutput) === 'table');
      expect(conn != null, true, 'dragged output auto-connected to the chosen node');
      expect(document.querySelector('.ff-suggest-popup'), null, 'popup closed');
    } finally {
      ((view as never as {flow?: {destroy?: () => void}}).flow)?.destroy?.();
      view.root.remove();
    }
  }, {timeout: 60000});
});
