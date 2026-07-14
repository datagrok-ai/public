/** Precise run-result invalidation: a graph edit marks only the affected
 *  downstream cone "Out of date" (not the whole canvas), and the autorun
 *  scheduler turns the accumulated dirty set into a debounced slice re-run.
 *
 *  Covers: the forward slice (`sliceDownFrom`), the per-edit invalidation
 *  semantics (`ExecutionController.applyGraphEdit`), the classified edit
 *  stream the editor emits (`onGraphEdited`), the live-boundary expansion
 *  that plans a partial re-run, and the `AutorunScheduler` debounce. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, createNode} from '../rete/node-factory';
import {FlowEditor, GraphEdit} from '../rete/flow-editor';
import {missingRequiredProps, nodeMissingRequirements} from '../rete/scheme';
import {sliceDownFrom} from '../compiler/graph-compiler';
import {emitScript} from '../compiler/script-emitter';
import {ExecutionController, expandToLiveBoundary} from '../execution/execution-controller';
import {NodeExecStatus} from '../execution/execution-state';
import {AutorunScheduler, isAutorunByDefault} from '../execution/autorun';
import {PropertyPanel} from '../panel/property-panel';
import {makeEditor, destroyEditor, addNode, until, TestEditor} from './test-utils';

const SETTINGS = {name: 'InvalidationFlow', description: '', tags: ['funcflow']};

const sleepMs = (ms: number): Promise<void> => new Promise((r) => setTimeout(r, ms));

function funcTypeName(name: string): string | null {
  return getRegisteredFuncs().find((f) => f.func.name === name)?.nodeTypeName ?? null;
}

/** A detached editor that records every classified edit (what `makeEditor`
 *  can't do — it passes no callbacks). */
function makeRecordingEditor(): {flow: FlowEditor; container: HTMLElement; edits: GraphEdit[]} {
  const container = ui.div([], {style: {width: '800px', height: '600px', position: 'absolute', left: '-10000px'}});
  document.body.appendChild(container);
  const edits: GraphEdit[] = [];
  const flow = new FlowEditor(container, {onGraphEdited: (edit) => edits.push(edit)});
  return {flow, container, edits};
}

const paramsEdits = (edits: GraphEdit[]): GraphEdit[] => edits.filter((x) => x.kind === 'params-changed');

/** Set the tab-global live-value registry for the duration of a test. */
function withRegistry(reg: Record<string, Record<string, unknown>>, body: () => void): void {
  const g = globalThis as {__ffFlowLive?: unknown};
  const prev = g.__ffFlowLive;
  g.__ffFlowLive = reg;
  try {
    body();
  } finally {
    g.__ffFlowLive = prev;
  }
}

/** A three-node chain: Constants/String → ToString → ToString. */
async function makeChain(e: TestEditor): Promise<{a: string; b: string; c: string}> {
  const a = await addNode(e.flow, 'Constants/String');
  const b = await addNode(e.flow, 'Utilities/ToString');
  const c = await addNode(e.flow, 'Utilities/ToString');
  await e.flow.addConnectionByKeys(a.id, 'value', b.id, 'value');
  await e.flow.addConnectionByKeys(b.id, 'text', c.id, 'value');
  return {a: a.id, b: b.id, c: c.id};
}

function seedCompleted(ctrl: ExecutionController, ids: string[]): void {
  for (const id of ids) ctrl.state.setNodeStatus(id, NodeExecStatus.completed, {});
}

function status(ctrl: ExecutionController, id: string): NodeExecStatus | undefined {
  return ctrl.state.getNodeState(id)?.status;
}

category('Flow: invalidation', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('sliceDownFrom returns the node and its transitive successors', async () => {
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      expect([...sliceDownFrom(e.flow, a)].sort().join(), [a, b, c].sort().join(), 'from the source: everything');
      const fromB = sliceDownFrom(e.flow, b);
      expect(fromB.has(a), false, 'upstream is NOT downstream');
      expect(fromB.has(b) && fromB.has(c), true, 'the node and its successor are');
      expect([...sliceDownFrom(e.flow, c)].join(), c, 'a sink is alone');
    } finally {
      destroyEditor(e);
    }
  });

  test('adding a node invalidates nothing', async () => {
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      const ctrl = new ExecutionController(e.flow);
      seedCompleted(ctrl, [a, b, c]);
      const lone = await addNode(e.flow, 'Constants/String', 500, 0);
      const affected = ctrl.applyGraphEdit({kind: 'node-added', nodeId: lone.id});
      expect(affected.size, 0, 'nothing to recompute');
      for (const id of [a, b, c])
        expect(status(ctrl, id), NodeExecStatus.completed, 'results survive an unrelated node');
    } finally {
      destroyEditor(e);
    }
  });

  test('a connection change invalidates only its target and downstream', async () => {
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      const ctrl = new ExecutionController(e.flow);
      seedCompleted(ctrl, [a, b, c]);
      withRegistry({[a]: {value: 1}, [b]: {text: 2}, [c]: {text: 3}}, () => {
        const affected = ctrl.applyGraphEdit({kind: 'connection-removed', sourceId: a, targetId: b});
        expect(affected.has(b) && affected.has(c), true, 'target and downstream invalidated');
        expect(affected.has(a), false, 'source is untouched');
        expect(status(ctrl, a), NodeExecStatus.completed, 'upstream result survives');
        expect(status(ctrl, b), NodeExecStatus.stale, 'target went stale');
        expect(status(ctrl, c), NodeExecStatus.stale, 'downstream went stale');
        const reg = (globalThis as {__ffFlowLive?: Record<string, unknown>}).__ffFlowLive!;
        expect(a in reg, true, 'upstream captured value kept');
        expect(b in reg || c in reg, false, 'invalidated captured values dropped');
      });
    } finally {
      destroyEditor(e);
    }
  });

  test('a parameter edit invalidates the node and its downstream', async () => {
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      const ctrl = new ExecutionController(e.flow);
      seedCompleted(ctrl, [a, b, c]);
      const affected = ctrl.applyGraphEdit({kind: 'params-changed', nodeId: b});
      expect([...affected].sort().join(), [b, c].sort().join(), 'the edited node and its successors');
      expect(status(ctrl, a), NodeExecStatus.completed, 'upstream result survives');
      expect(status(ctrl, b), NodeExecStatus.stale);
      expect(status(ctrl, c), NodeExecStatus.stale);
    } finally {
      destroyEditor(e);
    }
  });

  test('removing a node drops its state and nothing else', async () => {
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      const ctrl = new ExecutionController(e.flow);
      seedCompleted(ctrl, [a, b, c]);
      ctrl.applyGraphEdit({kind: 'node-removed', nodeId: c});
      // (`expect`'s second arg defaults to true — compare explicitly.)
      expect(status(ctrl, c) === undefined, true, 'the removed node is forgotten');
      expect(status(ctrl, a), NodeExecStatus.completed);
      expect(status(ctrl, b), NodeExecStatus.completed);
    } finally {
      destroyEditor(e);
    }
  });

  test('the editor emits classified edits for every result-affecting change', async () => {
    const container = ui.div([], {style: {width: '600px', height: '400px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    const edits: GraphEdit[] = [];
    const flow = new FlowEditor(container, {onGraphEdited: (edit) => edits.push(edit)});
    try {
      const a = await addNode(flow, 'Constants/String');
      const b = await addNode(flow, 'Utilities/ToString');
      expect(edits.filter((x) => x.kind === 'node-added').length, 2);

      await flow.addConnectionByKeys(a.id, 'value', b.id, 'value');
      const added = edits.find((x) => x.kind === 'connection-added');
      expect(added !== undefined, true, 'connection-added emitted');
      expect((added as {sourceId: string}).sourceId, a.id);
      expect((added as {targetId: string}).targetId, b.id);

      flow.notifyNodeParamsChanged(b.id);
      const params = edits[edits.length - 1];
      expect(params.kind, 'params-changed');
      expect((params as {nodeId: string}).nodeId, b.id);

      await flow.removeNode(b.id);
      expect(edits.some((x) => x.kind === 'connection-removed'), true,
        'deleting a node reports its connections first');
      expect(edits[edits.length - 1].kind, 'node-removed');

      await flow.clear();
      expect(edits[edits.length - 1].kind, 'cleared');
    } finally {
      try {
        flow.destroy();
      } finally {
        container.remove();
      }
    }
  });

  test('clicking a node repeatedly (panel rebuilds) is never a parameter edit', async () => {
    // The user-side regression: select node → context panel builds → its DG
    // inputs fire onValueChanged during init → params-changed → invalidation
    // (and, with autorun on, a rerun) on EVERY click. Three "clicks" with the
    // panel open must produce zero params-changed edits.
    const rec = makeRecordingEditor();
    try {
      const panel = new PropertyPanel(rec.flow);
      // Every node carries a NON-EMPTY stored value: `initInputValue` only
      // initializes non-empty values, and initializing a DG input is exactly
      // what fires the spurious onValueChanged this test protects against.
      const constant = await addNode(rec.flow, 'Constants/String', 0, 0); // utility textarea
      constant.properties['value'] = 'hello';
      const selTable = await addNode(rec.flow, 'Utilities/Select Table', 0, 120); // utility textarea prop
      selTable.properties['tableName'] = 'demog';
      const selColumn = await addNode(rec.flow, 'Utilities/Select Column', 0, 240); // column field (DG input)
      selColumn.properties['columnName'] = 'age';
      const output = await addNode(rec.flow, 'Outputs/Table Output', 0, 360); // textarea + combo
      output.properties['paramName'] = 'result';
      const intInput = await addNode(rec.flow, 'Inputs/Int Input', 0, 480); // number/toggle props
      intInput.properties['defaultValue'] = 5;
      const nodes = [constant, selTable, selColumn, output, intInput];
      const openFile = funcTypeName('OpenFile');
      if (openFile) {
        // The reported scenario: an OpenFile node with a file path, clicked
        // three times with the context panel open (DG forProperty input).
        const of = await addNode(rec.flow, openFile, 0, 600);
        of.inputValues['fullPath'] = 'System:AppData/demo.csv';
        nodes.push(of);
      }
      rec.edits.length = 0;

      for (const node of nodes) {
        for (let click = 0; click < 3; click++) {
          panel.showNode(node);
          await sleepMs(30); // let any async init onValueChanged land
        }
        expect(paramsEdits(rec.edits).length, 0,
          `opening the panel for '${node.label}' must not report an edit`);
      }
    } finally {
      try {
        rec.flow.destroy();
      } finally {
        rec.container.remove();
      }
    }
  });

  test('a real edit reports once; re-entering the same value reports nothing', async () => {
    const rec = makeRecordingEditor();
    try {
      const panel = new PropertyPanel(rec.flow);
      const node = await addNode(rec.flow, 'Utilities/Select Table');
      panel.showNode(node);
      await sleepMs(30);
      rec.edits.length = 0;

      const textarea = panel.root.querySelector('[data-param="tableName"] textarea') as HTMLTextAreaElement;
      expect(textarea !== null, true, 'the tableName editor is rendered');
      const type = (v: string): void => {
        textarea.value = v;
        textarea.dispatchEvent(new Event('input', {bubbles: true}));
      };

      type('demog');
      expect(paramsEdits(rec.edits).length, 1, 'a real change reports');
      expect(String(node.properties['tableName']), 'demog', 'the value was applied');
      type('demog');
      expect(paramsEdits(rec.edits).length, 1, 'the same value again reports nothing');
      type('cars');
      expect(paramsEdits(rec.edits).length, 2, 'a different value reports again');
    } finally {
      try {
        rec.flow.destroy();
      } finally {
        rec.container.remove();
      }
    }
  });

  test('pendingNodes: everything without a fresh result, closed downstream', async () => {
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      const ctrl = new ExecutionController(e.flow);
      expect([...ctrl.pendingNodes()].sort().join(), [a, b, c].sort().join(), 'never ran → everything pends');
      seedCompleted(ctrl, [a]);
      expect([...ctrl.pendingNodes()].sort().join(), [b, c].sort().join(), 'completed upstream is excluded');
      seedCompleted(ctrl, [b, c]);
      expect(ctrl.pendingNodes().size, 0, 'a fully fresh flow pends nothing');
      ctrl.state.markStale([b]);
      expect([...ctrl.pendingNodes()].sort().join(), [b, c].sort().join(), 'stale counts as not run');
    } finally {
      destroyEditor(e);
    }
  });

  test('expandToLiveBoundary grows upstream only past uncaptured values', async () => {
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      // b's output is captured → the slice for a dirty {c} stops at c.
      const withLive = expandToLiveBoundary(e.flow, [c], (nodeId) => nodeId === b);
      expect([...withLive].join(), c, 'captured boundary feeds the slice');
      // Nothing captured → the slice must include the whole ancestry.
      const noLive = expandToLiveBoundary(e.flow, [c], () => false);
      expect([...noLive].sort().join(), [a, b, c].sort().join(), 'uncaptured upstream is pulled in');
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: in-place isolation', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  /** Column names of the first dataframe summary in a node's captured outputs. */
  function previewColNames(ctrl: ExecutionController, nodeId: string): string[] | null {
    const outputs = ctrl.state.getNodeState(nodeId)?.outputs ?? {};
    for (const s of Object.values(outputs))
      if (s.type === 'dataframe' && Array.isArray(s.colNames)) return s.colNames as string[];
    return null;
  }

  test('instrumented emission snapshots dataframe inputs; clean emission does not', async () => {
    const addCol = funcTypeName('AddNewColumn');
    if (!addCol) {
      expect(true, true);
      return;
    }
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');
      const anc = await addNode(e.flow, addCol);
      await e.flow.addConnectionByKeys(src.id, 'table', anc.id, 'table');

      const inst = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'clone-1'});
      expect(inst.includes('function __ff_clone'), true, 'clone helper defined');
      const m = inst.match(/(\w+_table_in) = __ff_clone\(/);
      expect(m !== null, true, 'the table input is snapshot-cloned');
      const snap = m![1];
      expect(inst.includes(`table: ${snap}`), true, 'the call receives the snapshot');
      expect(inst.includes(`"table__pt": ${snap}`), true,
        'the pass-through/stash carry the snapshot — the mutated copy flows on, not the upstream instance');

      const clean = emitScript(e.flow, SETTINGS);
      expect(clean.includes('__ff_clone'), false, 'clean scripts keep the platform in-place idiom');
    } finally {
      destroyEditor(e);
    }
  });

  test('in-place transform runs on a clone: previews stay at-node, reruns are idempotent', async () => {
    // The reported scenario (OpenFile → in-place calc → …): the calc mutated
    // the very instance the upstream node had stashed, so its preview showed
    // downstream columns and an autorun slice re-run applied the calc twice.
    const addColType = funcTypeName('AddNewColumn');
    const addColFunc = DG.Func.find({name: 'AddNewColumn'})[0];
    if (!addColType || !addColFunc) {
      expect(true, true);
      return;
    }
    const nameParam = addColFunc.inputs.find((p) =>
      String(p.propertyType) === 'string' && p.name.toLowerCase() === 'name')?.name;
    const exprParam = addColFunc.inputs.find((p) =>
      p.name.toLowerCase() === 'expression')?.name;
    if (!nameParam || !exprParam) {
      expect(true, true);
      return;
    }

    const e = makeEditor();
    const df = DG.DataFrame.fromCsv('x\n1\n2\n');
    df.name = 'ffInplaceIsolation';
    const shellTable = grok.shell.addTable(df);
    try {
      const sel = await addNode(e.flow, 'Utilities/Select Table');
      sel.properties['tableName'] = 'ffInplaceIsolation';
      const anc = await addNode(e.flow, addColType, 300, 0);
      anc.inputValues[nameParam] = 'extra';
      anc.inputValues[exprParam] = '1';
      await e.flow.addConnectionByKeys(sel.id, 'table', anc.id, 'table');

      const ctrl = new ExecutionController(e.flow);
      const completedBoth = (): boolean => [sel.id, anc.id].every((id) =>
        ctrl.state.getNodeState(id)?.status === NodeExecStatus.completed);

      expect(ctrl.runAutorun(new Set(), SETTINGS), 'started', 'full run starts');
      expect(await until(completedBoth, 10000), true, 'the full run completed');

      const extraCols = (names: string[] | null): number =>
        (names ?? []).filter((n) => n.startsWith('extra')).length;

      expect(shellTable.columns.names().includes('extra'), false,
        'the open shell table was never mutated (the calc worked on a clone)');
      expect(extraCols(previewColNames(ctrl, sel.id)), 0,
        'the upstream preview shows the state at that node');
      expect(extraCols(previewColNames(ctrl, anc.id)), 1, 'the calc preview has its column');

      // The autorun path after a parameter edit: rerun only the calc slice.
      const affected = ctrl.applyGraphEdit({kind: 'params-changed', nodeId: anc.id});
      expect(ctrl.runAutorun(affected, SETTINGS), 'started', 'slice rerun starts');
      expect(await until(() =>
        ctrl.state.getNodeState(anc.id)?.status === NodeExecStatus.completed, 10000), true,
      'the slice rerun completed');

      expect(extraCols(previewColNames(ctrl, anc.id)), 1,
        'rerun is idempotent — the column is NOT added a second time');
      expect(extraCols(previewColNames(ctrl, sel.id)), 0,
        'the upstream captured value stayed pristine through the rerun');

      // The column output now previews the whole in-place table, scrolled to the
      // produced column — captured by __ff_col_summary (by-instance detection).
      const outputs = ctrl.state.getNodeState(anc.id)?.outputs ?? {};
      const colS = Object.values(outputs).find((s) => s.type === 'column');
      expect(colS != null, true, 'the calc has a column output summary');
      expect(colS!.scrollToColumn, 'extra', 'the preview scrolls to the produced column');
      expect(colS!.tableClone instanceof DG.DataFrame, true, 'the whole input table is captured');
      expect((colS!.tableClone as DG.DataFrame).columns.names().includes('extra'), true,
        'the captured table contains the produced column (it belongs to that table)');
    } finally {
      try {
        grok.shell.closeTable(shellTable);
      } catch {/* best effort */}
      destroyEditor(e);
    }
  });

  test('a single-input column output is emitted with in-place table detection', async () => {
    const addCol = funcTypeName('AddNewColumn');
    if (!addCol) {
      expect(true, true);
      return;
    }
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');
      const anc = await addNode(e.flow, addCol);
      await e.flow.addConnectionByKeys(src.id, 'table', anc.id, 'table');

      const inst = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'colsum-1'});
      expect(inst.includes('function __ff_col_summary'), true, 'the helper is defined');
      // The column output is summarized against the same snapshot the call
      // mutates (`<var>_table_in`), so an added column is detectable by instance.
      expect(/__ff_col_summary\(\w+, \w+_table_in, 'column'\)/.test(inst), true,
        'the column output is summarized against its snapshot input table');
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: autorun', () => {
  const edit = (nodeId: string): GraphEdit => ({kind: 'params-changed', nodeId});
  const sleep = (ms: number): Promise<void> => new Promise((r) => setTimeout(r, ms));

  test('debounce coalesces a burst of edits into one run with the union', async () => {
    const runs: Set<string>[] = [];
    const s = new AutorunScheduler((dirty) => {
      runs.push(dirty);
      return 'started';
    }, 20);
    s.toggle();
    s.onEdit(edit('a'), new Set(['a', 'b']));
    s.onEdit(edit('b'), new Set(['b']));
    s.onEdit(edit('c'), new Set(['c']));
    await sleep(80);
    expect(runs.length, 1, 'one run for the whole burst');
    expect([...runs[0]].sort().join(), 'a,b,c', 'dirty sets are unioned');
  });

  test('disabled scheduler never runs; node-added never schedules', async () => {
    let runs = 0;
    const s = new AutorunScheduler(() => {
      runs++;
      return 'started';
    }, 20);
    s.onEdit(edit('a'), new Set(['a'])); // disabled
    s.toggle();
    s.onEdit({kind: 'node-added', nodeId: 'x'}, new Set()); // non-invalidating
    await sleep(60);
    expect(runs, 0);
  });

  test('live-by-default nodes schedule runs even while the toggle is off', async () => {
    const runs: Array<{dirty: string[]; liveOnly: boolean}> = [];
    const s = new AutorunScheduler((dirty, liveOnly) => {
      runs.push({dirty: [...dirty].sort(), liveOnly});
      return 'started';
    }, 20, (id) => id === 'viewer');
    // The toggle stays OFF. A non-live edit schedules nothing…
    s.onEdit(edit('a'), new Set(['a']));
    await sleep(60);
    expect(runs.length, 0, 'non-live edit ignored while off');
    // …an edit whose affected cone touches a live node runs ONLY that node —
    // the rest of the cone must never run uninvited.
    s.onEdit(edit('a'), new Set(['a', 'viewer']));
    await sleep(60);
    expect(runs.length, 1, 'a live node admitted the run');
    expect(runs[0].dirty.join(), 'viewer', 'only the live node is handed over');
    expect(runs[0].liveOnly, true, 'flagged as a live-only run');
    // With the toggle ON the whole affected slice goes through as before.
    s.toggle();
    s.onEdit(edit('a'), new Set(['a', 'b']));
    await sleep(60);
    expect(runs[1].dirty.join(), 'a,b', 'toggle on → full affected slice');
    expect(runs[1].liveOnly, false, 'flagged as a normal autorun');
  });

  test('a dropped live node and a loaded flow schedule live runs (node-added + kickLive)', async () => {
    const runs: Array<{dirty: string[]; liveOnly: boolean}> = [];
    const s = new AutorunScheduler((dirty, liveOnly) => {
      runs.push({dirty: [...dirty].sort(), liveOnly});
      return 'started';
    }, 20, (id) => id.startsWith('live'));
    // Toggle OFF. A plain node drop never schedules…
    s.onEdit({kind: 'node-added', nodeId: 'plain'}, new Set());
    await sleep(60);
    expect(runs.length, 0, 'plain node-added ignored');
    // …a LIVE node drop does (a file dragged onto the canvas creates a ready
    // Open File — readiness is re-checked at fire time anyway).
    s.onEdit({kind: 'node-added', nodeId: 'live1'}, new Set());
    await sleep(60);
    expect(runs.length, 1, 'live node-added schedules');
    expect(runs[0].dirty.join(), 'live1');
    expect(runs[0].liveOnly, true);
    // Loading a flow kicks its live nodes only.
    s.kickLive(['a', 'live2', 'b']);
    await sleep(60);
    expect(runs.length, 2, 'kickLive schedules');
    expect(runs[1].dirty.join(), 'live2', 'only live nodes enter the set');
    // With the toggle ON, node-added keeps its old no-op semantics.
    s.toggle();
    s.onEdit({kind: 'node-added', nodeId: 'live3'}, new Set());
    await sleep(60);
    expect(runs.length, 2, 'toggle on → node-added still never schedules');
  });

  test('runLiveNodes executes only the ready live nodes, never the rest of the canvas', async () => {
    registerBuiltinNodes();
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      const ctrl = new ExecutionController(e.flow);
      // c's input is connected but no upstream value was ever captured →
      // inputs not satisfied → nothing runs (the old behavior ran EVERYTHING).
      expect(ctrl.runLiveNodes(new Set([c]), SETTINGS), 'skipped',
        'unsatisfied inputs never trigger a run');
      expect(ctrl.state.nodeStates.size, 0, 'no node was touched');

      // The const source runs from its parameters alone — and ONLY it runs.
      e.flow.getNodeById(a)!.properties['value'] = 'hello';
      expect(ctrl.runLiveNodes(new Set([a]), SETTINGS), 'started');
      const done = await until(() => status(ctrl, a) === NodeExecStatus.completed, 8000);
      expect(done, true, 'the live node ran');
      expect(status(ctrl, b) === undefined, true, 'downstream b stays untouched');
      expect(status(ctrl, c) === undefined, true, 'downstream c stays untouched');
    } finally {
      destroyEditor(e);
    }
  });

  test('isAutorunByDefault: Open File, Add New Column, and viewers are live', async () => {
    registerBuiltinNodes();
    registerAllFunctions();
    const viewer = createNode('Viewers/Scatter Plot');
    if (viewer) expect(isAutorunByDefault(viewer), true, 'viewers are live');
    for (const name of ['OpenFile', 'AddNewColumn']) {
      const typeName = funcTypeName(name);
      if (!typeName) continue; // function not present on this stand
      expect(isAutorunByDefault(createNode(typeName)!), true, `${name} is live`);
    }
    expect(isAutorunByDefault(createNode('Utilities/Info')!), false, 'ordinary utilities are not live');
    expect(isAutorunByDefault(createNode('Inputs/Table Input')!), false, 'inputs are not live');
  });

  test('busy postpones and keeps the dirty set; skipped waits for the next edit', async () => {
    const runs: Array<{dirty: string[]; outcome: string}> = [];
    let outcome: 'started' | 'busy' | 'skipped' = 'busy';
    const s = new AutorunScheduler((dirty) => {
      runs.push({dirty: [...dirty].sort(), outcome});
      return outcome;
    }, 20);
    s.toggle();
    s.onEdit(edit('a'), new Set(['a']));
    await sleep(30); // fires → busy → rescheduled
    outcome = 'started';
    await sleep(40); // retry succeeds with the same set
    expect(runs.length >= 2, true, 'busy outcome retried');
    expect(runs[runs.length - 1].dirty.join(), 'a', 'the set survived the busy attempt');

    outcome = 'skipped';
    s.onEdit(edit('b'), new Set(['b']));
    const afterSkip = runs.length + 1;
    await sleep(70);
    expect(runs.length, afterSkip, 'skipped is not retried on a timer');
    outcome = 'started';
    s.onEdit(edit('c'), new Set(['c']));
    await sleep(40);
    expect(runs[runs.length - 1].dirty.join(), 'b,c', 'the next edit re-schedules with the kept set');
  });

  test('kick: switching autorun on schedules the pending set without an edit', async () => {
    const runs: Set<string>[] = [];
    const s = new AutorunScheduler((dirty) => {
      runs.push(dirty);
      return 'started';
    }, 20);
    s.kick(new Set(['x']));
    await sleep(50);
    expect(runs.length, 0, 'kick is a no-op while disabled');
    s.toggle();
    s.kick(new Set(['a', 'b']));
    await sleep(60);
    expect(runs.length, 1, 'the pending set runs after the debounce');
    expect([...runs[0]].sort().join(), 'a,b', 'exactly the kicked set');
  });

  test('runAutorun skips when the run would prompt for script inputs', async () => {
    registerBuiltinNodes();
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input');
      const output = await addNode(e.flow, 'Outputs/Table Output');
      await e.flow.addConnectionByKeys(input.id, 'table', output.id, 'table');
      const ctrl = new ExecutionController(e.flow);
      expect(ctrl.runAutorun(new Set([output.id]), SETTINGS), 'skipped',
        'an input node in the run set would open a dialog — autorun must not');
    } finally {
      destroyEditor(e);
    }
  });

  test('runAutorun really executes: full run first, then only the edited slice', async () => {
    registerBuiltinNodes();
    const e = makeEditor();
    try {
      const {a, b, c} = await makeChain(e);
      const ctrl = new ExecutionController(e.flow);

      expect(ctrl.runAutorun(new Set(), SETTINGS), 'started', 'nothing captured yet → full run');
      const ranAll = await until(() =>
        [a, b, c].every((id) => ctrl.state.getNodeState(id)?.status === NodeExecStatus.completed), 8000);
      expect(ranAll, true, 'the full autorun completed');

      // Edit b → {b, c} invalidated; a keeps its result and captured value.
      const affected = ctrl.applyGraphEdit({kind: 'params-changed', nodeId: b});
      const aStateBefore = ctrl.state.getNodeState(a);
      const canSlice = ctrl.hasLiveValue(a, 'value');
      expect(ctrl.runAutorun(affected, SETTINGS), 'started', 'edited slice reruns');
      const ranSlice = await until(() =>
        [b, c].every((id) => ctrl.state.getNodeState(id)?.status === NodeExecStatus.completed), 8000);
      expect(ranSlice, true, 'the slice rerun completed');
      if (canSlice) {
        // preserveState: a's state object was never touched by the slice run.
        expect(ctrl.state.getNodeState(a) === aStateBefore, true, 'upstream was not re-run');
      }
    } finally {
      destroyEditor(e);
    }
  });

  test('hold suspends autorun (editor dialog open); release fires the backlog', async () => {
    const runs: Set<string>[] = [];
    const s = new AutorunScheduler((dirty) => {
      runs.push(dirty);
      return 'started';
    }, 20);
    s.toggle();
    s.hold();
    s.onEdit(edit('a'), new Set(['a']));
    await sleep(60);
    expect(runs.length, 0, 'held — nothing fires');
    s.hold(); // re-entrant: a second holder
    s.release();
    await sleep(60);
    expect(runs.length, 0, 'still one hold outstanding');
    s.release();
    await sleep(60);
    expect(runs.length, 1, 'the backlog fired after the last release');
    expect([...runs[0]].join(), 'a', 'with the accumulated set');

    // A hold placed while the debounce timer is already pending cancels it.
    s.onEdit(edit('b'), new Set(['b']));
    s.hold();
    await sleep(60);
    expect(runs.length, 1, 'pending timer canceled by hold');
    s.release();
    await sleep(60);
    expect(runs.length, 2, 'released — the kept set ran');
    expect([...runs[1]].join(), 'b');
  });

  test('turning autorun off cancels the pending run and clears the set', async () => {
    let runs = 0;
    const s = new AutorunScheduler(() => {
      runs++;
      return 'started';
    }, 20);
    s.toggle();
    s.onEdit(edit('a'), new Set(['a']));
    s.toggle(); // off before the debounce fires
    await sleep(60);
    expect(runs, 0);
  });
});

category('Flow: run readiness', () => {
  const SET = {name: 'ReadinessFlow', description: '', tags: ['funcflow']};

  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('custom nodes declare their required inputs and properties', async () => {
    const selCol = createNode('Utilities/Select Column')!;
    expect(selCol.requiredInputs.includes('table'), true, 'Select Column needs a table');
    expect(selCol.requiredProps.includes('columnName'), true, 'Select Column needs a column name');

    const selCols = createNode('Utilities/Select Columns')!;
    expect(selCols.requiredInputs.includes('table'), true, 'Select Columns needs a table');
    expect(selCols.requiredProps.includes('columnNames'), true, 'Select Columns needs column names');

    const selTable = createNode('Utilities/Select Table')!;
    expect(selTable.requiredProps.includes('tableName'), true, 'Select Table needs a table name');

    const atv = createNode('Utilities/Add Table View')!;
    expect(atv.requiredInputs.includes('table'), true, 'Add Table View needs a table');

    const scatter = createNode('Viewers/Scatter Plot')!;
    expect(scatter.requiredInputs.includes('table'), true, 'a plot needs a table');
  });

  test('missing requirements combine unset inputs and unset properties', async () => {
    const selTable = createNode('Utilities/Select Table')!;
    expect(missingRequiredProps(selTable).includes('tableName'), true, 'empty table name is missing');
    expect(nodeMissingRequirements(selTable, () => false).length > 0, true, 'the node is not ready');
    selTable.properties['tableName'] = 'demog';
    expect(missingRequiredProps(selTable).length, 0, 'a set table name satisfies the requirement');
    expect(nodeMissingRequirements(selTable, () => false).length, 0, 'the node is ready');
  });

  test('runnableNodes excludes an unready node and everything downstream of it', async () => {
    const e = makeEditor();
    try {
      // Unready Select Table (no name) → Add Table View (its own table requirement
      // is met by the connection, yet it is downstream of an unready node).
      const sel = await addNode(e.flow, 'Utilities/Select Table');
      const atv = await addNode(e.flow, 'Utilities/Add Table View', 300, 0);
      await e.flow.addConnectionByKeys(sel.id, 'table', atv.id, 'table');
      const konst = await addNode(e.flow, 'Constants/String', 0, 300);

      const ctrl = new ExecutionController(e.flow);
      let runnable = ctrl.runnableNodes();
      expect(runnable.has(sel.id), false, 'the unready Select Table is excluded');
      expect(runnable.has(atv.id), false, 'the node downstream of it is excluded too');
      expect(runnable.has(konst.id), true, 'an unrelated ready node still runs');

      // Filling the required property makes the whole chain runnable.
      sel.properties['tableName'] = 'demog';
      runnable = ctrl.runnableNodes();
      expect(runnable.has(sel.id) && runnable.has(atv.id) && runnable.has(konst.id), true,
        'once the requirement is set, the chain is runnable');
    } finally {
      destroyEditor(e);
    }
  });

  test('a plot without a table is excluded from the run set', async () => {
    const e = makeEditor();
    try {
      const scatter = await addNode(e.flow, 'Viewers/Scatter Plot');
      const konst = await addNode(e.flow, 'Constants/String', 0, 300);
      const ctrl = new ExecutionController(e.flow);
      const runnable = ctrl.runnableNodes();
      expect(runnable.has(scatter.id), false, 'the table-less plot is not run');
      expect(runnable.has(konst.id), true, 'the ready node is');
    } finally {
      destroyEditor(e);
    }
  });

  test('autorun skips entirely when nothing is ready', async () => {
    const e = makeEditor();
    try {
      await addNode(e.flow, 'Viewers/Scatter Plot'); // no table → unready
      const ctrl = new ExecutionController(e.flow);
      expect(ctrl.runAutorun(new Set(), SET), 'skipped', 'no runnable node → nothing autoruns');
    } finally {
      destroyEditor(e);
    }
  });

  test('autorun runs the ready part and never runs an unready node', async () => {
    const e = makeEditor();
    const df = DG.DataFrame.fromCsv('x\n1\n2\n');
    df.name = 'ffReadiness';
    const shellTable = grok.shell.addTable(df);
    try {
      const sel = await addNode(e.flow, 'Utilities/Select Table');
      sel.properties['tableName'] = 'ffReadiness';
      // Add Table View with NO table: without the readiness gate this would emit
      // `addTableView(undefined)` and fault at run time — the gate must keep it
      // out of the run entirely (it stays idle, never errors).
      const atv = await addNode(e.flow, 'Utilities/Add Table View', 300, 0);

      const ctrl = new ExecutionController(e.flow);
      expect(ctrl.runnableNodes().has(atv.id), false, 'the table-less Add Table View is not runnable');
      expect(ctrl.runAutorun(new Set(), SET), 'started', 'the ready Select Table triggers a run');
      expect(await until(() =>
        ctrl.state.getNodeState(sel.id)?.status === NodeExecStatus.completed, 10000), true,
      'the ready node ran to completion');
      const atvStatus = ctrl.state.getNodeState(atv.id)?.status;
      expect(atvStatus === NodeExecStatus.completed || atvStatus === NodeExecStatus.errored, false,
        'the unready node was never run — no error, no completion');
    } finally {
      try {
        grok.shell.closeTable(shellTable);
      } catch {/* best effort */}
      destroyEditor(e);
    }
  });
});
