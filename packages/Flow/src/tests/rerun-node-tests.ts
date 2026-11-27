/** "Rerun this node only": a single node re-executes from values captured by a
 *  prior run (upstream inputs resolve to `_ffLive(...)` registry reads). */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {emitScript} from '../compiler/script-emitter';
import {ExecutionController} from '../execution/execution-controller';
import {FlowEditor} from '../rete/flow-editor';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';
import * as ui from 'datagrok-api/ui';

const SETTINGS = {name: 'T', description: '', tags: []};

function funcTypeName(name: string): string | null {
  return getRegisteredFuncs().find((f) => f.func.name === name)?.nodeTypeName ?? null;
}

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

category('Flow: rerun node', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('a full instrumented run stashes each node\'s live outputs', async () => {
    const addCol = funcTypeName('AddNewColumn');
    if (!addCol) return;
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');
      const anc = await addNode(e.flow, addCol);
      await e.flow.addConnectionByKeys(src.id, 'table', anc.id, 'table');

      const full = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'r1'});
      expect(full.includes('__ff_stash('), true, 'stashes live values');
      expect(full.includes('"table__pt":'), true, 'passthrough table stashed');
      expect(full.includes(`__ff_stash(${JSON.stringify(src.id)}`), true, 'the table input is stashed too');
    } finally {
      destroyEditor(e);
    }
  });

  test('re-running one node reads upstream values from the registry, not a re-run', async () => {
    const addCol = funcTypeName('AddNewColumn');
    if (!addCol) return;
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');
      const anc = await addNode(e.flow, addCol);
      await e.flow.addConnectionByKeys(src.id, 'table', anc.id, 'table');

      const rerun = emitScript(e.flow, SETTINGS,
        {instrumented: true, runId: 'r2', onlyNodeIds: new Set([anc.id]), liveExternalInputs: true});
      expect(rerun.includes(`_ffLive(${JSON.stringify(src.id)}, "table")`), true,
        'the connected table resolves to a registry read');
      expect(rerun.includes('grok.functions.call(\'AddNewColumn\''), true, 'the node itself runs');
      expect(/^\/\/input:/m.test(rerun), false, 'upstream is not re-declared as a script input');
      expect(rerun.includes('function _ffLive'), true, 'the registry accessor is defined');
    } finally {
      destroyEditor(e);
    }
  });

  test('re-running a viewer wired to a passthrough reads the modified table', async () => {
    const addCol = funcTypeName('AddNewColumn');
    if (!addCol) return;
    const e = makeEditor();
    try {
      const src = await addNode(e.flow, 'Inputs/Table Input');
      const anc = await addNode(e.flow, addCol);
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      await e.flow.addConnectionByKeys(src.id, 'table', anc.id, 'table');
      await e.flow.addConnectionByKeys(anc.id, 'table__pt', viewer.id, 'table');

      const rerun = emitScript(e.flow, SETTINGS,
        {instrumented: true, runId: 'r3', onlyNodeIds: new Set([viewer.id]), liveExternalInputs: true});
      expect(rerun.includes(`_ffLive(${JSON.stringify(anc.id)}, "table__pt")`), true,
        'the viewer plots the passthrough table from the registry');
      expect(rerun.includes('.plot.fromType('), true, 'the viewer is created');
    } finally {
      destroyEditor(e);
    }
  });

  test('canRerunNode: needs required inputs AND captured connection values', async () => {
    const addCol = funcTypeName('AddNewColumn');
    if (!addCol) return;
    const e = makeEditor();
    try {
      const ctrl = new ExecutionController(e.flow);

      const lone = await addNode(e.flow, addCol);
      expect(ctrl.canRerunNode(lone.id), false, 'missing required table input');

      const src = await addNode(e.flow, 'Inputs/Table Input');
      const anc = await addNode(e.flow, addCol);
      await e.flow.addConnectionByKeys(src.id, 'table', anc.id, 'table');

      withRegistry({}, () => {
        expect(ctrl.canRerunNode(anc.id), false, 'no captured upstream value');
      });
      withRegistry({[src.id]: {table: {}}}, () => {
        expect(ctrl.canRerunNode(anc.id), false, 'blank required params (name, expression) block the rerun');
        anc.inputValues['name'] = 'x2';
        anc.inputValues['expression'] = '${x} * 2';
        expect(ctrl.canRerunNode(anc.id), true, 'ready to re-run');
      });
      withRegistry({[src.id]: {table: {}}}, () => {
        expect(ctrl.canRerunNode(src.id), false, 'input nodes are not rerunnable');
      });
    } finally {
      destroyEditor(e);
    }
  });

  test('Run from here: the node and its downstream only, other branches read from what was captured', async () => {
    const addCol = funcTypeName('AddNewColumn');
    const join = funcTypeName('JoinTables');
    if (!addCol || !join) return;
    const e = makeEditor();
    try {
      // two branches meeting at a join: a → fa ↘ join ↙ fb ← b
      const a = await addNode(e.flow, 'Inputs/Table Input');
      const b = await addNode(e.flow, 'Inputs/Table Input');
      const fa = await addNode(e.flow, addCol);
      const fb = await addNode(e.flow, addCol);
      const j = await addNode(e.flow, join);
      await e.flow.addConnectionByKeys(a.id, 'table', fa.id, 'table');
      await e.flow.addConnectionByKeys(b.id, 'table', fb.id, 'table');
      await e.flow.addConnectionByKeys(fa.id, 'table__pt', j.id, 'table1');
      await e.flow.addConnectionByKeys(fb.id, 'table__pt', j.id, 'table2');

      const ctrl = new ExecutionController(e.flow);
      const runs: Array<{onlyNodeIds?: Set<string>; liveExternalInputs?: boolean; preserveState?: boolean}> = [];
      (ctrl as unknown as {executeInstrumented: (...args: unknown[]) => void}).executeInstrumented =
        (_s, _d, opts) => {runs.push(opts as typeof runs[number]);};
      const everything = {[a.id]: {table: {}}, [b.id]: {table: {}}, [fa.id]: {table__pt: {}},
        [fb.id]: {table__pt: {}}, [j.id]: {result: {}}};

      withRegistry(everything, () => ctrl.runDownstream(b.id, SETTINGS));
      const ran = [...runs[0].onlyNodeIds!].sort();
      expect(ran.join(','), [b.id, fb.id, j.id].sort().join(','), 'the second input, its branch and the join');
      expect(runs[0].liveExternalInputs && runs[0].preserveState, true, 'the first branch is read, not rerun');

      withRegistry({[b.id]: {table: {}}}, () => ctrl.runDownstream(fb.id, SETTINGS));
      expect(runs[1].onlyNodeIds!.has(fa.id) && runs[1].onlyNodeIds!.has(a.id), true,
        'an upstream result never captured is computed first');
      expect(runs[1].onlyNodeIds!.has(b.id), false, 'a captured one is not');
    } finally {
      destroyEditor(e);
    }
  });

  test('the node context menu offers Run from here, which runs that node', async () => {
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    const asked: string[] = [];
    const flow = new FlowEditor(container, {onRunDownstream: (id) => asked.push(id)});
    try {
      const node = await addNode(flow, 'Inputs/String Input', 0, 0);
      await until(() => !!container.querySelector(`.ff-node[data-node-id="${node.id}"]`));
      const nodeEl = container.querySelector<HTMLElement>(`.ff-node[data-node-id="${node.id}"]`)!;
      nodeEl.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 10, clientY: 10}));
      const item = (): HTMLElement | undefined =>
        Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'Run from here');
      expect(await until(() => item() != null), true, 'the item is offered on an input node');
      item()!.click();
      expect(await until(() => asked.length > 0), true, 'clicking the item asks for a run');
      expect(asked.join(','), node.id, 'from that node');
    } finally {
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();
      flow.destroy();
      container.remove();
    }
  });
});
