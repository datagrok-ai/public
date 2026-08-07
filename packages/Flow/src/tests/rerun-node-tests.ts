/** "Rerun this node only": a single node re-executes from values captured by a
 *  prior run (upstream inputs resolve to `_ffLive(...)` registry reads). */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {emitScript} from '../compiler/script-emitter';
import {ExecutionController} from '../execution/execution-controller';
import {makeEditor, destroyEditor, addNode} from './test-utils';

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
});
