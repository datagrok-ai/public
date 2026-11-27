/** Table links made in a flow: the tables a link node passes on, its preview and the flow's outputs are the
 *  same linked instances; a rerun of the node replaces its link, deleting the node removes it. */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs} from '../rete/node-factory';
import {emitScript} from '../compiler/script-emitter';
import {ExecutionController} from '../execution/execution-controller';
import {NodeExecStatus} from '../execution/execution-state';
import {FlowNode} from '../rete/scheme';
import {makeEditor, destroyEditor, addNode, until, TestEditor} from './test-utils';

const SETTINGS = {name: 'LinkTest', description: '', tags: []};

function linkTablesType(): string | null {
  return getRegisteredFuncs().find((f) => f.func.nqName.toLowerCase() === 'core:linktables')?.nodeTypeName ?? null;
}

/** parents → LinkTables (selection to filter on `k`) → two table outputs. */
async function linkedFlow(e: TestEditor, typeName: string): Promise<{
  link: FlowNode; out1: FlowNode; out2: FlowNode; parents: DG.DataFrame; children: DG.DataFrame;
}> {
  const parents = DG.DataFrame.fromCsv('k,v\na,1\nb,2\nc,3\n');
  const children = DG.DataFrame.fromCsv('k,w\na,1\na,2\nb,3\nc,4\nc,5\nc,6\n');
  parents.name = 'ff link parents';
  children.name = 'ff link children';
  const in1 = await addNode(e.flow, 'Inputs/Table Input');
  in1.transientValue = parents;
  in1.properties['paramName'] = 'parents';
  const in2 = await addNode(e.flow, 'Inputs/Table Input');
  in2.transientValue = children;
  in2.properties['paramName'] = 'children';
  const link = await addNode(e.flow, typeName);
  link.inputValues['keys1'] = 'k';
  link.inputValues['keys2'] = 'k';
  link.inputValues['linkType'] = DG.SYNC_TYPE.SELECTION_TO_FILTER;
  await e.flow.addConnectionByKeys(in1.id, 'table', link.id, 'table1');
  await e.flow.addConnectionByKeys(in2.id, 'table', link.id, 'table2');
  const out1 = await addNode(e.flow, 'Outputs/Table Output');
  out1.properties['paramName'] = 'linkedParents';
  const out2 = await addNode(e.flow, 'Outputs/Table Output');
  out2.properties['paramName'] = 'linkedChildren';
  await e.flow.addConnectionByKeys(link.id, 'table1__pt', out1.id, 'table');
  await e.flow.addConnectionByKeys(link.id, 'table2__pt', out2.id, 'table');
  return {link, out1, out2, parents, children};
}

function dataframeSummary(ctrl: ExecutionController, nodeId: string, key?: string): DG.DataFrame | null {
  const outputs = ctrl.state.getNodeState(nodeId)?.outputs ?? {};
  const summary = key ? outputs[key] : Object.values(outputs).find((s) => s?.type === 'dataframe');
  return summary?.clone instanceof DG.DataFrame ? summary.clone : null;
}

const same = (a: DG.DataFrame | null | undefined, b: DG.DataFrame | null | undefined): boolean =>
  a != null && b != null && a.dart === b.dart;

category('Flow: table links', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('a flow run links through grok.data.linkTables and keeps the link; a script calls LinkTables', async () => {
    const typeName = linkTablesType();
    if (!typeName) return;
    const e = makeEditor();
    try {
      const {link} = await linkedFlow(e, typeName);
      const run = emitScript(e.flow, SETTINGS, {instrumented: true, runId: 'link1'});
      expect(run.includes('grok.data.linkTables('), true, 'the run links with the call that returns the link');
      expect(run.includes(`__ff_stash(${JSON.stringify(link.id)}, {`) && run.includes('"__link":'), true,
        'the link is kept with the node\'s values');
      const script = emitScript(e.flow, SETTINGS, {});
      expect(script.includes('grok.functions.call(\'LinkTables\''), true, 'a plain script calls the function');
    } finally {
      destroyEditor(e);
    }
  });

  test('outputs and previews are the linked instances; a rerun replaces the link, deleting the node removes it',
    async () => {
      const typeName = linkTablesType();
      if (!typeName) return;
      const e = makeEditor();
      const ctrl = new ExecutionController(e.flow);
      try {
        const {link, out1, out2} = await linkedFlow(e, typeName);
        const ids = [link.id, out1.id, out2.id];
        const done = (): boolean => ids.every((id) => ctrl.state.getNodeState(id)?.status === NodeExecStatus.completed);
        ctrl.runDownstream(link.id, SETTINGS);
        expect(await until(done, 15_000), true, 'the flow ran');

        const p = ctrl.liveValue(link.id, 'table1__pt') as DG.DataFrame;
        const c = ctrl.liveValue(link.id, 'table2__pt') as DG.DataFrame;
        expect(same(dataframeSummary(ctrl, out1.id), p) && same(dataframeSummary(ctrl, out2.id), c), true,
          'each output is the instance the link node passed on, not a copy');
        expect(same(dataframeSummary(ctrl, link.id, 'table1 (modified)'), p), true,
          'the link node\'s preview shows that instance too');
        p.selection.set(0, true);
        expect(await until(() => c.filter.trueCount === 2), true, 'selecting a parent filters its children');

        const oldP = p;
        const oldC = c;
        oldP.selection.setAll(false);
        await until(() => oldC.filter.trueCount === oldC.rowCount);
        ctrl.runDownstream(link.id, SETTINGS);
        expect(await until(() => ctrl.liveValue(link.id, 'table1__pt') !== oldP && done() && !ctrl.state.isRunning,
          15_000), true, 'the link node ran again');
        const newP = ctrl.liveValue(link.id, 'table1__pt') as DG.DataFrame;
        const newC = ctrl.liveValue(link.id, 'table2__pt') as DG.DataFrame;
        oldP.selection.set(2, true);
        await new Promise((r) => setTimeout(r, 200));
        expect(oldC.filter.trueCount, oldC.rowCount, 'the old link is gone: the old children no longer follow');
        newP.selection.set(2, true);
        expect(await until(() => newC.filter.trueCount === 3), true, 'the new tables are linked');

        newP.selection.setAll(false);
        await until(() => newC.filter.trueCount === newC.rowCount);
        ctrl.applyGraphEdit({kind: 'node-removed', nodeId: link.id});
        newP.selection.set(0, true);
        await new Promise((r) => setTimeout(r, 200));
        expect(newC.filter.trueCount, newC.rowCount, 'deleting the link node removed its link');
      } finally {
        ctrl.dispose();
        destroyEditor(e);
      }
    });
});
