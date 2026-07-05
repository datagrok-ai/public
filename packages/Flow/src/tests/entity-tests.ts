/** First-class Flow entity: the `.flow` body format, deterministic loading,
 *  the context-panel widget, and the live save → run → round-trip path
 *  (Script with language 'flow' executed through the package script handler). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {ensureFunctionsRegistered} from '../rete/node-factory';
import {flowScriptText, parseFlowBody, isFlowBody, FLOW_LANGUAGE} from '../serialization/flow-script-format';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {FlowEntityHandler} from '../entity/flow-entity-handler';
import {FuncFlowView} from '../funcflow-view';
import {makeEditor, destroyEditor, addNode, TestEditor} from './test-utils';

const SETTINGS = {scriptName: 'EntityFlow', scriptDescription: 'entity test', tags: ['funcflow']};

/** Int Input('a', default 3) → Value Output('result') — the smallest runnable flow. */
async function buildPassThroughFlow(e: TestEditor): Promise<void> {
  const input = await addNode(e.flow, 'Inputs/Int Input', 10, 20);
  input.properties['paramName'] = 'a';
  input.properties['defaultValue'] = 3;
  const output = await addNode(e.flow, 'Outputs/Value Output', 340, 20);
  output.properties['paramName'] = 'result';
  output.properties['outputType'] = 'int';
  await e.flow.addConnectionByKeys(input.id, 'value', output.id, 'value');
}

category('Flow: entity format', () => {
  before(async () => {
    ensureFunctionsRegistered();
  });

  test('flowScriptText emits header + json, parseFlowBody round-trips', async () => {
    const e = makeEditor();
    try {
      await buildPassThroughFlow(e);
      const text = flowScriptText(e.flow, SETTINGS);

      expect(text.includes('//name: EntityFlow'), true);
      expect(text.includes(`//language: ${FLOW_LANGUAGE}`), true);
      expect(text.includes('//input: int a = 3'), true);
      expect(text.includes('//output: int result'), true);
      expect(text.includes('//tags: funcflow, flow'), true); // flow tag auto-added

      const {header, doc} = parseFlowBody(text);
      expect(header.startsWith('//name: EntityFlow'), true);
      expect(doc.version, '2.0');
      expect(doc.nodes.length, 2);
      expect(doc.connections.length, 1);
      expect(doc.metadata?.settings?.scriptName, 'EntityFlow');

      // The parsed doc reloads into an editor with the same topology.
      await deserializeFlow(doc, e.flow);
      expect(e.flow.getNodeCount(), 2);
      expect(e.flow.getConnectionCount(), 1);
    } finally {
      destroyEditor(e);
    }
  });

  test('header carries input qualifiers and captions', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Int Input', 0, 0);
      input.properties['paramName'] = 'threshold';
      input.properties['min'] = 0;
      input.properties['max'] = 100;
      input.properties['caption'] = 'Cutoff';
      const text = flowScriptText(e.flow, SETTINGS);
      const line = text.split('\n').find((l) => l.startsWith('//input: int threshold'));
      expect(line != null, true);
      expect(line!.includes('min: 0'), true);
      expect(line!.includes('max: 100'), true);
      expect(line!.includes('caption: Cutoff'), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('isFlowBody rejects plain scripts and accepts flow bodies', async () => {
    expect(isFlowBody('//name: X\n//language: javascript\nlet a = 1;'), false);
    expect(isFlowBody(''), false);
    const e = makeEditor();
    try {
      await buildPassThroughFlow(e);
      expect(isFlowBody(flowScriptText(e.flow, SETTINGS)), true);
    } finally {
      destroyEditor(e);
    }
  });

  test('widget renders facts for a flow script entity', async () => {
    const e = makeEditor();
    try {
      await buildPassThroughFlow(e);
      const script = DG.Script.create(flowScriptText(e.flow, SETTINGS));
      const w = FlowEntityHandler.instance.widget(script);
      expect(w.root.textContent!.includes('2 nodes'), true);
      expect(w.root.textContent!.includes('Open editor'), true);
    } finally {
      destroyEditor(e);
    }
  });
});

category('Flow: entity loading', () => {
  test('loadFromJson right after construction is race-free', async () => {
    const view = new FuncFlowView();
    try {
      const e = makeEditor();
      let json: string;
      try {
        await buildPassThroughFlow(e);
        json = JSON.stringify(serializeFlow(e.flow, SETTINGS));
      } finally {
        destroyEditor(e);
      }
      // No waiting, no timers: loadFromJson awaits editor readiness itself.
      await view.loadFromJson(json);
      const flow = (view as any).flow;
      expect(flow.getNodeCount(), 2);
      expect(flow.getConnectionCount(), 1);
    } finally {
      ((view as any).flow)?.destroy?.();
      view.root.remove();
    }
  }, {timeout: 30000});

  test('compileToJs produces a runnable header on a detached editor', async () => {
    const e = makeEditor();
    try {
      await buildPassThroughFlow(e);
      const {doc} = parseFlowBody(flowScriptText(e.flow, SETTINGS));
      const js = await FlowEntityHandler.instance.compileToJs(doc);
      expect(js.includes('//language: javascript'), true);
      expect(js.includes('//input: int a = 3'), true);
      expect(js.includes('//output: int result'), true);
      expect(js.includes('result = a;'), true);
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});

  // Output nodes keep paramName as the assigned variable (no label camel-casing),
  // so the `//output:` header, the body assignment, and the copy key in
  // FlowEntityHandler.run() are one identifier — even for a non-'result' name.
  test('output paramName drives both the header and the body assignment', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Int Input', 10, 20);
      input.properties['paramName'] = 'a';
      input.properties['defaultValue'] = 3;
      const output = await addNode(e.flow, 'Outputs/Value Output', 340, 20);
      output.properties['paramName'] = 'netResult';
      output.properties['outputType'] = 'int';
      await e.flow.addConnectionByKeys(input.id, 'value', output.id, 'value');
      const {doc} = parseFlowBody(flowScriptText(e.flow, SETTINGS));
      const js = await FlowEntityHandler.instance.compileToJs(doc);
      expect(js.includes('//output: int netResult'), true);
      expect(js.includes('netResult = a;'), true);
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});
});

category('Flow: entity server round-trip', () => {
  before(async () => {
    ensureFunctionsRegistered();
  });

  test('save → find → run via script handler → delete', async () => {
    const e = makeEditor();
    let saved: DG.Script | null = null;
    try {
      await buildPassThroughFlow(e);
      const name = `FFEntityTest${Math.floor(Math.random() * 1e6)}`;
      const text = flowScriptText(e.flow, {...SETTINGS, scriptName: name});
      const script = DG.Script.create(text);
      expect(script.language as string, FLOW_LANGUAGE);

      saved = await grok.dapi.scripts.save(script);
      expect(saved != null, true);

      const found = await grok.dapi.scripts.find(saved!.id);
      expect(found != null, true);
      expect(found!.language as string, FLOW_LANGUAGE);
      const {doc} = parseFlowBody(found!.script);
      expect(doc.nodes.length, 2);
      expect(doc.metadata?.settings?.scriptName, name);
      expect(found!.inputs.length, 1);
      expect(found!.inputs[0].name, 'a');

      // Execute the entity like any function — resolves the 'flow' language
      // through the package script handler, compiles, runs, copies outputs.
      const result = await grok.functions.call(found!.nqName, {a: 5});
      expect(result, 5);

      // Default value path: prepare() with no explicit inputs.
      const call = (found! as DG.Func).prepare({a: 7});
      await call.call(undefined, undefined, {processed: true});
      expect(call.outputs['result'], 7);
    } finally {
      destroyEditor(e);
      if (saved != null)
        await grok.dapi.scripts.delete(saved).catch(() => {});
    }
  }, {timeout: 120000});
});

// The exact server contract the SpacePicker and Save As dialog rely on:
// root/subspace creation, subspace listing, and binding a flow into a space.
category('Flow: space binding', () => {
  before(async () => {
    ensureFunctionsRegistered();
  });

  test('create space/subspace, save flow into it, list children', async () => {
    const stamp = Math.floor(Math.random() * 1e6);
    const rootName = `FFSpace${stamp}`;
    const subName = `Sub${stamp}`;
    let root: DG.Project | null = null;
    let sub: DG.Project | null = null;
    let saved: DG.Script | null = null;
    const e = makeEditor();
    try {
      await buildPassThroughFlow(e);

      root = await grok.dapi.spaces.createRootSpace(rootName);
      expect(root != null, true);
      sub = await grok.dapi.spaces.id(root!.id).addSubspace(subName);
      expect(await grok.dapi.spaces.id(root!.id).subspaceExists(subName), true);
      // What the picker's lazy loader runs on expand.
      const children = await grok.dapi.spaces.id(root!.id).children.filter('Project', false).list();
      expect(children.some((c) => c.id === sub!.id), true);

      // What Save As does: save the entity, then bind it into the space.
      const name = `FFSpaceFlow${stamp}`;
      saved = await grok.dapi.scripts.save(DG.Script.create(flowScriptText(e.flow, {...SETTINGS, scriptName: name})));
      await grok.dapi.spaces.id(sub!.id).addEntity(saved!.id, false);
      const found = await grok.dapi.scripts.find(saved!.id);
      // Ownership moved: the namespace is rewritten to the subspace path.
      expect(found!.nqName.toLowerCase(), `${sub!.nqName}:${found!.name}`.toLowerCase());

      // A flow living in a space stays runnable through the script handler.
      const result = await grok.functions.call(found!.nqName, {a: 4});
      expect(result, 4);
    } finally {
      destroyEditor(e);
      if (saved != null)
        await grok.dapi.scripts.delete(saved).catch(() => {});
      if (sub != null)
        await grok.dapi.spaces.delete(sub).catch(() => {});
      if (root != null)
        await grok.dapi.spaces.delete(root).catch(() => {});
    }
  }, {timeout: 120000});
});
