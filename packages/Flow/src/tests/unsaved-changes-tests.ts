/** Unsaved changes (#4677): what counts as one, and what asks before they are lost. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {ensureFunctionsRegistered, ensureFuncNodeType} from '../rete/node-factory';
import {flowScriptText, parseFlowBody} from '../serialization/flow-script-format';
import {serializeFlow} from '../serialization/flow-serializer';
import {FlowSettings, FuncFlowDocument} from '../serialization/flow-schema';
import {FuncFlowView} from '../funcflow-view';
import {_package} from '../package';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

interface Probe {
  flow: FlowEditor;
  editorReady: Promise<void>;
  boundScript: DG.Script | null;
  flowSettings: FlowSettings;
  propertyPanel: {showNode(node: FlowNode): void};
  newFlow(): Promise<void>;
  openRecentFlow(f: DG.Script): Promise<void>;
  loadTemplate(file: string): Promise<void>;
  importCreationScript(script: string): Promise<void>;
}

interface OpenView {
  view: FuncFlowView;
  probe: Probe;
}

async function openView(): Promise<OpenView> {
  const view = new FuncFlowView();
  const probe = view as unknown as Probe;
  await probe.editorReady;
  return {view, probe};
}

function closeView(v: OpenView): void {
  for (const title of ['Unsaved changes', 'Save Flow'])
    dialog(title)?.close();
  v.probe.flow?.destroy();
  v.view.root.remove();
}

const dialog = (title: string): DG.Dialog | undefined => DG.Dialog.getOpenDialogs().find((d) => d.title === title);

async function answer(title: string, button: string): Promise<void> {
  if (!await until(() => dialog(title) != null, 10000))
    throw new Error(`The "${title}" dialog did not open`);
  const buttons = Array.from(dialog(title)!.root.querySelectorAll('button'));
  buttons.find((b) => b.textContent!.trim().toLowerCase() === button.toLowerCase())!.click();
  expect(await until(() => dialog(title) == null, 10000), true, `"${title}" closes on ${button}`);
}

/** A flow of a function with a description of its own, an input and an output. */
async function flowDocument(): Promise<FuncFlowDocument> {
  const e = makeEditor();
  try {
    await buildFlow(e.flow);
    return serializeFlow(e.flow, settings('FFUnsaved'));
  } finally {
    destroyEditor(e);
  }
}

async function buildFlow(flow: FlowEditor): Promise<void> {
  const func = DG.Func.find({name: 'AddNewColumn'})[0];
  await addNode(flow, ensureFuncNodeType(func), 0, 200);
  const input = await addNode(flow, 'Inputs/Int Input', 10, 20);
  input.properties['defaultValue'] = 3;
  const output = await addNode(flow, 'Outputs/Value Output', 340, 20);
  output.properties['outputType'] = 'int';
  await flow.addConnectionByKeys(input.id, 'value', output.id, 'value');
}

const settings = (name: string): FlowSettings => ({scriptName: name, scriptDescription: '', tags: ['funcflow']});

async function saveScript(): Promise<DG.Script> {
  const e = makeEditor();
  try {
    await buildFlow(e.flow);
    return await grok.dapi.scripts.save(DG.Script.create(
      flowScriptText(e.flow, settings(`FFUnsaved${Math.floor(Math.random() * 1e6)}`))));
  } finally {
    destroyEditor(e);
  }
}

async function forScript(script: DG.Script): Promise<OpenView> {
  const view = FuncFlowView.forScript(script);
  const probe = view as unknown as Probe;
  if (!await until(() => probe.flow?.getNodeCount() === 3 && view.name === script.friendlyName, 15000))
    throw new Error('The saved flow did not load');
  return {view, probe};
}

/** The default value the saved flow's input node has on the server. */
async function savedValue(script: DG.Script): Promise<number> {
  const {doc} = parseFlowBody((await grok.dapi.scripts.find(script.id)).script);
  return Number(doc.nodes.find((n) => n.typeName === 'Inputs/Int Input')!.properties['defaultValue']);
}

const isChanged = (v: OpenView): boolean => v.view.unsavedWork !== null;
const nodeOf = (v: OpenView, type: string): FlowNode => v.probe.flow.getNodes().find((n) => n.dgTypeName === type)!;

function edit(v: OpenView): void {
  const input = nodeOf(v, 'Inputs/Int Input');
  input.properties['defaultValue'] = Number(input.properties['defaultValue']) + 1;
  v.probe.flow.notifyNodeParamsChanged(input.id);
}

category('Flow: unsaved changes', () => {
  before(async () => {
    ensureFunctionsRegistered();
  });

  test('a loaded flow has none, whatever the user looks at; an edit is one until it is undone', async () => {
    const v = await openView();
    try {
      expect(isChanged(v), false, 'an empty canvas');
      await v.view.loadFromDoc(await flowDocument());
      expect(isChanged(v), false, 'a loaded flow');
      for (const node of v.probe.flow.getNodes()) {
        v.probe.propertyPanel.showNode(node);
        await new Promise((r) => setTimeout(r, 30));
        expect(isChanged(v), false, `the context panel of "${node.label}"`);
      }

      const added = await addNode(v.probe.flow, 'Constants/String', 100, 400);
      expect(v.view.unsavedWork!.changes.join(), 'Edited the flow');
      await v.probe.flow.removeNode(added.id);
      expect(isChanged(v), false, 'the flow is as it was loaded');

      const input = nodeOf(v, 'Inputs/Int Input');
      const {x, y} = input.pos;
      await v.probe.flow.area.translate(input.id, {x: x + 130, y: y + 60});
      expect(isChanged(v), true, 'a moved node');
      await v.probe.flow.area.translate(input.id, {x, y});
      expect(isChanged(v), false, 'a node moved back');

      edit(v);
      expect(isChanged(v), true, 'a parameter edit');
      await v.probe.flow.clear();
      expect(isChanged(v), false, 'an empty canvas has nothing to save');
    } finally {
      closeView(v);
    }
  }, {timeout: 60000});

  test('saving the unsaved work of a saved flow updates it; a cancelled save keeps the changes', async () => {
    const script = await saveScript();
    const v = await forScript(script);
    try {
      expect(isChanged(v), false, 'a flow opened from the platform');
      edit(v);
      const work = v.view.unsavedWork!;
      expect(work.saveLink, 'save the flow');
      expect(await work.save(), true);
      expect(isChanged(v), false, 'saved');
      expect(await savedValue(script), 4, 'the saved flow has the edit');

      await grok.dapi.scripts.delete(script);
      edit(v);
      const saved = v.view.unsavedWork!.save(); // gone from the server: asks for a name, like a flow never saved
      await answer('Save Flow', 'Cancel');
      expect(await saved, false, 'a cancelled save');
      expect(isChanged(v), true);
    } finally {
      closeView(v);
      await grok.dapi.scripts.delete(script).catch(() => {});
    }
  }, {timeout: 60000});

  test('a new flow asks before it replaces changes: cancel, save or discard', async () => {
    const script = await saveScript();
    const v = await forScript(script);
    try {
      edit(v);
      let replaced = v.probe.newFlow();
      await answer('Unsaved changes', 'Cancel');
      await replaced;
      expect(v.probe.flow.getNodeCount(), 3, 'CANCEL keeps the flow');
      expect(isChanged(v), true);

      replaced = v.probe.newFlow();
      await answer('Unsaved changes', 'Save');
      await replaced;
      expect(v.probe.flow.getNodeCount(), 0, 'SAVE goes on');
      expect(`${v.view.name}, ${v.probe.flowSettings.scriptName}`, 'Flow, My flow', 'the new flow has its own name');
      expect(await savedValue(script), 4, 'after it saved');

      await v.probe.openRecentFlow((await grok.dapi.scripts.find(script.id))!);
      expect(await until(() => v.probe.flow.getNodeCount() === 3 && !isChanged(v), 15000), true, 'the flow reopens');
      edit(v);
      replaced = v.probe.newFlow();
      await answer('Unsaved changes', 'Discard');
      await replaced;
      expect(v.probe.flow.getNodeCount(), 0, 'DISCARD goes on');
      expect(await savedValue(script), 4, 'and saves nothing');

      await v.probe.newFlow();
      expect(dialog('Unsaved changes') == null, true, 'nothing to lose, nothing asked');
    } finally {
      closeView(v);
      await grok.dapi.scripts.delete(script).catch(() => {});
    }
  }, {timeout: 60000});

  test('opening a flow, a template or a table\'s history asks too; a template is no version of it', async () => {
    const script = await saveScript();
    const v = await forScript(script);
    const template = 'Workflow Demo.flow';
    // a test bundle's package has no files on the server
    const json = JSON.stringify({...await flowDocument(), name: 'Template'});
    Object.defineProperty(_package, 'files', {value: {readAsText: async () => json}, configurable: true});
    try {
      edit(v);
      for (const replace of [() => v.probe.openRecentFlow(script), () => v.probe.loadTemplate(template),
        () => v.probe.importCreationScript('Demog = OpenFile("System:DemoFiles/demog.csv")')]) {
        const replaced = replace();
        await answer('Unsaved changes', 'Cancel');
        await replaced;
        expect(v.probe.flow.getNodeCount() === 3 && isChanged(v), true, 'CANCEL keeps the flow');
        expect(v.probe.boundScript?.id, script.id);
      }

      const replaced = v.probe.loadTemplate(template);
      await answer('Unsaved changes', 'Discard');
      await replaced;
      expect(v.view.name, 'Template', 'the template is on the canvas');
      expect(isChanged(v), false, 'an opened template has no changes');
      expect(v.probe.boundScript === null, true, 'a save of the template must not overwrite the flow opened before');
    } finally {
      delete (_package as {files?: unknown}).files;
      closeView(v);
      await grok.dapi.scripts.delete(script).catch(() => {});
    }
  }, {timeout: 60000});
});
