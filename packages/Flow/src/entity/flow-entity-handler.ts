/** JS-side behavior behind the Flow entity — a Script (language 'flow') holding the flow JSON document. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {ensureFunctionsRegistered} from '../rete/node-factory';
import {deserializeFlow} from '../serialization/flow-serializer';
import {parseFlowBody} from '../serialization/flow-script-format';
import {FuncFlowDocument} from '../serialization/flow-schema';
import {emitScript, ScriptSettings} from '../compiler/script-emitter';
import {FuncFlowView} from '../funcflow-view';

function settingsOf(doc: FuncFlowDocument): ScriptSettings {
  const s = doc.metadata?.settings;
  return {
    name: s?.scriptName ?? doc.name ?? 'Flow',
    description: s?.scriptDescription ?? doc.description ?? '',
    tags: s?.tags ?? [],
  };
}

export class FlowEntityHandler {
  static readonly instance = new FlowEntityHandler();

  /** Runs `fn` against a graph deserialized into an off-screen editor; destroyed afterwards. */
  async withDetachedEditor<T>(doc: FuncFlowDocument, fn: (flow: FlowEditor) => T | Promise<T>): Promise<T> {
    ensureFunctionsRegistered();
    const container = ui.div([], {
      style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'},
    });
    document.body.appendChild(container);
    const flow = new FlowEditor(container);
    try {
      await deserializeFlow(doc, flow);
      return await fn(flow);
    } finally {
      try {
        flow.destroy();
      } finally {
        container.remove();
      }
    }
  }

  compileToJs(doc: FuncFlowDocument): Promise<string> {
    return this.withDetachedEditor(doc, (flow) => emitScript(flow, settingsOf(doc)));
  }

  /** Executes a flow-script call: compile to the JS twin, run it, copy the outputs back onto the call. */
  async run(scriptCall: DG.FuncCall): Promise<void> {
    const script = scriptCall.func as DG.Script;
    const {doc} = parseFlowBody(script.script);
    const js = await this.compileToJs(doc);
    const params: Record<string, unknown> = {};
    for (const name of Object.keys(scriptCall.inputParams))
      params[name] = scriptCall.inputs[name];
    const jsCall = DG.Script.create(js).prepare(params);
    await jsCall.call(undefined, undefined, {processed: true});
    for (const name of Object.keys(scriptCall.outputParams))
      scriptCall.setParamValue(name, jsCall.outputs[name]);
  }

  editorView(script: DG.Script): DG.ViewBase {
    return FuncFlowView.forScript(script);
  }

  previewView(script: DG.Script): DG.ViewBase {
    return FuncFlowView.forScript(script);
  }

  widget(script: DG.Script): DG.Widget {
    const host = ui.divV([], 'ff-entity-widget');
    try {
      const {doc} = parseFlowBody(script.script);
      const inputs = doc.nodes.filter((n) => n.typeName.startsWith('Inputs/')).length;
      const outputs = doc.nodes.filter((n) => n.typeName.startsWith('Outputs/')).length;
      const facts: {[key: string]: string} = {
        'Steps': `${doc.nodes.length} nodes, ${doc.connections.length} links`,
        'Parameters': `${inputs} in, ${outputs} out`,
      };
      if (doc.author && doc.author !== 'unknown') facts['Author'] = doc.author;
      if (doc.modified) facts['Modified'] = doc.modified.substring(0, 10);
      if (doc.description) host.appendChild(ui.divText(doc.description));
      host.appendChild(ui.tableFromMap(facts));
      host.appendChild(ui.divH([
        ui.button('Open editor', () => grok.shell.addView(this.editorView(script) as DG.ViewBase)),
        ui.button('Run…', () => (script as DG.Func).prepare().edit()),
      ]));
    } catch (e) {
      host.appendChild(ui.divText(`Cannot read flow body: ${e instanceof Error ? e.message : e}`));
    }
    return DG.Widget.fromRoot(host);
  }
}
