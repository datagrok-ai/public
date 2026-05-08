/** Saves and loads FuncFlow documents. .ffjson v2 — Rete-native. */

import * as grok from 'datagrok-api/grok';

import {FlowEditor} from '../rete/flow-editor';
import {createNode} from '../rete/node-factory';
import {FlowSettings, FuncFlowDocument, FuncFlowConnection} from './flow-schema';

export function serializeFlow(flow: FlowEditor, settings: FlowSettings): FuncFlowDocument {
  const nodes = flow.getNodes().map((n): FuncFlowDocument['nodes'][number] => ({
    id: n.id,
    typeName: n.dgTypeName ?? '',
    label: n.label,
    pos: {x: n.pos.x, y: n.pos.y},
    properties: {...n.properties},
    inputValues: {...n.inputValues},
  }));

  const connections: FuncFlowConnection[] = flow.getConnections().map((c) => ({
    id: c.id,
    source: c.source,
    sourceOutput: String(c.sourceOutput),
    target: c.target,
    targetInput: String(c.targetInput),
  }));

  let author = 'unknown';
  try {author = grok.shell.user?.login ?? 'unknown';} catch { /* no shell */ }

  return {
    version: '2.0',
    name: settings.scriptName,
    description: settings.scriptDescription,
    author,
    created: new Date().toISOString(),
    modified: new Date().toISOString(),
    nodes,
    connections,
    metadata: {settings},
  };
}

/** Deserialize into the editor. Clears the editor first. Skips nodes whose
 *  registered type is not currently known (e.g. a DG.Func that disappeared).
 *  Connections referencing missing nodes are silently skipped. */
export async function deserializeFlow(doc: FuncFlowDocument, flow: FlowEditor): Promise<void> {
  await flow.clear();

  // Map old node ids → new node ids (Rete generates a fresh id on construction).
  const idMap = new Map<string, string>();

  for (const docNode of doc.nodes) {
    const node = createNode(docNode.typeName);
    if (!node) {
      console.warn(`FuncFlow: skipped unknown node type "${docNode.typeName}"`);
      continue;
    }
    node.label = docNode.label;
    node.properties = {...docNode.properties};
    node.inputValues = {...docNode.inputValues};
    await flow.addNodeAt(node, docNode.pos.x, docNode.pos.y);
    idMap.set(docNode.id, node.id);
  }

  for (const c of doc.connections) {
    const source = idMap.get(c.source);
    const target = idMap.get(c.target);
    if (!source || !target) continue;
    await flow.addConnectionByKeys(source, c.sourceOutput, target, c.targetInput);
  }
}

export function downloadFlow(doc: FuncFlowDocument): void {
  const json = JSON.stringify(doc, null, 2);
  const blob = new Blob([json], {type: 'application/json'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `${doc.name || 'flow'}.ffjson`;
  a.click();
  URL.revokeObjectURL(url);
}

export async function loadFlowFromFile(file: File): Promise<FuncFlowDocument> {
  const text = await file.text();
  const doc = JSON.parse(text) as FuncFlowDocument;
  if (doc.version !== '2.0')
    throw new Error(`Unsupported flow file version "${doc.version}"; expected 2.0`);
  return doc;
}

