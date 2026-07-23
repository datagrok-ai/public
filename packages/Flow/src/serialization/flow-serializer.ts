/** Saves and loads FuncFlow documents. .ffjson v2 — Rete-native. */

import * as grok from 'datagrok-api/grok';

import {FlowEditor} from '../rete/flow-editor';
import {createNode, ensureFunctionsRegistered} from '../rete/node-factory';
import {FlowSettings, FuncFlowDocument, FuncFlowConnection} from './flow-schema';

export function serializeFlow(flow: FlowEditor, settings: FlowSettings): FuncFlowDocument {
  const nodes = flow.getNodes().map((n): FuncFlowDocument['nodes'][number] => ({
    id: n.id,
    typeName: n.dgTypeName ?? '',
    label: n.label,
    description: n.description,
    collapsed: n.collapsed || undefined, // omit when false to keep saves tidy
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
    waypoints: c.waypoints?.map((w) => ({x: w.x, y: w.y})),
  }));

  const annotations = flow.getAnnotations().map((a) => a.toDoc());
  const groups = flow.getGroups().map((g) => g.toDoc());

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
    annotations,
    groups,
    metadata: {settings},
  };
}

/** Deserialize into the editor. Clears the editor first. Skips nodes whose
 *  registered type is not currently known (e.g. a DG.Func that disappeared).
 *  Connections referencing missing nodes are silently skipped. */
export async function deserializeFlow(doc: FuncFlowDocument, flow: FlowEditor): Promise<void> {
  // DG-func node factories exist only after catalog registration; make it
  // deterministic here so no load path can race the view's deferred timer.
  ensureFunctionsRegistered();
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
    // Migration: older saves carried the per-node annotation in
    // properties.description before we promoted it to a top-level field.
    node.description = docNode.description ?? String(docNode.properties?.description ?? '');
    delete (node.properties as Record<string, unknown>).description;
    node.collapsed = docNode.collapsed === true;
    await flow.addNodeAt(node, docNode.pos.x, docNode.pos.y);
    idMap.set(docNode.id, node.id);
  }

  for (const c of doc.connections) {
    const source = idMap.get(c.source);
    const target = idMap.get(c.target);
    if (!source || !target) continue;
    await flow.addConnectionByKeys(source, c.sourceOutput, target, c.targetInput);
    if (c.waypoints && c.waypoints.length > 0) {
      // Match by source/target/keys since the new connection has a fresh id.
      const newConn = flow.getConnections().find((nc) =>
        nc.source === source && nc.target === target &&
        String(nc.sourceOutput) === c.sourceOutput &&
        String(nc.targetInput) === c.targetInput,
      );
      if (newConn) newConn.waypoints = c.waypoints.map((w) => ({x: w.x, y: w.y}));
    }
  }

  if (doc.annotations) {
    for (const a of doc.annotations) flow.addAnnotation(a);
  }

  if (doc.groups) {
    for (const gd of doc.groups) {
      // Node ids remap on load — resolve members through idMap; a group whose
      // nodes all vanished (unknown types skipped above) is dropped.
      const memberIds = gd.memberIds
        .map((mid) => idMap.get(mid))
        .filter((mid): mid is string => mid !== undefined);
      if (memberIds.length === 0) continue;
      flow.createGroup(memberIds, {
        title: gd.title,
        description: gd.description,
        minimized: gd.minimized,
        pos: gd.pos,
      });
    }
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

