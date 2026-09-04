/** Saves and loads FuncFlow documents (the `.flow` JSON payload, v2). */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {createNode, ensureFunctionsRegistered} from '../rete/node-factory';
import {applyPropertyInputShape} from '../rete/nodes/input-nodes';
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

/** Clears the editor first; unknown node types and connections referencing missing nodes are skipped. */
export async function deserializeFlow(doc: FuncFlowDocument, flow: FlowEditor): Promise<void> {
  // DG-func factories exist only after catalog registration — no load path may race the view's deferred timer.
  ensureFunctionsRegistered();
  await flow.clear();

  // Old node id → new node id (Rete generates a fresh id on construction).
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
    // Migration: older saves carried the per-node annotation in properties.description.
    node.description = docNode.description ?? String(docNode.properties?.description ?? '');
    delete (node.properties as Record<string, unknown>).description;
    node.collapsed = docNode.collapsed === true;
    // A Property Input's live shape (type, socket, qualifier keys) derives from its properties.
    applyPropertyInputShape(node);
    await flow.addNodeAt(node, docNode.pos?.x ?? 0, docNode.pos?.y ?? 0);
    idMap.set(docNode.id, node.id);
  }

  for (const c of doc.connections) {
    const source = idMap.get(c.source);
    const target = idMap.get(c.target);
    if (!source || !target) continue;
    try {
      await flow.addConnectionByKeys(source, c.sourceOutput, target, c.targetInput);
    } catch (e) {
      // A slot key that no longer exists (renamed/removed func parameter) — skip the wire, keep loading.
      console.warn(`FuncFlow: skipped connection ${c.sourceOutput} → ${c.targetInput}: ${(e as Error).message}`);
      continue;
    }
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
      // Members resolve through idMap; a group whose nodes all vanished is dropped.
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
  DG.Utils.download(`${doc.name || 'flow'}.flow`, json, 'application/json');
}

/** Reads a `.flow` file: bare JSON or annotated script body. Inlined rather than
 *  reusing parseFlowBody — flow-script-format imports this module. */
export async function loadFlowFromFile(file: File): Promise<FuncFlowDocument> {
  const text = await file.text();
  const lines = text.split('\n');
  let i = 0;
  while (i < lines.length && (lines[i].trimStart().startsWith('//') || lines[i].trim() === ''))
    i++;
  const doc = JSON.parse(lines.slice(i).join('\n')) as FuncFlowDocument;
  if (doc.version !== '2.0')
    throw new Error(`Unsupported flow file version "${doc.version}"; expected 2.0`);
  return doc;
}

