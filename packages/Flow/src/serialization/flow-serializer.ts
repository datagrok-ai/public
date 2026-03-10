import {LGraph, LGraphNode} from 'litegraph.js';
import * as grok from 'datagrok-api/grok';
import {FuncFlowDocument, FuncFlowMetadata, FuncFlowNodeMeta, FlowSettings} from './flow-schema';

/** Serialize a graph to FuncFlow document format */
export function serializeFlow(graph: LGraph, settings: FlowSettings): FuncFlowDocument {
  const nodes = (graph as any)._nodes as LGraphNode[];
  const nodesMeta: Record<number, FuncFlowNodeMeta> = {};

  for (const node of nodes) {
    nodesMeta[node.id] = {
      dgFuncName: (node as any).dgFuncName || node.title,
      dgNodeType: (node as any).dgNodeType || 'func',
      paramName: node.properties?.['paramName'],
      defaultValue: node.properties?.['defaultValue'],
      description: node.properties?.['description'],
    };
  }

  let author = 'unknown';
  try {author = grok.shell.user?.login ?? 'unknown';} catch {/* ok */}

  return {
    version: '1.0',
    name: settings.scriptName,
    description: settings.scriptDescription,
    author,
    created: new Date().toISOString(),
    modified: new Date().toISOString(),
    graph: graph.serialize(),
    metadata: {nodes: nodesMeta, settings},
  };
}

/** Deserialize a FuncFlow document into a graph */
export function deserializeFlow(doc: FuncFlowDocument, graph: LGraph): void {
  graph.configure(doc.graph);

  // Restore FuncFlow-specific metadata
  const nodes = (graph as any)._nodes as LGraphNode[];
  for (const node of nodes) {
    const meta = doc.metadata?.nodes?.[node.id];
    if (meta) {
      (node as any).dgFuncName = meta.dgFuncName;
      (node as any).dgNodeType = meta.dgNodeType;
      if (meta.paramName !== undefined) node.properties['paramName'] = meta.paramName;
      if (meta.defaultValue !== undefined) node.properties['defaultValue'] = meta.defaultValue;
      if (meta.description !== undefined) node.properties['description'] = meta.description;
    }
  }
}

/** Download a flow as a .funcflow.json file */
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

/** Load a flow from a File object */
export async function loadFlowFromFile(file: File): Promise<FuncFlowDocument> {
  const text = await file.text();
  return JSON.parse(text) as FuncFlowDocument;
}
