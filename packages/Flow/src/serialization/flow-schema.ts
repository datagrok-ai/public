/** JSON schema for the .ffjson save format (Rete-based, version 2).
 *
 * **Breaking change from v1**: stores Rete nodes & connections directly,
 * no LiteGraph payload. Old .ffjson files will not load. */

export interface FuncFlowDocument {
  version: '2.0';
  name: string;
  description: string;
  author: string;
  created: string;
  modified: string;

  /** Flat node list. Each node carries its concrete type name + properties. */
  nodes: FuncFlowNode[];

  /** Connections by source/target node id and slot key. */
  connections: FuncFlowConnection[];

  metadata: FuncFlowMetadata;
}

export interface FuncFlowNode {
  id: string;
  /** The registered type name from `node-factory.ts` (e.g. "Inputs/Table Input"
   *  or "DG Functions/Transform/MyFunc"). */
  typeName: string;
  /** Human-friendly label shown on the node title bar. */
  label: string;
  /** Canvas position. */
  pos: {x: number; y: number};
  /** Free-form node properties (paramName, defaultValue, etc.). */
  properties: Record<string, unknown>;
  /** Hardcoded values for unconnected primitive func inputs. */
  inputValues: Record<string, unknown>;
}

export interface FuncFlowConnection {
  id: string;
  source: string;
  sourceOutput: string;
  target: string;
  targetInput: string;
}

export interface FuncFlowMetadata {
  settings: FlowSettings;
}

export interface FlowSettings {
  scriptName: string;
  scriptDescription: string;
  tags: string[];
}
