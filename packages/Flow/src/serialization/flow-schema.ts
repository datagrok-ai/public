/** JSON schema for the .flow save format (Rete-based, version 2).
 *  Breaking change from v1: stores Rete nodes & connections directly — old .flow files will not load. */

export interface FuncFlowDocument {
  version: '2.0';
  name: string;
  description: string;
  author: string;
  created: string;
  modified: string;

  nodes: FuncFlowNode[];

  connections: FuncFlowConnection[];

  /** Purely visual. Optional for back-compat. */
  annotations?: FuncFlowAnnotation[];

  /** Visual only — `nodes`/`connections` stay flat. Optional for back-compat. */
  groups?: FuncFlowGroup[];

  /** Output-tab layouts keyed by paramName (node ids remap on load); excluded from dirty tracking. */
  outputViews?: {[paramName: string]: {layout: string}};

  /** Bound dashboard project — re-publishing updates it instead of creating a new one. */
  dashboard?: {projectId: string};

  metadata: FuncFlowMetadata;
}

export interface FuncFlowAnnotation {
  id: string;
  pos: {x: number; y: number};
  size: {w: number; h: number};
  text: string;
  color: string;
  /** Title font size in px; absent in older saves = the default (13). */
  fontSize?: number;
}

export interface FuncFlowGroup {
  id: string;
  title: string;
  description: string;
  /** Node ids — remapped through the loader's idMap. */
  memberIds: string[];
  minimized: boolean;
  /** Card anchor (canvas coords); the expanded frame derives from members. */
  pos: {x: number; y: number};
}

export interface FuncFlowNode {
  id: string;
  /** The registered type name from `node-factory.ts`. */
  typeName: string;
  label: string;
  /** Optional — the deserializer copies it from `properties.description` in older saves. */
  description?: string;
  /** Optional — absent in older saves; treated as false. */
  collapsed?: boolean;
  pos: {x: number; y: number};
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
  /** Optional routing waypoints (canvas coords). */
  waypoints?: Array<{x: number; y: number}>;
}

export interface FuncFlowMetadata {
  settings: FlowSettings;
}

export interface FlowSettings {
  scriptName: string;
  scriptDescription: string;
  tags: string[];
  /** Ribbon autorun toggle — saved with the flow so it reopens live. */
  autorun?: boolean;
}
