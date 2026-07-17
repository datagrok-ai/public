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

  /** Workflow annotations — purely visual, not part of the executable graph.
   *  Optional for back-compat with .ffjson files that pre-date this field. */
  annotations?: FuncFlowAnnotation[];

  /** Node groups — collapsible titled frames around member node sets. Visual
   *  only: the `nodes`/`connections` above stay flat; a minimized group just
   *  renders its members as one card. Optional for back-compat. */
  groups?: FuncFlowGroup[];

  /** Saved TableView layouts of the output-view tabs, keyed by the output's
   *  paramName (node ids remap on load — never key by them). Optional for
   *  back-compat; excluded from dirty tracking (viewState strings are not
   *  canonical across serializations). */
  outputViews?: {[paramName: string]: {layout: string}};

  /** The dashboard project this flow publishes into — re-publishing updates
   *  that project instead of creating a new one per save. Optional; excluded
   *  from dirty tracking. */
  dashboard?: {projectId: string};

  metadata: FuncFlowMetadata;
}

export interface FuncFlowAnnotation {
  id: string;
  pos: {x: number; y: number};
  size: {w: number; h: number};
  text: string;
  color: string;
}

export interface FuncFlowGroup {
  id: string;
  title: string;
  description: string;
  /** Node ids (of the `nodes` array) — remapped through the loader's idMap. */
  memberIds: string[];
  minimized: boolean;
  /** Card anchor (canvas coords); the expanded frame derives from members. */
  pos: {x: number; y: number};
}

export interface FuncFlowNode {
  id: string;
  /** The registered type name from `node-factory.ts` (e.g. "Inputs/Table Input"
   *  or "DG Functions/Transform/MyFunc"). */
  typeName: string;
  /** Human-friendly label shown on the node title bar. */
  label: string;
  /** KNIME-style annotation rendered under the title. Optional — older saves
   *  that pre-date this field will have it copied from `properties.description`
   *  by the deserializer. */
  description?: string;
  /** Whether the node was collapsed (title-bar-only render) when saved.
   *  Optional — absent in older saves; treated as false. */
  collapsed?: boolean;
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
  /** Optional routing waypoints (canvas-coord points the line bends through). */
  waypoints?: Array<{x: number; y: number}>;
}

export interface FuncFlowMetadata {
  settings: FlowSettings;
}

export interface FlowSettings {
  scriptName: string;
  scriptDescription: string;
  tags: string[];
}
