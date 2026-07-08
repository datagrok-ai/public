/** Scheme types and base node class for FuncFlow's Rete editor. */
import {ClassicPreset, GetSchemes} from 'rete';
import * as DG from 'datagrok-api/dg';
import {TypedSocket} from './sockets';

export type DgNodeType = 'input' | 'output' | 'utility' | 'func';

/** The narrow callback surface a `FlowEditor` exposes to the React node
 *  components it renders. Stamped onto every node that enters an editor's data
 *  layer (`FlowNode.editorBridge`), so a component always talks to the editor
 *  that owns it — several editors can coexist on a page (file previews, the
 *  creation-script dialog, detached compile editors). */
export interface FlowEditorBridge {
  toggleCollapsed(id: string): void;
  isSocketConnected(nodeId: string, side: 'input' | 'output', key: string): boolean;
}

/** Base class for every node we put on the canvas.
 *
 * Extends `ClassicPreset.Node` with FuncFlow-specific metadata (Datagrok
 * function reference, role, pass-through bookkeeping, free-form properties).
 *
 * The Rete `inputs` / `outputs` / `controls` collections are managed by the
 * superclass; everything else lives directly on the instance. */
export class FlowNode extends ClassicPreset.Node<
  Record<string, TypedSocket>,
  Record<string, TypedSocket>,
  Record<string, ClassicPreset.Control>
> {
  /** What kind of node this is — drives compiler behavior. */
  dgNodeType: DgNodeType = 'func';

  /** For input nodes: the DG type emitted on the single output slot.
   *  For output nodes: the DG type expected on the single input slot. */
  dgOutputType?: string;

  /** For func nodes: the underlying DG.Func reference. */
  dgFunc?: DG.Func;

  /** Qualified name (`Pkg:funcName`) used in `grok.functions.call(...)`. */
  dgFuncName?: string;

  /** DG role string, used for theming and categorization. */
  dgRole?: string | null;

  /** Number of pass-through outputs at the start of `outputs` (func nodes). */
  passthroughCount = 0;

  /** Hardcoded values for unconnected primitive inputs. Keyed by input name.
   *  Used by func nodes for `_input_<name>` defaults shown in the property
   *  panel — Rete `Control` objects render *inside* the node, but we want
   *  these in the side panel only, so they go here instead. */
  inputValues: Record<string, any> = {};

  /** Free-form node properties (paramName, defaultValue, qualifiers). */
  properties: Record<string, any> = {};

  /** KNIME-style annotation rendered below the title. For input/output nodes
   *  this also becomes the `[description]` suffix in the generated `//input:` /
   *  `//output:` line. Empty string means no annotation. */
  description: string = '';

  /** Registered type name from `node-factory.ts`, stamped at construction
   *  via `createNode()`. Needed for serialization. */
  dgTypeName?: string;

  /** When true the node renders as a single title bar with no body.
   *  Toggled by the caret in the title bar (the status dot is display-only). */
  collapsed = false;

  /** Short, plain-language run status shown under the node title
   *  (e.g. "Running…", "Done · 1,204 × 8", "Error"). Set by
   *  `ExecutionVisualizer`; empty when idle. */
  statusText = '';

  /** Input keys that must be satisfied (connected, or filled in the panel) for
   *  the node to do anything — the structural inputs (a table, a column).
   *  Populated by `FuncNode`/output nodes; drives the "Needs input" hint. */
  requiredInputs: string[] = [];

  /** Property keys that must carry a non-empty value for the node to run — the
   *  node's non-socket requirements (Select Column's column name, Select Table's
   *  table name). The panel-property analogue of {@link requiredInputs}; also
   *  drives the "Needs input" hint and the run gate. */
  requiredProps: string[] = [];

  /** Per-slot descriptions (from the DG.Func param `description`/`caption`),
   *  keyed by input/output slot key — shown as hover tooltips on the node's
   *  sockets and in the context panel. Empty entries are omitted. */
  inputDescriptions: Record<string, string> = {};
  outputDescriptions: Record<string, string> = {};

  /** Source package of the underlying function (`''` for core / built-ins).
   *  Shown in the context panel so a vague function name is disambiguated. */
  dgPackageName = '';

  /** Visual position — kept in sync with AreaPlugin's NodeView for
   *  serialization. Updated by `FlowEditor` on `nodetranslated`. */
  pos: {x: number; y: number} = {x: 0, y: 0};

  /** Back-reference to the owning editor's callback surface, stamped by
   *  `FlowEditor` when the node enters its data layer. Runtime-only — the
   *  serializer picks fields explicitly, so this never reaches `.ffjson`. */
  editorBridge?: FlowEditorBridge;

  /** Human-friendly title — `label` from the Rete superclass is what we
   *  render, so this is just an alias for symmetry with the LiteGraph world. */
  get title(): string {
    return this.label;
  }
  set title(v: string) {
    this.label = v;
  }

  /** Whether a hardcoded value is recorded for an input slot. */
  hasInputValue(name: string): boolean {
    return Object.prototype.hasOwnProperty.call(this.inputValues, name);
  }
}

/** Connection type. Generic params left at `ClassicPreset.Node` (rather than
 *  narrowed to `FlowNode`) so it satisfies the connection-plugin's
 *  `ClassicScheme` constraint without TS variance complaints. */
export class FlowConnection extends ClassicPreset.Connection<ClassicPreset.Node, ClassicPreset.Node> {
  isPseudo?: boolean;
  /** Optional routing points. The React Connection component chains
   *  classicConnectionPath through start → waypoints[…] → end. Stored in
   *  canvas coords; the editor exposes `addWaypoint`/`removeWaypoint` to
   *  mutate and trigger a re-render. */
  waypoints?: Array<{x: number; y: number}>;
}

export type FlowScheme = GetSchemes<FlowNode, FlowConnection>;

/** Execution-ordering ports added to *every* node by `createNode`. A connection
 *  exec-out → exec-in is a pure topological dependency ("run the source, and
 *  therefore everything before it, before the target"): the topological sort
 *  treats it like any other edge, but it carries no data — the compiler ignores
 *  it and it never produces a variable. Lets users express run-order explicitly
 *  instead of relying on vertical node position. */
export const ORDER_SOCKET_TYPE = 'order';
export const EXEC_IN_KEY = '__exec_in';
export const EXEC_OUT_KEY = '__exec_out';

/** Whether a port key is one of the execution-ordering ports. */
export function isExecKey(key: string): boolean {
  return key === EXEC_IN_KEY || key === EXEC_OUT_KEY;
}

/** Whether a func node is the platform `SetVar`. Flow treats SetVar and Value
 *  Output as the same concept: both register their value in the run context
 *  under their name AND declare a script output of that name — they compile to
 *  the same thing (see script-emitter / creation-script-emitter). */
export function isSetVarNode(node: FlowNode): boolean {
  return (node.dgFunc?.name?.toLowerCase() ?? '') === 'setvar';
}

/** The labels of a node's required inputs that are neither connected nor filled
 *  with a value — i.e. what the user still has to provide. Pure (no editor
 *  dependency): `isConnected(key)` is supplied by the caller. Drives the
 *  "Needs input" hint on the node and lightweight pre-run validation. */
export function missingRequiredInputs(node: FlowNode, isConnected: (key: string) => boolean): string[] {
  const missing: string[] = [];
  for (const key of node.requiredInputs) {
    if (isConnected(key)) continue;
    const v = node.inputValues[key];
    if (v !== undefined && v !== null && String(v).trim() !== '') continue;
    const input = (node.inputs as Record<string, {label?: string} | undefined>)[key];
    missing.push(input?.label ?? key);
  }
  return missing;
}

/** The labels of a node's {@link FlowNode.requiredProps} left empty (undefined /
 *  null / blank) — the panel-property requirements the user still has to fill
 *  (a Select Column's column name, a Select Table's table name). */
export function missingRequiredProps(node: FlowNode): string[] {
  const missing: string[] = [];
  for (const key of node.requiredProps) {
    const v = node.properties[key];
    if (v === undefined || v === null || String(v).trim() === '') missing.push(key);
  }
  return missing;
}

/** Everything the user still has to provide before the node can run: required
 *  inputs (sockets) not wired/filled AND required properties left empty. Drives
 *  the node's "Needs input" hint and the run gate — a node with any missing
 *  requirement (a plot with no table, a Select Column with no column) is not
 *  run, and neither is anything downstream of it. */
export function nodeMissingRequirements(node: FlowNode, isConnected: (key: string) => boolean): string[] {
  return [...missingRequiredInputs(node, isConnected), ...missingRequiredProps(node)];
}
