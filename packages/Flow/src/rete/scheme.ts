/** Scheme types and base node class for FuncFlow's Rete editor. */
import {ClassicPreset, GetSchemes} from 'rete';
import * as DG from 'datagrok-api/dg';
import {TypedSocket} from './sockets';

export type DgNodeType = 'input' | 'output' | 'utility' | 'func';

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
   *  Toggled by clicking the status circle. */
  collapsed = false;

  /** Visual position — kept in sync with AreaPlugin's NodeView for
   *  serialization. Updated by `FlowEditor` on `nodetranslated`. */
  pos: {x: number; y: number} = {x: 0, y: 0};

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
