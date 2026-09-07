/** Scheme types and base node class for FuncFlow's Rete editor. */
import {ClassicPreset, GetSchemes} from 'rete';
import * as DG from 'datagrok-api/dg';
import {TypedSocket} from './sockets';
import {isPrimitiveSlotType} from '../types/type-map';
import type {FuncWrapper} from '../utils/func-input-overrides';

export type DgNodeType = 'input' | 'output' | 'utility' | 'func';

/** The callback surface a `FlowEditor` exposes to the React node components — stamped per node
 *  so a component always talks to the editor that owns it (several editors can coexist on a page). */
export interface FlowEditorBridge {
  toggleCollapsed(id: string): void;
  isSocketConnected(nodeId: string, side: 'input' | 'output', key: string): boolean;
  /** Same contract as the property panel's `paramsChanged`. */
  notifyParamsChanged(nodeId: string): void;
  /** Pop the "Shown inputs" checkbox menu (the ⋯ hidden-inputs indicator). */
  showShownInputsMenu(nodeId: string, event: MouseEvent): void;
  /** Flip the in-node viewer/widget preview (`properties['inlinePreview']`). */
  toggleInlinePreview(nodeId: string): void;
  /** The live viewer/widget root captured for this node — shown by the in-node preview. */
  getInlinePreviewContent(nodeId: string): HTMLElement | null;
  /** Bind/refresh the screen-space preview portal tracking a node's in-card container. */
  syncInlinePreview(nodeId: string, host: HTMLElement): void;
  /** Tear the node's preview portal down (toggle-off, collapse, unmount). */
  releaseInlinePreview(nodeId: string): void;
  /** Current canvas zoom — the in-node sketcher exists only at native zoom (1)
   *  and folds away when the user zooms. */
  getZoom(): number;
  /** Snap the canvas to native zoom and pan so a `w`×`h` box anchored at the
   *  node's top-left is fully visible — what expanding the in-node sketcher does. */
  focusNodeForEditing(nodeId: string, w: number, h: number): void;
  /** True while a run in progress still has this node ahead of it — the empty
   *  preview shows a loader instead of the "run the flow" placeholder. */
  isInlinePreviewPending(nodeId: string): boolean;
}

/** Base class for every canvas node — `ClassicPreset.Node` plus FuncFlow metadata. */
export class FlowNode extends ClassicPreset.Node<
  Record<string, TypedSocket>,
  Record<string, TypedSocket>,
  Record<string, ClassicPreset.Control>
> {
  /** What kind of node this is — drives compiler behavior. */
  dgNodeType: DgNodeType = 'func';

  /** Input nodes: the DG type emitted on the output slot; output nodes: the type expected on the input slot. */
  dgOutputType?: string;

  dgFunc?: DG.Func;

  /** Qualified name (`Pkg:funcName`) used in `grok.functions.call(...)`. */
  dgFuncName?: string;

  dgRole?: string | null;

  /** Number of pass-through outputs at the start of `outputs` (func nodes). */
  passthroughCount = 0;

  /** Hardcoded values for unconnected primitive inputs, keyed by input name — kept here
   *  (not as Rete `Control`s, which render inside the node) so they're edited in the side panel only. */
  inputValues: Record<string, any> = {};

  /** Free-form node properties (paramName, defaultValue, qualifiers). */
  properties: Record<string, any> = {};

  /** Annotation rendered below the title; for input/output nodes also the `[description]`
   *  suffix in the generated `//input:` / `//output:` line. */
  description: string = '';

  /** Registered type name stamped by `createNode()`; needed for serialization. */
  dgTypeName?: string;

  collapsed = false;

  /** Run status shown under the title (e.g. "Done · 1,204 × 8"); empty when idle. */
  statusText = '';

  /** Input keys hidden from the node body but still data-carrying (slot, seeded value,
   *  compilation, import/emit untouched); a connected hidden socket still renders. */
  hiddenInputs: ReadonlySet<string> = new Set();

  /** Real output keys hidden the same way; a connected one still renders. */
  hiddenOutputs: ReadonlySet<string> = new Set();

  /** Extra readiness check (`FUNC_NODE_VALIDATORS`); returns labels of what's missing;
   *  must be synchronous (runs on every render and in the run gate). */
  extraValidator?: (node: FlowNode) => string[];

  /** When set, the node's input slots are the wrapper's exposed inputs, not the function's own;
   *  the compiler folds their resolved expressions into the real call args via `mapInputs`. */
  funcWrapper?: FuncWrapper;

  /** Input keys that must be connected or filled for the node to run; drives the "Needs input" hint. */
  requiredInputs: string[] = [];

  /** Property keys that must carry a non-empty value — the panel-property analogue of {@link requiredInputs}. */
  requiredProps: string[] = [];

  /** Per-slot descriptions keyed by slot key, shown as socket tooltips and in the context panel. */
  inputDescriptions: Record<string, string> = {};
  outputDescriptions: Record<string, string> = {};

  /** Source package of the underlying function (`''` for core / built-ins). */
  dgPackageName = '';

  /** Visual position, kept in sync with AreaPlugin's NodeView for serialization. */
  pos: {x: number; y: number} = {x: 0, y: 0};

  /** Back-reference to the owning editor; runtime-only (the serializer picks fields explicitly). */
  editorBridge?: FlowEditorBridge;

  /** Runtime-only companion to `properties['defaultValue']` for values a string can't rebuild
   *  (a picked DataFrame, a FileInfo); never serialized — re-resolved from the workspace on reload. */
  transientValue?: unknown;

  get title(): string {
    return this.label;
  }
  set title(v: string) {
    this.label = v;
  }

  hasInputValue(name: string): boolean {
    return Object.prototype.hasOwnProperty.call(this.inputValues, name);
  }

  /** Type-based default for {@link inputSlotShown}: primitive-typed inputs start hidden — edited in the panel. */
  inputSlotShownByDefault(name: string): boolean {
    const slot = (this.inputs as Record<string, {socket?: TypedSocket} | undefined>)[name];
    return !isPrimitiveSlotType(slot?.socket?.dgType ?? '');
  }

  /** Whether an input's socket row renders on the card; per-node deviations live in
   *  `properties['shownSlots']` so save/load, duplicate, and copy-paste carry them for free. */
  inputSlotShown(name: string): boolean {
    const overrides = this.properties['shownSlots'] as Record<string, boolean> | undefined;
    return overrides?.[name] ?? this.inputSlotShownByDefault(name);
  }
}

/** Whether a socket row is absent from the rendered card — unless the slot is connected,
 *  in which case it always renders so the wire keeps its endpoint. */
export function hiddenSocketRow(node: FlowNode, side: 'input' | 'output', key: string,
  isConnected: (side: 'input' | 'output', key: string) => boolean): boolean {
  if (isConnected(side, key)) return false;
  const base = side === 'input' ? key : key.endsWith('__pt') ? key.slice(0, -'__pt'.length) : null;
  if (base === null) return node.hiddenOutputs.has(key);
  return node.hiddenInputs.has(base) || !node.inputSlotShown(base);
}

/** Generic params left at `ClassicPreset.Node` so the connection-plugin's `ClassicScheme`
 *  constraint is satisfied without TS variance complaints. */
export class FlowConnection extends ClassicPreset.Connection<ClassicPreset.Node, ClassicPreset.Node> {
  isPseudo?: boolean;
  /** Optional routing points in canvas coords; the connection path chains start → waypoints → end. */
  waypoints?: Array<{x: number; y: number}>;
}

export type FlowScheme = GetSchemes<FlowNode, FlowConnection>;

/** Execution-ordering ports added to every node: a pure topological dependency
 *  that carries no data — the compiler ignores it. */
export const ORDER_SOCKET_TYPE = 'order';
export const EXEC_IN_KEY = '__exec_in';
export const EXEC_OUT_KEY = '__exec_out';

export function isExecKey(key: string): boolean {
  return key === EXEC_IN_KEY || key === EXEC_OUT_KEY;
}

/** Flow treats SetVar and Value Output as the same concept — both register their value
 *  in the run context and declare a script output, compiling to the same thing. */
export function isSetVarNode(node: FlowNode): boolean {
  return (node.dgFunc?.name?.toLowerCase() ?? '') === 'setvar';
}

/** In-node preview of a viewer/widget output. Both keys live in `node.properties`
 *  so save/load, duplicate, and copy-paste carry the choice for free. */
export const INLINE_PREVIEW_PROP = 'inlinePreview';
export const INLINE_PREVIEW_SIZE_PROP = 'inlinePreviewSize';
export const INLINE_PREVIEW_DEFAULT_SIZE = 300;

/** `dataset` key stamped on a live root while the in-node preview hosts it — the
 *  bottom output panel renders a note instead of stealing the element. */
export const INLINE_HOSTED_DATA_KEY = 'ffInlineHosted';

/** Whether the node can show an in-node preview: it produces a viewer, a widget,
 *  or graphics (a manual viewer node, or a function declaring such an output —
 *  e.g. Chem's Gasteiger Partial Charges script). */
export function supportsInlinePreview(node: FlowNode): boolean {
  if (node.properties['viewerType'] != null) return true;
  for (const [key, out] of Object.entries(
    node.outputs as Record<string, {socket?: TypedSocket} | undefined>)) {
    if (isExecKey(key) || key.endsWith('__pt')) continue;
    const t = out?.socket?.dgType;
    if (t === 'viewer' || t === 'widget' || t === 'graphics') return true;
  }
  return false;
}

export function inlinePreviewEnabled(node: FlowNode): boolean {
  return node.properties[INLINE_PREVIEW_PROP] === true && supportsInlinePreview(node);
}

/** Size of the expanded in-node sketcher — it only exists at native zoom (1),
 *  so layout px equal screen px. */
export const INLINE_SKETCHER_SIZE = 500;

/** Default in-card box of the HELM editor on a Helm Input node's body. */
export const INLINE_HELM_WIDTH = 300;
export const INLINE_HELM_HEIGHT = 250;

/** User-resized helm box, `{width, height}` in `node.properties` (serializes). */
export const EDITOR_BOX_SIZE_PROP = 'editorSize';

/** A Molecule-tagged string input hosts the compact-preview / expand-in-place
 *  sketcher on the node body (the panel's Value row keeps the standard editor).
 *  Reads the `semType` qualifier, not the node type, so a String Input the user
 *  tags `Molecule` gets it too — same routing as `inputValueProperty`. */
export function hostsInlineSketcher(node: FlowNode): boolean {
  return node.dgNodeType === 'input' && node.dgOutputType === 'string' &&
    String(node.properties['semType'] ?? '').trim() === 'Molecule';
}

/** A Macromolecule-tagged string input hosts the HELM editor in a resizable
 *  in-card box on the node body. */
export function hostsHelmEditor(node: FlowNode): boolean {
  return node.dgNodeType === 'input' && node.dgOutputType === 'string' &&
    String(node.properties['semType'] ?? '').trim() === 'Macromolecule';
}

/** The helm box; malformed/missing stored values fall back to the 300×250 default. */
export function editorBoxSize(node: FlowNode): {width: number; height: number} {
  const s = node.properties[EDITOR_BOX_SIZE_PROP] as {width?: unknown; height?: unknown} | undefined;
  const dim = (v: unknown, d: number): number =>
    typeof v === 'number' && isFinite(v) && v > 0 ? Math.round(v) : d;
  return {width: dim(s?.width, INLINE_HELM_WIDTH), height: dim(s?.height, INLINE_HELM_HEIGHT)};
}


/** Preview container size; malformed/missing values fall back to the 300×300 default. */
export function inlinePreviewSize(node: FlowNode): {width: number; height: number} {
  const s = node.properties[INLINE_PREVIEW_SIZE_PROP] as {width?: unknown; height?: unknown} | undefined;
  const dim = (v: unknown): number =>
    typeof v === 'number' && isFinite(v) && v > 0 ? Math.round(v) : INLINE_PREVIEW_DEFAULT_SIZE;
  return {width: dim(s?.width), height: dim(s?.height)};
}

/** Labels of required inputs neither connected nor filled; pure — `isConnected` is supplied by the caller. */
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

export function missingRequiredProps(node: FlowNode): string[] {
  const missing: string[] = [];
  for (const key of node.requiredProps) {
    const v = node.properties[key];
    if (v === undefined || v === null || String(v).trim() === '') missing.push(key);
  }
  return missing;
}

/** Everything the user still has to provide before the node can run —
 *  drives the "Needs input" hint and the run gate. */
export function nodeMissingRequirements(node: FlowNode, isConnected: (key: string) => boolean): string[] {
  return [
    ...missingRequiredInputs(node, isConnected),
    ...missingRequiredProps(node),
    // Stamped by FuncNode so the node model never depends on the panel's editors.
    ...(node.extraValidator?.(node) ?? []),
  ];
}
