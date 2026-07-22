/** FlowEditor — owns the NodeEditor + AreaPlugin + ConnectionPlugin pipeline.
 *
 * Replaces LiteGraph's `CanvasController` and `GraphManager`. The host view
 * supplies a container; we build the editor inside it and expose just enough
 * surface for the rest of the package (compiler, view, property panel) to
 * operate on a Rete editor without knowing what's underneath. */

import {NodeEditor} from 'rete';
import {AreaExtensions, AreaPlugin} from 'rete-area-plugin';
import {
  ClassicFlow,
  ConnectionPlugin,
  SocketData,
} from 'rete-connection-plugin';
import {Presets as ReactPresets, ReactArea2D, ReactPlugin} from 'rete-react-plugin';
import {
  HistoryActions, HistoryExtensions, HistoryPlugin,
  Presets as HistoryPresets,
} from 'rete-history-plugin';
import {getDOMSocketPosition} from 'rete-render-utils';
import {createRoot} from 'react-dom/client';
import * as DG from 'datagrok-api/dg';

import {
  FlowConnection, FlowEditorBridge, FlowNode, FlowScheme, isExecKey, isSetVarNode,
  EXEC_IN_KEY, EXEC_OUT_KEY,
} from './scheme';
import {TypedSocket} from './sockets';
import {FlowConnectionComponent, FlowNodeComponent, FlowSocketComponent} from './node-component';
import {getSlotColor, getSlotLetter} from '../types/type-map';
import {tid, setTid} from '../utils/test-ids';
import {FlowAnnotation, AnnotationDoc, ANNOTATION_COLORS} from './annotation';
import {
  FlowGroup, GroupDoc, GROUP_TITLE_H, GROUP_PAD, GROUP_DOT_TOP, GROUP_DOT_STEP,
} from './node-group';
import {computeLayers, layoutGraph, LayoutEdge} from './graph-layout';

/** A classified graph edit — tells listeners *what* changed, so run results
 *  can be invalidated precisely (only downstream of the change) instead of
 *  wholesale. Cosmetic changes (node moves, annotations, titles) do not emit
 *  one. `params-changed` is reported by the property panel via
 *  {@link FlowEditor.notifyNodeParamsChanged}. */
export type GraphEdit =
  | {kind: 'node-added'; nodeId: string}
  | {kind: 'node-removed'; nodeId: string}
  | {kind: 'connection-added'; sourceId: string; targetId: string}
  | {kind: 'connection-removed'; sourceId: string; targetId: string}
  | {kind: 'params-changed'; nodeId: string}
  | {kind: 'cleared'};

export interface FlowEditorCallbacks {
  onNodeSelected?: (node: FlowNode) => void;
  onNodeDeselected?: (node: FlowNode) => void;
  /** Fired after selection changes that never go through `nodepicked` — the
   *  marquee (whose release is swallowed before it can bubble, see
   *  `installRectSelect`), Ctrl+A / Ctrl+Shift+A, the pointerup modifier
   *  semantics (toggle-off / remove / collapse), and programmatic
   *  select/unselect. Hosts that track "what is selected now" (the suggestion
   *  pane) listen here; per-node callbacks above stay click-driven. */
  onSelectionChanged?: () => void;
  onGraphChanged?: () => void;
  /** Fired with the classified edit for every change that can affect run
   *  results — drives precise invalidation and autorun. Fires alongside (not
   *  instead of) `onGraphChanged`, which remains the coarse "refresh UI" hook
   *  and also covers cosmetic changes (annotations). */
  onGraphEdited?: (edit: GraphEdit) => void;
  /** Run the slice up to this node and preview its output ("inspect anywhere").
   *  Wired from the node's right-click menu in addition to the output-port menu. */
  onPreviewNode?: (nodeId: string) => void;
  /** Re-run just this node using values captured from a prior run (no upstream
   *  re-run). Offered in the node menu only when `canRerunNode` returns true. */
  onRerunNode?: (nodeId: string) => void;
  /** Whether the "Rerun this node only" menu item should be shown for a node. */
  canRerunNode?: (nodeId: string) => boolean;
}

export type ConnectionStatus = 'idle' | 'active' | 'completed' | 'errored' | 'stale';

/** A copyable snapshot of a node set: the node payloads plus every connection
 *  whose BOTH endpoints are inside the set (data, pass-through, and order edges
 *  alike). Positions are the originals — materializing applies an offset.
 *  Deep-copied at snapshot time, so later edits to the originals never leak
 *  into a paste. */
interface GraphClip {
  nodes: Array<{
    id: string; typeName: string; label: string; description: string;
    collapsed: boolean; pos: {x: number; y: number};
    properties: Record<string, any>; inputValues: Record<string, any>;
  }>;
  connections: Array<{source: string; sourceOutput: string; target: string; targetInput: string}>;
}

export class FlowEditor {
  readonly editor = new NodeEditor<FlowScheme>();
  readonly area: AreaPlugin<FlowScheme>;
  readonly connection = new ConnectionPlugin<FlowScheme>();
  readonly render: ReactPlugin<FlowScheme, ReactArea2D<FlowScheme>>;
  readonly history = new HistoryPlugin<FlowScheme, HistoryActions<FlowScheme>>();
  readonly container: HTMLElement;
  /** Absolutely-positioned wrapper filling the host: a flex row of
   *  [canvasEl | Outputs strip]. */
  private readonly canvasWrap!: HTMLElement;
  /** Inner element the AreaPlugin mounts on — the actual canvas viewport. */
  readonly canvasEl!: HTMLElement;

  private selector = AreaExtensions.selector();
  /** Tracks which node is under the cursor at the moment of pointerdown.
   *  Read by `accumulating.active()` to decide whether to *preserve* an
   *  existing multi-selection when one of the already-selected nodes is
   *  clicked — required for KNIME/Figma-style group-drag. */
  private lastPointerDownNodeId: string | null = null;
  /** Mouse/pointer button last pressed (0 = primary, 2 = secondary).
   *  Several plugins (notably rete-connection-plugin) don't filter by button;
   *  right-clicking an output socket would otherwise start a fake connection
   *  drag that ends with `created:false`, accidentally triggering the
   *  suggestion menu. Handlers that should be left-click-only consult this. */
  private lastPointerButton = 0;
  /** The node most recently reported to the host as selected (`nodepicked` /
   *  chip click) — i.e. what the context panel currently shows. Re-picking it
   *  while it's still selected is a no-op and must not re-fire the host
   *  callbacks (panel/suggestion rebuilds on every click or grab). */
  private lastPickedId: string | null = null;
  /** Snapshot of the last pointerdown (modifiers, position, whether the node
   *  under the cursor was already selected), taken by
   *  `installPointerDownTracker` in the capture phase — before rete's
   *  `nodepicked` fires. Node clicks follow the platform's `selectRows`
   *  modifier convention (d4 `viewer_utils.dart`): plain click selects
   *  exclusively, Shift adds, Ctrl toggles, Ctrl+Shift removes. Rete's
   *  `nodepicked` can only ever ADD, so `accumulating` admits any modifier
   *  and the removal half runs on a clean release in the pointerup tracker. */
  private lastPointerDownWasSelected = false;
  private lastPointerDownModifier = false;
  private lastPointerDownPos = {x: 0, y: 0};
  private accumulating = {
    active: (): boolean => {
      if (this.lastPointerDownModifier) return true;
      const id = this.lastPointerDownNodeId;
      if (!id) return false;
      const node = this.editor.getNode(id) as {selected?: boolean} | undefined;
      return node?.selected === true;
    },
  };
  /** Returned by `AreaExtensions.selectableNodes` — gives us programmatic
   *  select/unselect on top of the click-to-select that the extension wires
   *  up automatically. Used by `selectNode` (e.g. for auto-select on
   *  run-complete) and the rectangle-select tool. */
  private selectableApi!: {
    select: (nodeId: string, accumulate: boolean) => Promise<void>;
    unselect: (nodeId: string) => Promise<void>;
  };
  private callbacks: FlowEditorCallbacks;
  private keydownHandler: ((e: KeyboardEvent) => void) | null = null;
  private pointerDownTracker: ((e: PointerEvent) => void) | null = null;
  private pointerUpTracker: ((e: PointerEvent) => void) | null = null;
  /** Per-connection status (for execution coloring). */
  private connectionStatuses = new Map<string, ConnectionStatus>();

  /** Workflow annotations — colored frames behind the graph (KNIME pattern).
   *  Owned by the editor (not by Rete), persisted alongside the graph. */
  private annotations = new Map<string, FlowAnnotation>();

  /** Node groups — collapsible frames around member node sets (see
   *  `node-group.ts`). Editor-level like annotations: the graph stays flat. */
  private groups = new Map<string, FlowGroup>();
  /** Wire-endpoint subscriptions for sockets on groupable (non-output) nodes,
   *  keyed by node id — while a node hides inside a minimized group, the
   *  editor pushes card-edge anchors through these instead of DOM positions
   *  (same pattern as `chipSocketSubs`). */
  private groupSocketSubs = new Map<string, Set<{
    side: 'input' | 'output'; key: string; cb: (pos: {x: number; y: number}) => void;
  }>>();
  /** Group ids with a frame refit already scheduled for the next frame. */
  private groupRefitScheduled = new Set<string>();

  /** Snap-to-grid step in canvas units; alignment guides override grid when within threshold. */
  private readonly gridSize = 20;
  /** Distance (canvas units) within which a node edge/center snaps to another node's edge/center. */
  private readonly alignThreshold = 6;
  private guideOverlay: HTMLElement | null = null;
  private vGuide: HTMLElement | null = null;
  private hGuide: HTMLElement | null = null;

  /** Bottom-right overview minimap (screen-space overlay; not part of the
   *  transformed canvas). `null` until `installMinimap`. */
  private minimapEl: HTMLElement | null = null;
  private minimapSvg: SVGSVGElement | null = null;
  private minimapRedrawScheduled = false;
  /** Minimap inner drawing area in px (SVG viewport). */
  private readonly minimapW = 200;
  private readonly minimapH = 130;

  /** The Outputs strip — a thin column OUTSIDE the canvas viewport that hosts
   *  every output node as a screen-space chip (our own DOM, not a rete node
   *  view — the canvas view is hidden). The nodes stay real graph citizens
   *  (data model, serialization, compiler untouched); only their visual form
   *  is the chip. `null` until `installOutputStrip`. */
  private outputStripEl: HTMLElement | null = null;
  /** Container inside the strip holding the chips. */
  private stripChipsEl: HTMLElement | null = null;
  private stripResizeObserver: ResizeObserver | null = null;
  private stripSyncScheduled = false;
  /** Whether the pending strip sync must rebuild the chip DOM (graph or
   *  selection changed) or only refresh wire endpoints (pan/zoom). */
  private stripRenderPending = false;
  /** Wire-endpoint subscriptions for sockets on output nodes, keyed by node id
   *  (see the chip-aware socketPositionWatcher in the constructor). */
  private chipSocketSubs = new Map<string, Set<(pos: {x: number; y: number}) => void>>();

  /** Suggestion-menu drag state. Set on `connectionpick` for an output
   *  socket; cleared on `connectiondrop`. If the drop didn't create a
   *  connection AND wasn't on a target socket, the suggestion popup opens. */
  private dragOutSource: {nodeId: string; outputKey: string; dgType: string} | null = null;

  /** Input-side drag state (dragging out of an input socket, or the tail of an
   *  existing connection). Drives the reverse drop-on-node shortcut — dropping
   *  on a node body connects from that node's compatible output (a real output
   *  wins over a passthrough) — and, on empty canvas, the reverse suggestion
   *  menu ("what produces this?"). */
  private dragInSource: {nodeId: string; inputKey: string; dgType: string} | null = null;

  constructor(container: HTMLElement, callbacks: FlowEditorCallbacks = {}) {
    this.callbacks = callbacks;
    this.container = container;
    // The editor splits into [canvas | Outputs strip]: the area plugin mounts
    // on an inner element, so the strip column is OUTSIDE the canvas viewport —
    // pan, zoom-to-fit, and drops can never put graph content behind it. The
    // pair lives in an absolutely-positioned wrapper (not host-level flex): the
    // host keeps its own normal-flow children (the view's start panel).
    this.canvasWrap = document.createElement('div');
    this.canvasWrap.className = 'ff-canvas-wrap';
    container.appendChild(this.canvasWrap);
    this.canvasEl = document.createElement('div');
    this.canvasEl.className = 'ff-canvas';
    setTid(this.canvasEl, 'canvas-viewport');
    this.canvasWrap.appendChild(this.canvasEl);
    this.area = new AreaPlugin<FlowScheme>(this.canvasEl);
    this.render = new ReactPlugin<FlowScheme, ReactArea2D<FlowScheme>>({createRoot});

    // Output nodes have NO canvas view — their visible form is a screen-space
    // chip inside the Outputs strip. A wire into one must still end somewhere,
    // so the DOM-measuring watcher is wrapped: sockets on output nodes resolve
    // analytically to the canvas' right edge at the chip's row (the wire runs
    // to the edge and visually plugs into the adjacent strip chip), refreshed
    // on pan/zoom/reorder via `notifyChipSockets`.
    const domWatcher = getDOMSocketPosition({
      // No arrow markers anymore — direction comes from the dash-flow CSS
      // animation. The line just needs to land on the dot edge, so a small
      // symmetric offset that puts both endpoints just inside the socket dot
      // (radius 4.5 px) keeps everything visually attached.
      offset: (pos, _id, side) => ({
        x: pos.x + (side === 'output' ? 2 : -2),
        y: pos.y,
      }),
    });
    const chipAwareWatcher = {
      attach: (scope: never) => (domWatcher as {attach(s: never): void}).attach(scope),
      listen: (nodeId: string, side: 'input' | 'output', key: string,
        onChange: (pos: {x: number; y: number}) => void): (() => void) => {
        if (this.editor.getNode(nodeId)?.dgNodeType === 'output')
          return this.listenChipSocket(nodeId, onChange);
        // Group-aware: while the node hides inside a minimized group, the
        // wire endpoint is the group card's edge, not the (display:none) DOM.
        return this.listenGroupableSocket(nodeId, side, key, onChange, domWatcher as unknown as {
          listen(n: string, s: string, k: string, cb: (p: {x: number; y: number}) => void): () => void;
        });
      },
    };

    this.render.addPreset(ReactPresets.classic.setup({
      socketPositionWatcher: chipAwareWatcher as never,
      customize: {
        node: () => FlowNodeComponent as never,
        socket: () => FlowSocketComponent as never,
        connection: () => FlowConnectionComponent as never,
      },
    }));

    this.installTypeValidation();

    this.editor.use(this.area);
    // Casts: Rete's `Scope.use` does a structural variance check that gets
    // pessimistic with our narrowed schemes. Runtime contracts are exact.
    this.area.use(this.connection as never);
    this.area.use(this.render as never);

    this.selectableApi = AreaExtensions.selectableNodes(this.area, this.selector, {
      accumulating: this.accumulating,
    });
    AreaExtensions.simpleNodesOrder(this.area);
    AreaExtensions.restrictor(this.area, {
      scaling: () => ({min: 0.2, max: 2.5}),
    });

    // Undo/redo: history-plugin tracks add/remove/drag of nodes & connections.
    this.history.addPreset(HistoryPresets.classic.setup());
    this.area.use(this.history as never);
    HistoryExtensions.keyboard(this.history); // Ctrl+Z / Ctrl+Shift+Z / Ctrl+Y

    this.installPointerDownTracker();
    this.wireEvents();
    this.installGuideOverlay();
    this.installContextMenu();
    this.installKeyboardShortcuts();
    this.installDoubleClickToFit();
    this.installSuggestionMenu();
    this.installRectSelect();
    this.installHoverDocs();
    this.installWaypointInteractions();
    this.installMinimap();
    this.installOutputStrip();
  }

  /** Narrow callback surface for the React node components, stamped onto every
   *  node this editor owns (`FlowNode.editorBridge` — see the `nodecreate` pipe
   *  in `wireEvents`). Resolving it from the node instead of a page-level
   *  global keeps each component bound to its own editor: several editors
   *  coexist on a page (file previews, the creation-script dialog, detached
   *  compile editors), and a global bridge that any construction rebinds and
   *  any `destroy()` deletes broke collapse toggling and collapsed-socket
   *  rendering in whichever editor didn't own it last. */
  private readonly bridge: FlowEditorBridge = {
    toggleCollapsed: (id) => void this.toggleCollapsed(id),
    isSocketConnected: (nodeId, side, key) => this.isSocketConnected(nodeId, side, key),
  };

  /** Configure ClassicFlow to reject incompatible socket connections at pick
   *  time, before any connection ever enters the editor's data layer. */
  private installTypeValidation(): void {
    this.connection.addPreset(() =>
      new ClassicFlow<FlowScheme, never[]>({
        canMakeConnection: (from, to) => this.canConnect(from, to),
        makeConnection: (from, to) => {
          const [out, inp] = from.side === 'output' ? [from, to] : [to, from];
          if (out.side !== 'output' || inp.side !== 'input') return false;
          const outNode = this.editor.getNode(out.nodeId);
          const inNode = this.editor.getNode(inp.nodeId);
          if (!outNode || !inNode) return false;
          void this.editor.addConnection(new FlowConnection(
            outNode as never, out.key, inNode as never, inp.key,
          ));
          return true;
        },
      }),
    );
  }

  /** When a connection lands on a ValueOutput node and the source slot has a
   *  meaningful type, copy that type into the output node's `outputType`. */
  private maybeAutoTypeValueOutput(connection: FlowScheme['Connection']): void {
    // Execution-ordering edges carry no data type — never derive an output type from one.
    if (isExecKey(String(connection.targetInput)) || isExecKey(String(connection.sourceOutput))) return;
    const targetNode = this.editor.getNode(connection.target) as FlowNode | undefined;
    // Match by registered type, not label — titles are user-editable.
    if (!targetNode || targetNode.dgTypeName !== 'Outputs/Value Output') return;
    const sourceNode = this.editor.getNode(connection.source) as FlowNode | undefined;
    if (!sourceNode) return;
    const sourceOutput = sourceNode.outputs[String(connection.sourceOutput)] as
      {socket: TypedSocket} | undefined;
    if (!sourceOutput) return;
    const detected = sourceOutput.socket.dgType;
    if (detected && detected !== 'dynamic' && detected !== 'object') {
      targetNode.properties['outputType'] = detected;
      void this.area.update('node', targetNode.id);
    }
  }

  /** Collapsed nodes render socket DOM only for *connected* sockets (see
   *  node-component.tsx). A connection created or removed while an endpoint is
   *  collapsed changes which sockets must exist, so re-render those nodes.
   *  Without this, a connection added to an already-collapsed node (creation-
   *  script import, .ffjson load) has no socket element to attach to and stays
   *  invisible until the node is expanded and collapsed again. */
  private refreshCollapsedEndpoints(conn: FlowScheme['Connection']): void {
    for (const id of [conn.source, conn.target]) {
      const node = this.editor.getNode(id);
      if (node?.collapsed) void this.area.update('node', id);
    }
  }

  private canConnect(from: SocketData, to: SocketData): boolean {
    if (from.nodeId === to.nodeId) return false;
    if (from.side === to.side) return false;
    const [out, inp] = from.side === 'output' ? [from, to] : [to, from];
    const outNode = this.editor.getNode(out.nodeId);
    const inNode = this.editor.getNode(inp.nodeId);
    if (!outNode || !inNode) return false;
    const outSocket = outNode.outputs[out.key]?.socket as TypedSocket | undefined;
    const inSocket = inNode.inputs[inp.key]?.socket as TypedSocket | undefined;
    if (!outSocket || !inSocket) return false;
    return outSocket.isCompatibleWith(inSocket);
  }

  private wireEvents(): void {
    this.editor.addPipe((context) => {
      // Stamp the owning editor's bridge BEFORE the node ever renders, so the
      // React node component always talks back to this editor (not a global).
      if (context.type === 'nodecreate')
        context.data.editorBridge = this.bridge;
      // Stamp `_color` on every new connection BEFORE the area-plugin emits
      // 'render', so the React Connection component picks up the right color
      // on its very first render.
      if (context.type === 'connectioncreate')
        this.decorateConnection(context.data);
      if (context.type === 'noderemoved') {
        this.chipSocketSubs.delete(context.data.id);
        this.groupSocketSubs.delete(context.data.id);
        this.handleGroupMemberRemoved(context.data.id);
      }
      if (context.type === 'connectioncreated')
        this.maybeAutoTypeValueOutput(context.data);
      if (context.type === 'connectionremoved')
        this.connectionStatuses.delete(context.data.id);
      if (context.type === 'connectioncreated' || context.type === 'connectionremoved') {
        this.refreshCollapsedEndpoints(context.data);
        // A wire into/out of a minimized group changes its boundary dot rows.
        for (const id of [context.data.source, context.data.target]) {
          const g = this.minimizedGroupOf(id);
          if (g) this.refreshGroupCard(g);
        }
      }
      if (
        context.type === 'nodecreated' || context.type === 'noderemoved' ||
        context.type === 'connectioncreated' || context.type === 'connectionremoved' ||
        context.type === 'cleared'
      ) {
        this.callbacks.onGraphChanged?.();
        this.callbacks.onGraphEdited?.(this.classifyEdit(context));
        this.scheduleMinimapRedraw();
        this.scheduleStripSync(true);
      }
      return context;
    });

    this.area.addPipe((context) => {
      if (context.type === 'nodepicked') {
        const node = this.editor.getNode(context.data.id);
        if (node) {
          // Re-picking the node that is ALREADY the current object (click it
          // again, grab it to drag) changes nothing — don't make the host
          // rebuild its panels. `lastPointerDownWasSelected` is the state
          // snapshot from BEFORE rete's add-only pick, so a click that
          // re-selects after a deselect-all still fires.
          const samePick = node.id === this.lastPickedId && this.lastPointerDownWasSelected;
          if (this.lastPickedId && this.lastPickedId !== node.id) {
            const prev = this.editor.getNode(this.lastPickedId);
            if (prev) this.callbacks.onNodeDeselected?.(prev);
          }
          this.lastPickedId = node.id;
          if (!samePick) this.callbacks.onNodeSelected?.(node);
          this.refreshChipSelection(); // chip selected-state may have changed
        }
      }
      // Intercept the during-drag translate intent: snap position and show guides.
      // Mutating `data.position` in-place propagates to the actual translate.
      // Only snap the *picked* node — when several nodes are selected and one
      // is dragged, the selectable extension translates the rest by the same
      // delta to preserve relative offsets. Snapping each follower would
      // recompute their positions independently and break the group geometry.
      if (context.type === 'nodetranslate') {
        const data = context.data as {id: string; position: {x: number; y: number}};
        if (this.selector.isPicked({id: data.id, label: 'node'})) {
          const snap = this.computeSnap(data.id, data.position);
          data.position.x = snap.x;
          data.position.y = snap.y;
          this.showGuides(snap.guideX, snap.guideY);
        }
      }
      if (context.type === 'nodetranslated') {
        const node = this.editor.getNode(context.data.id);
        if (node) node.pos = {...context.data.position};
        // An expanded group's frame hugs its members — refit when one moves.
        const g = this.groupOf(context.data.id);
        if (g && !g.minimized) this.scheduleGroupRefit(g);
        this.scheduleMinimapRedraw();
      }
      if (context.type === 'nodedragged') this.hideGuides();
      // Keep the CSS dot-grid background aligned to the area transform — the
      // grid is screen-space, the canvas content lives in a transformed
      // space, so we rescale and shift the bg whenever pan/zoom changes.
      if (context.type === 'translated' || context.type === 'zoomed' ||
          context.type === 'render') {
        this.updateGridTransform();
        this.scheduleMinimapRedraw();
        this.scheduleStripSync(false); // wire endpoints track the transform
      }
      // Tag each rendered connection wrapper with its id + status for the
      // CSS-driven execution-state animations (`[data-status="active"]` etc).
      if (context.type === 'rendered' && (context.data as {type?: string}).type === 'connection')
        this.tagConnectionElement(context.data as {element: HTMLElement; payload: FlowConnection});
      return context;
    });
  }

  /** Map a rete editor event to the classified {@link GraphEdit} handed to
   *  `onGraphEdited`. Only called for the five event types listed in
   *  `wireEvents` — anything else would be a programming error. */
  private classifyEdit(context:
    | {type: 'nodecreated' | 'noderemoved'; data: {id: string}}
    | {type: 'connectioncreated' | 'connectionremoved'; data: {source: string; target: string}}
    | {type: 'cleared'},
  ): GraphEdit {
    switch (context.type) {
    case 'nodecreated': return {kind: 'node-added', nodeId: context.data.id};
    case 'noderemoved': return {kind: 'node-removed', nodeId: context.data.id};
    case 'connectioncreated':
      return {kind: 'connection-added', sourceId: context.data.source, targetId: context.data.target};
    case 'connectionremoved':
      return {kind: 'connection-removed', sourceId: context.data.source, targetId: context.data.target};
    default: return {kind: 'cleared'};
    }
  }

  /** Report that a node's parameters (its `inputValues` / `properties`) were
   *  edited — the property panel calls this so run results downstream of the
   *  node can be invalidated. Cosmetic edits (title, description) must NOT be
   *  reported. */
  notifyNodeParamsChanged(nodeId: string): void {
    this.callbacks.onGraphEdited?.({kind: 'params-changed', nodeId});
  }

  /** Update the canvas dot-grid background to track the AreaPlugin transform.
   *  Without this, panning would drift the dots out of alignment with snapped
   *  nodes. Called from translate/zoom events in the area pipe. */
  private updateGridTransform(): void {
    const t = this.area.area.transform;
    const size = 20 * t.k;
    this.container.style.backgroundSize = `${size}px ${size}px`;
    this.container.style.backgroundPosition = `${t.x}px ${t.y}px`;
  }

  // ---------- snap + alignment guides ----------

  /** Build a screen-space overlay layer for the alignment guides. The two
   *  `<div>` rules act as 1px dashed crosshair lines that we move on demand. */
  private installGuideOverlay(): void {
    if (this.container.style.position === '') this.container.style.position = 'relative';
    const overlay = document.createElement('div');
    overlay.className = 'ff-guide-overlay';
    setTid(overlay, 'guide-overlay');
    this.guideOverlay = overlay;
    const v = document.createElement('div');
    v.className = 'ff-guide-v';
    overlay.appendChild(v);
    this.vGuide = v;
    const h = document.createElement('div');
    h.className = 'ff-guide-h';
    overlay.appendChild(h);
    this.hGuide = h;
    this.container.appendChild(overlay);
  }

  /** Measure a node's rendered size in canvas units (DOM / scale). */
  private measureNode(id: string): {w: number; h: number} {
    // `nodeViews` is a public Map<id, {element}> on AreaPlugin.
    const views = (this.area as unknown as {nodeViews: Map<string, {element: HTMLElement}>}).nodeViews;
    const view = views?.get(id);
    if (!view) return {w: 220, h: 80};
    const r = view.element.getBoundingClientRect();
    const k = this.area.area.transform.k || 1;
    return {w: r.width / k, h: r.height / k};
  }

  /** Given the dragged node's id and its proposed position, decide where it
   *  should actually land. Alignment with another node's edge/center wins
   *  over the grid; otherwise the position snaps to the grid. */
  private computeSnap(
    draggedId: string, pos: {x: number; y: number},
  ): {x: number; y: number; guideX: number | null; guideY: number | null} {
    const dragged = this.measureNode(draggedId);
    const dxFrom = [pos.x, pos.x + dragged.w / 2, pos.x + dragged.w]; // left, cx, right
    const dyFrom = [pos.y, pos.y + dragged.h / 2, pos.y + dragged.h]; // top, cy, bottom

    let bestX: {delta: number; guide: number} | null = null;
    let bestY: {delta: number; guide: number} | null = null;

    for (const other of this.editor.getNodes()) {
      // Strip-pinned rows sit at viewport-dependent positions — never a
      // meaningful alignment target for canvas nodes. Neither are nodes
      // hidden inside a minimized group.
      if (other.id === draggedId || other.dgNodeType === 'output' ||
          this.minimizedGroupOf(other.id)) continue;
      const sz = this.measureNode(other.id);
      const ox = [other.pos.x, other.pos.x + sz.w / 2, other.pos.x + sz.w];
      const oy = [other.pos.y, other.pos.y + sz.h / 2, other.pos.y + sz.h];
      for (const f of dxFrom) for (const t of ox) {
        const d = t - f;
        if (Math.abs(d) <= this.alignThreshold && (!bestX || Math.abs(d) < Math.abs(bestX.delta)))
          bestX = {delta: d, guide: t};
      }
      for (const f of dyFrom) for (const t of oy) {
        const d = t - f;
        if (Math.abs(d) <= this.alignThreshold && (!bestY || Math.abs(d) < Math.abs(bestY.delta)))
          bestY = {delta: d, guide: t};
      }
    }

    const x = bestX ? pos.x + bestX.delta : Math.round(pos.x / this.gridSize) * this.gridSize;
    const y = bestY ? pos.y + bestY.delta : Math.round(pos.y / this.gridSize) * this.gridSize;
    return {x, y, guideX: bestX?.guide ?? null, guideY: bestY?.guide ?? null};
  }

  /** Position the dashed guide lines (canvas-space coords → screen-space px). */
  private showGuides(canvasX: number | null, canvasY: number | null): void {
    if (!this.vGuide || !this.hGuide) return;
    const t = this.area.area.transform;
    if (canvasX !== null) {
      this.vGuide.style.left = `${t.x + canvasX * t.k}px`;
      this.vGuide.style.display = 'block';
    } else this.vGuide.style.display = 'none';
    if (canvasY !== null) {
      this.hGuide.style.top = `${t.y + canvasY * t.k}px`;
      this.hGuide.style.display = 'block';
    } else this.hGuide.style.display = 'none';
  }

  private hideGuides(): void {
    if (this.vGuide) this.vGuide.style.display = 'none';
    if (this.hGuide) this.hGuide.style.display = 'none';
  }

  // ---------- minimap ----------

  /** Build the bottom-right overview minimap: an SVG that draws every node as a
   *  small rect plus the current viewport rectangle. Click/drag inside it pans
   *  the canvas; the header button minimizes it to a title bar. Lives directly
   *  in `this.container` (screen-space), so it doesn't pan/zoom with the graph. */
  private installMinimap(): void {
    const el = document.createElement('div');
    el.className = 'ff-minimap';
    el.dataset.collapsed = 'false';
    setTid(el, 'minimap');

    const header = document.createElement('div');
    header.className = 'ff-minimap-header';
    setTid(header, 'minimap-header');
    const title = document.createElement('span');
    title.className = 'ff-minimap-title';
    title.textContent = 'Overview';
    const toggle = document.createElement('span');
    toggle.className = 'ff-minimap-toggle';
    setTid(toggle, 'minimap-toggle');
    toggle.title = 'Minimize';
    toggle.textContent = '▾';
    header.appendChild(title);
    header.appendChild(toggle);

    const body = document.createElement('div');
    body.className = 'ff-minimap-body';
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('class', 'ff-minimap-svg');
    svg.setAttribute('width', String(this.minimapW));
    svg.setAttribute('height', String(this.minimapH));
    svg.setAttribute('viewBox', `0 0 ${this.minimapW} ${this.minimapH}`);
    body.appendChild(svg);

    el.appendChild(header);
    el.appendChild(body);
    // Inside the canvas viewport (not the host) — the host's right column is
    // the Outputs strip, which the minimap must not cover.
    this.canvasEl.appendChild(el);
    this.minimapEl = el;
    this.minimapSvg = svg;

    // Clicking anywhere on the header minimizes/restores. stopPropagation so it
    // never pans the canvas. (The chevron is just a visual affordance — the
    // click bubbles to the header handler.)
    header.addEventListener('pointerdown', (e) => e.stopPropagation());
    header.addEventListener('click', (e) => {
      e.stopPropagation();
      this.toggleMinimapCollapsed();
    });

    this.installMinimapNavigation(body);
    this.scheduleMinimapRedraw();
  }

  /** Click/drag in the minimap body → pan so the clicked graph point centers in
   *  the viewport. We map minimap px → canvas coords using the same fit the
   *  draw step computed (stored on the svg via data-* for reuse). */
  private installMinimapNavigation(body: HTMLElement): void {
    const panToEvent = (e: PointerEvent): void => {
      const svg = this.minimapSvg;
      if (!svg) return;
      const fit = this.readMinimapFit();
      if (!fit) return;
      const rect = svg.getBoundingClientRect();
      const mmx = e.clientX - rect.left;
      const mmy = e.clientY - rect.top;
      // minimap px → canvas coords (inverse of the draw transform)
      const cx = (mmx - fit.offsetX) / fit.scale + fit.minX;
      const cy = (mmy - fit.offsetY) / fit.scale + fit.minY;
      const cont = this.container.getBoundingClientRect();
      const k = this.area.area.transform.k;
      void this.area.area.translate(cont.width / 2 - cx * k, cont.height / 2 - cy * k);
    };

    body.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      e.preventDefault();
      e.stopPropagation();
      panToEvent(e);
      body.setPointerCapture(e.pointerId);
      const onMove = (ev: PointerEvent): void => panToEvent(ev);
      const onUp = (): void => {
        body.removeEventListener('pointermove', onMove);
        body.removeEventListener('pointerup', onUp);
        body.removeEventListener('pointercancel', onUp);
      };
      body.addEventListener('pointermove', onMove);
      body.addEventListener('pointerup', onUp);
      body.addEventListener('pointercancel', onUp);
    });
  }

  /** Collapse the minimap to its header bar, or restore it. Public so hosts can
   *  set the initial state (e.g. collapsed inside a preview dialog). */
  setMinimapCollapsed(collapsed: boolean): void {
    const el = this.minimapEl;
    if (!el) return;
    el.dataset.collapsed = collapsed ? 'true' : 'false';
    const toggle = el.querySelector<HTMLElement>('.ff-minimap-toggle');
    if (toggle) {
      toggle.textContent = collapsed ? '▸' : '▾';
      toggle.title = collapsed ? 'Expand' : 'Minimize';
    }
    if (!collapsed) this.scheduleMinimapRedraw();
  }

  private toggleMinimapCollapsed(): void {
    if (this.minimapEl) this.setMinimapCollapsed(this.minimapEl.dataset.collapsed !== 'true');
  }

  /** Read the fit transform stashed on the svg by the last redraw. */
  private readMinimapFit(): {scale: number; offsetX: number; offsetY: number; minX: number; minY: number} | null {
    const svg = this.minimapSvg;
    if (!svg || svg.dataset.scale === undefined) return null;
    return {
      scale: parseFloat(svg.dataset.scale),
      offsetX: parseFloat(svg.dataset.offsetX!),
      offsetY: parseFloat(svg.dataset.offsetY!),
      minX: parseFloat(svg.dataset.minX!),
      minY: parseFloat(svg.dataset.minY!),
    };
  }

  /** Coalesce minimap redraws to one per frame — the area pipe fires many
   *  translate/render events during a single drag. */
  /** Re-evaluate the overview after a graph edit (visibility + redraw). */
  refreshMinimap(): void {
    this.scheduleMinimapRedraw();
  }

  private scheduleMinimapRedraw(): void {
    if (!this.minimapEl) return;
    // Nothing to overview on an empty canvas — hide the panel entirely (it
    // reappears the moment the first node lands). Done before the collapsed/
    // scheduled early-returns so an empty canvas hides regardless of either.
    this.minimapEl.style.display = this.editor.getNodes().length === 0 ? 'none' : '';
    if (this.minimapRedrawScheduled) return;
    if (this.minimapEl.dataset.collapsed === 'true') return;
    this.minimapRedrawScheduled = true;
    requestAnimationFrame(() => {
      this.minimapRedrawScheduled = false;
      this.redrawMinimap();
    });
  }

  private redrawMinimap(): void {
    const svg = this.minimapSvg;
    if (!svg) return;
    while (svg.firstChild) svg.removeChild(svg.firstChild);

    const nodes = this.editor.getNodes();
    const pad = 8;
    if (nodes.length === 0) {
      delete svg.dataset.scale;
      return;
    }

    // Graph bounds in canvas coords.
    let minX = Infinity; let minY = Infinity; let maxX = -Infinity; let maxY = -Infinity;
    const boxes: Array<{x: number; y: number; w: number; h: number; color: string}> = [];
    for (const node of nodes) {
      // Strip-pinned rows track the viewport, not the graph — including them
      // would smear the overview bounds on every pan. Members hidden inside a
      // minimized group are drawn as their group's card below.
      if (node.dgNodeType === 'output' || this.minimizedGroupOf(node.id)) continue;
      const sz = this.measureNode(node.id);
      const color = (node as unknown as {color?: string}).color ?? '#90a4ae';
      boxes.push({x: node.pos.x, y: node.pos.y, w: sz.w, h: sz.h, color});
      minX = Math.min(minX, node.pos.x);
      minY = Math.min(minY, node.pos.y);
      maxX = Math.max(maxX, node.pos.x + sz.w);
      maxY = Math.max(maxY, node.pos.y + sz.h);
    }
    for (const g of this.groups.values()) {
      if (!g.minimized) continue;
      const w = g.element.offsetWidth || 180; const h = g.element.offsetHeight || 40;
      boxes.push({x: g.pos.x, y: g.pos.y, w, h, color: '#607d8b'});
      minX = Math.min(minX, g.pos.x);
      minY = Math.min(minY, g.pos.y);
      maxX = Math.max(maxX, g.pos.x + w);
      maxY = Math.max(maxY, g.pos.y + h);
    }
    if (!Number.isFinite(minX)) {
      delete svg.dataset.scale;
      return;
    }

    const graphW = Math.max(1, maxX - minX);
    const graphH = Math.max(1, maxY - minY);
    const scale = Math.min((this.minimapW - 2 * pad) / graphW, (this.minimapH - 2 * pad) / graphH);
    const offsetX = (this.minimapW - graphW * scale) / 2 - minX * scale;
    const offsetY = (this.minimapH - graphH * scale) / 2 - minY * scale;
    svg.dataset.scale = String(scale);
    svg.dataset.offsetX = String(offsetX);
    svg.dataset.offsetY = String(offsetY);
    svg.dataset.minX = String(minX);
    svg.dataset.minY = String(minY);

    const toMapX = (cx: number): number => cx * scale + offsetX;
    const toMapY = (cy: number): number => cy * scale + offsetY;

    for (const b of boxes) {
      const r = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      r.setAttribute('x', String(toMapX(b.x)));
      r.setAttribute('y', String(toMapY(b.y)));
      r.setAttribute('width', String(Math.max(2, b.w * scale)));
      r.setAttribute('height', String(Math.max(2, b.h * scale)));
      r.setAttribute('rx', '1.5');
      r.setAttribute('fill', b.color);
      r.setAttribute('class', 'ff-minimap-node');
      svg.appendChild(r);
    }

    // Viewport rectangle (canvas coords visible through the container).
    const t = this.area.area.transform;
    const cont = this.container.getBoundingClientRect();
    const vx = -t.x / t.k; const vy = -t.y / t.k;
    const vw = cont.width / t.k; const vh = cont.height / t.k;
    const vp = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    vp.setAttribute('x', String(toMapX(vx)));
    vp.setAttribute('y', String(toMapY(vy)));
    vp.setAttribute('width', String(Math.max(0, vw * scale)));
    vp.setAttribute('height', String(Math.max(0, vh * scale)));
    vp.setAttribute('class', 'ff-minimap-viewport');
    svg.appendChild(vp);
  }

  // ---------- output strip ----------

  /** Chip layout constants (screen px — chips are plain DOM inside the strip
   *  column, so zoom never touches them; keep in sync with the
   *  `.ff-output-strip-chips` / `.ff-output-row` CSS). `GAP`/`H` also drive
   *  the analytic wire-endpoint math in {@link chipSocketPos}. */
  private static readonly STRIP_CHIP_GAP = 6;
  private static readonly STRIP_CHIP_H = 24;

  /** Build the strip column: a thin flex sibling to the RIGHT of the canvas
   *  viewport — graph content can never pan or fit behind it. Hosts the chips
   *  (one per output node) above a vertical "Outputs" label. Chip interaction
   *  is delegated here: click selects the node (property panel, Delete key),
   *  right-click opens the node context menu; drops are hit-tested by rect in
   *  `handleOutputDrop` (chips carry `.ff-node`, so the drop-on-node branch
   *  binds a chip's free input first). */
  private installOutputStrip(): void {
    const strip = document.createElement('div');
    strip.className = 'ff-output-strip';
    setTid(strip, 'output-strip');
    strip.title = 'Flow outputs — drag any output socket here to publish its value';

    const chips = document.createElement('div');
    chips.className = 'ff-output-strip-chips';
    strip.appendChild(chips);
    this.stripChipsEl = chips;

    const header = document.createElement('div');
    header.className = 'ff-output-strip-header';
    header.textContent = 'Outputs';
    strip.appendChild(header);

    strip.addEventListener('click', (ev) => {
      const chip = (ev.target as HTMLElement | null)?.closest('[data-node-id]') as HTMLElement | null;
      const id = chip?.dataset.nodeId;
      if (!id) return;
      // Re-clicking the chip that is already the sole-selected current object
      // changes nothing — don't re-fire the host callbacks (panel rebuilds).
      const accumulate = ev.ctrlKey || ev.metaKey;
      const node = this.editor.getNode(id) as {selected?: boolean} | undefined;
      if (node?.selected && this.lastPickedId === id &&
          (accumulate || this.getSelectedNodeIds().length === 1)) return;
      void this.selectNode(id, accumulate);
    });
    strip.addEventListener('contextmenu', (ev) => {
      const chip = (ev.target as HTMLElement | null)?.closest('[data-node-id]') as HTMLElement | null;
      const node = chip?.dataset.nodeId ? this.editor.getNode(chip.dataset.nodeId) : undefined;
      if (!node) return;
      ev.preventDefault();
      ev.stopPropagation();
      this.showNodeContextMenu(ev, node);
    });

    this.canvasWrap.appendChild(strip); // after canvasEl → right column
    this.outputStripEl = strip;

    this.stripResizeObserver = new ResizeObserver(() => this.scheduleStripSync(false));
    this.stripResizeObserver.observe(this.canvasEl);
    this.scheduleStripSync(true);
  }

  /** Coalesce strip work to one pass per event-loop tick (a pan emits many
   *  transform events per pointermove). `render: true` also rebuilds the chip
   *  DOM (graph / selection / params changed); `false` only refreshes the wire
   *  endpoints (pan/zoom/resize). Microtask — not rAF — so endpoints update
   *  within the same frame and wires never visibly lag. */
  private scheduleStripSync(render: boolean): void {
    if (render) this.stripRenderPending = true;
    if (!this.outputStripEl || this.stripSyncScheduled) return;
    this.stripSyncScheduled = true;
    queueMicrotask(() => {
      this.stripSyncScheduled = false;
      const doRender = this.stripRenderPending;
      this.stripRenderPending = false;
      this.syncOutputStrip(doRender);
    });
  }

  private syncOutputStrip(render: boolean): void {
    const strip = this.outputStripEl;
    const chips = this.stripChipsEl;
    if (!strip || !chips) return;
    const rows = this.outputNodes();
    strip.dataset.empty = rows.length === 0 ? 'true' : 'false';
    if (render) {
      chips.textContent = '';
      for (const n of rows) chips.appendChild(this.buildChip(n));
    }
    this.notifyChipSockets();
  }

  private outputNodes(): FlowNode[] {
    return this.editor.getNodes().filter((n) => n.dgNodeType === 'output');
  }

  /** Update every chip's `data-selected` IN PLACE. Deliberately not a rebuild:
   *  selection changes fire mid-click-gesture (pointerup), and replacing the
   *  pressed chip element there would keep the browser from dispatching its
   *  `click` — the very event that selects the node. Deferred a microtask so
   *  the (async) selectable-extension calls have landed on `node.selected`. */
  private refreshChipSelection(): void {
    queueMicrotask(() => {
      const chips = this.stripChipsEl;
      if (!chips) return;
      for (const el of Array.from(chips.children) as HTMLElement[]) {
        const node = el.dataset.nodeId ? this.editor.getNode(el.dataset.nodeId) : undefined;
        el.dataset.selected = (node as {selected?: boolean} | undefined)?.selected ? 'true' : 'false';
      }
    });
  }

  /** One chip: [socket dot | type letter], fixed 40×24, screen-space. Carries
   *  the same identity attributes as a canvas node card (`.ff-node`,
   *  `data-node-id`, `data-node-type-name`, `data-selected`, the
   *  `socket-input` test-id), so guides, connect hints, the drop-on-node
   *  branch, and the tests address chips exactly like nodes. */
  private buildChip(node: FlowNode): HTMLElement {
    const inputKeys = Object.keys(node.inputs).filter((k) => !isExecKey(k));
    const boundKey = inputKeys.find((k) => this.isSocketConnected(node.id, 'input', k));
    const paramName = String(node.properties['paramName'] ?? '') || node.label;
    const typeText = String(node.properties['outputType'] ?? node.dgOutputType ?? 'dynamic');
    const dgStatus = (node as unknown as {dgStatus?: string}).dgStatus ?? 'idle';
    const statusText = (node as unknown as {statusText?: string}).statusText ?? '';

    const el = document.createElement('div');
    el.className = 'ff-node ff-output-row';
    setTid(el, 'node');
    el.dataset.nodeId = node.id;
    el.dataset.nodeType = 'output';
    el.dataset.nodeTypeName = node.dgTypeName ?? '';
    el.dataset.nodeLabel = node.label;
    el.dataset.selected = (node as {selected?: boolean}).selected ? 'true' : 'false';
    el.dataset.status = dgStatus;
    el.dataset.bound = boundKey ? 'true' : 'false';

    let source: string | null = null;
    if (boundKey) {
      const src = this.getInputSource(node.id, boundKey);
      if (src) {
        const slot = src.node.outputs[src.outputKey] as {label?: string} | undefined;
        // A pass-through's label is just '→' — name it after its input instead.
        source = `${src.node.label} › ${src.outputKey.endsWith('__pt') ?
          `${src.outputKey.slice(0, -'__pt'.length)} (pass-through)` :
          (slot?.label ?? src.outputKey)}`;
      }
    }
    el.title = `Flow output "${paramName}" (${typeText})\n` +
      (source ? `← ${source}` : 'Not connected yet — wire any output socket into it') +
      (statusText ? `\n${statusText}` : '');

    for (const key of inputKeys) {
      const wrap = document.createElement('span');
      wrap.className = 'ff-output-row-socket';
      wrap.dataset.testid = tid('socket-input', key);
      const dot = document.createElement('span');
      dot.className = 'ff-socket ff-output-chip-socket';
      const slot = node.inputs[key] as {socket: TypedSocket} | undefined;
      dot.style.setProperty('--socket-color', getSlotColor(slot?.socket.dgType ?? 'dynamic'));
      wrap.appendChild(dot);
      el.appendChild(wrap);
    }

    const letter = document.createElement('span');
    letter.className = 'ff-output-row-letter';
    letter.dataset.testid = tid('output-row-letter');
    letter.style.color = getSlotColor(typeText);
    letter.textContent = getSlotLetter(typeText);
    el.appendChild(letter);
    return el;
  }

  /** Canvas-coord wire endpoint for any socket on an output node: the canvas'
   *  right edge, at the vertical center of that node's chip. Chips are
   *  vertically centered as a group (flex `justify-content: center`), so chip
   *  i's center sits at the strip's middle, offset by its index's distance
   *  from the group middle. The wire runs to the edge and visually plugs into
   *  the adjacent strip chip. */
  private chipSocketPos(nodeId: string): {x: number; y: number} {
    const rows = this.outputNodes();
    const i = Math.max(0, rows.findIndex((n) => n.id === nodeId));
    const t = this.area.area.transform;
    const yPx = this.canvasEl.clientHeight / 2 +
      (i - (rows.length - 1) / 2) * (FlowEditor.STRIP_CHIP_H + FlowEditor.STRIP_CHIP_GAP);
    return {x: (this.canvasEl.clientWidth - t.x) / t.k, y: (yPx - t.y) / t.k};
  }

  /** `socketPositionWatcher.listen` for sockets on output nodes (see the
   *  constructor): emit the analytic chip endpoint now and on every
   *  {@link notifyChipSockets}. */
  private listenChipSocket(nodeId: string, onChange: (pos: {x: number; y: number}) => void): () => void {
    let subs = this.chipSocketSubs.get(nodeId);
    if (!subs) {
      subs = new Set();
      this.chipSocketSubs.set(nodeId, subs);
    }
    subs.add(onChange);
    onChange(this.chipSocketPos(nodeId));
    return () => void this.chipSocketSubs.get(nodeId)?.delete(onChange);
  }

  /** Push fresh endpoints to every subscribed wire — called from the strip
   *  sync (transform changes, chip reorders, canvas resize). */
  private notifyChipSockets(): void {
    for (const [nodeId, subs] of this.chipSocketSubs) {
      if (subs.size === 0) continue;
      const pos = this.chipSocketPos(nodeId);
      for (const cb of subs) cb(pos);
    }
  }

  /** Whether a client-space point lies inside the strip column. Used by the
   *  drop handler. */
  private stripContains(clientX: number, clientY: number): boolean {
    const strip = this.outputStripEl;
    if (!strip) return false;
    const r = strip.getBoundingClientRect();
    return clientX >= r.left && clientX <= r.right && clientY >= r.top && clientY <= r.bottom;
  }

  /** An output-socket drag dropped on the strip: publish that value as a flow
   *  output. Creates the matching output node (Table Output for a dataframe,
   *  Value Output otherwise — its declared type auto-set by
   *  `maybeAutoTypeValueOutput` on connect), names it after the source slot,
   *  and wires it up; the `nodecreated` strip sync renders the chip. */
  private async bindOutputToStrip(src: {nodeId: string; outputKey: string; dgType: string}): Promise<void> {
    const {createNode} = await import('./node-factory');
    const isTable = src.dgType === 'dataframe';
    const node = createNode(isTable ? 'Outputs/Table Output' : 'Outputs/Value Output');
    if (!node) return;
    node.properties['paramName'] = this.uniqueOutputParamName(src.outputKey);
    await this.editor.addNode(node);
    await this.addConnectionByKeys(src.nodeId, src.outputKey, node.id, isTable ? 'table' : 'value');
  }

  /** A script-identifier param name derived from the source slot key (`__pt`
   *  suffix stripped), made unique among the existing output nodes' names. */
  private uniqueOutputParamName(outputKey: string): string {
    let base = outputKey.replace(/__pt$/, '').replace(/[^a-zA-Z0-9_]/g, '');
    if (!/^[a-zA-Z_]/.test(base)) base = 'result';
    const taken = new Set(this.editor.getNodes()
      .filter((n) => n.dgNodeType === 'output')
      .map((n) => String(n.properties['paramName'] ?? '')));
    if (!taken.has(base)) return base;
    for (let i = 2; ; i++)
      if (!taken.has(`${base}${i}`)) return `${base}${i}`;
  }

  // ---------- drag-output-to-empty suggestion menu ----------

  /** Hook into the connection plugin's own `connectionpick` / `connectiondrop`
   *  signals — these fire reliably regardless of how the plugin captures
   *  pointer events. (Earlier we tried watching raw DOM pointerup, but the
   *  connection plugin captures the gesture first and the bubbling never
   *  reached our handler.) `connectiondrop` carries `created: boolean` so
   *  we know exactly when a real connection didn't happen. */
  private installSuggestionMenu(): void {
    let lastPointer = {x: 0, y: 0};
    const trackPointer = (e: PointerEvent): void => {
      lastPointer = {x: e.clientX, y: e.clientY};
    };
    // Always track — cheap, and lets us fall back without the listener
    // dance per-pick.
    window.addEventListener('pointermove', trackPointer, true);
    // Safety net: a pick that ends without a `connectiondrop` (e.g. Esc) still
    // releases the pointer — clear the compatibility hints then. Idempotent.
    window.addEventListener('pointerup', () => this.endConnectHints(), true);

    this.connection.addPipe((context) => {
      const c = context as {type: string; data: any};
      if (c.type === 'connectionpick') {
        // Right-click also fires connectionpick (the plugin doesn't filter
        // by button). Without this guard, a right-click on an output socket
        // would arm the suggestion menu, then drop with no connection,
        // opening it accidentally.
        if (this.lastPointerButton !== 0) {
          this.dragOutSource = null;
          this.dragInSource = null;
          return context;
        }
        const sock = c.data.socket as {nodeId: string; key: string; side: 'input' | 'output'};
        // Dim the canvas and light up only the sockets/nodes this pick can
        // legally connect to (compatible type, opposite side).
        this.beginConnectHints(sock.nodeId, sock.key, sock.side);
        const node = this.editor.getNode(sock.nodeId);
        if (sock.side === 'output') {
          this.dragInSource = null;
          const slot = node?.outputs[sock.key] as {socket: TypedSocket} | undefined;
          if (node && slot)
            this.dragOutSource = {nodeId: sock.nodeId, outputKey: sock.key, dgType: slot.socket.dgType};
        } else {
          // Input-side pick (a fresh drag out of an input, or the tail of an
          // existing connection) — arms the reverse drop-on-node shortcut
          // (exec-in included: its body drop connects the node's exec-out).
          this.dragOutSource = null;
          const slot = node?.inputs[sock.key] as {socket: TypedSocket} | undefined;
          this.dragInSource = (node && slot) ?
            {nodeId: sock.nodeId, inputKey: sock.key, dgType: slot.socket.dgType} : null;
        }
      }
      if (c.type === 'connectiondrop') {
        this.endConnectHints();
        const data = c.data as {created: boolean; socket: {nodeId: string} | null};
        const srcOut = this.dragOutSource;
        const srcIn = this.dragInSource;
        this.dragOutSource = null;
        this.dragInSource = null;
        if (!srcOut && !srcIn) return context;
        // A real connection happened, or the user dropped on a socket the
        // plugin handled (accepted, or rejected on type mismatch) — either way
        // they aimed at a specific socket, so don't second-guess them.
        if (data.created || data.socket) return context;
        // Dropped without hitting a socket: if it landed on a node body,
        // connect to its one obvious counterpart slot (no need to hit the tiny
        // pin). Empty-canvas drops open the suggestion menu for the matching
        // direction — consumers for an output drag, producers for an input drag.
        if (srcOut) void this.handleOutputDrop(srcOut, lastPointer.x, lastPointer.y);
        else if (srcIn) void this.handleInputDrop(srcIn, lastPointer.x, lastPointer.y);
      }
      return context;
    });
  }

  /** The rendered DOM element for a node (or null if not painted). */
  private nodeEl(nodeId: string): HTMLElement | null {
    return this.container.querySelector(`.ff-node[data-node-id="${CSS.escape(nodeId)}"]`);
  }

  /** Drop of an output-drag that missed every socket: connect to a node's sole
   *  compatible free input if it landed on one, else open the suggestion menu. */
  private async handleOutputDrop(
    src: {nodeId: string; outputKey: string; dgType: string}, x: number, y: number,
  ): Promise<void> {
    // `elementsFromPoint` (not `elementFromPoint`) so a transient overlay on top
    // of the node doesn't hide it.
    const stack = document.elementsFromPoint(x, y) as HTMLElement[];
    let nodeEl: HTMLElement | null = null;
    for (const el of stack) {
      const n = el.closest?.('.ff-node') as HTMLElement | null;
      if (n) {nodeEl = n; break;}
    }
    const targetNodeId = nodeEl?.dataset.nodeId;
    if (targetNodeId && targetNodeId !== src.nodeId) {
      // An order drag connects straight to the target's exec-in — every node
      // has one and it accepts many predecessors, so a body drop is
      // unambiguous (no aiming at the small square); duplicates are skipped.
      const key = isExecKey(src.outputKey) ?
        (this.hasConnection(src.nodeId, src.outputKey, targetNodeId, EXEC_IN_KEY) ? null : EXEC_IN_KEY) :
        this.soleCompatibleInput(src.nodeId, src.outputKey, targetNodeId);
      // Dropped on a node: connect to its one obvious input, or do nothing when
      // it has zero / several candidates (don't guess, don't pop the menu).
      if (key) await this.addConnectionByKeys(src.nodeId, src.outputKey, targetNodeId, key);
      return;
    }
    // No suggestion menu for order drags — nothing "produces" or "consumes"
    // an order signal; an empty-canvas drop is simply a no-op.
    if (isExecKey(src.outputKey)) return;
    // Dropped on the Outputs strip backdrop (between rows): publish the value
    // as a new flow output. Drops ON an existing row were handled by the
    // node-body branch above (they bind that row's free input).
    if (this.stripContains(x, y)) {
      await this.bindOutputToStrip(src);
      return;
    }
    await this.openSuggestionMenu(x, y, src);
  }

  /** Whether this exact connection already exists. */
  private hasConnection(source: string, sourceOutput: string, target: string, targetInput: string): boolean {
    return this.editor.getConnections().some((c) => c.source === source &&
      String(c.sourceOutput) === sourceOutput && c.target === target && String(c.targetInput) === targetInput);
  }

  /** Drop of an input-drag that missed every socket: if it landed on another
   *  node's body, connect from that node's one obvious output (a real output
   *  wins over a passthrough — see `soleCompatibleOutput`; ambiguity aborts).
   *  On empty canvas: open the reverse suggestion menu ("what produces this?"). */
  private async handleInputDrop(
    src: {nodeId: string; inputKey: string; dgType: string}, x: number, y: number,
  ): Promise<void> {
    const stack = document.elementsFromPoint(x, y) as HTMLElement[];
    let nodeEl: HTMLElement | null = null;
    for (const el of stack) {
      const n = el.closest?.('.ff-node') as HTMLElement | null;
      if (n) {
        nodeEl = n;
        break;
      }
    }
    const sourceNodeId = nodeEl?.dataset.nodeId;
    if (sourceNodeId && sourceNodeId !== src.nodeId) {
      // Order drag out of an exec-in: the dropped-on node becomes the
      // predecessor via its exec-out (mirror of the output-drop shortcut).
      const key = isExecKey(src.inputKey) ?
        (this.hasConnection(sourceNodeId, EXEC_OUT_KEY, src.nodeId, src.inputKey) ? null : EXEC_OUT_KEY) :
        this.soleCompatibleOutput(src.nodeId, src.inputKey, sourceNodeId);
      if (key) await this.addConnectionByKeys(sourceNodeId, key, src.nodeId, src.inputKey);
      return;
    }
    if (!sourceNodeId && !isExecKey(src.inputKey)) await this.openReverseSuggestionMenu(x, y, src);
  }

  /** The reverse suggestion menu: an input drag dropped on empty canvas offers
   *  every node type with a compatible output — or pass-through — (real
   *  producers first), creates the chosen one at the drop point, and wires its
   *  first compatible output (real over pass-through) into the dragged input. */
  private async openReverseSuggestionMenu(
    clientX: number, clientY: number,
    target: {nodeId: string; inputKey: string; dgType: string},
  ): Promise<void> {
    const {findNodeTypesProducingOutput, createNode} = await import('./node-factory');
    const nodes = this.editor.getNodes();
    const candidates = findNodeTypesProducingOutput(target.dgType, {
      sourcePackageName: this.editor.getNode(target.nodeId)?.dgPackageName,
      graphPackageNames: nodes.map((n) => n.dgPackageName).filter(Boolean),
      graphFuncNames: nodes.map((n) => n.dgFunc?.name ?? '').filter(Boolean),
    });
    if (candidates.length === 0) return;

    const choice = await this.promptSuggestion(clientX, clientY, candidates);
    if (!choice) return;

    const node = createNode(choice);
    if (!node) return;
    const {x, y} = this.screenToCanvas(clientX, clientY);
    await this.addNodeAt(node, x, y);

    // Auto-connect: the new node's first compatible output drives the dragged
    // input — a real output wins over a pass-through.
    const targetSocket = (this.editor.getNode(target.nodeId)?.inputs[target.inputKey] as
      {socket: TypedSocket} | undefined)?.socket;
    if (!targetSocket) return;
    let realKey: string | null = null;
    let ptKey: string | null = null;
    for (const [key, out] of Object.entries(node.outputs) as Array<[string, {socket: TypedSocket} | undefined]>) {
      if (!out || isExecKey(key)) continue;
      if (!out.socket.isCompatibleWith(targetSocket)) continue;
      if (key.endsWith('__pt')) {
        if (!ptKey) ptKey = key;
      } else {
        realKey = key;
        break;
      }
    }
    const outKey = realKey ?? ptKey;
    if (outKey) await this.addConnectionByKeys(node.id, outKey, target.nodeId, target.inputKey);
  }

  /** The output on `sourceNodeId` that can drive the dragged input: the sole
   *  compatible **real** output wins; only when no real output is compatible
   *  does the sole compatible **pass-through** qualify. Zero or several
   *  candidates in the winning group → null (don't guess — aim at a pin).
   *  Execution-order ports are ignored. Already-wired outputs stay eligible
   *  (an output legitimately feeds many consumers). Drives the reverse
   *  drop-on-node shortcut (drop an input drag anywhere on a producer node). */
  soleCompatibleOutput(dragNodeId: string, dragInputKey: string, sourceNodeId: string): string | null {
    const inSocket = (this.editor.getNode(dragNodeId)?.inputs[dragInputKey] as
      {socket: TypedSocket} | undefined)?.socket;
    const sourceNode = this.editor.getNode(sourceNodeId);
    if (!inSocket || !sourceNode) return null;
    const real: string[] = [];
    const passthrough: string[] = [];
    for (const [key, out] of Object.entries(sourceNode.outputs) as
      Array<[string, {socket: TypedSocket} | undefined]>) {
      if (!out || isExecKey(key)) continue;
      if (!out.socket.isCompatibleWith(inSocket)) continue;
      (key.endsWith('__pt') ? passthrough : real).push(key);
    }
    if (real.length === 1) return real[0];
    if (real.length === 0 && passthrough.length === 1) return passthrough[0];
    return null;
  }

  /** The single input on `targetNodeId` that the source output can drive AND is
   *  not already wired — or null if there are zero or several such inputs.
   *  Execution-order ports are ignored (they're aimed at deliberately). Drives
   *  the drop-on-node shortcut (drop anywhere on a node with one obvious input). */
  soleCompatibleInput(srcNodeId: string, srcOutputKey: string, targetNodeId: string): string | null {
    const srcSocket = (this.editor.getNode(srcNodeId)?.outputs[srcOutputKey] as
      {socket: TypedSocket} | undefined)?.socket;
    const targetNode = this.editor.getNode(targetNodeId);
    if (!srcSocket || !targetNode) return null;
    const candidates: string[] = [];
    for (const [key, input] of Object.entries(targetNode.inputs) as Array<[string, {socket: TypedSocket} | undefined]>) {
      if (!input || isExecKey(key)) continue;
      if (!srcSocket.isCompatibleWith(input.socket)) continue;
      if (this.isInputConnected(targetNodeId, key)) continue;
      candidates.push(key);
    }
    return candidates.length === 1 ? candidates[0] : null;
  }

  /** Dim the canvas and highlight only the opposite-side sockets (and their
   *  nodes) that a pick from `srcKey` could legally connect to. Cleared by
   *  `endConnectHints` on drop / pointer release. Symmetric: an output pick
   *  lights compatible inputs, an input pick (existing-connection tail) lights
   *  compatible outputs. */
  private beginConnectHints(srcNodeId: string, srcKey: string, srcSide: 'input' | 'output'): void {
    this.endConnectHints();
    const srcNode = this.editor.getNode(srcNodeId);
    const srcSlot = (srcSide === 'output' ? srcNode?.outputs[srcKey] : srcNode?.inputs[srcKey]) as
      {socket: TypedSocket} | undefined;
    const srcSocket = srcSlot?.socket;
    if (!srcSocket) return;
    this.container.classList.add('ff-connecting');
    // An order-port drag keeps every node's (normally hover-only) exec squares
    // visible for the whole gesture — the drag has visible targets and the
    // source square can't vanish mid-drag.
    if (isExecKey(srcKey))
      this.container.classList.add('ff-connecting-order');
    // A data-output drag can always land on the Outputs strip — light it up.
    if (srcSide === 'output' && !isExecKey(srcKey))
      this.outputStripEl?.classList.add('ff-strip-droptarget');
    this.nodeEl(srcNodeId)?.classList.add('ff-node-source');
    const targetSide: 'input' | 'output' = srcSide === 'output' ? 'input' : 'output';
    if (isExecKey(srcKey)) {
      // An order drag: every other node is a legal run-order neighbor — light
      // its opposite exec square (the wrapper carries the same compat class the
      // data-socket rows use, so the green glow rule applies as-is).
      const targetTid = targetSide === 'input' ? tid('exec-in') : tid('exec-out');
      for (const node of this.editor.getNodes()) {
        if (node.id === srcNodeId) continue;
        const nodeEl = this.nodeEl(node.id);
        if (!nodeEl) continue;
        nodeEl.querySelector(`[data-testid="${targetTid}"]`)?.classList.add('ff-socket-compat');
        nodeEl.classList.add('ff-node-compat');
      }
      return;
    }
    for (const node of this.editor.getNodes()) {
      if (node.id === srcNodeId) continue;
      const nodeEl = this.nodeEl(node.id);
      if (!nodeEl) continue;
      const slots = (targetSide === 'input' ? node.inputs : node.outputs) as
        Record<string, {socket: TypedSocket} | undefined>;
      let nodeCompat = false;
      for (const [key, slot] of Object.entries(slots)) {
        if (!slot || isExecKey(key)) continue;
        const ok = srcSide === 'output' ?
          srcSocket.isCompatibleWith(slot.socket) :
          slot.socket.isCompatibleWith(srcSocket);
        if (!ok) continue;
        nodeCompat = true;
        nodeEl.querySelector(`[data-testid="${tid('socket-' + targetSide, key)}"]`)
          ?.classList.add('ff-socket-compat');
      }
      if (nodeCompat) nodeEl.classList.add('ff-node-compat');
    }
  }

  /** Remove every connect-mode hint class. Idempotent. */
  private endConnectHints(): void {
    this.outputStripEl?.classList.remove('ff-strip-droptarget');
    if (!this.container.classList.contains('ff-connecting')) return;
    this.container.classList.remove('ff-connecting', 'ff-connecting-order');
    this.container.querySelectorAll('.ff-node-source, .ff-node-compat, .ff-socket-compat')
      .forEach((el) => el.classList.remove('ff-node-source', 'ff-node-compat', 'ff-socket-compat'));
  }

  private async openSuggestionMenu(
    clientX: number, clientY: number,
    source: {nodeId: string; outputKey: string; dgType: string},
  ): Promise<void> {
    const {findNodeTypesAcceptingInput, createNode} = await import('./node-factory');
    // Canvas context for the ranking heuristics: the science the drag came
    // from (source node's package), what's already on the canvas (packages →
    // domain fallback), and which functions the user already reached for.
    const nodes = this.editor.getNodes();
    const candidates = findNodeTypesAcceptingInput(source.dgType, {
      sourcePackageName: this.editor.getNode(source.nodeId)?.dgPackageName,
      graphPackageNames: nodes.map((n) => n.dgPackageName).filter(Boolean),
      graphFuncNames: nodes.map((n) => n.dgFunc?.name ?? '').filter(Boolean),
    });
    if (candidates.length === 0) return;

    const choice = await this.promptSuggestion(clientX, clientY, candidates);
    if (!choice) return;

    const node = createNode(choice);
    if (!node) return;
    const {x, y} = this.screenToCanvas(clientX, clientY);
    await this.addNodeAt(node, x, y);

    // Auto-connect: pick the first input on the new node whose socket the
    // dragged source can drive.
    const sourceSocket = (this.editor.getNode(source.nodeId)?.outputs[source.outputKey] as
      {socket: TypedSocket} | undefined)?.socket;
    if (!sourceSocket) return;
    let connectedKey: string | null = null;
    for (const [key, input] of Object.entries(node.inputs) as Array<[string, {socket: TypedSocket} | undefined]>) {
      if (!input) continue;
      if (sourceSocket.isCompatibleWith(input.socket)) {connectedKey = key; break;}
    }
    if (connectedKey)
      await this.addConnectionByKeys(source.nodeId, source.outputKey, node.id, connectedKey);
  }

  // ---------- hover docs ----------

  private hoverDocsEl: HTMLElement | null = null;
  private hoverDocsTimer: number | null = null;
  private hoverDocsNodeId: string | null = null;

  /** KNIME-style hover popup: a card next to the node with description,
   *  type, and (for func nodes) the input/output list. Shows after a short
   *  delay so flicking through the canvas doesn't spam. Stays open while
   *  the cursor is on the popup itself. */
  private installHoverDocs(): void {
    const popup = document.createElement('div');
    popup.className = 'ff-hover-docs';
    setTid(popup, 'hover-docs');
    popup.style.display = 'none';
    document.body.appendChild(popup);
    this.hoverDocsEl = popup;

    popup.addEventListener('mouseenter', () => {
      if (this.hoverDocsTimer !== null) {
        clearTimeout(this.hoverDocsTimer);
        this.hoverDocsTimer = null;
      }
    });
    popup.addEventListener('mouseleave', () => this.hideHoverDocs());

    // Hover trigger is the status circle in the title bar — narrow target so
    // the docs don't pop up every time the cursor passes over a node.
    this.container.addEventListener('mouseover', (ev) => {
      const target = ev.target as HTMLElement | null;
      if (!target?.classList?.contains('ff-node-status')) return;
      const nodeEl = target.closest('.ff-node') as HTMLElement | null;
      if (!nodeEl) return;
      const nodeId = nodeEl.dataset.nodeId;
      if (!nodeId || nodeId === this.hoverDocsNodeId) return;
      this.scheduleHoverDocs(nodeId, nodeEl);
    });

    this.container.addEventListener('mouseout', (ev) => {
      const fromEl = ev.target as HTMLElement | null;
      if (!fromEl?.classList?.contains('ff-node-status')) return;
      const toEl = ev.relatedTarget as HTMLElement | null;
      // Cursor moved onto the popup → keep it open.
      if (toEl && popup.contains(toEl)) return;
      this.hideHoverDocs();
    });
  }

  private scheduleHoverDocs(nodeId: string, nodeEl: HTMLElement): void {
    if (this.hoverDocsTimer !== null) clearTimeout(this.hoverDocsTimer);
    this.hoverDocsNodeId = nodeId;
    this.hoverDocsTimer = window.setTimeout(() => {
      this.showHoverDocs(nodeId, nodeEl);
    }, 500);
  }

  private showHoverDocs(nodeId: string, nodeEl: HTMLElement): void {
    const popup = this.hoverDocsEl;
    if (!popup) return;
    const node = this.editor.getNode(nodeId);
    if (!node) return;

    popup.innerHTML = '';

    const title = document.createElement('div');
    title.className = 'ff-hover-title';
    title.textContent = node.label;
    popup.appendChild(title);

    const badge = document.createElement('div');
    badge.className = 'ff-hover-badge';
    badge.textContent = node.dgNodeType ?? 'function';
    popup.appendChild(badge);

    const description = node.description ||
      (node.dgFunc?.description ?? '');
    if (description) {
      const desc = document.createElement('div');
      desc.className = 'ff-hover-desc';
      desc.textContent = description;
      popup.appendChild(desc);
    }

    if (node.dgFunc) {
      const fn = node.dgFunc;
      this.appendHoverParamSection(popup, 'Inputs', fn.inputs);
      this.appendHoverParamSection(popup, 'Outputs', fn.outputs);
    }

    // Position next to the node — right side preferred, left if no room.
    popup.style.display = 'block';
    const nodeRect = nodeEl.getBoundingClientRect();
    const popupRect = popup.getBoundingClientRect();
    const vw = window.innerWidth, vh = window.innerHeight;
    let left = nodeRect.right + 12;
    if (left + popupRect.width > vw - 8)
      left = nodeRect.left - popupRect.width - 12;
    if (left < 8) left = 8;
    let top = nodeRect.top;
    if (top + popupRect.height > vh - 8)
      top = vh - popupRect.height - 8;
    if (top < 8) top = 8;
    popup.style.left = `${left}px`;
    popup.style.top = `${top}px`;
  }

  private appendHoverParamSection(
    popup: HTMLElement, title: string, params: Array<DG.Property>,
  ): void {
    if (params.length === 0) return;
    const sec = document.createElement('div');
    sec.className = 'ff-hover-section';
    const h = document.createElement('div');
    h.className = 'ff-hover-section-title';
    h.textContent = title;
    sec.appendChild(h);
    for (const p of params) {
      const row = document.createElement('div');
      row.className = 'ff-hover-row';
      row.textContent = `${p.name}: ${String(p.propertyType)}`;
      if (p.description) row.title = p.description;
      sec.appendChild(row);
    }
    popup.appendChild(sec);
  }

  private hideHoverDocs(): void {
    if (this.hoverDocsTimer !== null) {
      clearTimeout(this.hoverDocsTimer);
      this.hoverDocsTimer = null;
    }
    this.hoverDocsNodeId = null;
    if (this.hoverDocsEl) this.hoverDocsEl.style.display = 'none';
  }

  // ---------- workflow annotations ----------

  /** Create a new annotation, mount its element in the transformed canvas
   *  layer (so it pans/zooms with the graph), and wire interactions. */
  addAnnotation(opts: Partial<AnnotationDoc> = {}): FlowAnnotation {
    const ann = new FlowAnnotation(opts);
    this.annotations.set(ann.id, ann);
    const content = this.area.area.content;
    content.add(ann.element);
    // Send to back: insert at the start of the holder's children list. With
    // `simpleNodesOrder` driving picked nodes to the end, this guarantees
    // every node and every existing connection paints on top of new
    // annotations (later DOM children paint over earlier ones in absolute-
    // positioned siblings, and z-index alone wasn't enough across stacking
    // contexts established by transforms).
    const firstChild = content.holder.firstChild;
    if (firstChild && firstChild !== ann.element)
      void content.reorder(ann.element, firstChild);
    this.installAnnotationInteractions(ann);
    this.callbacks.onGraphChanged?.();
    return ann;
  }

  removeAnnotation(id: string): void {
    const ann = this.annotations.get(id);
    if (!ann) return;
    this.annotations.delete(id);
    this.area.area.content.remove(ann.element);
    this.callbacks.onGraphChanged?.();
  }

  getAnnotations(): FlowAnnotation[] {
    return Array.from(this.annotations.values());
  }

  /** Everything an annotation drag carries along: nodes whose CENTER sits
   *  inside the annotation rect (strip-pinned output rows excluded), smaller
   *  annotations fully inside it, and the waypoints of connections linking two
   *  carried nodes (so routed wires travel with their endpoints). Computed at
   *  drag START — a stateless "capture": a node dragged out of the frame simply
   *  isn't inside at the next grab, so nothing has to be remembered or saved. */
  private annotationCargo(ann: FlowAnnotation): {
    nodes: Array<{id: string; start: {x: number; y: number}}>;
    annotations: Array<{ann: FlowAnnotation; start: {x: number; y: number}}>;
    groups: Array<{g: FlowGroup; start: {x: number; y: number}}>;
    waypoints: Array<{wp: {x: number; y: number}; start: {x: number; y: number}}>;
    connIds: string[];
  } {
    const x1 = ann.pos.x, y1 = ann.pos.y;
    const x2 = x1 + ann.size.w, y2 = y1 + ann.size.h;
    const nodes: Array<{id: string; start: {x: number; y: number}}> = [];
    const carried = new Set<string>();
    for (const node of this.editor.getNodes()) {
      // Members of a minimized group travel with their CARD (below), never
      // individually — moving them without the card would tear the group.
      if (node.dgNodeType === 'output' || this.minimizedGroupOf(node.id)) continue;
      const sz = this.measureNode(node.id);
      const cx = node.pos.x + sz.w / 2, cy = node.pos.y + sz.h / 2;
      if (cx >= x1 && cx <= x2 && cy >= y1 && cy <= y2) {
        nodes.push({id: node.id, start: {...node.pos}});
        carried.add(node.id);
      }
    }
    // Minimized group cards, carried card + hidden members together.
    const groups: Array<{g: FlowGroup; start: {x: number; y: number}}> = [];
    for (const g of this.groups.values()) {
      if (!g.minimized) continue;
      const w = g.element.offsetWidth || 180, h = g.element.offsetHeight || 40;
      const cx = g.pos.x + w / 2, cy = g.pos.y + h / 2;
      if (cx >= x1 && cx <= x2 && cy >= y1 && cy <= y2) {
        groups.push({g, start: {...g.pos}});
        for (const id of g.memberIds) {
          const n = this.editor.getNode(id);
          if (n) {
            nodes.push({id, start: {...n.pos}});
            carried.add(id);
          }
        }
      }
    }
    const annotations: Array<{ann: FlowAnnotation; start: {x: number; y: number}}> = [];
    for (const other of this.annotations.values()) {
      if (other === ann) continue;
      if (other.pos.x >= x1 && other.pos.y >= y1 &&
          other.pos.x + other.size.w <= x2 && other.pos.y + other.size.h <= y2)
        annotations.push({ann: other, start: {...other.pos}});
    }
    const waypoints: Array<{wp: {x: number; y: number}; start: {x: number; y: number}}> = [];
    const connIds: string[] = [];
    for (const c of this.editor.getConnections() as FlowConnection[]) {
      if (!c.waypoints || !carried.has(c.source) || !carried.has(c.target)) continue;
      connIds.push(c.id);
      for (const wp of c.waypoints) waypoints.push({wp, start: {...wp}});
    }
    return {nodes, annotations, groups, waypoints, connIds};
  }

  /** Drag-to-move on the body, drag-to-resize on the corner handle, custom
   *  contextmenu (Color · Delete), inline contenteditable for the title.
   *  Pointer deltas are divided by zoom so visible movement matches cursor. */
  private installAnnotationInteractions(ann: FlowAnnotation): void {
    const el = ann.element;
    const handle = ann.resizeHandle;
    const title = ann.titleEl;

    // ---- contextmenu: color palette + delete ----
    el.addEventListener('contextmenu', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      const menu = DG.Menu.popup();
      const colorMenu = menu.group('Color');
      for (const c of ANNOTATION_COLORS) {
        colorMenu.item(c.name, () => {
          ann.color = c.bg;
          ann.applyColor();
        });
      }
      colorMenu.endGroup()
        .separator()
        .item('Delete', () => this.removeAnnotation(ann.id))
        .show({causedBy: ev});
    });

    // ---- title editing: stopPropagation so AreaPlugin doesn't pan ----
    title.addEventListener('pointerdown', (ev) => ev.stopPropagation());

    // ---- drag-to-move (body, not title, not handle) ----
    el.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      const target = ev.target as HTMLElement | null;
      if (target && (target === title || title.contains(target))) return;
      if (target === handle) return;
      ev.preventDefault();
      ev.stopPropagation();
      const startPos = {...ann.pos};
      const startClient = {x: ev.clientX, y: ev.clientY};
      // The frame carries its contents: whatever is inside NOW moves with it.
      const cargo = this.annotationCargo(ann);
      // Synthetic pointers (tests) aren't active — capture is best-effort.
      try {el.setPointerCapture(ev.pointerId);} catch { /* no active pointer */ }
      const onMove = (e: PointerEvent): void => {
        const k = this.area.area.transform.k || 1;
        const dx = (e.clientX - startClient.x) / k;
        const dy = (e.clientY - startClient.y) / k;
        ann.pos.x = startPos.x + dx;
        ann.pos.y = startPos.y + dy;
        ann.applyPos();
        // Carried nodes are never "picked", so the snap interception in the
        // nodetranslate pipe skips them — group geometry stays intact.
        for (const n of cargo.nodes)
          void this.area.translate(n.id, {x: n.start.x + dx, y: n.start.y + dy});
        for (const a of cargo.annotations) {
          a.ann.pos.x = a.start.x + dx;
          a.ann.pos.y = a.start.y + dy;
          a.ann.applyPos();
        }
        for (const gr of cargo.groups) {
          gr.g.pos = {x: gr.start.x + dx, y: gr.start.y + dy};
          gr.g.applyCardPos();
          this.notifyGroupSockets(gr.g);
        }
        for (const w of cargo.waypoints) {
          w.wp.x = w.start.x + dx;
          w.wp.y = w.start.y + dy;
        }
        for (const id of cargo.connIds) void this.area.update('connection', id);
      };
      const onUp = (): void => {
        el.removeEventListener('pointermove', onMove);
        el.removeEventListener('pointerup', onUp);
        el.removeEventListener('pointercancel', onUp);
      };
      el.addEventListener('pointermove', onMove);
      el.addEventListener('pointerup', onUp);
      el.addEventListener('pointercancel', onUp);
    });

    // ---- drag-to-resize (bottom-right handle) ----
    handle.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      ev.preventDefault();
      ev.stopPropagation();
      const startSize = {...ann.size};
      const startClient = {x: ev.clientX, y: ev.clientY};
      try {handle.setPointerCapture(ev.pointerId);} catch { /* no active pointer */ }
      // Tiny floor only so the resize handle stays grabbable; the user
      // explicitly wanted no real lower bound.
      const minSize = 8;
      const onMove = (e: PointerEvent): void => {
        const k = this.area.area.transform.k || 1;
        ann.size.w = Math.max(minSize, startSize.w + (e.clientX - startClient.x) / k);
        ann.size.h = Math.max(minSize, startSize.h + (e.clientY - startClient.y) / k);
        ann.applySize();
      };
      const onUp = (): void => {
        handle.removeEventListener('pointermove', onMove);
        handle.removeEventListener('pointerup', onUp);
        handle.removeEventListener('pointercancel', onUp);
      };
      handle.addEventListener('pointermove', onMove);
      handle.addEventListener('pointerup', onUp);
      handle.addEventListener('pointercancel', onUp);
    });
  }

  // ---------- node groups ----------

  /** Create a group around the given nodes. Output nodes (strip-pinned) and
   *  nodes already in a group are filtered out; an empty remainder aborts.
   *  `opts` carries deserialized state (title, minimized, card pos). */
  createGroup(memberIds: string[], opts: Partial<GroupDoc> = {}): FlowGroup | null {
    const ids = memberIds.filter((id) => {
      const n = this.editor.getNode(id);
      return !!n && n.dgNodeType !== 'output' && !this.groupOf(id);
    });
    if (ids.length === 0) return null;
    const g = new FlowGroup({...opts, memberIds: ids});
    this.groups.set(g.id, g);
    const content = this.area.area.content;
    content.add(g.element);
    // Behind the nodes and wires, like annotations (see addAnnotation).
    const firstChild = content.holder.firstChild;
    if (firstChild && firstChild !== g.element)
      void content.reorder(g.element, firstChild);
    this.installGroupInteractions(g);
    if (g.minimized) {
      // Load path: the card lands at its saved pos with members hidden.
      this.setGroupHidden(g, true);
      g.applyMode();
      this.refreshGroupCard(g);
    } else {
      this.fitGroupFrame(g);
      // Member sizes settle when React mounts the views — refit shortly after
      // (fresh .ffjson loads create groups before the first paint).
      this.scheduleGroupRefit(g);
      setTimeout(() => {
        if (this.groups.has(g.id) && !g.minimized) this.fitGroupFrame(g);
      }, 300);
    }
    this.scheduleMinimapRedraw();
    this.callbacks.onGraphChanged?.();
    return g;
  }

  /** Group the current multi-selection (Ctrl+G / node context menu). */
  async createGroupFromSelection(): Promise<FlowGroup | null> {
    const ids = this.getSelectedNodeIds().filter((id) =>
      this.editor.getNode(id)?.dgNodeType !== 'output' && !this.groupOf(id));
    if (ids.length < 2) return null;
    const g = this.createGroup(ids);
    if (g) await this.unselectAllNodes();
    return g;
  }

  /** Dissolve a group — members stay on the canvas exactly where they are. */
  ungroup(id: string): void {
    const g = this.groups.get(id);
    if (!g) return;
    if (g.minimized) {
      g.minimized = false;
      this.setGroupHidden(g, false);
      for (const mid of g.memberIds) void this.area.update('node', mid);
    }
    this.groups.delete(id);
    this.area.area.content.remove(g.element);
    this.scheduleMinimapRedraw();
    this.callbacks.onGraphChanged?.();
  }

  /** Delete the group AND every member node (clearly-labeled menu item). */
  async deleteGroupWithNodes(id: string): Promise<void> {
    const g = this.groups.get(id);
    if (!g) return;
    const ids = Array.from(g.memberIds);
    this.ungroup(id);
    await this.removeNodes(ids);
  }

  getGroups(): FlowGroup[] {
    return Array.from(this.groups.values());
  }

  getGroupById(id: string): FlowGroup | undefined {
    return this.groups.get(id);
  }

  /** The group a node belongs to (a node is in at most one group). */
  groupOf(nodeId: string): FlowGroup | undefined {
    for (const g of this.groups.values())
      if (g.memberIds.has(nodeId)) return g;
    return undefined;
  }

  private minimizedGroupOf(nodeId: string): FlowGroup | undefined {
    const g = this.groupOf(nodeId);
    return g?.minimized ? g : undefined;
  }

  /** Pull one node out of its group (node context menu). The node stays put;
   *  a group left empty dissolves. */
  removeFromGroup(nodeId: string): void {
    const g = this.groupOf(nodeId);
    if (!g) return;
    g.memberIds.delete(nodeId);
    if (g.minimized) {
      this.setNodeHidden(nodeId, false);
      void this.area.update('node', nodeId);
    }
    if (g.memberIds.size === 0) {
      this.ungroup(g.id);
      return;
    }
    if (g.minimized) this.refreshGroupCard(g);
    else this.scheduleGroupRefit(g);
    this.callbacks.onGraphChanged?.();
  }

  async toggleGroupMinimized(id: string): Promise<void> {
    const g = this.groups.get(id);
    if (!g) return;
    if (g.minimized) await this.maximizeGroup(id);
    else await this.minimizeGroup(id);
  }

  /** Collapse the frame into a card at the frame's top-left: members and
   *  internal wires hide, boundary wires re-anchor to the card's edge dots. */
  async minimizeGroup(id: string): Promise<void> {
    const g = this.groups.get(id);
    if (!g || g.minimized) return;
    g.minimized = true;
    // Invisible nodes can't stay selected.
    for (const mid of g.memberIds) {
      const n = this.editor.getNode(mid) as {selected?: boolean} | undefined;
      if (n?.selected) await this.selectableApi.unselect(mid);
    }
    this.setGroupHidden(g, true);
    g.applyMode();
    this.refreshGroupCard(g);
    this.scheduleMinimapRedraw();
    this.callbacks.onSelectionChanged?.();
    this.callbacks.onGraphChanged?.();
  }

  /** Expand the card back into the frame. Members reappear where they are —
   *  card drags translated them live, so contents track the card's travels. */
  async maximizeGroup(id: string): Promise<void> {
    const g = this.groups.get(id);
    if (!g || !g.minimized) return;
    g.minimized = false;
    this.setGroupHidden(g, false);
    g.applyMode();
    this.fitGroupFrame(g);
    // Re-render members so the DOM socket watcher re-measures and re-emits
    // endpoints — releases the card anchors the boundary wires were glued to.
    for (const mid of g.memberIds) void this.area.update('node', mid);
    this.scheduleMinimapRedraw();
    this.callbacks.onGraphChanged?.();
  }

  /** Hide/show one node's canvas view (the wrapper element persists across
   *  React re-renders, so the class survives status updates). */
  private setNodeHidden(nodeId: string, hidden: boolean): void {
    const views = (this.area as unknown as {nodeViews: Map<string, {element: HTMLElement}>}).nodeViews;
    views?.get(nodeId)?.element.classList.toggle('ff-group-hidden', hidden);
  }

  /** The rendered wrapper element of a connection. */
  private connectionViewEl(connId: string): HTMLElement | null {
    const views = (this.area as unknown as {
      connectionViews?: Map<string, {element: HTMLElement}>;
    }).connectionViews;
    return views?.get(connId)?.element ??
      this.container.querySelector(`[data-connection-id="${CSS.escape(connId)}"]`);
  }

  /** Hide/show a group's member views and fully-internal connections. */
  private setGroupHidden(g: FlowGroup, hidden: boolean): void {
    for (const id of g.memberIds) this.setNodeHidden(id, hidden);
    for (const c of this.editor.getConnections()) {
      if (!g.memberIds.has(c.source) || !g.memberIds.has(c.target)) continue;
      this.connectionViewEl(c.id)?.classList.toggle('ff-group-hidden', hidden);
    }
  }

  /** Distinct member sockets with at least one connection crossing the group
   *  boundary, in stable connection order — each gets a dot row on the card. */
  private groupBoundarySockets(g: FlowGroup, side: 'input' | 'output'): Array<{nodeId: string; key: string}> {
    const result: Array<{nodeId: string; key: string}> = [];
    const seen = new Set<string>();
    for (const c of this.editor.getConnections()) {
      const memberEnd = side === 'input' ? c.target : c.source;
      const otherEnd = side === 'input' ? c.source : c.target;
      if (!g.memberIds.has(memberEnd) || g.memberIds.has(otherEnd)) continue;
      const key = side === 'input' ? String(c.targetInput) : String(c.sourceOutput);
      const dedupe = `${memberEnd}:${key}`;
      if (seen.has(dedupe)) continue;
      seen.add(dedupe);
      result.push({nodeId: memberEnd, key});
    }
    return result;
  }

  /** Canvas-coord wire anchor for a hidden member's socket: the matching dot
   *  on the minimized card's edge (inputs left, outputs right, stacked by the
   *  boundary row order — same math as `FlowGroup.renderDots`). */
  private groupSocketAnchor(
    g: FlowGroup, nodeId: string, side: 'input' | 'output', key: string,
  ): {x: number; y: number} {
    const rows = this.groupBoundarySockets(g, side);
    const i = Math.max(0, rows.findIndex((s) => s.nodeId === nodeId && s.key === key));
    const w = g.element.offsetWidth || 180;
    return {
      x: side === 'input' ? g.pos.x : g.pos.x + w,
      y: g.pos.y + GROUP_DOT_TOP + i * GROUP_DOT_STEP,
    };
  }

  /** `socketPositionWatcher.listen` for canvas nodes (see the constructor):
   *  forwards DOM-measured positions while the node is visible, and swallows
   *  them in favor of card-edge anchors while it hides in a minimized group
   *  (a display:none view measures at 0,0). */
  private listenGroupableSocket(
    nodeId: string, side: 'input' | 'output', key: string,
    onChange: (pos: {x: number; y: number}) => void,
    domWatcher: {listen(n: string, s: string, k: string, cb: (p: {x: number; y: number}) => void): () => void},
  ): () => void {
    let subs = this.groupSocketSubs.get(nodeId);
    if (!subs) {
      subs = new Set();
      this.groupSocketSubs.set(nodeId, subs);
    }
    const entry = {side, key, cb: onChange};
    subs.add(entry);
    const unsub = domWatcher.listen(nodeId, side, key, (pos) => {
      if (!this.minimizedGroupOf(nodeId)) onChange(pos);
    });
    const g = this.minimizedGroupOf(nodeId);
    if (g) onChange(this.groupSocketAnchor(g, nodeId, side, key));
    return () => {
      this.groupSocketSubs.get(nodeId)?.delete(entry);
      unsub();
    };
  }

  /** Push fresh card-edge anchors to every subscribed wire endpoint of a
   *  minimized group's members (card drags, dot-row changes). */
  private notifyGroupSockets(g: FlowGroup): void {
    if (!g.minimized) return;
    for (const id of g.memberIds) {
      const subs = this.groupSocketSubs.get(id);
      if (!subs) continue;
      for (const e of subs) e.cb(this.groupSocketAnchor(g, id, e.side, e.key));
    }
  }

  /** Re-sync a minimized card: boundary dots, wire anchors, internal-wire
   *  hiding (idempotent), aggregate status. */
  private refreshGroupCard(g: FlowGroup): void {
    if (!g.minimized) return;
    g.renderDots(
      this.groupBoundarySockets(g, 'input').length,
      this.groupBoundarySockets(g, 'output').length,
    );
    this.setGroupHidden(g, true);
    this.notifyGroupSockets(g);
    this.refreshGroupStatus(g);
  }

  /** Aggregate member run status → the card's title-bar dot. */
  private refreshGroupStatus(g: FlowGroup): void {
    const statuses = Array.from(g.memberIds)
      .map((id) => (this.editor.getNode(id) as {dgStatus?: string} | undefined)?.dgStatus ?? 'idle');
    let status = 'idle';
    if (statuses.some((s) => s === 'running')) status = 'running';
    else if (statuses.some((s) => s === 'errored')) status = 'errored';
    else if (statuses.some((s) => s === 'stale')) status = 'stale';
    else if (statuses.length > 0 && statuses.every((s) => s === 'completed')) status = 'completed';
    g.setStatus(status);
  }

  /** Size the expanded frame to the member bounding box (+ title bar and
   *  padding), keeping `g.pos` at the frame's top-left so minimizing collapses
   *  the group in place. */
  private fitGroupFrame(g: FlowGroup): void {
    if (g.minimized) return;
    let minX = Infinity; let minY = Infinity; let maxX = -Infinity; let maxY = -Infinity;
    for (const id of g.memberIds) {
      const n = this.editor.getNode(id);
      if (!n) continue;
      const sz = this.measureNode(id);
      minX = Math.min(minX, n.pos.x);
      minY = Math.min(minY, n.pos.y);
      maxX = Math.max(maxX, n.pos.x + sz.w);
      maxY = Math.max(maxY, n.pos.y + sz.h);
    }
    if (!Number.isFinite(minX)) return;
    g.pos = {x: minX - GROUP_PAD, y: minY - GROUP_PAD - GROUP_TITLE_H};
    g.frameSize = {
      w: maxX - minX + 2 * GROUP_PAD,
      h: maxY - minY + 2 * GROUP_PAD + GROUP_TITLE_H,
    };
    g.applyFrame();
  }

  /** Coalesce frame refits to one per animation frame (member drags emit many
   *  translates per pointermove). */
  private scheduleGroupRefit(g: FlowGroup): void {
    if (this.groupRefitScheduled.has(g.id)) return;
    this.groupRefitScheduled.add(g.id);
    requestAnimationFrame(() => {
      this.groupRefitScheduled.delete(g.id);
      if (this.groups.has(g.id) && !g.minimized) this.fitGroupFrame(g);
    });
  }

  /** A member node was removed from the editor — shrink (or dissolve) its group. */
  private handleGroupMemberRemoved(nodeId: string): void {
    const g = this.groupOf(nodeId);
    if (!g) return;
    g.memberIds.delete(nodeId);
    if (g.memberIds.size === 0) {
      this.ungroup(g.id);
      return;
    }
    if (g.minimized) this.refreshGroupCard(g);
    else this.scheduleGroupRefit(g);
  }

  /** Re-arrange just a group's members with the layered layout, anchored at
   *  the current bounding box's top-left (the rest of the canvas stays put). */
  async tidyGroup(id: string): Promise<void> {
    const g = this.groups.get(id);
    if (!g || g.minimized) return;
    const members = Array.from(g.memberIds)
      .map((mid) => this.editor.getNode(mid))
      .filter((n): n is FlowNode => !!n);
    if (members.length === 0) return;
    const inGroup = new Set(members.map((n) => n.id));
    const edges: LayoutEdge[] = [];
    for (const c of this.editor.getConnections()) {
      if (inGroup.has(c.source) && inGroup.has(c.target))
        edges.push({source: this.editor.getNode(c.source)!, target: this.editor.getNode(c.target)!});
    }
    const oldMinX = Math.min(...members.map((n) => n.pos.x));
    const oldMinY = Math.min(...members.map((n) => n.pos.y));
    layoutGraph(members, edges, computeLayers(members, edges));
    const newMinX = Math.min(...members.map((n) => n.pos.x));
    const newMinY = Math.min(...members.map((n) => n.pos.y));
    const dx = oldMinX - newMinX; const dy = oldMinY - newMinY;
    for (const n of members) await this.area.translate(n.id, {x: n.pos.x + dx, y: n.pos.y + dy});
    this.fitGroupFrame(g);
  }

  /** Drag (frame or card) carries the members; caret toggles minimized;
   *  double-click on the card maximizes; contextmenu offers group actions.
   *  Title/description are inline-editable (pointerdown stops propagation so
   *  the drag never starts on them). */
  private installGroupInteractions(g: FlowGroup): void {
    const el = g.element;

    el.addEventListener('contextmenu', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      (ev as MouseEvent & {_ffHandled?: boolean})._ffHandled = true;
      this.showGroupContextMenu(ev, g);
    });

    // While a title/description edit is in progress, keep pointer gestures out
    // of the drag/pan machinery; otherwise the title bar is the drag handle.
    for (const editable of [g.titleEl, g.descEl])
      editable.addEventListener('pointerdown', (ev) => {
        if (editable.isContentEditable) ev.stopPropagation();
      });
    g.descEl.addEventListener('blur', () => {
      delete g.descEl.dataset.editing;
      g.descEl.contentEditable = 'false';
      g.applyMode();
      this.callbacks.onGraphChanged?.();
    });
    g.titleEl.addEventListener('blur', () => this.callbacks.onGraphChanged?.());

    g.caretEl.addEventListener('pointerdown', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
    });
    g.caretEl.addEventListener('click', (ev) => {
      ev.stopPropagation();
      void this.toggleGroupMinimized(g.id);
    });

    // Double-press on the card maximizes. Detected on pointerdown ourselves:
    // the drag handler preventDefaults pointerdown, which suppresses the
    // browser's compatibility `dblclick` for that pointer.
    let lastDown = {t: 0, x: 0, y: 0};

    el.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      const target = ev.target as HTMLElement | null;
      if (target && (target === g.caretEl || target.closest('[contenteditable="true"]'))) return;
      ev.preventDefault();
      ev.stopPropagation();
      const now = Date.now();
      const isDouble = now - lastDown.t < 400 &&
        Math.abs(ev.clientX - lastDown.x) < 5 && Math.abs(ev.clientY - lastDown.y) < 5;
      lastDown = {t: now, x: ev.clientX, y: ev.clientY};
      if (isDouble) {
        // Double-press on the title renames, on the description edits it,
        // anywhere else on a minimized card maximizes.
        if (target && (target === g.titleEl || g.titleEl.contains(target))) {
          g.startTitleEdit();
          return;
        }
        if (target && (target === g.descEl || g.descEl.contains(target))) {
          g.startDescEdit();
          return;
        }
        if (g.minimized) {
          void this.toggleGroupMinimized(g.id);
          return;
        }
      }
      const startClient = {x: ev.clientX, y: ev.clientY};
      const startPos = {...g.pos};
      const members = Array.from(g.memberIds)
        .map((id) => this.editor.getNode(id))
        .filter((n): n is FlowNode => !!n)
        .map((n) => ({id: n.id, start: {...n.pos}}));
      // Waypoints of internal connections travel too.
      const waypoints: Array<{wp: {x: number; y: number}; start: {x: number; y: number}}> = [];
      const connIds: string[] = [];
      for (const c of this.editor.getConnections() as FlowConnection[]) {
        if (!c.waypoints || !g.memberIds.has(c.source) || !g.memberIds.has(c.target)) continue;
        connIds.push(c.id);
        for (const wp of c.waypoints) waypoints.push({wp, start: {...wp}});
      }
      try {el.setPointerCapture(ev.pointerId);} catch { /* synthetic pointer */ }
      const onMove = (e: PointerEvent): void => {
        const k = this.area.area.transform.k || 1;
        const dx = (e.clientX - startClient.x) / k;
        const dy = (e.clientY - startClient.y) / k;
        if (g.minimized) {
          g.pos = {x: startPos.x + dx, y: startPos.y + dy};
          g.applyCardPos();
        }
        // Members follow (hidden ones too — contents must track the card).
        // The expanded frame follows via the nodetranslated → refit path.
        for (const m of members)
          void this.area.translate(m.id, {x: m.start.x + dx, y: m.start.y + dy});
        for (const w of waypoints) {
          w.wp.x = w.start.x + dx;
          w.wp.y = w.start.y + dy;
        }
        if (!g.minimized)
          for (const id of connIds) void this.area.update('connection', id);
        if (g.minimized) this.notifyGroupSockets(g);
      };
      const onUp = (): void => {
        el.removeEventListener('pointermove', onMove);
        el.removeEventListener('pointerup', onUp);
        el.removeEventListener('pointercancel', onUp);
      };
      el.addEventListener('pointermove', onMove);
      el.addEventListener('pointerup', onUp);
      el.addEventListener('pointercancel', onUp);
    });
  }

  private showGroupContextMenu(event: MouseEvent, g: FlowGroup): void {
    const menu = DG.Menu.popup()
      .item(g.minimized ? 'Maximize' : 'Minimize', () => void this.toggleGroupMinimized(g.id));
    if (!g.minimized)
      menu.item('Tidy layout', () => void this.tidyGroup(g.id));
    menu.item('Edit description', () => g.startDescEdit());
    menu.separator()
      .item('Ungroup', () => this.ungroup(g.id))
      .item('Delete group and nodes', () => void this.deleteGroupWithNodes(g.id))
      .show({causedBy: event});
  }

  // ---------- connection waypoints ----------

  /** Right-click → "Add waypoint here" inserts a routing point on the
   *  connection. The React component chains classicConnectionPath segments
   *  through `start → waypoints → end`. */
  addWaypoint(conn: FlowConnection, at: {x: number; y: number}): void {
    if (!conn.waypoints) conn.waypoints = [];
    // Insert in the position closest to where the user clicked: pick the
    // segment whose midpoint is nearest, and insert after that endpoint.
    const insertIndex = this.bestWaypointInsertIndex(conn, at);
    conn.waypoints.splice(insertIndex, 0, {x: at.x, y: at.y});
    void this.area.update('connection', conn.id);
  }

  removeWaypoint(connId: string, index: number): void {
    const conn = this.editor.getConnections().find((c) => c.id === connId) as FlowConnection | undefined;
    if (!conn || !conn.waypoints) return;
    conn.waypoints.splice(index, 1);
    if (conn.waypoints.length === 0) delete conn.waypoints;
    void this.area.update('connection', conn.id);
  }

  /** Pick the segment whose endpoints sandwich the click best. Walk the
   *  start→…→end polyline, and for each segment compute distance from the
   *  click to the segment midpoint; insert the new waypoint after the
   *  start of the closest segment. Keeps the path geometry monotonic. */
  private bestWaypointInsertIndex(conn: FlowConnection, at: {x: number; y: number}): number {
    const waypoints = conn.waypoints ?? [];
    if (waypoints.length === 0) return 0;
    const sourceNode = this.editor.getNode(conn.source);
    const targetNode = this.editor.getNode(conn.target);
    if (!sourceNode || !targetNode) return waypoints.length;
    // We don't have exact socket positions here, so approximate with node
    // centers — good enough for picking which segment to split.
    const sSize = this.measureNode(sourceNode.id);
    const tSize = this.measureNode(targetNode.id);
    const start = {x: sourceNode.pos.x + sSize.w, y: sourceNode.pos.y + sSize.h / 2};
    const end = {x: targetNode.pos.x, y: targetNode.pos.y + tSize.h / 2};
    const points = [start, ...waypoints, end];
    let bestIdx = 0;
    let bestDist = Infinity;
    for (let i = 0; i < points.length - 1; i++) {
      const mx = (points[i].x + points[i + 1].x) / 2;
      const my = (points[i].y + points[i + 1].y) / 2;
      const d = (mx - at.x) ** 2 + (my - at.y) ** 2;
      if (d < bestDist) {bestDist = d; bestIdx = i;}
    }
    return bestIdx; // insert at this index (between point i and point i+1)
  }

  /** Wire pointerdown/contextmenu on every rendered waypoint circle. We use
   *  delegation on the canvas container since the React component re-creates
   *  circles whenever waypoints change. */
  private installWaypointInteractions(): void {
    this.container.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      const target = ev.target as Element | null;
      if (!target || !(target as HTMLElement).classList?.contains('ff-waypoint')) return;
      ev.preventDefault();
      ev.stopPropagation();
      const connId = (target as HTMLElement).dataset.connectionId!;
      const wpIdx = parseInt((target as HTMLElement).dataset.waypointIndex!, 10);
      const conn = this.editor.getConnections().find((c) => c.id === connId) as FlowConnection | undefined;
      if (!conn || !conn.waypoints) return;
      const wp = conn.waypoints[wpIdx];
      if (!wp) return;
      const startWP = {...wp};
      const startClient = {x: ev.clientX, y: ev.clientY};
      const onMove = (e: PointerEvent): void => {
        const k = this.area.area.transform.k || 1;
        wp.x = startWP.x + (e.clientX - startClient.x) / k;
        wp.y = startWP.y + (e.clientY - startClient.y) / k;
        void this.area.update('connection', conn.id);
      };
      const onUp = (): void => {
        window.removeEventListener('pointermove', onMove, true);
        window.removeEventListener('pointerup', onUp, true);
      };
      window.addEventListener('pointermove', onMove, true);
      window.addEventListener('pointerup', onUp, true);
    }, true);

    this.container.addEventListener('contextmenu', (ev) => {
      const target = ev.target as Element | null;
      if (!target || !(target as HTMLElement).classList?.contains('ff-waypoint')) return;
      ev.preventDefault();
      ev.stopPropagation();
      const connId = (target as HTMLElement).dataset.connectionId!;
      const wpIdx = parseInt((target as HTMLElement).dataset.waypointIndex!, 10);
      DG.Menu.popup()
        .item('Delete waypoint', () => this.removeWaypoint(connId, wpIdx))
        .show({causedBy: ev});
    }, true);
  }

  // ---------- shift+drag rectangle multi-select ----------

  /** Shift+drag on empty canvas → draw a marquee over the nodes. Mirrors the
   *  platform's area select (d4 `areaSelector` + `selectRows`): Shift+drag
   *  ADDS every node whose bounding box intersects the rectangle to the
   *  selection, Ctrl+Shift+drag REMOVES them (Ctrl read at mouse-up). The
   *  existing selection is never replaced — a plain empty-canvas click still
   *  clears it. */
  private installRectSelect(): void {
    let startClient: {x: number; y: number} | null = null;
    let rectEl: HTMLElement | null = null;

    const updateRect = (x: number, y: number): void => {
      if (!startClient || !rectEl) return;
      const containerRect = this.container.getBoundingClientRect();
      const left = Math.min(startClient.x, x) - containerRect.left;
      const top = Math.min(startClient.y, y) - containerRect.top;
      rectEl.style.left = `${left}px`;
      rectEl.style.top = `${top}px`;
      rectEl.style.width = `${Math.abs(startClient.x - x)}px`;
      rectEl.style.height = `${Math.abs(startClient.y - y)}px`;
    };

    const onMove = (e: PointerEvent): void => updateRect(e.clientX, e.clientY);
    const onUp = (e: PointerEvent): void => {
      // The AreaPlugin keeps an always-on window pointerup listener and treats
      // a release with few prior moves as an empty-canvas click → unselectAll.
      // This release ends OUR marquee (its pointerdown never reached the
      // area) — don't let it clear what the marquee just selected.
      e.stopImmediatePropagation();
      window.removeEventListener('pointermove', onMove, true);
      window.removeEventListener('pointerup', onUp, true);
      const sc = startClient;
      startClient = null;
      if (rectEl) {rectEl.remove(); rectEl = null;}
      if (!sc) return;
      void this.completeRectSelect(sc, {x: e.clientX, y: e.clientY}, e.ctrlKey || e.metaKey);
    };

    // Capture phase so we beat the AreaPlugin's pan handler — the user's
    // Shift+drag must produce a rectangle, not a canvas pan.
    this.container.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      if (!ev.shiftKey) return;
      const target = ev.target as HTMLElement | null;
      if (target?.closest('.ff-node, .ff-socket, .ff-minimap')) return;
      ev.preventDefault();
      ev.stopPropagation();

      startClient = {x: ev.clientX, y: ev.clientY};
      rectEl = document.createElement('div');
      rectEl.className = 'ff-rect-select';
      this.guideOverlay?.appendChild(rectEl);
      updateRect(ev.clientX, ev.clientY);

      window.addEventListener('pointermove', onMove, true);
      window.addEventListener('pointerup', onUp, true);
    }, true);
  }

  /** Hit-test every node's canvas-space bounding box against the marquee
   *  (also in canvas space). Adds the hits to the selection, or removes them
   *  when `remove` (Ctrl held at mouse-up). */
  private async completeRectSelect(
    startClient: {x: number; y: number},
    endClient: {x: number; y: number},
    remove: boolean,
  ): Promise<void> {
    const a = this.screenToCanvas(startClient.x, startClient.y);
    const b = this.screenToCanvas(endClient.x, endClient.y);
    const rx1 = Math.min(a.x, b.x), ry1 = Math.min(a.y, b.y);
    const rx2 = Math.max(a.x, b.x), ry2 = Math.max(a.y, b.y);
    if (rx2 - rx1 < 3 || ry2 - ry1 < 3) return; // ignore fat-finger clicks

    let touched = false;
    for (const node of this.editor.getNodes()) {
      if (this.minimizedGroupOf(node.id)) continue; // invisible — not selectable
      const sz = this.measureNode(node.id);
      const nx1 = node.pos.x, ny1 = node.pos.y;
      const nx2 = nx1 + sz.w, ny2 = ny1 + sz.h;
      if (!(rx1 < nx2 && nx1 < rx2 && ry1 < ny2 && ny1 < ry2)) continue;
      touched = true;
      if (remove) await this.selectableApi.unselect(node.id);
      else await this.selectableApi.select(node.id, true);
    }
    if (touched) {
      this.callbacks.onSelectionChanged?.();
      this.refreshChipSelection();
    }
  }

  /** Build a transient floating popup with a search input and a scrollable
   *  list of candidates. Resolves with the chosen typeName (or null on
   *  dismiss / Escape / click-outside). Keyboard nav: Up/Down/Enter. */
  private promptSuggestion(
    clientX: number, clientY: number,
    candidates: Array<{typeName: string; label: string; isBuiltin: boolean}>,
  ): Promise<string | null> {
    return new Promise((resolve) => {
      let resolved = false;
      const close = (val: string | null): void => {
        if (resolved) return;
        resolved = true;
        document.removeEventListener('mousedown', onDocMouseDown, true);
        document.removeEventListener('keydown', onKeyDown, true);
        popup.remove();
        resolve(val);
      };

      const popup = document.createElement('div');
      popup.className = 'ff-suggest-popup';
      setTid(popup, 'suggest-popup');
      popup.style.left = `${clientX}px`;
      popup.style.top = `${clientY}px`;

      const search = document.createElement('input');
      search.type = 'text';
      search.placeholder = 'Add node…';
      search.className = 'ff-suggest-search';
      setTid(search, 'suggest-search');
      popup.appendChild(search);

      const list = document.createElement('div');
      list.className = 'ff-suggest-list';
      setTid(list, 'suggest-list');
      popup.appendChild(list);

      let filtered = candidates;
      let activeIdx = 0;

      const renderList = (): void => {
        list.innerHTML = '';
        filtered.forEach((c, i) => {
          const row = document.createElement('div');
          row.className = 'ff-suggest-item' + (i === activeIdx ? ' ff-suggest-item-active' : '');
          row.textContent = c.label;
          row.dataset.testid = tid('suggest-item', c.typeName);
          row.dataset.nodeTypeName = c.typeName;
          if (c.isBuiltin) row.classList.add('ff-suggest-item-builtin');
          row.addEventListener('mouseenter', () => {
            activeIdx = i;
            for (const el of Array.from(list.children))
              (el as HTMLElement).classList.remove('ff-suggest-item-active');
            row.classList.add('ff-suggest-item-active');
          });
          row.addEventListener('mousedown', (ev) => {
            if (ev.button !== 0) return;
            ev.preventDefault();
            close(c.typeName);
          });
          list.appendChild(row);
        });
      };

      search.addEventListener('input', () => {
        const q = search.value.toLowerCase().trim();
        filtered = q === '' ? candidates :
          candidates.filter((c) => c.label.toLowerCase().includes(q) || c.typeName.toLowerCase().includes(q));
        activeIdx = 0;
        renderList();
      });

      const onKeyDown = (ev: KeyboardEvent): void => {
        if (ev.key === 'Escape') {ev.preventDefault(); close(null);}
        else if (ev.key === 'Enter') {
          ev.preventDefault();
          const c = filtered[activeIdx];
          if (c) close(c.typeName);
        } else if (ev.key === 'ArrowDown') {
          ev.preventDefault();
          activeIdx = Math.min(filtered.length - 1, activeIdx + 1);
          renderList();
          (list.children[activeIdx] as HTMLElement | undefined)?.scrollIntoView({block: 'nearest'});
        } else if (ev.key === 'ArrowUp') {
          ev.preventDefault();
          activeIdx = Math.max(0, activeIdx - 1);
          renderList();
          (list.children[activeIdx] as HTMLElement | undefined)?.scrollIntoView({block: 'nearest'});
        }
      };

      const onDocMouseDown = (ev: MouseEvent): void => {
        if (!popup.contains(ev.target as Node)) close(null);
      };

      document.body.appendChild(popup);
      // Clamp to viewport.
      const r = popup.getBoundingClientRect();
      const vw = window.innerWidth, vh = window.innerHeight;
      if (r.right > vw) popup.style.left = `${Math.max(8, vw - r.width - 8)}px`;
      if (r.bottom > vh) popup.style.top = `${Math.max(8, vh - r.height - 8)}px`;

      renderList();
      search.focus();

      document.addEventListener('mousedown', onDocMouseDown, true);
      document.addEventListener('keydown', onKeyDown, true);
    });
  }

  // ---------- connection styling ----------

  /** Look up the type color for a connection's source slot. */
  private connectionColor(conn: FlowConnection): string {
    const sourceNode = this.editor.getNode(conn.source);
    const sourceSlot = sourceNode?.outputs[String(conn.sourceOutput)] as
      {socket: TypedSocket} | undefined;
    return sourceSlot ? getSlotColor(sourceSlot.socket.dgType) : '#8892a0';
  }

  /** Stamp `_color` and a stable element id on a freshly-created connection
   *  so the React `<FlowConnectionComponent>` can paint it the right color
   *  without us touching the DOM. We also tag the connection wrapper element
   *  with `data-connection-id` so status-driven CSS in `funcflow.css` works. */
  private decorateConnection(conn: FlowConnection): void {
    (conn as FlowConnection & {_color?: string})._color = this.connectionColor(conn);
    // Tag the wrapper element after AreaPlugin mounts it. The 'rendered' signal
    // we listen to in wireEvents handles connection status; we only need to
    // stamp data-connection-id once when the area emits the render.
  }

  /** Stamp the wrapper with data attributes used by status-driven CSS. */
  private tagConnectionElement(data: {element: HTMLElement; payload: FlowConnection}): void {
    data.element.dataset.connectionId = data.payload.id;
    data.element.dataset.status = this.connectionStatuses.get(data.payload.id) ?? 'idle';
    // Execution-ordering edges render dashed/gray (CSS keys off data-order).
    data.element.dataset.order = isExecKey(String(data.payload.sourceOutput)) ? 'true' : 'false';
    // A connection internal to a minimized group that mounts late (load path,
    // wires created while collapsed) must come up hidden.
    const g = this.minimizedGroupOf(data.payload.source);
    if (g && g.memberIds.has(data.payload.target))
      data.element.classList.add('ff-group-hidden');
  }

  /** Set the status of a connection (drives the data-flow animation). */
  setConnectionStatus(connectionId: string, status: ConnectionStatus): void {
    this.connectionStatuses.set(connectionId, status);
    const el = this.container.querySelector<HTMLElement>(`[data-connection-id="${connectionId}"]`);
    if (el) el.dataset.status = status;
  }

  /** Reset all connections to the idle styling (used between runs). */
  resetConnectionStatuses(): void {
    for (const id of this.connectionStatuses.keys())
      this.setConnectionStatus(id, 'idle');
    this.connectionStatuses.clear();
  }

  /** Show (or clear, when `text` is null) a small data-count label at a
   *  connection's midpoint — the row/value count flowing through it after a run.
   *  Stuffed into the payload as `_count` and re-rendered, mirroring `_color`. */
  setConnectionLabel(connectionId: string, text: string | null): void {
    const conn = this.editor.getConnections().find((c) => c.id === connectionId) as
      (FlowConnection & {_count?: string}) | undefined;
    if (!conn) return;
    conn._count = text ?? undefined;
    void this.area.update('connection', connectionId);
  }

  /** Drop every wire's count label (between/after runs, or on edit). */
  clearConnectionLabels(): void {
    for (const c of this.editor.getConnections() as Array<FlowConnection & {_count?: string}>) {
      if (c._count !== undefined) {
        c._count = undefined;
        void this.area.update('connection', c.id);
      }
    }
  }

  // ---------- context menu + delete key ----------

  /** Right-click on a node or connection opens a `DG.Menu` popup with the
   *  appropriate actions. The platform menu handles positioning, dismissal,
   *  styling, and z-index for us. */
  private installContextMenu(): void {
    this.area.addPipe((context) => {
      if (context.type !== 'contextmenu') return context;
      // The DOM event bubbles from connection-wrapper → container, and the
      // area-plugin's emit-on-bubble fires on EACH listener — so a single
      // right-click on a connection produces two `contextmenu` signals: one
      // with `context: connection`, then one with `context: 'root'`. Without
      // this guard, the second (root) menu would clobber the first.
      const data = context.data as {
        event: MouseEvent & {_ffHandled?: boolean};
        context: 'root' | FlowNode | FlowConnection;
      };
      data.event.preventDefault();
      if (data.context === 'root') {
        if (data.event._ffHandled) return context; // already shown specific menu
        data.event._ffHandled = true;
        this.showRootContextMenu(data.event);
        return context;
      }
      data.event._ffHandled = true;
      // FlowNode has `inputs` (a Record of Inputs); FlowConnection has `source`/`target` ids.
      if ((data.context as FlowNode).inputs !== undefined)
        this.showNodeContextMenu(data.event, data.context as FlowNode);
      else
        this.showConnectionContextMenu(data.event, data.context as FlowConnection);
      return context;
    });
  }

  private showNodeContextMenu(event: MouseEvent, node: FlowNode): void {
    const menu = DG.Menu.popup();
    let hasRunItem = false;
    if (this.callbacks.onPreviewNode) {
      menu.item('Run up to here & preview', () => this.callbacks.onPreviewNode!(node.id));
      hasRunItem = true;
    }
    if (this.callbacks.onRerunNode && this.callbacks.canRerunNode?.(node.id)) {
      menu.item('Rerun this node only', () => this.callbacks.onRerunNode!(node.id));
      hasRunItem = true;
    }
    if (hasRunItem) menu.separator();
    const sel = this.getSelectedNodeIds();
    const inSelection = (node as {selected?: boolean}).selected === true;
    menu
      .item(node.collapsed ? 'Expand' : 'Collapse', () => void this.toggleCollapsed(node.id))
      .item('Duplicate', () => {
        // Right-clicking a node that is part of a multi-selection duplicates
        // the whole selection (with its internal connections).
        void this.duplicateNodes(inSelection && sel.length > 1 ? sel : [node.id]);
      });
    // Grouping: a multi-selection of ungrouped canvas nodes can become a
    // group; a grouped node offers the way out.
    const groupable = inSelection && sel.length > 1 && sel.filter((id) =>
      this.editor.getNode(id)?.dgNodeType !== 'output' && !this.groupOf(id)).length > 1;
    if (groupable)
      menu.item('Group selected', () => void this.createGroupFromSelection());
    if (this.groupOf(node.id))
      menu.item('Remove from group', () => this.removeFromGroup(node.id));
    menu
      .separator()
      .item('Delete', () => void this.removeNode(node.id))
      .show({causedBy: event});
  }

  private showConnectionContextMenu(event: MouseEvent, conn: FlowConnection): void {
    const canvasPt = this.screenToCanvas(event.clientX, event.clientY);
    DG.Menu.popup()
      .item('Add waypoint here', () => this.addWaypoint(conn, canvasPt))
      .separator()
      .item('Delete connection', () => void this.editor.removeConnection(conn.id))
      .show({causedBy: event});
  }

  /** Right-click on empty canvas (or canvas background between nodes). */
  private showRootContextMenu(event: MouseEvent): void {
    const canvasPt = this.screenToCanvas(event.clientX, event.clientY);
    DG.Menu.popup()
      .item('Add annotation here', () => void this.addAnnotation({
        pos: {x: canvasPt.x - 120, y: canvasPt.y - 70},
      }))
      .show({causedBy: event});
  }

  /** Capture-phase listeners implementing the selectRows click semantics.
   *  pointerdown snapshots the node under the cursor, its selection state,
   *  and the modifiers — `accumulating.active()` is called synchronously
   *  inside the selectable extension's `nodepicked` handler, so the snapshot
   *  must exist before that fires. pointerup applies the removals (Ctrl
   *  toggle-off, Ctrl+Shift remove, plain-click collapse) that rete's
   *  add-only `nodepicked` can't express. */
  private installPointerDownTracker(): void {
    this.pointerDownTracker = (ev: PointerEvent): void => {
      this.lastPointerButton = ev.button;
      const target = ev.target as HTMLElement | null;
      const nodeEl = target?.closest('.ff-node') as HTMLElement | null;
      const id = nodeEl?.dataset.nodeId ?? null;
      this.lastPointerDownNodeId = target?.closest('.ff-socket') ? null : id;
      const node = id ? this.editor.getNode(id) as {selected?: boolean} | undefined : undefined;
      this.lastPointerDownWasSelected = node?.selected === true;
      this.lastPointerDownModifier = ev.ctrlKey || ev.metaKey || ev.shiftKey;
      this.lastPointerDownPos = {x: ev.clientX, y: ev.clientY};
    };
    window.addEventListener('pointerdown', this.pointerDownTracker, true);

    // The removal half of the selectRows semantics. Rete's `nodepicked` only
    // ever ADDS, so Ctrl-toggle-off, Ctrl+Shift-remove, and the collapse of a
    // multi-selection on a plain click all run here, on a clean release (a
    // click, not a drag — a drag of a selected node must keep the group).
    this.pointerUpTracker = (ev: PointerEvent): void => {
      // Any release can end in a selection change (incl. the area extension's
      // own empty-canvas unselect-all, which bypasses our callbacks) — refresh
      // the chips' selected state after the handlers have run. In-place
      // attribute update, NEVER a chip rebuild: replacing the pressed element
      // mid-gesture would keep the browser from ever dispatching its `click`,
      // killing chip selection.
      if (this.container.contains(ev.target as Node)) this.refreshChipSelection();
      const id = this.lastPointerDownNodeId;
      if (ev.button !== 0 || id == null) return;
      if (Math.abs(ev.clientX - this.lastPointerDownPos.x) > 4 ||
          Math.abs(ev.clientY - this.lastPointerDownPos.y) > 4) return;
      const node = this.editor.getNode(id);
      if (!node) return;
      const ctrl = ev.ctrlKey || ev.metaKey;
      if (ctrl && (ev.shiftKey || this.lastPointerDownWasSelected)) {
        // Ctrl+Shift+click removes; Ctrl+click on a selected node toggles it off.
        void this.selectableApi.unselect(id);
        this.callbacks.onNodeDeselected?.(node);
        this.callbacks.onSelectionChanged?.();
      }
      else if (!ctrl && !ev.shiftKey && this.lastPointerDownWasSelected &&
               this.getSelectedNodeIds().length > 1) {
        void this.selectableApi.select(id, false); // plain click → exclusive
        this.callbacks.onSelectionChanged?.();
      }
    };
    window.addEventListener('pointerup', this.pointerUpTracker, true);

    // Block right- or middle-click pointer events on sockets and node bodies
    // from reaching the rete plugins — neither filters by button:
    // - the connection plugin's socket pointerdown starts a fake
    //   pseudoconnection drag that follows the cursor until pointerup;
    // - the node view ignores non-left pointerdowns, so the event bubbles to
    //   the area, whose selectable extension counts pointerdown→pointerup with
    //   <4 moves as an empty-canvas click and unselects ALL — clearing the very
    //   multi-selection the context menu's "Duplicate" is about to act on. The
    //   area's pointerup listener sits on `window`, so BOTH halves of the
    //   gesture must be swallowed here (capture phase, before the bubble path).
    // Pan, rect-select, and the `contextmenu` DOM event (a separate event — the
    // node menu still opens) are unaffected.
    const guardNonPrimary = (ev: PointerEvent): void => {
      if (ev.button === 0) return;
      const target = ev.target as HTMLElement | null;
      if (target?.closest('.ff-socket-row-input, .ff-socket-row-output, .ff-socket, .ff-node, ' +
          '.ff-group, .ff-annotation'))
        ev.stopPropagation();
    };
    this.container.addEventListener('pointerdown', guardNonPrimary, true);
    this.container.addEventListener('pointerup', guardNonPrimary, true);
  }

  private installKeyboardShortcuts(): void {
    this.keydownHandler = (e: KeyboardEvent) => {
      // Ignore key events while typing in form controls.
      const target = e.target as HTMLElement | null;
      const tag = target?.tagName ?? '';
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT' ||
          target?.isContentEditable) return;

      if (e.key === 'Delete' || e.key === 'Backspace') {
        const selectedIds = this.getSelectedNodeIds();
        if (selectedIds.length > 0) {
          e.preventDefault();
          void this.removeNodes(selectedIds);
        }
      }

      // Platform selection keys (scatterplot navigation.dart): Ctrl+A selects
      // every node, Ctrl+Shift+A deselects all.
      if ((e.key === 'a' || e.key === 'A') && (e.ctrlKey || e.metaKey) &&
          this.container.isConnected) {
        e.preventDefault();
        if (e.shiftKey)
          void this.unselectAllNodes();
        else
          // Members hidden inside minimized groups are invisible — selecting
          // them would arm Delete/copy on nodes the user can't see.
          for (const n of this.editor.getNodes())
            if (!this.minimizedGroupOf(n.id)) void this.selectableApi.select(n.id, true);
        this.callbacks.onSelectionChanged?.();
        this.refreshChipSelection();
      }

      // Ctrl+G groups the selection; Ctrl+Shift+G ungroups every group any
      // selected node belongs to.
      if ((e.key === 'g' || e.key === 'G') && (e.ctrlKey || e.metaKey) &&
          this.container.isConnected) {
        e.preventDefault();
        if (e.shiftKey) {
          const seen = new Set<string>();
          for (const id of this.getSelectedNodeIds()) {
            const g = this.groupOf(id);
            if (g) seen.add(g.id);
          }
          for (const gid of seen) this.ungroup(gid);
        }
        else
          void this.createGroupFromSelection();
      }

      if (e.key === 'Escape' && this.container.isConnected &&
          this.getSelectedNodeIds().length > 0)
        void this.unselectAllNodes();

      // Copy / paste nodes. A live text selection means the user is copying
      // text — leave the event to the browser.
      if ((e.key === 'c' || e.key === 'C') && (e.ctrlKey || e.metaKey) && !e.shiftKey &&
          this.container.isConnected && !document.getSelection()?.toString())
        this.copySelection();

      if ((e.key === 'v' || e.key === 'V') && (e.ctrlKey || e.metaKey) && !e.shiftKey &&
          this.container.isConnected && this.clipboard) {
        e.preventDefault();
        void this.pasteClipboard();
      }
    };
    window.addEventListener('keydown', this.keydownHandler);
  }

  getSelectedNodeIds(): string[] {
    return this.editor.getNodes().filter((n) => (n as {selected?: boolean}).selected === true).map((n) => n.id);
  }

  /** Double-click on empty canvas → zoom to fit all nodes. Clicks that hit a
   *  node, socket, control, or context menu are ignored. */
  private installDoubleClickToFit(): void {
    this.container.addEventListener('dblclick', (e) => {
      const target = e.target as HTMLElement | null;
      if (!target) return;
      if (target.closest(
        '.ff-node, .ff-socket, .ff-control, .ff-minimap, .ff-group, .d4-menu-item, input, textarea, select'))
        return;
      e.preventDefault();
      void this.zoomToFit();
    });
  }

  // ---------- undo / redo ----------

  async undo(): Promise<void> {await this.history.undo();}
  async redo(): Promise<void> {await this.history.redo();}

  // ---------- node operations ----------

  /** Programmatic selection — fires the same callback chain a click would.
   *  `accumulate=true` keeps existing selection (additive), false clears
   *  first. Used by auto-select-on-run-complete and the rectangle tool. */
  async selectNode(nodeId: string, accumulate = false): Promise<void> {
    const node = this.editor.getNode(nodeId);
    if (!node) return;
    await this.selectableApi.select(nodeId, accumulate);
    // Deliberately NOT deduped: programmatic selection (run-complete
    // auto-select) must re-fire even for an already-selected node — the host
    // re-shows the panel with fresh execution state. Pointer paths dedupe
    // at their source (nodepicked, the chip click handler).
    this.lastPickedId = nodeId;
    this.callbacks.onNodeSelected?.(node);
    this.callbacks.onSelectionChanged?.();
    this.refreshChipSelection();
  }

  async unselectAllNodes(): Promise<void> {
    await this.selector.unselectAll();
    this.callbacks.onSelectionChanged?.();
    this.refreshChipSelection();
  }

  /** Toggle a node's collapsed flag and re-render. */
  async toggleCollapsed(nodeId: string): Promise<void> {
    const node = this.editor.getNode(nodeId);
    if (!node) return;
    node.collapsed = !node.collapsed;
    await this.area.update('node', nodeId);
  }

  /** Remove a node and any connections touching it. */
  async removeNode(nodeId: string): Promise<void> {
    const conns = this.editor.getConnections().filter(
      (c) => c.source === nodeId || c.target === nodeId,
    );
    for (const c of conns) await this.editor.removeConnection(c.id);
    await this.editor.removeNode(nodeId);
  }

  /** Batch-remove multiple nodes and their connections. Safer than calling
   *  `removeNode` in a loop without awaiting — when two selected nodes share
   *  a connection, the second call would otherwise try to remove an already-
   *  gone connection. We dedupe connections first, then remove sequentially. */
  async removeNodes(ids: string[]): Promise<void> {
    if (ids.length === 0) return;
    const idSet = new Set(ids);
    const seen = new Set<string>();
    const connIds: string[] = [];
    for (const c of this.editor.getConnections()) {
      if (!seen.has(c.id) && (idSet.has(c.source) || idSet.has(c.target))) {
        seen.add(c.id);
        connIds.push(c.id);
      }
    }
    for (const cid of connIds) {
      try { await this.editor.removeConnection(cid); } catch { /* already gone */ }
    }
    for (const nid of ids) {
      try { await this.editor.removeNode(nid); } catch { /* already gone */ }
    }
  }

  // ---------- duplicate / copy / paste ----------

  private clipboard: GraphClip | null = null;
  /** How many times the current clipboard was pasted — each paste fans out
   *  further so repeated Ctrl+V doesn't stack copies on the same spot. */
  private pasteCount = 0;

  private snapshotNodes(ids: string[]): GraphClip {
    const idSet = new Set(ids);
    const nodes = this.editor.getNodes()
      .filter((n) => idSet.has(n.id) && n.dgTypeName != null)
      .map((n) => ({
        id: n.id, typeName: n.dgTypeName!, label: n.label, description: n.description,
        collapsed: n.collapsed, pos: {...n.pos},
        properties: JSON.parse(JSON.stringify(n.properties)),
        inputValues: JSON.parse(JSON.stringify(n.inputValues)),
      }));
    const kept = new Set(nodes.map((n) => n.id));
    const connections = this.editor.getConnections()
      .filter((c) => kept.has(c.source) && kept.has(c.target))
      .map((c) => ({source: c.source, sourceOutput: String(c.sourceOutput),
        target: c.target, targetInput: String(c.targetInput)}));
    return {nodes, connections};
  }

  /** Output paramNames and SetVar variableNames share one namespace (the
   *  validator flags duplicates as errors) — a copy must land with a unique
   *  name instead of instantly invalidating the graph. */
  private dedupeVariableName(fresh: FlowNode): void {
    const taken = new Set<string>();
    for (const n of this.editor.getNodes()) {
      const p = n.properties?.paramName;
      if (typeof p === 'string' && p !== '') taken.add(p);
      const v = n.inputValues?.variableName;
      if (typeof v === 'string' && v !== '') taken.add(v);
    }
    const bump = (name: string): string => {
      if (!taken.has(name)) return name;
      const m = name.match(/^(.*?)(\d+)$/);
      const base = m ? m[1] : name;
      let i = m ? parseInt(m[2], 10) + 1 : 2;
      while (taken.has(`${base}${i}`)) i++;
      return `${base}${i}`;
    };
    const paramName = fresh.properties?.paramName;
    if (typeof paramName === 'string' && paramName !== '')
      fresh.properties.paramName = bump(paramName);
    const varName = fresh.inputValues?.variableName;
    if (isSetVarNode(fresh) && typeof varName === 'string' && varName !== '')
      fresh.inputValues.variableName = bump(varName);
  }

  /** Instantiate a clip's nodes (offset from their recorded positions) and the
   *  connections among them; the copies become the new selection so they can
   *  be dragged as a group right away. */
  private async materializeClip(clip: GraphClip, offset: number): Promise<FlowNode[]> {
    // Lazy require to avoid a circular import.
    const {createNode} = await import('./node-factory');
    const idMap = new Map<string, FlowNode>();
    for (const snap of clip.nodes) {
      const fresh = createNode(snap.typeName);
      if (!fresh) continue;
      fresh.label = snap.label;
      fresh.description = snap.description;
      fresh.collapsed = snap.collapsed;
      fresh.properties = JSON.parse(JSON.stringify(snap.properties));
      fresh.inputValues = JSON.parse(JSON.stringify(snap.inputValues));
      this.dedupeVariableName(fresh);
      await this.editor.addNode(fresh);
      fresh.pos = {x: snap.pos.x + offset, y: snap.pos.y + offset};
      await this.area.translate(fresh.id, fresh.pos);
      idMap.set(snap.id, fresh);
    }
    for (const c of clip.connections) {
      const source = idMap.get(c.source);
      const target = idMap.get(c.target);
      if (source && target)
        await this.addConnectionByKeys(source.id, c.sourceOutput, target.id, c.targetInput);
    }
    const created = [...idMap.values()];
    if (created.length > 0) {
      await this.selector.unselectAll();
      for (const n of created) await this.selectableApi.select(n.id, true);
      this.lastPickedId = created[0].id;
      this.callbacks.onSelectionChanged?.();
      this.refreshChipSelection();
    }
    return created;
  }

  /** Duplicate the given nodes next to the originals. Connections whose both
   *  endpoints are duplicated are duplicated too, and the copies become the
   *  selection (so they're immediately movable as a group). */
  async duplicateNodes(ids: string[]): Promise<FlowNode[]> {
    return this.materializeClip(this.snapshotNodes(ids), 30);
  }

  /** Snapshot the selected nodes into the editor clipboard (Ctrl+C). Returns
   *  how many nodes were copied (0 = nothing selected, clipboard untouched). */
  copySelection(): number {
    const clip = this.snapshotNodes(this.getSelectedNodeIds());
    if (clip.nodes.length === 0) return 0;
    this.clipboard = clip;
    this.pasteCount = 0;
    return clip.nodes.length;
  }

  /** Materialize the clipboard (Ctrl+V). Each repeated paste offsets further. */
  async pasteClipboard(): Promise<FlowNode[]> {
    if (!this.clipboard) return [];
    this.pasteCount++;
    return this.materializeClip(this.clipboard, 30 * this.pasteCount);
  }

  // ---------- public API consumed by the rest of the package ----------

  async addNodeAtCenter(node: FlowNode): Promise<FlowNode> {
    await this.editor.addNode(node);
    const spot = this.findFreeSpot(this.viewportCenter(), node.id);
    node.pos = spot;
    await this.area.translate(node.id, spot);
    // Pan only when the chosen spot isn't already in view — chasing every new
    // node re-centers the viewport and walks the EARLIER nodes off-screen
    // (by the third added node the first one sat behind the toolbox).
    if (!this.isSpotVisible(spot))
      await this.panToNode(node.id);
    return node;
  }

  /** Whether a node placed at `spot` (canvas coords, assumed ~220×140) would
   *  be fully inside the current viewport. */
  private isSpotVisible(spot: {x: number; y: number}, w = 220, h = 140): boolean {
    const t = this.area.area.transform;
    const rect = this.area.container.getBoundingClientRect();
    const x1 = spot.x * t.k + t.x;
    const y1 = spot.y * t.k + t.y;
    const x2 = (spot.x + w) * t.k + t.x;
    const y2 = (spot.y + h) * t.k + t.y;
    return x1 >= 0 && y1 >= 0 && x2 <= rect.width && y2 <= rect.height;
  }

  /** A position at/near `start` where a node of roughly `w`×`h` overlaps no
   *  existing canvas node — a freshly added node must never bury the previous
   *  one (that hides the very sockets the user is about to wire). Prefers
   *  moving right (data flows left→right), then down, ring by ring. */
  private findFreeSpot(
    start: {x: number; y: number}, skipId: string, w = 220, h = 140,
  ): {x: number; y: number} {
    const others = this.editor.getNodes()
      .filter((n) => n.id !== skipId && n.dgNodeType !== 'output' && !this.minimizedGroupOf(n.id))
      .map((n) => {
        const sz = this.measureNode(n.id);
        return {x: n.pos.x, y: n.pos.y, w: sz.w, h: sz.h};
      });
    const margin = 30;
    const free = (x: number, y: number): boolean => others.every((o) =>
      x + w + margin <= o.x || x >= o.x + o.w + margin ||
      y + h + margin <= o.y || y >= o.y + o.h + margin);
    if (free(start.x, start.y)) return start;
    for (let ring = 1; ring <= 8; ring++) {
      const dx = ring * 260;
      const dy = ring * 180;
      const candidates = [
        {x: start.x + dx, y: start.y}, {x: start.x, y: start.y + dy},
        {x: start.x + dx, y: start.y + dy}, {x: start.x - dx, y: start.y},
        {x: start.x, y: start.y - dy}, {x: start.x - dx, y: start.y + dy},
        {x: start.x + dx, y: start.y - dy}, {x: start.x - dx, y: start.y - dy},
      ];
      for (const c of candidates) {
        if (free(c.x, c.y)) return c;
      }
    }
    return start;
  }

  async addNodeAt(node: FlowNode, x: number, y: number): Promise<FlowNode> {
    await this.editor.addNode(node);
    node.pos = {x, y};
    await this.area.translate(node.id, {x, y});
    return node;
  }

  /** Pan the canvas so the given node sits at the viewport center. Zoom is
   *  preserved — this is purely a translate. Wait one rAF so the node is
   *  rendered and measurable before we read its size. */
  async panToNode(id: string): Promise<void> {
    const node = this.editor.getNode(id);
    if (!node) return;
    await new Promise<void>((r) => requestAnimationFrame(() => r()));
    const sz = this.measureNode(id);
    const cx = node.pos.x + sz.w / 2;
    const cy = node.pos.y + sz.h / 2;
    const rect = this.area.container.getBoundingClientRect();
    const k = this.area.area.transform.k;
    await this.area.area.translate(rect.width / 2 - cx * k, rect.height / 2 - cy * k);
  }

  async addConnectionByKeys(
    sourceId: string, sourceKey: string,
    targetId: string, targetKey: string,
  ): Promise<boolean> {
    const source = this.editor.getNode(sourceId);
    const target = this.editor.getNode(targetId);
    if (!source || !target) return false;
    return this.editor.addConnection(new FlowConnection(
      source as never, sourceKey, target as never, targetKey,
    ));
  }

  async translate(nodeId: string, x: number, y: number): Promise<void> {
    const node = this.editor.getNode(nodeId);
    if (!node) return;
    node.pos = {x, y};
    await this.area.translate(nodeId, {x, y});
  }

  getNodes(): FlowNode[] {return this.editor.getNodes();}
  getConnections(): FlowConnection[] {return this.editor.getConnections();}
  getNodeById(id: string): FlowNode | undefined {return this.editor.getNode(id);}
  getNodeCount(): number {return this.editor.getNodes().length;}
  getConnectionCount(): number {return this.editor.getConnections().length;}

  getInputSource(nodeId: string, inputKey: string): {node: FlowNode; outputKey: string} | undefined {
    for (const c of this.editor.getConnections()) {
      if (c.target === nodeId && c.targetInput === inputKey) {
        const src = this.editor.getNode(c.source);
        if (src) return {node: src, outputKey: String(c.sourceOutput)};
      }
    }
    return undefined;
  }

  isInputConnected(nodeId: string, inputKey: string): boolean {
    return this.getInputSource(nodeId, inputKey) !== undefined;
  }

  /** Whether a node's slot has at least one connection touching it. */
  isSocketConnected(nodeId: string, side: 'input' | 'output', key: string): boolean {
    for (const c of this.editor.getConnections()) {
      if (side === 'input' && c.target === nodeId && String(c.targetInput) === key) return true;
      if (side === 'output' && c.source === nodeId && String(c.sourceOutput) === key) return true;
    }
    return false;
  }

  /** Find every connection going *out* of a given node output. */
  getOutgoingConnections(nodeId: string, outputKey?: string): FlowConnection[] {
    return this.editor.getConnections().filter((c) =>
      c.source === nodeId && (outputKey === undefined || String(c.sourceOutput) === outputKey),
    );
  }

  async updateNode(nodeId: string): Promise<void> {
    await this.area.update('node', nodeId);
    // An output node's visible form is its strip chip — re-render it too
    // (param renames, declared-type changes, run status).
    if (this.editor.getNode(nodeId)?.dgNodeType === 'output') this.scheduleStripSync(true);
    // Run-status changes arrive here (ExecutionVisualizer) — keep the node's
    // group card dot in sync.
    const g = this.groupOf(nodeId);
    if (g) this.refreshGroupStatus(g);
  }

  // ---------- viewport ----------

  private viewportCenter(): {x: number; y: number} {
    const t = this.area.area.transform;
    const rect = this.area.container.getBoundingClientRect();
    return {
      x: (rect.width / 2 - t.x) / t.k,
      y: (rect.height / 2 - t.y) / t.k,
    };
  }

  /** Convert a (clientX, clientY) point — e.g. from a pointer event — into
   *  canvas-space coords matching `node.pos`. Used by drop handlers and the
   *  drag-out suggestion menu. */
  screenToCanvas(clientX: number, clientY: number): {x: number; y: number} {
    const r = this.container.getBoundingClientRect();
    const t = this.area.area.transform;
    return {
      x: (clientX - r.left - t.x) / t.k,
      y: (clientY - r.top - t.y) / t.k,
    };
  }

  zoomIn(): void {void this.area.area.zoom(this.area.area.transform.k * 1.2);}
  zoomOut(): void {void this.area.area.zoom(this.area.area.transform.k * 0.8);}

  async zoomToFit(): Promise<void> {
    // Fit the graph proper — strip-pinned output rows follow the viewport, so
    // including them would chase a moving target (and they're always visible
    // anyway). Fall back to everything when only output rows exist.
    const nodes = this.editor.getNodes();
    const inner = nodes.filter((n) => n.dgNodeType !== 'output');
    if (!Array.from(this.groups.values()).some((g) => g.minimized)) {
      await AreaExtensions.zoomAt(this.area, inner.length > 0 ? inner : nodes, {scale: 0.9});
      return;
    }
    // Manual fit: zoomAt measures node views, and hidden members measure at
    // (0,0) — compute the bounds from visible nodes + minimized cards instead.
    let minX = Infinity; let minY = Infinity; let maxX = -Infinity; let maxY = -Infinity;
    const grow = (x: number, y: number, w: number, h: number): void => {
      minX = Math.min(minX, x); minY = Math.min(minY, y);
      maxX = Math.max(maxX, x + w); maxY = Math.max(maxY, y + h);
    };
    for (const n of inner) {
      if (this.minimizedGroupOf(n.id)) continue;
      const sz = this.measureNode(n.id);
      grow(n.pos.x, n.pos.y, sz.w, sz.h);
    }
    for (const g of this.groups.values()) {
      if (!g.minimized) continue;
      grow(g.pos.x, g.pos.y, g.element.offsetWidth || 180, g.element.offsetHeight || 40);
    }
    if (!Number.isFinite(minX)) return;
    const rect = this.canvasEl.getBoundingClientRect();
    const gw = Math.max(1, maxX - minX); const gh = Math.max(1, maxY - minY);
    const k = Math.min(2.5, Math.max(0.2, Math.min(rect.width / gw, rect.height / gh) * 0.9));
    await this.area.area.zoom(k);
    await this.area.area.translate(
      rect.width / 2 - (minX + gw / 2) * k, rect.height / 2 - (minY + gh / 2) * k);
  }

  /** Re-arrange the whole graph with the layered/banded layout used by the
   *  creation-script importer (`rete/graph-layout.ts`): layers from the
   *  connection structure (every edge points right), one band per disjoint
   *  path, producer paths above the paths that consume them. Repositions every
   *  node and zooms to fit. */
  async autoLayout(): Promise<void> {
    // A whole-canvas re-layout repositions every node — expand minimized
    // groups first so their (hidden) members don't land under a stale card.
    for (const g of this.groups.values())
      if (g.minimized) await this.maximizeGroup(g.id);
    // Output rows are strip-pinned — lay out the graph proper without them
    // (their translates would be canceled by the pin guard anyway).
    const nodes = this.editor.getNodes().filter((n) => n.dgNodeType !== 'output');
    if (nodes.length === 0) return;
    const byId = new Map(nodes.map((n) => [n.id, n]));
    const edges: LayoutEdge[] = [];
    for (const c of this.editor.getConnections()) {
      // Order edges ARE included. computeLayers runs fresh here (unlike the
      // importer's stale incremental layer map), so an order edge is just
      // another forward dependency — it places the "after" node further right,
      // letting explicit run-order shape the layout left-to-right.
      const source = byId.get(c.source);
      const target = byId.get(c.target);
      if (source && target) edges.push({source, target});
    }
    layoutGraph(nodes, edges, computeLayers(nodes, edges));
    for (const node of nodes) await this.area.translate(node.id, {x: node.pos.x, y: node.pos.y});
    await this.zoomToFit();
  }

  // ---------- lifecycle ----------

  async clear(): Promise<void> {
    this.connectionStatuses.clear();
    for (const ann of Array.from(this.annotations.values()))
      this.removeAnnotation(ann.id);
    for (const g of Array.from(this.groups.values()))
      this.ungroup(g.id);
    await this.editor.clear();
  }

  destroy(): void {
    if (this.keydownHandler) window.removeEventListener('keydown', this.keydownHandler);
    if (this.pointerDownTracker)
      window.removeEventListener('pointerdown', this.pointerDownTracker, true);
    if (this.pointerUpTracker)
      window.removeEventListener('pointerup', this.pointerUpTracker, true);
    if (this.hoverDocsEl) this.hoverDocsEl.remove();
    if (this.minimapEl) this.minimapEl.remove();
    // Null these so a still-pending rAF redraw after teardown is a no-op.
    this.minimapEl = null;
    this.minimapSvg = null;
    this.stripResizeObserver?.disconnect();
    this.stripResizeObserver = null;
    if (this.outputStripEl) this.outputStripEl.remove();
    this.outputStripEl = null;
    this.stripChipsEl = null;
    this.chipSocketSubs.clear();
    this.groupSocketSubs.clear();
    this.groups.clear();
    this.area.destroy();
    this.canvasWrap.remove();
  }
}
