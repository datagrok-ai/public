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

import {FlowConnection, FlowEditorBridge, FlowNode, FlowScheme, isExecKey} from './scheme';
import {TypedSocket} from './sockets';
import {FlowConnectionComponent, FlowNodeComponent, FlowSocketComponent} from './node-component';
import {getSlotColor} from '../types/type-map';
import {tid, setTid} from '../utils/test-ids';
import {FlowAnnotation, AnnotationDoc, ANNOTATION_COLORS} from './annotation';
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

export class FlowEditor {
  readonly editor = new NodeEditor<FlowScheme>();
  readonly area: AreaPlugin<FlowScheme>;
  readonly connection = new ConnectionPlugin<FlowScheme>();
  readonly render: ReactPlugin<FlowScheme, ReactArea2D<FlowScheme>>;
  readonly history = new HistoryPlugin<FlowScheme, HistoryActions<FlowScheme>>();
  readonly container: HTMLElement;

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
  /** Underlying ctrl-held accumulator from rete; we wrap it to also accumulate
   *  when the click landed on an already-selected node. */
  private ctrlAccumulating = AreaExtensions.accumulateOnCtrl();
  private accumulating = {
    active: (): boolean => {
      if (this.ctrlAccumulating.active()) return true;
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
  /** Per-connection status (for execution coloring). */
  private connectionStatuses = new Map<string, ConnectionStatus>();

  /** Workflow annotations — colored frames behind the graph (KNIME pattern).
   *  Owned by the editor (not by Rete), persisted alongside the graph. */
  private annotations = new Map<string, FlowAnnotation>();

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
    this.area = new AreaPlugin<FlowScheme>(container);
    this.render = new ReactPlugin<FlowScheme, ReactArea2D<FlowScheme>>({createRoot});

    this.render.addPreset(ReactPresets.classic.setup({
      // No arrow markers anymore — direction comes from the dash-flow CSS
      // animation. The line just needs to land on the dot edge, so a small
      // symmetric offset that puts both endpoints just inside the socket dot
      // (radius 4.5 px) keeps everything visually attached.
      socketPositionWatcher: getDOMSocketPosition({
        offset: (pos, _id, side) => ({
          x: pos.x + (side === 'output' ? 2 : -2),
          y: pos.y,
        }),
      }) as never,
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
      if (context.type === 'connectioncreated')
        this.maybeAutoTypeValueOutput(context.data);
      if (context.type === 'connectionremoved')
        this.connectionStatuses.delete(context.data.id);
      if (context.type === 'connectioncreated' || context.type === 'connectionremoved')
        this.refreshCollapsedEndpoints(context.data);
      if (
        context.type === 'nodecreated' || context.type === 'noderemoved' ||
        context.type === 'connectioncreated' || context.type === 'connectionremoved' ||
        context.type === 'cleared'
      ) {
        this.callbacks.onGraphChanged?.();
        this.callbacks.onGraphEdited?.(this.classifyEdit(context));
        this.scheduleMinimapRedraw();
      }
      return context;
    });

    let lastPickedId: string | null = null;
    this.area.addPipe((context) => {
      if (context.type === 'nodepicked') {
        const node = this.editor.getNode(context.data.id);
        if (node) {
          if (lastPickedId && lastPickedId !== node.id) {
            const prev = this.editor.getNode(lastPickedId);
            if (prev) this.callbacks.onNodeDeselected?.(prev);
          }
          lastPickedId = node.id;
          this.callbacks.onNodeSelected?.(node);
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
      if (other.id === draggedId) continue;
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
    this.container.appendChild(el);
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
      const sz = this.measureNode(node.id);
      const color = (node as unknown as {color?: string}).color ?? '#90a4ae';
      boxes.push({x: node.pos.x, y: node.pos.y, w: sz.w, h: sz.h, color});
      minX = Math.min(minX, node.pos.x);
      minY = Math.min(minY, node.pos.y);
      maxX = Math.max(maxX, node.pos.x + sz.w);
      maxY = Math.max(maxY, node.pos.y + sz.h);
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
          // existing connection) — arms the reverse drop-on-node shortcut.
          this.dragOutSource = null;
          const slot = node?.inputs[sock.key] as {socket: TypedSocket} | undefined;
          this.dragInSource = (node && slot && !isExecKey(sock.key)) ?
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
      const key = this.soleCompatibleInput(src.nodeId, src.outputKey, targetNodeId);
      // Dropped on a node: connect to its one obvious input, or do nothing when
      // it has zero / several candidates (don't guess, don't pop the menu).
      if (key) await this.addConnectionByKeys(src.nodeId, src.outputKey, targetNodeId, key);
      return;
    }
    await this.openSuggestionMenu(x, y, src);
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
      const key = this.soleCompatibleOutput(src.nodeId, src.inputKey, sourceNodeId);
      if (key) await this.addConnectionByKeys(sourceNodeId, key, src.nodeId, src.inputKey);
      return;
    }
    if (!sourceNodeId) await this.openReverseSuggestionMenu(x, y, src);
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
    this.nodeEl(srcNodeId)?.classList.add('ff-node-source');
    const targetSide: 'input' | 'output' = srcSide === 'output' ? 'input' : 'output';
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
    if (!this.container.classList.contains('ff-connecting')) return;
    this.container.classList.remove('ff-connecting');
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
      el.setPointerCapture(ev.pointerId);
      const onMove = (e: PointerEvent): void => {
        const k = this.area.area.transform.k || 1;
        ann.pos.x = startPos.x + (e.clientX - startClient.x) / k;
        ann.pos.y = startPos.y + (e.clientY - startClient.y) / k;
        ann.applyPos();
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
      handle.setPointerCapture(ev.pointerId);
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

  // ---------- ctrl+drag rectangle multi-select ----------

  /** Ctrl+drag (or Cmd+drag on macOS) on empty canvas → draw a marquee and
   *  select every node whose bounding box intersects it. Hold Shift to add
   *  to the existing selection instead of replacing it. */
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
      window.removeEventListener('pointermove', onMove, true);
      window.removeEventListener('pointerup', onUp, true);
      const sc = startClient;
      startClient = null;
      if (rectEl) {rectEl.remove(); rectEl = null;}
      if (!sc) return;
      void this.completeRectSelect(sc, {x: e.clientX, y: e.clientY}, e.shiftKey);
    };

    // Capture phase so we beat the AreaPlugin's pan handler — the user's
    // Ctrl+drag must produce a rectangle, not a canvas pan.
    this.container.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      if (!ev.ctrlKey && !ev.metaKey) return;
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
   *  (also in canvas space). Replaces the existing selection unless `additive`. */
  private async completeRectSelect(
    startClient: {x: number; y: number},
    endClient: {x: number; y: number},
    additive: boolean,
  ): Promise<void> {
    const a = this.screenToCanvas(startClient.x, startClient.y);
    const b = this.screenToCanvas(endClient.x, endClient.y);
    const rx1 = Math.min(a.x, b.x), ry1 = Math.min(a.y, b.y);
    const rx2 = Math.max(a.x, b.x), ry2 = Math.max(a.y, b.y);
    if (rx2 - rx1 < 3 || ry2 - ry1 < 3) return; // ignore fat-finger clicks

    if (!additive) await this.selector.unselectAll();

    for (const node of this.editor.getNodes()) {
      const sz = this.measureNode(node.id);
      const nx1 = node.pos.x, ny1 = node.pos.y;
      const nx2 = nx1 + sz.w, ny2 = ny1 + sz.h;
      const intersects = rx1 < nx2 && nx1 < rx2 && ry1 < ny2 && ny1 < ry2;
      if (intersects) await this.selectableApi.select(node.id, true);
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
    menu
      .item(node.collapsed ? 'Expand' : 'Collapse', () => void this.toggleCollapsed(node.id))
      .item('Duplicate', () => void this.duplicateNode(node.id))
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

  /** Capture-phase listener that records which canvas node was under the
   *  cursor at pointerdown. Read by `accumulating.active()` to preserve a
   *  multi-selection when the click lands on an already-selected node — the
   *  selectable extension calls `accumulating.active()` synchronously inside
   *  its `nodepicked` handler, so the value must be set before that fires. */
  private installPointerDownTracker(): void {
    this.pointerDownTracker = (ev: PointerEvent): void => {
      this.lastPointerButton = ev.button;
      const target = ev.target as HTMLElement | null;
      const nodeEl = target?.closest('.ff-node') as HTMLElement | null;
      this.lastPointerDownNodeId = nodeEl?.dataset.nodeId ?? null;
    };
    window.addEventListener('pointerdown', this.pointerDownTracker, true);

    // Block right- or middle-click pointerdown from reaching the connection
    // plugin's socket handler. The plugin doesn't filter by button; without
    // this guard, a right-click on an output socket starts a fake
    // pseudoconnection drag that follows the cursor until pointerup.
    // We only stop the event when it landed inside a socket — pan, rect-
    // select, and contextmenu emission elsewhere are unaffected.
    this.container.addEventListener('pointerdown', (ev) => {
      if (ev.button === 0) return;
      const target = ev.target as HTMLElement | null;
      if (target?.closest('.ff-socket-row-input, .ff-socket-row-output, .ff-socket'))
        ev.stopPropagation();
    }, true);
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
    };
    window.addEventListener('keydown', this.keydownHandler);
  }

  private getSelectedNodeIds(): string[] {
    return this.editor.getNodes().filter((n) => (n as {selected?: boolean}).selected === true).map((n) => n.id);
  }

  /** Double-click on empty canvas → zoom to fit all nodes. Clicks that hit a
   *  node, socket, control, or context menu are ignored. */
  private installDoubleClickToFit(): void {
    this.container.addEventListener('dblclick', (e) => {
      const target = e.target as HTMLElement | null;
      if (!target) return;
      if (target.closest('.ff-node, .ff-socket, .ff-control, .ff-minimap, .d4-menu-item, input, textarea, select'))
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
    this.callbacks.onNodeSelected?.(node);
  }

  async unselectAllNodes(): Promise<void> {
    await this.selector.unselectAll();
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

  /** Duplicate a node next to the original (no connections copied). */
  async duplicateNode(nodeId: string): Promise<FlowNode | null> {
    const original = this.editor.getNode(nodeId);
    if (!original || !original.dgTypeName) return null;
    // Lazy require to avoid a circular import.
    const {createNode} = await import('./node-factory');
    const fresh = createNode(original.dgTypeName);
    if (!fresh) return null;
    fresh.label = original.label;
    fresh.properties = JSON.parse(JSON.stringify(original.properties));
    fresh.inputValues = JSON.parse(JSON.stringify(original.inputValues));
    await this.editor.addNode(fresh);
    const offset = 30;
    fresh.pos = {x: original.pos.x + offset, y: original.pos.y + offset};
    await this.area.translate(fresh.id, fresh.pos);
    return fresh;
  }

  // ---------- public API consumed by the rest of the package ----------

  async addNodeAtCenter(node: FlowNode): Promise<FlowNode> {
    await this.editor.addNode(node);
    const center = this.viewportCenter();
    node.pos = center;
    await this.area.translate(node.id, center);
    await this.panToNode(node.id);
    return node;
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
    await AreaExtensions.zoomAt(this.area, this.editor.getNodes(), {scale: 0.9});
  }

  /** Re-arrange the whole graph with the layered/banded layout used by the
   *  creation-script importer (`rete/graph-layout.ts`): layers from the
   *  connection structure (every edge points right), one band per disjoint
   *  path, producer paths above the paths that consume them. Repositions every
   *  node and zooms to fit. */
  async autoLayout(): Promise<void> {
    const nodes = this.editor.getNodes();
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
    await this.editor.clear();
  }

  destroy(): void {
    if (this.keydownHandler) window.removeEventListener('keydown', this.keydownHandler);
    if (this.pointerDownTracker)
      window.removeEventListener('pointerdown', this.pointerDownTracker, true);
    if (this.hoverDocsEl) this.hoverDocsEl.remove();
    if (this.minimapEl) this.minimapEl.remove();
    // Null these so a still-pending rAF redraw after teardown is a no-op.
    this.minimapEl = null;
    this.minimapSvg = null;
    this.area.destroy();
  }
}
