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

import {FlowConnection, FlowNode, FlowScheme} from './scheme';
import {TypedSocket} from './sockets';
import {FlowConnectionComponent, FlowNodeComponent, FlowSocketComponent} from './node-component';
import {getSlotColor} from '../types/type-map';

export interface FlowEditorCallbacks {
  onNodeSelected?: (node: FlowNode) => void;
  onNodeDeselected?: (node: FlowNode) => void;
  onGraphChanged?: () => void;
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
  private accumulating = AreaExtensions.accumulateOnCtrl();
  private callbacks: FlowEditorCallbacks;
  private keydownHandler: ((e: KeyboardEvent) => void) | null = null;
  /** Per-connection status (for execution coloring). */
  private connectionStatuses = new Map<string, ConnectionStatus>();

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

    AreaExtensions.selectableNodes(this.area, this.selector, {
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

    this.wireEvents();
    this.applyInitialBackground(container);
    this.installContextMenu();
    this.installKeyboardShortcuts();
    this.installDoubleClickToFit();

    // Expose a narrow callback surface to React components that need to talk
    // back into the editor.
    (window as unknown as {__ff_editor: {
      toggleCollapsed(id: string): void;
      isSocketConnected(nodeId: string, side: 'input' | 'output', key: string): boolean;
    }}).__ff_editor = {
      toggleCollapsed: (id) => void this.toggleCollapsed(id),
      isSocketConnected: (nodeId, side, key) => this.isSocketConnected(nodeId, side, key),
    };
  }

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
    const targetNode = this.editor.getNode(connection.target) as FlowNode | undefined;
    if (!targetNode || targetNode.label !== 'Value Output') return;
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
      // Stamp `_color` on every new connection BEFORE the area-plugin emits
      // 'render', so the React Connection component picks up the right color
      // on its very first render.
      if (context.type === 'connectioncreate')
        this.decorateConnection(context.data);
      if (context.type === 'connectioncreated')
        this.maybeAutoTypeValueOutput(context.data);
      if (context.type === 'connectionremoved')
        this.connectionStatuses.delete(context.data.id);
      if (
        context.type === 'nodecreated' || context.type === 'noderemoved' ||
        context.type === 'connectioncreated' || context.type === 'connectionremoved' ||
        context.type === 'cleared'
      )
        this.callbacks.onGraphChanged?.();
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
      if (context.type === 'nodetranslated') {
        const node = this.editor.getNode(context.data.id);
        if (node) node.pos = {...context.data.position};
      }
      // Tag each rendered connection wrapper with its id + status for the
      // CSS-driven execution-state animations (`[data-status="active"]` etc).
      if (context.type === 'rendered' && (context.data as {type?: string}).type === 'connection')
        this.tagConnectionElement(context.data as {element: HTMLElement; payload: FlowConnection});
      return context;
    });
  }

  private applyInitialBackground(container: HTMLElement): void {
    container.style.background = `#ebedf2 url('${this.makeDotGridDataUrl()}')`;
  }

  private makeDotGridDataUrl(): string {
    const size = 20;
    const half = size / 2;
    const c = document.createElement('canvas');
    c.width = size;
    c.height = size;
    const ctx = c.getContext('2d')!;
    ctx.fillStyle = '#ebedf2';
    ctx.fillRect(0, 0, size, size);
    ctx.fillStyle = '#8b8b8b';
    ctx.fillRect(half - 1, half - 1, 2, 2);
    return c.toDataURL();
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

  // ---------- context menu + delete key ----------

  /** Right-click on a node or connection opens a `DG.Menu` popup with the
   *  appropriate actions. The platform menu handles positioning, dismissal,
   *  styling, and z-index for us. */
  private installContextMenu(): void {
    this.area.addPipe((context) => {
      if (context.type === 'contextmenu') {
        const data = context.data as {event: MouseEvent; context: 'root' | FlowNode | FlowConnection};
        data.event.preventDefault();
        if (data.context === 'root') return context;
        // FlowNode has `inputs` (a Record of Inputs); FlowConnection has `source`/`target` ids.
        if ((data.context as FlowNode).inputs !== undefined)
          this.showNodeContextMenu(data.event, data.context as FlowNode);
        else
          this.showConnectionContextMenu(data.event, data.context as FlowConnection);
      }
      return context;
    });
  }

  private showNodeContextMenu(event: MouseEvent, node: FlowNode): void {
    DG.Menu.popup()
      .item(node.collapsed ? 'Expand' : 'Collapse', () => void this.toggleCollapsed(node.id))
      .item('Duplicate', () => void this.duplicateNode(node.id))
      .separator()
      .item('Delete', () => void this.removeNode(node.id))
      .show({causedBy: event});
  }

  private showConnectionContextMenu(event: MouseEvent, conn: FlowConnection): void {
    DG.Menu.popup()
      .item('Delete connection', () => void this.editor.removeConnection(conn.id))
      .show({causedBy: event});
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
          for (const id of selectedIds) void this.removeNode(id);
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
      if (target.closest('.ff-node, .ff-socket, .ff-control, .d4-menu-item, input, textarea, select'))
        return;
      e.preventDefault();
      void this.zoomToFit();
    });
  }

  // ---------- undo / redo ----------

  async undo(): Promise<void> {await this.history.undo();}
  async redo(): Promise<void> {await this.history.redo();}

  // ---------- node operations ----------

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
    return node;
  }

  async addNodeAt(node: FlowNode, x: number, y: number): Promise<FlowNode> {
    await this.editor.addNode(node);
    node.pos = {x, y};
    await this.area.translate(node.id, {x, y});
    return node;
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

  zoomIn(): void {void this.area.area.zoom(this.area.area.transform.k * 1.2);}
  zoomOut(): void {void this.area.area.zoom(this.area.area.transform.k * 0.8);}

  async zoomToFit(): Promise<void> {
    await AreaExtensions.zoomAt(this.area, this.editor.getNodes(), {scale: 0.9});
  }

  // ---------- lifecycle ----------

  async clear(): Promise<void> {
    this.connectionStatuses.clear();
    await this.editor.clear();
  }

  destroy(): void {
    if (this.keydownHandler) window.removeEventListener('keydown', this.keydownHandler);
    this.area.destroy();
    delete (window as unknown as {__ff_editor?: unknown}).__ff_editor;
  }
}
