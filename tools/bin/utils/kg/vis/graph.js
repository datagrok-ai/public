// The canvas: cosmos.gl over the compact arrays, the selection, the label overlay, and the events the
// page listens to. Everything cosmos knows is a dense index; `view.back` turns it into a graph index.

const SPACE = 8192;
/** Below this many visible nodes the forces are those of a small diagram, not of a hundred-thousand-node cloud. */
const SMALL = 300;
/** Fit never zooms past this: labels stay readable and points stay points. */
const MAX_ZOOM = 5;
/** Mean on-screen motion below which the layout counts as settled and is paused. */
const SETTLED_PX_PER_S = 0.6;
/** While the layout moves, labels are re-placed this often, not every frame. */
const LABEL_MS = 250;
/** How much faster than d3-zoom's default a wheel notch zooms. */
const WHEEL_GAIN = 4;
const AMPLIFIED = Symbol('kgv-wheel');

/** cosmos forces are per pair: what spreads 100k nodes into a readable cloud flings 7 nodes off the screen. */
function forcesFor(n) {
  if (n < SMALL)
    return {simulationRepulsion: 0.4, simulationGravity: 0.5, simulationCenter: 0.5, simulationLinkSpring: 1, simulationLinkDistance: 12,
      simulationFriction: 0.8, simulationDecay: 500, simulationCluster: 0};
  if (n < 5000)
    return {simulationRepulsion: 2, simulationGravity: 0.2, simulationCenter: 0.2, simulationLinkSpring: 0.8, simulationLinkDistance: 20,
      simulationFriction: 0.85, simulationDecay: 900, simulationCluster: 0.02};
  return {simulationRepulsion: 8, simulationGravity: 0.08, simulationCenter: 0, simulationLinkSpring: 0.6, simulationLinkDistance: 30,
    simulationFriction: 0.85, simulationDecay: 1500, simulationCluster: 0.02};
}

export class GraphView {
  constructor(container, labels, handlers) {
    this.container = container;
    this.labelsLayer = labels;
    this.handlers = handlers;
    this.view = null;
    this.selected = null;
    this.hovered = null;
    this.labelsOn = true;
    this.labelOrder = [];
    this.userMoved = false;
    const {Graph} = window.Cosmos;
    this.graph = new Graph(container, {
      spaceSize: SPACE,
      backgroundColor: getComputedStyle(document.documentElement).getPropertyValue('--canvas').trim() || '#f4f4f2',
      pointSizeScale: 1.6,
      scalePointsOnZoom: false,
      linkWidth: 0.5,
      linkWidthScale: 1,
      curvedLinks: false,
      linkArrows: false,
      renderHoveredPointRing: true,
      hoveredPointCursor: 'pointer',
      pointGreyoutOpacity: 0.12,
      linkGreyoutOpacity: 0.04,
      simulationRepulsion: 8,
      simulationGravity: 0.08,
      simulationCenter: 0,
      simulationLinkSpring: 0.6,
      simulationLinkDistance: 30,
      simulationFriction: 0.85,
      simulationDecay: 3000,
      simulationCluster: 0.02,
      fitViewOnInit: true,
      fitViewDelay: 400,
      fitViewPadding: 0.08,
      onClick: (i) => this.onClick(i),
      onPointMouseOver: (i, pos, event) => this.onHover(i, event),
      onPointMouseOut: () => this.onHover(undefined),
      onZoom: () => this.updateLabels(),
      onSimulationEnd: () => { if (!this.userMoved) this.fit(); },
    });
    container.addEventListener('dblclick', () => {
      if (this.hovered !== null && this.hovered !== undefined) handlers.onExpand(this.view.back[this.hovered]);
    });
    for (const type of ['mousedown', 'wheel', 'touchstart']) container.addEventListener(type, () => { this.userMoved = true; }, {passive: true});
    // cosmos zooms through d3-zoom on its canvas: a wheel over a label never reaches it, and a notch is timid.
    // Every wheel over the graph area is replayed on the canvas with the delta scaled up.
    for (const layer of [container, labels]) layer.addEventListener('wheel', (event) => {
      if (event[AMPLIFIED]) return;
      const canvas = container.querySelector('canvas');
      if (!canvas) return;
      event.preventDefault();
      event.stopPropagation();
      const replay = new WheelEvent('wheel', {bubbles: true, cancelable: true, clientX: event.clientX, clientY: event.clientY,
        deltaX: event.deltaX, deltaY: event.deltaY * WHEEL_GAIN, deltaMode: event.deltaMode, ctrlKey: event.ctrlKey, shiftKey: event.shiftKey});
      replay[AMPLIFIED] = true;
      canvas.dispatchEvent(replay);
      this.userMoved = true;
    }, {capture: true, passive: false});
    this.shown = new Set();
    let lastLabels = 0;
    this.tick = (now) => {
      if (this.graph.isSimulationRunning && now - lastLabels > LABEL_MS) {
        lastLabels = now;
        this.updateLabels();
      }
      requestAnimationFrame(this.tick);
    };
    requestAnimationFrame(this.tick);
  }

  async ready() {
    await this.graph.ready;
  }

  /** New arrays; positions of nodes still visible come from `positions`, which the caller refreshed from us.
   * `resimulate` restarts the layout and follows it with the camera, as after a preset switch. */
  setData(view, resimulate = false) {
    const first = !this.view;
    this.view = view;
    this.selected = null;
    this.hovered = null;
    // setConfig replaces the whole config and a config change after the data setters drops the pending positions
    this.graph.setConfigPartial({highlightedPointIndices: [], outlinedPointIndices: [], ...forcesFor(view.nodeCount)});
    this.graph.setPointPositions(view.positions);
    this.graph.setPointColors(view.colors);
    this.graph.setPointSizes(view.sizes);
    this.graph.setLinks(view.links);
    this.graph.setLinkColors(view.linkColors);
    this.graph.setPointClusters(view.clusters);
    this.graph.render();
    if (first || resimulate) {
      this.userMoved = false;
      this.graph.start(1);
      this.followSettling();
    }
    else if (this.graph.isSimulationRunning) this.graph.start(0.6);
    this.shown = new Set();
    this.pickLabelNodes();
    this.updateLabels();
  }

  /** Current positions written back into the graph-wide array so a re-filter keeps the layout. */
  savePositions(positions) {
    if (!this.view) return;
    const current = this.graph.getPointPositions();
    if (!current || current.length < this.view.nodeCount * 2) return;
    for (let k = 0; k < this.view.nodeCount; k++) {
      const i = this.view.back[k];
      positions[i * 2] = current[k * 2];
      positions[i * 2 + 1] = current[k * 2 + 1];
    }
  }

  onClick(dense) {
    if (dense === undefined || dense === null) {
      this.select(null);
      this.handlers.onSelect(null);
      return;
    }
    this.select(this.view.back[dense]);
    this.handlers.onSelect(this.view.back[dense]);
  }

  onHover(dense, event) {
    this.hovered = dense ?? null;
    this.handlers.onHover(dense === undefined || dense === null ? null : this.view.back[dense], event);
    this.updateLabels();
  }

  /** Selects a graph index (null clears): the node is outlined, its neighbours highlighted, the rest dimmed. */
  select(index) {
    this.selected = index;
    if (index === null || !this.view || this.view.map[index] < 0) {
      this.graph.setConfigPartial({highlightedPointIndices: [], outlinedPointIndices: []});
      this.updateLabels();
      return;
    }
    const dense = this.view.map[index];
    const around = this.graph.getNeighboringPointIndices(dense) ?? [];
    this.graph.setConfigPartial({highlightedPointIndices: [dense, ...around], outlinedPointIndices: [dense]});
    this.updateLabels();
  }

  /** Lights a set of graph indices (a Cypher answer); null restores the plain view. */
  highlight(indices) {
    if (!this.view) return;
    if (!indices) {
      this.graph.setConfigPartial({highlightedPointIndices: [], outlinedPointIndices: []});
      return;
    }
    const dense = [];
    for (const i of indices)
      if (this.view.map[i] >= 0) dense.push(this.view.map[i]);
    this.graph.setConfigPartial({highlightedPointIndices: dense, outlinedPointIndices: dense.length <= 500 ? dense : []});
    if (dense.length) this.graph.fitViewByPointIndices(dense, 600, 0.2);
  }

  isVisible(index) {
    return !!this.view && this.view.map[index] >= 0;
  }

  zoomTo(index) {
    if (!this.isVisible(index)) return;
    this.graph.zoomToPointByIndex(this.view.map[index], 600, 3, true);
  }

  fit(duration = 600) {
    this.graph.fitView(duration, 0.08);
    // a handful of nodes would otherwise fill the screen at a zoom where every point is a coin
    setTimeout(() => { if ((this.graph.getZoomLevel() ?? 1) > MAX_ZOOM) this.graph.setZoomLevel(MAX_ZOOM, 200); }, duration + 20);
  }

  /** Watches the layout settle: the camera fits twice early on, and the simulation is paused as soon as the
   * nodes stop moving on screen (or after a hard cap), because cosmos's decay counts frames and would keep
   * everything creeping for half a minute. */
  followSettling() {
    clearInterval(this.follow);
    const small = this.view.nodeCount < SMALL;
    let ticks = 0;
    let last = null;
    this.follow = setInterval(() => {
      ticks++;
      const now = this.graph.getPointPositions();
      const zoom = this.graph.getZoomLevel() ?? 1;
      let moved = 0;
      if (last && now.length === last.length) {
        for (let i = 0; i < now.length; i += 2) moved += Math.hypot(now[i] - last[i], now[i + 1] - last[i + 1]);
        moved = moved / (now.length / 2) * zoom * 2;
      }
      last = now;
      const settled = ticks >= 3 && moved < SETTLED_PX_PER_S;
      if (this.userMoved || !this.graph.isSimulationRunning || settled || ticks > (small ? 8 : 24)) {
        clearInterval(this.follow);
        if (this.graph.isSimulationRunning) this.graph.pause();
        if (!this.userMoved) this.fit();
        this.handlers.onSimulation?.();
        return;
      }
      if (!small && (ticks === 1 || ticks === 4)) this.fit(450);
    }, 500);
  }

  toggleSimulation() {
    if (this.graph.isSimulationRunning) this.graph.pause();
    else this.graph.start(0.6);
    return this.graph.isSimulationRunning;
  }

  get simulating() {
    return this.graph.isSimulationRunning;
  }

  /** The visible nodes with the highest degree get labels; how many depends on the zoom. */
  pickLabelNodes() {
    const {back, nodeCount} = this.view;
    const degree = this.handlers.degree;
    const order = Array.from({length: nodeCount}, (_, k) => k).sort((a, b) => degree[back[b]] - degree[back[a]]);
    this.labelOrder = order.slice(0, 400);
  }

  updateLabels() {
    if (!this.view || !this.labelsOn) {
      this.labelsLayer.replaceChildren();
      return;
    }
    const zoom = this.graph.getZoomLevel() ?? 1;
    const budget = Math.min(this.labelOrder.length, Math.round(25 + zoom * 40));
    const wanted = new Set(this.labelOrder.slice(0, budget));
    if (this.selected !== null && this.view.map[this.selected] >= 0) wanted.add(this.view.map[this.selected]);
    if (this.hovered !== null) wanted.add(this.hovered);
    const isolated = this.handlers.isolated();
    if (isolated !== null && this.view.map[isolated] >= 0) wanted.add(this.view.map[isolated]);
    const dense = [...wanted];
    this.graph.trackPointPositionsByIndices(dense);
    const positions = this.graph.getTrackedPointPositionsMap();
    const width = this.container.clientWidth, height = this.container.clientHeight;
    const names = this.handlers.names;
    const nodes = [];
    const placed = [];
    // a label already on screen keeps its place over a newcomer of equal rank, so settling does not flicker
    const priority = (k) => this.view.back[k] === this.selected ? 4 : this.view.back[k] === isolated ? 3 : k === this.hovered ? 2 : this.shown.has(k) ? 1 : 0;
    for (const k of dense.sort((a, b) => priority(b) - priority(a))) {
      const p = positions.get(k);
      if (!p) continue;
      const [x, y] = this.graph.spaceToScreenPosition(p);
      if (x < -20 || y < -20 || x > width + 20 || y > height + 20) continue;
      // a label that would sit on another is dropped: high degree, the selection and the hover win
      const w = Math.min(260, names[this.view.back[k]].length * 6.5 + 24), h = 16;
      if (placed.some((r) => x < r.x + r.w && x + w > r.x && y - h / 2 < r.y + r.h && y + h / 2 > r.y)) continue;
      placed.push({x, y: y - h / 2, w, h});
      nodes.push({k, x, y});
    }
    const frag = document.createDocumentFragment();
    for (const {k, x, y} of nodes) {
      const index = this.view.back[k];
      const el = document.createElement('div');
      el.className = 'kgv-label' + (index === this.selected ? ' kgv-label-selected' : '') + (index === isolated ? ' kgv-label-isolated' : '');
      const text = document.createElement('span');
      text.textContent = names[index];
      const iso = document.createElement('span');
      iso.className = 'kgv-iso';
      iso.textContent = index === isolated ? '↩' : '⤢';
      iso.title = index === isolated ? 'Back to the previous view' : 'Show this node alone';
      iso.addEventListener('click', (ev) => { ev.stopPropagation(); this.handlers.onIsolate(index); });
      el.addEventListener('click', () => { this.select(index); this.handlers.onSelect(index); });
      el.append(text, iso);
      el.style.left = `${x}px`;
      el.style.top = `${y}px`;
      frag.append(el);
    }
    this.labelsLayer.replaceChildren(frag);
    this.shown = new Set(nodes.map((n) => n.k));
  }

  setLabels(on) {
    this.labelsOn = on;
    this.updateLabels();
  }
}
