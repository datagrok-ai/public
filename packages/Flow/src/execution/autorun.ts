/** Autorun: debounced re-execution of the flow after every result-affecting
 *  edit.
 *
 *  The scheduler consumes the classified `GraphEdit` stream (after the
 *  execution controller invalidated the affected downstream cone) and
 *  accumulates the invalidated node ids. One debounce interval after the last
 *  edit it asks the host to run — handing over the accumulated dirty set so
 *  the host can re-run only the affected slice (live-boundary expansion in
 *  `ExecutionController.runAutorun`), not the whole graph.
 *
 *  Edits that cannot affect results (adding a disconnected node, removing an
 *  isolated one) schedule nothing. A run already in progress postpones the
 *  attempt — the dirty set keeps accumulating and the next interval retries.
 *
 *  Some nodes are LIVE BY DEFAULT (`AUTORUN_FUNC_NAMES` / `AUTORUN_NODE_TYPES`
 *  below): edits touching them schedule a run even while the ribbon autorun
 *  toggle is off, so an Open File loads on picking a file and a viewer
 *  refreshes on rewiring without ever pressing Run. */

import {GraphEdit} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';

/** Delay between the last result-affecting edit and the automatic run. Long
 *  enough to swallow a burst (typing a value, dragging several wires), short
 *  enough to feel live. */
export const AUTORUN_DEBOUNCE_MS = 1000;

// ---- Live-by-default nodes -------------------------------------------------
// These rerun automatically after an edit that touches them (their own params
// or connections, or an upstream change that invalidates them) even while the
// ribbon autorun toggle is OFF. Modify the two lists below to change what's
// live.

/** DG function simple names (case-insensitive) that are live by default. */
export const AUTORUN_FUNC_NAMES: string[] = [
  'OpenFile',
  'AddNewColumn',
];

/** Registered node-type names (`FlowNode.dgTypeName`) that are live by
 *  default; a trailing `*` matches a prefix — `'Viewers/*'` = every viewer. */
export const AUTORUN_NODE_TYPES: string[] = [
  'Viewers/*',
];

/** Whether edits touching this node trigger an automatic (debounced) rerun of
 *  the affected slice even when the global autorun toggle is off. */
export function isAutorunByDefault(node: FlowNode): boolean {
  const fn = node.dgFunc?.name?.toLowerCase();
  if (fn && AUTORUN_FUNC_NAMES.some((n) => n.toLowerCase() === fn)) return true;
  const typeName = node.dgTypeName ?? '';
  for (const pat of AUTORUN_NODE_TYPES) {
    if (pat.endsWith('*') ? typeName.startsWith(pat.slice(0, -1)) : typeName === pat)
      return true;
  }
  // Viewer nodes not registered under Viewers/ (constructed directly in code).
  return node.properties['viewerType'] != null;
}

export class AutorunScheduler {
  /** Whether autorun is on for EVERY node — toggled from the ribbon; off by
   *  default. Live-by-default nodes (see {@link isAutorunByDefault}) schedule
   *  runs regardless of this flag. */
  enabled = false;

  /** Invalidated node ids accumulated since the last successful run. */
  private dirty = new Set<string>();
  private timer: ReturnType<typeof setTimeout> | null = null;
  /** While > 0, firing is suspended (edits still accumulate) — see {@link hold}. */
  private holds = 0;

  /** @param run attempt a run for the accumulated dirty set. `liveOnly` is
   *  true when the global toggle is off (the set holds ONLY live-by-default
   *  node ids — the host must run just those, gated on satisfied inputs, via
   *  `ExecutionController.runLiveNodes`; with the toggle on it's the normal
   *  affected-slice `runAutorun`). Outcomes:
   *  - `'started'` — a run began; the set is consumed;
   *  - `'busy'`    — a run is still in progress; keep the set, retry after
   *                  another interval;
   *  - `'skipped'` — the graph can't autorun right now (validation errors,
   *                  a parameter dialog would be needed); keep the set but
   *                  don't poll — the next edit reschedules anyway.
   *  @param debounceMs override for tests.
   *  @param isLiveNode which node ids are live by default (host wires this to
   *  `isAutorunByDefault` over its graph); with the global toggle off, an edit
   *  schedules a run only when it affects at least one live node. */
  constructor(
    private readonly run: (dirty: Set<string>, liveOnly: boolean) => 'started' | 'busy' | 'skipped',
    private readonly debounceMs = AUTORUN_DEBOUNCE_MS,
    private readonly isLiveNode: (nodeId: string) => boolean = () => false,
  ) {}

  /** Flip the mode. Turning it off drops any pending run and dirty state. */
  toggle(): boolean {
    this.enabled = !this.enabled;
    if (!this.enabled) this.reset();
    return this.enabled;
  }

  /** Feed one classified edit and the node ids it invalidated (from
   *  `ExecutionController.applyGraphEdit`). Schedules a debounced run when the
   *  edit can affect results. */
  onEdit(edit: GraphEdit, affected: Set<string>): void {
    if (edit.kind === 'cleared') {
      this.reset();
      return;
    }
    // A removed node's connection-removal events have already been fed
    // separately.
    if (edit.kind === 'node-removed') return;
    if (edit.kind === 'node-added') {
      // A fresh node changes nothing for the graph until it's wired — but a
      // LIVE node dropped ready-to-run (a file dragged from the files tree
      // creates an Open File with the path pre-set) must run at once.
      // Readiness is re-checked at fire time, so a bare drop schedules and
      // then quietly no-ops. The toggle-on path keeps ignoring adds — a full
      // autorun would run far more than the dropped node.
      if (this.enabled || !this.isLiveNode(edit.nodeId)) return;
      this.dirty.add(edit.nodeId);
      this.schedule();
      return;
    }
    if (!this.enabled) {
      // Toggle off → only the live-by-default nodes the edit touched enter the
      // set (NOT the whole affected cone): the fired run executes exactly
      // those, so nothing else on the canvas ever runs uninvited.
      const live = [...affected].filter((id) => this.isLiveNode(id));
      if (live.length === 0) return;
      for (const id of live) this.dirty.add(id);
    }
    else {
      for (const id of affected) this.dirty.add(id);
    }
    this.schedule();
  }

  /** Schedule a debounced run for an externally computed dirty set — e.g. the
   *  moment autorun is switched on, everything without a fresh result
   *  (`ExecutionController.pendingNodes`) is kicked off without waiting for an
   *  edit. No-op while disabled. */
  kick(dirty: Set<string>): void {
    if (!this.enabled) return;
    for (const id of dirty) this.dirty.add(id);
    this.schedule();
  }

  /** Schedule the live-by-default nodes among `ids` without an edit — called
   *  when a flow is loaded, so its Open Files (and any other ready live node)
   *  run at once. Nodes whose inputs aren't satisfied are dropped at fire
   *  time by `runLiveNodes`, so kicking a whole freshly loaded graph is safe:
   *  a viewer with no captured upstream simply doesn't run. */
  kickLive(ids: Iterable<string>): void {
    let any = false;
    for (const id of ids) {
      if (!this.isLiveNode(id)) continue;
      this.dirty.add(id);
      any = true;
    }
    if (any) this.schedule();
  }

  /** Suspend firing while a modal interaction is in progress (edits still
   *  accumulate). Concretely: a function-editor dialog intercepts the global
   *  `d4-before-run-action` event, which fires for EVERY client funccall — an
   *  autorun kicking in mid-dialog would run the same function, get its call
   *  canceled by the dialog's interceptor, and resolve the dialog round-trip
   *  early with the wrong funccall (the user's OK then writes nothing back).
   *  Re-entrant: every `hold` needs a `release`. */
  hold(): void {
    this.holds++;
    if (this.timer !== null) {
      clearTimeout(this.timer);
      this.timer = null;
    }
  }

  /** Undo one {@link hold}; when the last one lifts, anything accumulated in
   *  the meantime is scheduled. */
  release(): void {
    this.holds = Math.max(0, this.holds - 1);
    // Anything in `dirty` was admitted by `onEdit`/`kick` (toggle on, or a
    // live-node edit) — reschedule it regardless of the current toggle state.
    if (this.holds === 0 && this.dirty.size > 0) this.schedule();
  }

  /** Cancel any pending run and forget the accumulated dirty set. */
  reset(): void {
    if (this.timer !== null) {
      clearTimeout(this.timer);
      this.timer = null;
    }
    this.dirty.clear();
  }

  private schedule(): void {
    if (this.holds > 0) return; // suspended — release() reschedules
    if (this.timer !== null) clearTimeout(this.timer);
    this.timer = setTimeout(() => {
      this.timer = null;
      // No `enabled` check: a pending set exists only if `onEdit` admitted it
      // (toggle on, or a live-node edit); toggling off resets and cancels.
      if (this.holds > 0) return;
      const outcome = this.run(new Set(this.dirty), !this.enabled);
      if (outcome === 'started')
        this.dirty.clear();
      else if (outcome === 'busy')
        this.schedule(); // a run is in progress — keep the set, retry later
      // 'skipped': keep the set; the next edit schedules again.
    }, this.debounceMs);
  }
}
