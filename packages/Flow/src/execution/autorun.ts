/** Autorun: debounced re-execution of the affected slice after every
 *  result-affecting edit, fed by the classified `GraphEdit` stream. */

import {GraphEdit} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {safeGet} from '../utils/dart-proxy-utils';

/** Delay between the last result-affecting edit and the automatic run. */
export const AUTORUN_DEBOUNCE_MS = 1000;

// Live-by-default nodes rerun even while the ribbon autorun toggle is OFF: the
// two lists below plus any function declaring `meta.autorun: true` on itself.
// Edit the lists for functions you don't own; prefer the meta for ones you do.

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

/** Whether edits touching this node rerun it even when the autorun toggle is off. */
export function isAutorunByDefault(node: FlowNode): boolean {
  const fn = node.dgFunc?.name?.toLowerCase();
  if (fn && AUTORUN_FUNC_NAMES.some((n) => n.toLowerCase() === fn)) return true;
  try {
    const meta = safeGet(node.dgFunc?.options, 'autorun');
    if (meta === true || String(meta).toLowerCase() === 'true') return true;
  } catch {/* options can throw on odd Dart proxies — fall through */}
  const typeName = node.dgTypeName ?? '';
  for (const pat of AUTORUN_NODE_TYPES) {
    if (pat.endsWith('*') ? typeName.startsWith(pat.slice(0, -1)) : typeName === pat)
      return true;
  }
  // Viewer nodes not registered under Viewers/ (constructed directly in code).
  return node.properties['viewerType'] != null;
}

export class AutorunScheduler {
  /** Autorun for EVERY node — toggled from the ribbon; live-by-default nodes
   *  schedule runs regardless of this flag. */
  enabled = false;

  private dirty = new Set<string>();
  private timer: ReturnType<typeof setTimeout> | null = null;
  private holds = 0;

  /** `run` outcomes: 'started' consumes the dirty set; 'busy' keeps it and
   *  retries after another interval; 'skipped' keeps it until the next edit.
   *  `liveOnly` is true when the toggle is off (the set holds only
   *  live-by-default node ids — run just those via `runLiveNodes`). */
  constructor(
    private readonly run: (dirty: Set<string>, liveOnly: boolean) => 'started' | 'busy' | 'skipped',
    private readonly debounceMs = AUTORUN_DEBOUNCE_MS,
    private readonly isLiveNode: (nodeId: string) => boolean = () => false,
  ) {}

  toggle(): boolean {
    this.enabled = !this.enabled;
    if (!this.enabled) this.reset();
    return this.enabled;
  }

  /** Feed one classified edit and the node ids it invalidated. */
  onEdit(edit: GraphEdit, affected: Set<string>): void {
    if (edit.kind === 'cleared') {
      this.reset();
      return;
    }
    // A removed node's connection-removal events have already been fed separately.
    if (edit.kind === 'node-removed') return;
    if (edit.kind === 'node-added') {
      // A live node dropped ready-to-run (a file dragged from the files tree)
      // must run at once; readiness is re-checked at fire time, so a bare
      // toolbox drop schedules and then quietly no-ops.
      if (this.enabled || !this.isLiveNode(edit.nodeId)) return;
      this.dirty.add(edit.nodeId);
      this.schedule();
      return;
    }
    if (!this.enabled) {
      // Only the live nodes the edit touched (never the whole cone) enter the
      // set, so nothing else on the canvas ever runs uninvited.
      const live = [...affected].filter((id) => this.isLiveNode(id));
      if (live.length === 0) return;
      for (const id of live) this.dirty.add(id);
    }
    else {
      for (const id of affected) this.dirty.add(id);
    }
    this.schedule();
  }

  /** Schedule an externally computed dirty set (e.g. on toggle-on); no-op while disabled. */
  kick(dirty: Set<string>): void {
    if (!this.enabled) return;
    for (const id of dirty) this.dirty.add(id);
    this.schedule();
  }

  /** Schedule the live-by-default nodes among `ids` without an edit (flow load). */
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
   *  accumulate) — an autorun mid-dialog gets hijacked by the function editor's
   *  `d4-before-run-action` interceptor. Re-entrant: every `hold` needs a `release`. */
  hold(): void {
    this.holds++;
    if (this.timer !== null) {
      clearTimeout(this.timer);
      this.timer = null;
    }
  }

  /** Undo one {@link hold}; when the last one lifts, the backlog is scheduled. */
  release(): void {
    this.holds = Math.max(0, this.holds - 1);
    if (this.holds === 0 && this.dirty.size > 0) this.schedule();
  }

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
      // No `enabled` check: a pending set exists only if `onEdit` admitted it.
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
