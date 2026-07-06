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
 *  attempt — the dirty set keeps accumulating and the next interval retries. */

import {GraphEdit} from '../rete/flow-editor';

/** Delay between the last result-affecting edit and the automatic run. Long
 *  enough to swallow a burst (typing a value, dragging several wires), short
 *  enough to feel live. */
export const AUTORUN_DEBOUNCE_MS = 800;

export class AutorunScheduler {
  /** Whether autorun is on — toggled from the ribbon; off by default. */
  enabled = false;

  /** Invalidated node ids accumulated since the last successful run. */
  private dirty = new Set<string>();
  private timer: ReturnType<typeof setTimeout> | null = null;

  /** @param run attempt a run for the accumulated dirty set:
   *  - `'started'` — a run began; the set is consumed;
   *  - `'busy'`    — a run is still in progress; keep the set, retry after
   *                  another interval;
   *  - `'skipped'` — the graph can't autorun right now (validation errors,
   *                  a parameter dialog would be needed); keep the set but
   *                  don't poll — the next edit reschedules anyway.
   *  @param debounceMs override for tests. */
  constructor(
    private readonly run: (dirty: Set<string>) => 'started' | 'busy' | 'skipped',
    private readonly debounceMs = AUTORUN_DEBOUNCE_MS,
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
    if (!this.enabled) return;
    if (edit.kind === 'cleared') {
      this.reset();
      return;
    }
    // A fresh node changes nothing until it's wired; a removed node's
    // connection-removal events have already been fed separately.
    if (edit.kind === 'node-added' || edit.kind === 'node-removed') return;
    for (const id of affected) this.dirty.add(id);
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

  /** Cancel any pending run and forget the accumulated dirty set. */
  reset(): void {
    if (this.timer !== null) {
      clearTimeout(this.timer);
      this.timer = null;
    }
    this.dirty.clear();
  }

  private schedule(): void {
    if (this.timer !== null) clearTimeout(this.timer);
    this.timer = setTimeout(() => {
      this.timer = null;
      if (!this.enabled) return;
      const outcome = this.run(new Set(this.dirty));
      if (outcome === 'started')
        this.dirty.clear();
      else if (outcome === 'busy')
        this.schedule(); // a run is in progress — keep the set, retry later
      // 'skipped': keep the set; the next edit schedules again.
    }, this.debounceMs);
  }
}
