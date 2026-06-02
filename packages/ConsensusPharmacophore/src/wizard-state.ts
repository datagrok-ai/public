/*
 * Wizard step state machine for the Consensus Pharmacophore pipeline.
 *
 * Holds the per-step status (pending | active | running | done | stale | error)
 * for the 5-stage pipeline. Subscribers (the wizard rail in `wizard-shell.ts`)
 * react to status emissions and repaint individual rail nodes.
 *
 * Inspired by:
 *  - `packages/Compute2/src/components/TreeWizard/types.ts` — 9-state Status
 *    taxonomy. Simplified here to 5 user-facing states + an `active` marker.
 *  - `libraries/compute-utils/old-views/src/pipeline-view.ts` — uses RxJS
 *    BehaviorSubject per step. We use a single Subject + per-step lookup
 *    instead — fewer allocations, easier to debug.
 *  - `js-api/src/ui/wizard.ts` — `WizardPage.onActivated` hook pattern
 *    (we use a separate per-step content builder, not an onActivated callback,
 *    because state changes from outside the wizard shell also need to refresh).
 */

import {Subject} from 'rxjs';

/** Stable numeric IDs for the 5 pipeline steps. */
export enum StepId {
  Fetch    = 0,
  Align    = 1,
  Pocket   = 2,
  Features = 3,
  Consensus= 4,
}

export const STEP_COUNT = 5;

/** Human-readable labels for the rail, one per step. */
export const STEP_LABELS: ReadonlyArray<string> = [
  'Fetch PDBs',
  'Align',
  'Pocket',
  'Features',
  'Consensus',
];

/** Status taxonomy borrowed from TreeWizard, condensed to 5 user states. */
export type StepStatus = 'pending' | 'active' | 'running' | 'done' | 'stale' | 'error';

/** Emitted by `WizardState.changes` whenever any step's status changes. */
export interface StepSnapshot {
  id: StepId;
  status: StepStatus;
}

/**
 * Per-step state holder + change notifier. The wizard shell builds a rail
 * element for each step and subscribes to `changes` to repaint it. The shell
 * also queries the current status of a step via `get(step)` when computing
 * navigation button enable/disable rules.
 *
 * Cascade semantics:
 *  - `markStaleFrom(step)` — user edited inputs in `step`; downstream steps
 *    that had `done` results become `stale` (results still on screen, but
 *    flagged). Downstream steps that were `pending` stay pending.
 *  - `reset()` — wipe all state; first step becomes `active`. Used by the
 *    "Reset" ribbon button.
 */
export class WizardState {
  private readonly _statuses: StepStatus[] = Array(STEP_COUNT).fill('pending');

  /** Stream of status changes — one emission per `set()` call. */
  readonly changes = new Subject<StepSnapshot>();

  /** True after the first step has ever completed. Used by the shell to
   *  decide whether Mol* shows the empty placeholder or the real viewer. */
  anyDone = false;

  constructor() {
    // Step 1 starts as `active` so the rail's first node shows the blue ring
    // instead of the grey clock icon.
    this._statuses[StepId.Fetch] = 'active';
  }

  get(step: StepId): StepStatus {
    return this._statuses[step];
  }

  set(step: StepId, status: StepStatus): void {
    this._statuses[step] = status;
    this.changes.next({id: step, status});
    if (status === 'done') this.anyDone = true;
  }

  /** Called by the wizard shell when the orchestrator begins running a step. */
  markRunning(step: StepId): void {
    this.set(step, 'running');
  }

  /** Called by the wizard shell when the orchestrator completes a step. */
  markDone(step: StepId): void {
    this.set(step, 'done');
  }

  /** Called by the wizard shell when the orchestrator throws on a step. */
  markError(step: StepId): void {
    this.set(step, 'error');
  }

  /**
   * User edited inputs upstream — cascade staleness to all downstream steps.
   * `fromStep` itself transitions to `active`. Downstream `done`/`stale`
   * steps become `stale` (results still rendered, but a yellow dot signals
   * the user to re-run). Downstream `pending`/`error` steps are left alone.
   *
   * After this call, the wizard shell typically also calls
   * `app.invalidatePreviewCache()` to drop the orchestrator's per-fingerprint
   * cache, and shows the amber Mol* overlay + stale banner.
   */
  markStaleFrom(fromStep: StepId): void {
    this.set(fromStep, 'active');
    for (let s = fromStep + 1; s < STEP_COUNT; s++) {
      const cur = this._statuses[s];
      if (cur === 'done' || cur === 'stale')
        this.set(s as StepId, 'stale');
      // pending / error / running / active downstream: untouched
    }
  }

  /**
   * Navigate to `step` without mutating its status if it already represents
   * a real result (done/stale/error). Only pending steps flip to active.
   * Used by the rail click handler when the user revisits an earlier step.
   */
  activate(step: StepId): void {
    const cur = this._statuses[step];
    if (cur === 'pending')
      this.set(step, 'active');
    // Always re-emit so the shell can re-render the rail (active highlight
    // moves to the new step).
    this.changes.next({id: step, status: this._statuses[step]});
  }

  /** Hard reset — used by the "Reset" ribbon button. */
  reset(): void {
    for (let s = 0; s < STEP_COUNT; s++)
      this._statuses[s] = 'pending';
    this._statuses[StepId.Fetch] = 'active';
    this.anyDone = false;
    for (let s = 0; s < STEP_COUNT; s++)
      this.changes.next({id: s as StepId, status: this._statuses[s]});
  }

  /**
   * Whether at least one downstream step is `stale`. Drives the visibility
   * of the amber "Re-run to refresh" CTA in the wizard footer and ribbon.
   */
  hasStaleDownstream(fromStep: StepId): boolean {
    for (let s = fromStep + 1; s < STEP_COUNT; s++) {
      if (this._statuses[s] === 'stale') return true;
    }
    return false;
  }

  /**
   * Like `markStaleFrom`, but does NOT mutate the from-step's status — only
   * flips downstream `done` steps to `stale`. Used by the per-option stale
   * cascade so that, e.g., changing a Pocket option marks Step 3+ stale
   * without forcing Step 3 to become `active` (the user may be elsewhere).
   *
   * `step` is the EARLIEST step whose result is invalidated by the change.
   * Steps before it are untouched; steps from `step` onwards (inclusive)
   * that are currently `done` become `stale`.
   */
  markStaleStartingFrom(step: StepId): void {
    for (let s = step; s < STEP_COUNT; s++) {
      if (this._statuses[s] === 'done')
        this.set(s as StepId, 'stale');
    }
  }

  /**
   * Inverse of `markStaleFrom`: revert `stale` steps to `done` for which the
   * caller asserts cached data still exists. `hasCachedResult[step]` must be
   * true exactly when the orchestrator still has a usable result for that
   * step under the current fingerprint.
   *
   * Used by the wizard's option-revert path: if the user edits an option and
   * then edits it back to the value the last run used, the displayed data is
   * still valid — flip the stale banner off and revert step statuses.
   */
  clearStaleIfCached(hasCachedResult: ReadonlyArray<boolean>): void {
    for (let s = 0; s < STEP_COUNT; s++) {
      if (this._statuses[s] === 'stale' && hasCachedResult[s])
        this.set(s as StepId, 'done');
    }
  }

  /** True if any step is currently `stale`. */
  hasAnyStale(): boolean {
    for (let s = 0; s < STEP_COUNT; s++) {
      if (this._statuses[s] === 'stale') return true;
    }
    return false;
  }
}
