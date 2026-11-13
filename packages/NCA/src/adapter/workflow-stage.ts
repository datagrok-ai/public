/**
 * Workflow stage — the ONE derivation of where an analysis stands, and the ONE wording of it
 * (nca-workflow-stages F-S1.5 / AD-9, research P4; D1, D7, D16, D18, D19).
 *
 * Every reader — the in-view stage strip and "Next: …" hint, the hub Analyses badge, the
 * Browse tree, the studies-list summary, the study-grid Warnings line and the D14 dialog —
 * calls these functions. None re-derives a stage or re-words a badge, so the surfaces
 * cannot drift (the same single-source principle as the engine-staleness predicate).
 *
 * The stage is DERIVED from the analysis state (D1), never stored as truth: the entry's
 * `stageCache` is a copy of {@link deriveStage} output for the hub, written only together with
 * the state it derives from (AD-8), and read only through {@link formatStageBadge} together with
 * the coordinator's `isStale` overlay.
 *
 * Pure: no `grok` / `ui` imports.
 */

import type {StageCache, WorkflowStage} from './study-schema';
/** Flags-cell line for a row whose cached numbers were computed by a DIFFERENT (or unknown)
 *  sci-comp version than the one installed (v14 `engineVersion`, F-N.6). The SAME wording as
 *  the overview's per-analysis stale badge (`markDrift`), placed at the ROW where the number
 *  is — a cache that predates a numbers-changing engine release must not read as current
 *  (thesis #4). Cleared by the analysis's next open (lazy recompute re-stamps the version).
 *  Lives here, with every other stage / badge wording (moved from `study-grid.ts`, which
 *  re-exports it — that removed a study-grid ↔ workflow-stage import cycle). */
export const STALE_ENGINE_FLAG = 'stale (refreshes on open)';

/** What {@link deriveStage} needs about one review unit. */
export interface StageUnit {
  /** `review` = the unit needs an analyst decision (Flagged / Not estimated / No exposure);
   *  `not-needed` = Accepted; `not-applicable` = decided at another unit (SF-3). */
  readonly need: 'review' | 'not-needed' | 'not-applicable';
  /** True when the unit has a decision that is still valid (fingerprint match). */
  readonly decided: boolean;
}

export interface StageInput {
  /** Fingerprint of the last Confirm, or `null` when never confirmed. */
  readonly confirmedFingerprint: string | null;
  /** Fingerprint of the recipe as it is now (the live Setup inputs). */
  readonly currentFingerprint: string;
  /** True when results computed from the confirmed recipe are held. */
  readonly hasResults: boolean;
  readonly units: ReadonlyArray<StageUnit>;
}

export interface StageResult extends StageCache {
  /** True when a Confirm exists but Setup was edited since (P6 keep-and-mark): Review and
   *  Report are locked until the next Confirm (the button then reads `Apply`), and `Discard`
   *  (formerly `Revert`) is offered. */
  readonly setupChanged: boolean;
}

/**
 * Derive the stage (AC-S1.5.1):
 *
 * | Confirmed (and unchanged since)? | Pending units | Stage      |
 * |----------------------------------|---------------|------------|
 * | no                               | any           | `setup`    |
 * | yes                              | > 0           | `review`   |
 * | yes                              | 0             | `complete` |
 *
 * "Confirmed" = a Confirm exists, results are held, AND the current recipe fingerprint equals
 * the confirmed one — a Setup edit that has not been re-confirmed puts the analysis back to
 * `setup` (P6: the hub shows `Setup`, and a reopen lands on Setup with Discard available). A
 * toggle-and-untoggle is not an edit: the fingerprint is the same (AC-S1.4.2).
 */
export function deriveStage(input: StageInput): StageResult {
  let pending = 0;
  let notApplicable = 0;
  for (const u of input.units) {
    if (u.need === 'not-applicable') notApplicable++;
    else if (u.need === 'review' && !u.decided) pending++;
  }
  const hasConfirm = input.confirmedFingerprint !== null;
  const setupChanged = hasConfirm && input.confirmedFingerprint !== input.currentFingerprint;
  const confirmed = hasConfirm && !setupChanged && input.hasResults;
  const stage: WorkflowStage = !confirmed ? 'setup' : (pending > 0 ? 'review' : 'complete');
  return {stage, pending, notApplicable, setupChanged};
}

/** Stage-strip / badge label of a stage (D7, D16). */
export const STAGE_LABEL: Readonly<Record<WorkflowStage, string>> = {
  setup: 'Setup',
  review: 'Review λz',
  complete: 'Report',
};

/** The remaining-count phrase shown beside `Review λz` (D16): `4 left`. The part the views tint
 *  amber (intake §1 colour row). */
export function formatPendingCount(pending: number): string {
  return `${pending} left`;
}

/**
 * The hub badge for a cached stage (D19): `Setup` / `Review λz · N left` / `Complete`. `Report`
 * folds into `Complete` (Report has no action of its own, D18).
 *
 * **Staleness is an overlay, never a replacement** (R1/F5): a stale entry reads e.g.
 * `Complete · stale (refreshes on open)` — never a bare `Complete` (AC-S1.5.2). The wording is
 * the shipped `STALE_ENGINE_FLAG`, so the analyst learns one phrase for "these numbers are not
 * current". A `null` cache (never saved with a derived stage) reads `Setup`.
 */
export function formatStageBadge(cache: StageCache | null | undefined, isStale: boolean): string {
  const base = cache === null || cache === undefined || cache.stage === 'setup' ? 'Setup' :
    cache.stage === 'review' ? `${STAGE_LABEL.review} · ${formatPendingCount(cache.pending)}` :
      'Complete';
  return isStale ? `${base} · ${STALE_ENGINE_FLAG}` : base;
}

/** The study-grid Warnings line for an analysis added with unreviewed profiles (D14/D19):
 *  `N not reviewed`. `''` when nothing is pending. */
export function formatNotReviewed(cache: StageCache | null | undefined): string {
  const n = cache?.stage === 'review' ? cache.pending : 0;
  return n > 0 ? `${n} not reviewed` : '';
}

/**
 * The studies-list Analyses summary (D19): `3 analyses · 1 in review`, plus `· 1 in setup` when
 * any analysis is not yet confirmed. The count phrase reuses the list's own singular/plural form.
 * `pipelineFlags` (the study's {@link pipelineFlagCount} total) adds `· N pipeline flags` — an
 * overlay, never a replacement: a study whose analyses are all Complete but carry pipeline flags
 * must not read as plainly done (nca-pipeline-save-rework AC-9.4a).
 */
export function formatStudyAnalysesSummary(
  caches: ReadonlyArray<StageCache | null | undefined>, pipelineFlags: number = 0,
): string {
  const n = caches.length;
  let inReview = 0;
  let inSetup = 0;
  for (const c of caches) {
    // A missing cache (an entry saved before v15, never re-saved) is an UNKNOWN stage — not
    // Setup: a study entry always has results (AC-S1.7.2), so claiming "in setup" would be false.
    if (c?.stage === 'review') inReview++;
    else if (c?.stage === 'setup') inSetup++;
  }
  const parts = [`${n} ${n === 1 ? 'analysis' : 'analyses'}`];
  if (inReview > 0) parts.push(`${inReview} in review`);
  if (inSetup > 0) parts.push(`${inSetup} in setup`);
  const flags = formatPipelineFlags(pipelineFlags);
  if (flags !== '') parts.push(flags);
  return parts.join(' · ');
}

// ─── Pipeline wordings (nca-pipeline-save-rework, approved design 2026-10-15) ──────────
// One home for every wording a pipeline-saved entry adds, so the study grid, the studies list and
// the analysis view say the same thing (the single-source principle of the engine-staleness row).

/** Profiles with at least one reason of the pipeline's OWN (`PipelineFlag.pipelineOnly`) — a
 *  reason the app's acceptance check does not raise. Absent flags ≡ 0. */
export function pipelineFlagCount(
  flags: ReadonlyArray<{readonly pipelineOnly?: ReadonlyArray<string>}> | undefined,
): number {
  return (flags ?? []).filter((f) => (f.pipelineOnly ?? []).length > 0).length;
}

/** `N pipeline flags` / `1 pipeline flag`; `''` for none. */
export function formatPipelineFlags(n: number): string {
  return n <= 0 ? '' : `${n} pipeline flag${n === 1 ? '' : 's'}`;
}

/** The study-grid line of a SPARSE analysis (AC-5.8b): its study-tier numbers are aggregates of the
 *  per-animal NCA (1–2 samples per animal — Cmax a mean of single samples, λz-dependent values
 *  blank), not the composite estimate. Measured on the mouse fixture at build step 0. */
export const PER_ANIMAL_SPARSE_NOTE = 'per-animal values, not the composite';

/** The line of a pipeline entry that superseded `label` when the sampling design changed. */
export function formatSupersedes(label: string): string {
  return `supersedes '${label}'`;
}

/** The line of a pipeline entry left in place, superseded by `label`. */
export function formatSupersededBy(label: string): string {
  return `superseded by '${label}'`;
}

/** The one-line next-action hint under the stage strip (D16) — the action itself, no "Next:"
 *  prefix (user-directed 2026-09-29). ONE message at a time (user-directed 2026-09-30): while there
 *  is a next action it is shown ALONE; `In study '<label>'` appears only when there is nothing left
 *  to do — a Complete analysis on Report that is in a study — as the visible answer to "did Add to
 *  study work?". `studyLabel` = the study the analysis belongs to (null when standalone / never
 *  saved). The stage itself is untouched (D18: Complete = nothing left to review, never "added").
 *
 *  `shown` = the stage currently SHOWN. A Complete analysis shown anywhere but Report says where to
 *  go next AND how to revisit a single profile (user-directed 2026-09-30). On Report outside a
 *  study it says what Report offers — grouping for aggregated results, and adding to a study. */
export function formatNextHint(
  s: StageResult, studyLabel: string | null = null, shown: 'setup' | 'review' | 'report' = 'report',
): string {
  if (s.setupChanged) return 'Apply changes'; // the button reads Apply then (user-directed 2026-10-13)
  if (s.stage === 'setup') return 'Complete the setup and click Confirm';
  if (s.stage === 'review') return `Review ${s.pending} ${s.pending === 1 ? 'profile' : 'profiles'}`;
  if (shown !== 'report')
    return 'All profiles reviewed — go to Report, or click a row or its Review cell to revisit a profile';
  if (studyLabel !== null) return `In study '${studyLabel}'`;
  return 'Group the data to see aggregated results. Add the results to a study.';
}
