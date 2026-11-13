/**
 * Stage strip — the one-glance answer to "where am I, what is done, what is next"
 * (nca-workflow-stages F-S3.1 / D16). ONE row (user-directed 2026-09-30, superseding the approved
 * two-line design): the chain as a native `ui.h1`, the next-action hint right-aligned beside it,
 * on a `--warm-grey-1` band.
 *
 *   ✓ Setup  ›  ● Review λz 7 left  ›  Report                      Review 7 profiles →
 *
 * **Variant A (user-directed 2026-09-29):** a compact custom chain — done stages carry ✓, the
 * current stage is bold (no underline, no severity fill), future stages are muted, and the
 * remaining count is amber TEXT in the Flagged verdict hue — exactly the charts' Flagged colour, no
 * fill (user-directed 2026-09-29). No stage numbers, no underline. The next-action hint has no
 * "Next:" prefix; it is a link only for `Review N profiles →` (and the study name in `In study X`).
 *
 * The stage + count come ONLY from `deriveStage` (via the caller's `StageResult`) and the wording
 * ONLY from `workflow-stage.ts` — this module never derives either (AC-S3.1.1).
 *
 * Leaf module (imports adapter + status-palette only), so it never joins the views/ import cycle.
 */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {STAGE_LABEL, formatPendingCount, formatNextHint} from '../adapter';
import type {StageResult} from '../adapter';
import {STATUS_FLAGGED} from './status-palette';

/** The three stages the analyst can SHOW (D7). `complete` is a derived state, shown as Report. */
export type ViewStage = 'setup' | 'review' | 'report';

export interface WorkflowStripModel {
  /** The derived state (D1). */
  readonly stage: StageResult;
  /** The stage currently shown. */
  readonly viewStage: ViewStage;
  /** The study the analysis belongs to, or null (standalone / unsaved) — named in the hint. */
  readonly studyLabel?: string | null;
  /** Confirm is enabled (no failing Setup term, something to compute). On Setup this makes the
   *  Review stage a link whose click IS Confirm (user-directed 2026-10-13). */
  readonly confirmable?: boolean;
}

export interface WorkflowStrip {
  readonly root: HTMLElement;
  update(model: WorkflowStripModel): void;
}

/** Which stages the analyst may enter from the current derived state. Setup is always
 *  reachable; Review and Report need confirmed, unchanged results (P6 lock). Report is reachable
 *  with profiles still to review — the D14 soft gate. */
export function reachableStages(s: StageResult): ReadonlySet<ViewStage> {
  const open = s.stage !== 'setup';
  return new Set<ViewStage>(open ? ['setup', 'review', 'report'] : ['setup']);
}

/** Build the strip. `onNavigate` fires when a reachable, non-current stage is clicked (on Setup
 *  with `confirmable`, Review is reachable and its click means Confirm — the caller decides);
 *  `onStartReview` makes the Review-stage hint an action ("Review N profiles →", user-directed
 *  2026-09-30): it jumps to the first profile still to review. `onOpenStudy` makes the study NAME in
 *  the hint a link to that study (user-directed 2026-09-30). Everything else is plain text. */
export function createWorkflowStrip(
  onNavigate: (s: ViewStage) => void, onStartReview?: () => void, onOpenStudy?: () => void,
): WorkflowStrip {
  // A native `ui.h1` (user-directed 2026-09-30: the chain reads at h1 size) — the platform's own
  // heading size, not a restated pixel value; margin + weight are reset in `.nca-stage-chain`.
  const chain = ui.h1('', 'nca-stage-chain');
  const hint = ui.divText('', 'nca-stage-hint');
  // ONE row: the chain, then the hint beside it (user-directed 2026-09-30).
  const root = ui.divH([chain, hint], 'nca-stage-strip');
  // Tooltip only while the hint is an action (user-directed 2026-09-30); bound once, read live.
  let hintTip: string | null = null;
  ui.tooltip.bind(hint, () => hintTip);

  const item = (
    key: ViewStage, label: string, mark: string, model: WorkflowStripModel,
    reachable: ReadonlySet<ViewStage>, count: string | null,
  ): HTMLElement => {
    const current = model.viewStage === key;
    const canGo = !current && reachable.has(key);
    const text = ui.span([mark === '' ? label : `${mark} ${label}`], 'nca-stage-text');
    const parts: HTMLElement[] = [text];
    if (count !== null) {
      // The remaining count as amber text — the ONE coloured element (approved intake row), the SAME
      // hue as the Flagged profiles on the charts (status-palette), no background fill.
      const pill = ui.span([count], 'nca-stage-count');
      pill.style.color = DG.Color.toHtml(STATUS_FLAGGED);
      parts.push(pill);
    }
    const el = ui.span(parts, 'nca-stage-item');
    el.classList.toggle('nca-stage-current', current);
    el.classList.toggle('nca-stage-link', canGo);
    el.classList.toggle('nca-stage-muted', !current && !canGo);
    if (canGo) el.onclick = () => onNavigate(key);
    return el;
  };

  const update = (model: WorkflowStripModel): void => {
    const s = model.stage;
    // On Setup, Review is reachable THROUGH Confirm whenever Confirm can run — the caller's
    // `onNavigate('review')` confirms first (user-directed 2026-10-13). Report still needs results.
    const base = reachableStages(s);
    const reachable: ReadonlySet<ViewStage> = s.stage === 'setup' && model.confirmable === true ?
      new Set<ViewStage>([...base, 'review']) : base;
    const setupDone = s.stage !== 'setup';
    const reviewDone = s.stage === 'complete';
    const pending = s.stage === 'review' ? formatPendingCount(s.pending) : null;
    ui.empty(chain);
    // A CSS-drawn chevron (no glyph, no icon widget): the "›" glyph sat below the text midline, and a
    // FontAwesome icon would hit `Widget.kill` when `ui.empty(chain)` re-renders the strip.
    const sep = (): HTMLElement => ui.span([], 'nca-stage-sep');
    chain.append(
      // No stage numbers (user-directed 2026-09-29): ✓ = done, ● = shown, nothing otherwise.
      item('setup', STAGE_LABEL.setup, setupDone ? '✓' : (model.viewStage === 'setup' ? '●' : ''),
        model, reachable, null), sep(),
      item('review', STAGE_LABEL.review, reviewDone ? '✓' : (model.viewStage === 'review' ? '●' : ''),
        model, reachable, pending), sep(),
      // Report is DONE once the analysis is in a study (user-directed 2026-09-30) — the ✓ is the
      // visible "it was added"; a Setup edit since then un-does it until the next Confirm.
      item('report', STAGE_LABEL.complete,
        model.studyLabel != null && !s.setupChanged ? '✓' : (model.viewStage === 'report' ? '●' : ''),
        model, reachable, null),
    );
    // While profiles are left to review (and Setup is unchanged), the hint IS the way in: a link
    // that starts the review at the first pending profile. Otherwise it is plain text.
    const studyLabel = model.studyLabel ?? null;
    const text = formatNextHint(s, studyLabel, model.viewStage);
    const review = !s.setupChanged && s.stage === 'review' && onStartReview !== undefined;
    hintTip = review ? 'Click to start the review' : null;
    hint.classList.toggle('nca-stage-hint-action', review);
    hint.onclick = review ? () => onStartReview!() : null;
    const full = review ? `${text} →` : text;
    // The study name is a link to the study; the rest of the line keeps its own behaviour.
    const prefix = studyLabel !== null ? `In study '${studyLabel}'` : '';
    if (prefix !== '' && onOpenStudy !== undefined && full.startsWith(prefix)) {
      const link = ui.span([studyLabel!], 'nca-stage-hint-link');
      link.onclick = (e: MouseEvent) => {
        e.stopPropagation(); // not the whole-hint "start review" action
        onOpenStudy();
      };
      ui.tooltip.bind(link, () => `Open study '${studyLabel}'`);
      hint.replaceChildren('In study ', link, full.slice(prefix.length));
    } else hint.replaceChildren(full);
  };

  return {root, update};
}
