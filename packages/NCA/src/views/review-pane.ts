/**
 * Context-panel `Review` pane (nca-workflow-stages F-S4.2 / F-S4.3 / F-S4.4; D11, D12, D15;
 * approved design 2026-09-26).
 *
 * One pane, first in the λz readout accordion, where the decision and the move to the next
 * profile sit together (D15): the unit's state, ONE untitled radio choice — `Adjust λz points…`
 * (the third remedy, the default where it applies), then its allowed actions (`Accept as is` /
 * `Exclude λz` / `Acknowledge`), hidden once a decision is saved (user-directed 2026-10-12) —
 * and `Reset`, an INLINE reason picker (criterion reasons
 * pre-checked from the profile's own violated flags, judgement reasons keyed by
 * (action, category), `Other` free text; Save dimmed until ≥ 1 reason — D12), and
 * `Next to review → (k of N)`. With ≥ 2 selected pending units that share the same allowed
 * actions, it offers ONE decision for all of them (bulk, RG-4 resolved "same action class"):
 * criterion reasons = the intersection of their violated criteria; each unit still gets its own
 * decision, fingerprint and audit entry (the caller's job).
 *
 * Read-only on Report. Pure rendering: every controller call goes through {@link ReviewPaneArgs}
 * callbacks, so this leaf module imports only the adapter + ui (never the views/ cycle).
 */

import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';

import {
  judgementReasons, criterionReasons, reasonsValid, formatReasons, REVIEW_ACTION_LABEL, REVIEW_STATE_LABEL,
  OTHER_REASON_CODE, formatOverrideTimestamp, sharedFlagCodes,
} from '../adapter';
import type {
  ReviewUnit, ReviewAction, ReviewReason, AcceptanceThresholds,
  ReviewCategory, ReasonDef,
} from '../adapter';

export interface ReviewPaneArgs {
  /** The unit shown (the current row's profile, or a sparse composite). */
  readonly unit: ReviewUnit;
  /** Other SELECTED pending units (bulk candidates), excluding `unit`. */
  readonly selected: ReadonlyArray<ReviewUnit>;
  /** The pending units in the problem-first order `Next to review` walks. */
  readonly pendingOrder: ReadonlyArray<string>;
  readonly thresholds: AcceptanceThresholds;
  /** Review λz = editable; Report = read-only (D20). */
  readonly editable: boolean;
  /** The pipeline run's own reasons for this unit, as one line (nca-pipeline-save-rework F10), or
   *  undefined when it has none. Shown muted in every state, so the scientist sees a CRO incident
   *  on a profile the λz check accepted, and can cite it while choosing a reason. */
  readonly pipelineLine?: string;
  onDecide(unitKeys: ReadonlyArray<string>, action: ReviewAction, reasons: ReadonlyArray<ReviewReason>): void;
  onReset(unitKey: string): void;
  onNext(): void;
}

const STATE_WORD = REVIEW_STATE_LABEL; // the ONE state wording (review-reasons.ts)

function muted(text: string): HTMLElement {
  const el = ui.divText(text);
  el.style.color = 'var(--grey-5, #7D8086)';
  el.style.whiteSpace = 'normal';
  return el;
}

/** The category's plain-language "why this needs a decision" line, or null for a Flagged profile —
 *  its violated criteria already head the readout (user-directed 2026-10-12). */
function whyLine(category: string): string | null {
  if (category === 'no-exposure') return 'No measurable concentration — all BLQ / zero / missing.';
  if (category === 'not-estimated') return 'λz could not be fitted, so t½, AUC∞ and CL are not estimated.';
  return null;
}

/** The third remedy, beside the decisions: re-selecting the terminal points in the λz editor. */
const ADJUST_OPTION = 'Adjust λz points in the λz editor'; // "click a point" rides its tooltip

/** Per-option tooltips of the review choice (the former buttons' tooltips). */
const OPTION_TIP: Readonly<Record<string, string>> = {
  [ADJUST_OPTION]: 'Re-select the terminal points: click a point in the λz editor to force it into or out of the fit',
  [REVIEW_ACTION_LABEL.accept]: 'Keep the numbers as they are, with a reason',
  [REVIEW_ACTION_LABEL.exclude]:
    'λz-dependent values (t½, AUC∞, CL, Vz, MRT, Vss, …) become not reportable; Cmax, Tmax, AUClast stay',
  [REVIEW_ACTION_LABEL.acknowledge]: 'Record that you have looked at it',
};

/** Adjusting λz points applies to a per-subject Flagged / Not estimated profile only — a
 *  no-exposure profile has no points to fit, and a sparse composite is not edited in the λz editor. */
function canAdjust(u: ReviewUnit): boolean {
  return u.kind === 'profile' && (u.category === 'flagged' || u.category === 'not-estimated');
}

/** The review choice (user-directed 2026-10-12): ONE untitled radio — `Adjust λz points…` (the
 *  default, where it applies), then the unit's decisions. Picking a decision opens its reason form;
 *  once a decision is saved the choice is hidden (only "To review" shows it). */
function reviewChoice(u: ReviewUnit, value: string | null,
  onPick: (v: string | null) => void): DG.InputBase<string | null> {
  const items = [...(canAdjust(u) ? [ADJUST_OPTION] : []), ...u.allowedActions.map((a) => REVIEW_ACTION_LABEL[a])];
  const radio = ui.input.radio('', {items, nullable: true, ...(value !== null ? {value} : {}),
    onValueChanged: onPick});
  radio.root.classList.add('nca-review-choice');
  for (const label of Array.from(radio.root.querySelectorAll<HTMLElement>('.ui-input-editor label'))) {
    const tip = OPTION_TIP[(label.textContent ?? '').trim()];
    if (tip !== undefined) ui.tooltip.bind(label, tip);
  }
  return radio;
}

/** The decision a choice label stands for, or null for `Adjust λz points…`. */
function actionOf(u: ReviewUnit, label: string | null): ReviewAction | null {
  return u.allowedActions.find((a) => REVIEW_ACTION_LABEL[a] === label) ?? null;
}


/** Same allowed-action SET (order-insensitive). */
function sameActions(a: ReadonlyArray<ReviewAction>, b: ReadonlyArray<ReviewAction>): boolean {
  return a.length === b.length && a.every((x) => b.includes(x));
}

/** Build the pane content into `host` (replacing it). */
export function renderReviewPane(host: HTMLElement, args: ReviewPaneArgs): void {
  ui.empty(host);
  const u = args.unit;
  const body = ui.divV([], {style: {gap: '6px', width: '100%'}});
  host.appendChild(body);
  // The pipeline run's own reasons lead the pane in every state (approved design 2026-10-15).
  if (args.pipelineLine !== undefined && args.pipelineLine !== '') body.appendChild(muted(args.pipelineLine));

  if (u.need === 'not-applicable') {
    body.appendChild(muted('Not applicable — reviewed at the group composite.'));
    if (u.blqOnly) {
      const times = [...new Set(u.blqTimes)].join(', ');
      body.appendChild(muted(`This animal's sample(s) are BLQ at ${times === '' ? 'every sampled time' : times}.`));
    }
    return;
  }
  if (u.need === 'not-needed') {
    body.appendChild(muted('Not needed — the λz fit meets every acceptance criterion.'));
    appendNext(body, args);
    return;
  }

  if (u.decision !== null) {
    const d = u.decision;
    // No review choice once a decision is saved (user-directed 2026-10-12): the state below says
    // what was decided; Reset re-renders "To review", whose choice starts on its top option.
    // Reset is an icon beside the decision title, not a button row (user-directed 2026-09-30).
    const title = ui.divH([ui.divText(STATE_WORD[d.action], {style: {fontWeight: '600'}})],
      {style: {gap: '8px', alignItems: 'center'}});
    if (args.editable) {
      const reset = ui.iconFA('undo', () => args.onReset(u.key),
        'Reset — return this profile to "To review" (the decision is kept in the history)');
      reset.style.color = 'var(--blue-1, #2083D5)'; // the pane's action blue (the buttons' colour)
      title.appendChild(reset);
    }
    body.appendChild(title);
    body.appendChild(muted(`Reasons: ${formatReasons(d.reasons)}`));
    body.appendChild(muted(`${d.user === '' ? 'unknown user' : d.user} · ${formatOverrideTimestamp(d.at)}`));
    appendNext(body, args);
    return;
  }

  // To review.
  const whyText = whyLine(u.category);
  if (whyText !== null) body.appendChild(muted(whyText));
  if (!args.editable) {
    body.appendChild(muted('To review — decide on the Review λz stage.'));
    return;
  }
  // Same action set AND same category (review F7): the judgement reasons are keyed by
  // (action, category) — AC-S1.6.0 — so a mixed Not-estimated / No-exposure batch would apply one
  // category's reasons (e.g. "true finding (BLQ)") to a profile of the other.
  const bulk = args.selected.filter((s) => s.need === 'review' && s.decision === null &&
    s.category === u.category && sameActions(s.allowedActions, u.allowedActions));
  const targets = [u, ...bulk];
  if (bulk.length > 0)
    body.appendChild(ui.divText(`Decide for ${targets.length} selected profiles`, {style: {fontWeight: '600'}}));
  // (the bulk line stays above the choice: it says WHO the decision is for, open form or not)
  const picker = ui.div([]);
  // The choice stays visible while a reason form is open: picking another decision swaps the
  // form, `Adjust λz points…` (or the form's Cancel) closes it (user-directed 2026-10-12).
  const initial = canAdjust(u) ? ADJUST_OPTION : null;
  const choice = reviewChoice(u, initial, (v) => {
    const action = actionOf(u, v);
    if (action === null) {
      ui.empty(picker);
      return;
    }
    openPicker(picker, action, targets, args, () => {choice.value = initial;});
  });
  body.append(choice.root, picker);
  appendNext(body, args);
}

/** The per-action heading of the reason form (user-directed 2026-09-30). */
const REASON_TITLE: Readonly<Record<ReviewAction, string>> = {
  accept: 'Accept reason', exclude: 'Exclude reason', acknowledge: 'Acknowledge reason',
};

/** The inline reason picker for one action (D12): criterion reasons (pre-checked), judgement
 *  reasons for (action, category), and Other — Save dimmed until ≥ 1 reason. ONE `ui.form` holding
 *  the heading (`Accept reason` / `Exclude reason` …), the inputs and the Save / Cancel buttons;
 *  short labels, the full reason on each tooltip (user-directed 2026-09-30). */
function openPicker(
  host: HTMLElement, action: ReviewAction, targets: ReadonlyArray<ReviewUnit>, args: ReviewPaneArgs,
  onCancel: () => void,
): void {
  ui.empty(host);
  const category = targets[0].category as ReviewCategory;
  const crit: ReadonlyArray<ReasonDef> = category === 'flagged' ?
    criterionReasons(sharedFlagCodes(targets), args.thresholds) : [];
  const judge = judgementReasons(action, category);
  const boxes: {def: ReasonDef; input: DG.InputBase<boolean | null>}[] = [];
  const save = ui.button('Save', () => {
    const reasons = collect();
    if (!reasonsValid(reasons)) return;
    args.onDecide(targets.map((t) => t.key), action, reasons);
  });
  const other = ui.input.string('Other', {value: '', tooltipText: 'Your own reason, in free text'});
  const collect = (): ReviewReason[] => {
    const out: ReviewReason[] = boxes.filter((b) => b.input.value === true)
      .map((b) => ({code: b.def.code, text: b.def.text}));
    const t = (other.value ?? '').trim();
    if (t !== '') out.push({code: OTHER_REASON_CODE, text: t});
    return out;
  };
  const refreshSave = (): void => {
    const ok = reasonsValid(collect());
    save.classList.toggle('ui-btn-ok', ok);
    save.style.opacity = ok ? '' : '0.5';
  };
  const addBox = (def: ReasonDef, checked: boolean): void => {
    const input = ui.input.bool(def.label, {
      value: checked, tooltipText: `Set "${def.text}" as a reason`, onValueChanged: () => refreshSave(),
    });
    boxes.push({def, input});
  };
  for (const d of crit) addBox(d, true);
  for (const d of judge) addBox(d, false);
  other.onChanged.subscribe(() => refreshSave());
  const cancel = ui.button('Cancel', () => {
    ui.empty(host);
    onCancel();
  });
  const form = ui.form([...boxes.map((b) => b.input), other]);
  form.prepend(ui.divText(REASON_TITLE[action], {style: {fontWeight: '600'}}));
  ui.forms.addButtons(form, [save, cancel]);
  // Right-aligned inside the form's own button row (user-directed 2026-09-30).
  const btnRow = form.querySelector<HTMLElement>('.ui-input-buttons > .ui-input-editor');
  if (btnRow !== null) btnRow.style.justifyContent = 'flex-end';
  host.append(form);
  refreshSave();
}

function appendNext(body: HTMLElement, args: ReviewPaneArgs): void {
  const n = args.pendingOrder.length;
  if (n === 0) return; // nothing left: the strip hint says so (user-directed 2026-09-30)
  const k = args.pendingOrder.indexOf(args.unit.key) + 1;
  const label = k > 0 ? `Next to review → (${k} of ${n})` : `Next to review → (${n} left)`;
  // A link (user-directed 2026-10-12, reversing the 2026-09-30 button: with the decisions now a
  // radio there is no button row for it to match). Wrapped so it keeps its own width.
  const next = ui.link(label, () => args.onNext(), 'Go to the next profile still to review');
  next.classList.add('nca-review-next');
  const row = ui.divH([next]);
  body.appendChild(row);
  // Left-aligned with the review choice's radio circles (user-directed 2026-10-12). The circles'
  // indent is DG's own input padding, so it is measured, not restated in CSS.
  ui.onSizeChanged(body).subscribe(() => {
    const circle = body.querySelector<HTMLElement>('.nca-review-choice input[type="radio"]');
    if (circle === null) return; // no choice (e.g. "Not needed") — nothing to align to
    const dx = circle.getBoundingClientRect().left - row.getBoundingClientRect().left;
    next.style.marginLeft = `${Math.max(0, dx)}px`;
  });
}
