/**
 * nca-workflow-stages — the QA-round UI refinements (manual QA on local, 2026-09-29 / 2026-09-30).
 *
 * Pins what the manual QA settled, headless (no view, no storage):
 *  - `formatNextHint` — ONE message at a time; `In study '<label>'` only when nothing is left to do;
 *  - the stage strip — the Review hint is the action link, the study NAME is a link (and never
 *    triggers the review), Report carries ✓ once the analysis is in a study;
 *  - the Review pane — the untitled radio choice (Adjust λz points — dense profiles only, the default —
 *    then the decisions; hidden once a decision is saved, 2026-10-12), the reason FORM a decision
 *    opens, the Reset icon, the `Next to review` link (a button until 2026-10-12), no "Nothing left" line;
 *  - reason labels — short picker labels, the FULL text is what a decision persists.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {formatNextHint} from '../adapter/workflow-stage';
import type {StageResult} from '../adapter/workflow-stage';
import {judgementReasons, criterionReasons} from '../adapter/review-reasons';
import {DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED} from '../adapter/defaults';
import type {AcceptanceThresholds} from '../adapter/acceptance';
import type {ReviewUnit} from '../adapter/review-units';
import type {ReviewAction, ReviewReason} from '../adapter/study-schema';
import {createWorkflowStrip} from '../views/workflow-strip';
import {renderReviewPane} from '../views/review-pane';
import {excludedCellTipSpec} from '../views/results-view';
import type {ReviewPaneArgs} from '../views/review-pane';

const T: AcceptanceThresholds = {
  lambdaAcceptanceRSquared: DEFAULT_ACCEPTANCE_RSQUARED,
  lambdaMinPoints: DEFAULT_RULES.lambdaZ.minPoints,
  extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
};

const st = (stage: StageResult['stage'], pending = 0, setupChanged = false): StageResult =>
  ({stage, pending, notApplicable: 0, setupChanged});

category('Workflow stages: next-step hint (QA round)', () => {
  test('one message at a time — the next action alone, never prefixed by the study', async () => {
    expect(formatNextHint(st('setup'), null, 'setup'), 'Complete the setup and click Confirm');
    expect(formatNextHint(st('setup', 0, true), 'S', 'setup'), 'Apply changes');
    expect(formatNextHint(st('review', 1), 'S', 'review'), 'Review 1 profile');
    expect(formatNextHint(st('review', 5), null, 'report'), 'Review 5 profiles');
  });
  test('Complete: revisit guidance off Report; the study name only when nothing is left to do', async () => {
    expect(formatNextHint(st('complete'), null, 'review'),
      'All profiles reviewed — go to Report, or click a row or its Review cell to revisit a profile');
    expect(formatNextHint(st('complete'), 'S', 'setup'),
      'All profiles reviewed — go to Report, or click a row or its Review cell to revisit a profile');
    expect(formatNextHint(st('complete'), null, 'report'),
      'Group the data to see aggregated results. Add the results to a study.');
    expect(formatNextHint(st('complete'), 'Rat PK', 'report'), 'In study \'Rat PK\'');
  });
});

category('Workflow stages: stage strip (QA round)', () => {
  const setup = (): {
    strip: ReturnType<typeof createWorkflowStrip>; calls: string[];
  } => {
    const calls: string[] = [];
    const strip = createWorkflowStrip((s) => calls.push(`nav:${s}`), () => calls.push('review'),
      () => calls.push('study'));
    return {strip, calls};
  };
  const hintOf = (root: HTMLElement): HTMLElement => root.querySelector('.nca-stage-hint') as HTMLElement;
  const items = (root: HTMLElement): string[] =>
    Array.from(root.querySelectorAll('.nca-stage-item')).map((e) => e.textContent ?? '');

  test('BUG-17: the strip keeps its height when the stage marks change (no line-height: normal)', async () => {
    const {strip} = setup();
    document.body.appendChild(strip.root); // mounted: computed style + layout need a live node
    try {
      strip.update({stage: st('setup'), viewStage: 'setup', studyLabel: null});
      const chain = strip.root.querySelector('.nca-stage-chain') as HTMLElement;
      // Under `normal` the ✓ / ● fallback-font metrics size the line, so the strip jumped at Confirm.
      expect(getComputedStyle(chain).lineHeight !== 'normal', true, 'the chain line-height is fixed');
      const before = strip.root.offsetHeight;
      strip.update({stage: st('review', 3), viewStage: 'review', studyLabel: null}); // ✓ Setup · ● Review
      expect(strip.root.offsetHeight, before, 'Setup → Review: same strip height');
    } finally {
      strip.root.remove();
    }
  });
  test('Setup: Review is a link exactly while Confirm can run, and its click navigates (= Confirm)', async () => {
    const {strip, calls} = setup();
    const review = (): HTMLElement => strip.root.querySelectorAll<HTMLElement>('.nca-stage-item')[1];
    strip.update({stage: st('setup'), viewStage: 'setup', studyLabel: null, confirmable: false});
    expect(review().classList.contains('nca-stage-link'), false, 'Confirm dimmed → Review not clickable');
    strip.update({stage: st('setup'), viewStage: 'setup', studyLabel: null, confirmable: true});
    expect(review().classList.contains('nca-stage-link'), true, 'Confirm ready → Review is a link');
    expect(strip.root.querySelectorAll('.nca-stage-item')[2].classList.contains('nca-stage-link'), false,
      'Report still needs results');
    review().click();
    expect(calls.join(), 'nav:review');
  });
  test('Review with profiles left: the whole hint is the action link and starts the review', async () => {
    const {strip, calls} = setup();
    strip.update({stage: st('review', 5), viewStage: 'report', studyLabel: null});
    const hint = hintOf(strip.root);
    expect(hint.textContent, 'Review 5 profiles →');
    expect(hint.classList.contains('nca-stage-hint-action'), true);
    hint.click();
    expect(calls.join(), 'review');
  });
  test('a Setup edit since Confirm: plain text, no action', async () => {
    const {strip, calls} = setup();
    strip.update({stage: st('setup', 0, true), viewStage: 'setup', studyLabel: null});
    const hint = hintOf(strip.root);
    expect(hint.textContent, 'Apply changes');
    expect(hint.classList.contains('nca-stage-hint-action'), false);
    hint.click();
    expect(calls.length, 0);
  });
  test('Complete on Report outside a study: plain text, no link element', async () => {
    const {strip} = setup();
    strip.update({stage: st('complete'), viewStage: 'report', studyLabel: null});
    const hint = hintOf(strip.root);
    expect(hint.querySelector('.nca-stage-hint-link'), null);
    expect(hint.classList.contains('nca-stage-hint-action'), false);
    expect(items(strip.root)[2].startsWith('✓'), false, 'Report is not done before it is in a study');
  });
  test('in a study: Report ✓, the study NAME is the link and opens the study (not the review)', async () => {
    const {strip, calls} = setup();
    strip.update({stage: st('complete'), viewStage: 'report', studyLabel: 'Rat PK'});
    expect(items(strip.root)[2].startsWith('✓'), true);
    const hint = hintOf(strip.root);
    expect(hint.textContent, 'In study Rat PK');
    const link = hint.querySelector('.nca-stage-hint-link') as HTMLElement;
    expect(link.textContent, 'Rat PK');
    link.click();
    expect(calls.join(), 'study');
  });
  test('in a study with profiles left: only the review action shows, and a Setup edit drops Report ✓', async () => {
    const {strip} = setup();
    strip.update({stage: st('review', 1), viewStage: 'review', studyLabel: 'Rat PK'});
    expect(hintOf(strip.root).textContent, 'Review 1 profile →');
    strip.update({stage: st('setup', 0, true), viewStage: 'setup', studyLabel: 'Rat PK'});
    expect(items(strip.root)[2].startsWith('✓'), false);
  });
});

/** A synthetic review unit — only the fields the pane reads. */
function unit(over: Partial<ReviewUnit>): ReviewUnit {
  return {
    key: '1', kind: 'profile', need: 'review', category: 'flagged', decision: null, staleDecision: null,
    allowedActions: ['accept', 'exclude'], flagCodes: ['SPAN_LT', 'NPTS_LE'], blqOnly: false, blqTimes: [],
    ...over,
  } as ReviewUnit;
}

interface PaneCalls {
  decide: {keys: ReadonlyArray<string>; action: ReviewAction; reasons: ReadonlyArray<ReviewReason>}[];
  reset: string[];
  next: number;
}

function renderPane(u: ReviewUnit, over: Partial<ReviewPaneArgs> = {}): {host: HTMLElement; calls: PaneCalls} {
  const calls: PaneCalls = {decide: [], reset: [], next: 0};
  const host = document.createElement('div');
  document.body.appendChild(host); // mounted: ui.form sizes its labels on a live node
  renderReviewPane(host, {
    unit: u, selected: [], pendingOrder: ['1', '2'], thresholds: T, editable: true,
    onDecide: (keys, action, reasons) => calls.decide.push({keys, action, reasons}),
    onReset: (k) => calls.reset.push(k),
    onNext: () => calls.next++,
    ...over,
  });
  return {host, calls};
}

const ADJUST = 'Adjust λz points in the λz editor';
const buttonByText = (host: HTMLElement, text: string): HTMLElement | undefined =>
  Array.from(host.querySelectorAll('button')).find((b) => (b.textContent ?? '').trim() === text) as
    HTMLElement | undefined;
/** The review choice's option labels, in order. */
const choiceLabels = (host: HTMLElement): string[] =>
  Array.from(host.querySelectorAll('.nca-review-choice .ui-input-editor label'))
    .map((l) => (l.textContent ?? '').trim());
/** The checked radio's label text, or null. */
const checkedChoice = (host: HTMLElement): string | null => {
  const labels = Array.from(host.querySelectorAll<HTMLElement>('.nca-review-choice .ui-input-editor label'));
  const radios = Array.from(host.querySelectorAll<HTMLInputElement>('.nca-review-choice input[type="radio"]'));
  const i = radios.findIndex((r) => r.checked);
  return i < 0 ? null : (labels[i]?.textContent ?? '').trim();
};
/** Pick a review-choice option, as a click on its label does. */
const pick = (host: HTMLElement, text: string): void => {
  const label = Array.from(host.querySelectorAll<HTMLElement>('.nca-review-choice .ui-input-editor label'))
    .find((l) => (l.textContent ?? '').trim() === text);
  label!.click();
};

category('Workflow stages: Review pane (QA round)', () => {
  test('the choice: Adjust first and selected (dense Flagged / Not estimated only), then the decisions', async () => {
    const opts = (u: ReviewUnit): {labels: string[]; checked: string | null} => {
      const {host} = renderPane(u);
      const out = {labels: choiceLabels(host), checked: checkedChoice(host)};
      host.remove();
      return out;
    };
    const flagged = opts(unit({}));
    expect(flagged.labels.join('|'), [ADJUST, 'Accept as is', 'Exclude λz'].join('|'), 'flagged profile');
    expect(flagged.checked, ADJUST, 'Adjust is the default');
    expect(opts(unit({category: 'not-estimated', allowedActions: ['acknowledge'], flagCodes: []})).labels[0], ADJUST);
    const noExp = opts(unit({category: 'no-exposure', allowedActions: ['acknowledge'], flagCodes: []}));
    expect(noExp.labels.includes(ADJUST), false, 'no points to fit');
    expect(noExp.checked === null, true, 'nothing pre-selected without Adjust');
    expect(opts(unit({kind: 'composite'})).labels.includes(ADJUST), false, 'a composite is not edited');
  });
  test('a Flagged profile has no "At least one criterion" line — its reasons head the readout', async () => {
    const {host} = renderPane(unit({}));
    expect((host.textContent ?? '').includes('At least one'), false);
    host.remove();
  });
  test('a decision opens ONE form titled by the action; the choice stays, Cancel returns it to Adjust', async () => {
    const {host} = renderPane(unit({}));
    pick(host, 'Exclude λz');
    expect(choiceLabels(host).length, 3, 'the choice stays visible while the form is open');
    const form = host.querySelector('.ui-form') as HTMLElement;
    expect(form !== null, true);
    expect((form.textContent ?? '').includes('Exclude reason'), true, 'title inside the form');
    expect(buttonByText(form, 'Save') !== undefined && buttonByText(form, 'Cancel') !== undefined, true,
      'Save / Cancel inside the form');
    const labels = Array.from(form.querySelectorAll('label')).map((l) => (l.textContent ?? '').trim());
    expect(labels.includes('Terminal phase not reached'), true, 'short label');
    expect(labels.includes('Other'), true, 'Other without the ellipsis');
    buttonByText(form, 'Cancel')!.click();
    expect(checkedChoice(host), ADJUST, 'Cancel returns the choice to Adjust');
    expect(host.querySelector('.ui-form'), null);
    pick(host, 'Accept as is');
    pick(host, ADJUST);
    expect(host.querySelector('.ui-form'), null, 'picking Adjust closes the form');
    host.remove();
  });
  test('Save persists the FULL reason text (short labels are display only)', async () => {
    const {host, calls} = renderPane(unit({}));
    pick(host, 'Accept as is');
    buttonByText(host, 'Save')!.click(); // the criterion reasons are pre-checked
    expect(calls.decide.length, 1);
    expect(calls.decide[0].action, 'accept');
    const expected = criterionReasons(['SPAN_LT', 'NPTS_LE'], T).map((r) => r.text).join('|');
    expect(calls.decide[0].reasons.map((r) => r.text).join('|'), expected);
    host.remove();
  });
  test('a decided unit: Reset is an icon beside the title (editable only); Next is a link', async () => {
    const decided = unit({decision: {action: 'exclude', reasons: [{code: 'X', text: 'x'}], user: 'a',
      at: '2026-09-30T12:00:00Z', fingerprint: 'f'}});
    const {host, calls} = renderPane(decided);
    expect(host.querySelector('.nca-review-choice'), null, 'no review choice once a decision is saved');
    const icon = host.querySelector('.fa-undo') as HTMLElement;
    expect(icon !== null, true);
    icon.click();
    expect(calls.reset.join(), '1');
    // (`expect(x, undefined)` would default the expected value to true — compare explicitly)
    expect(buttonByText(host, 'Reset') === undefined, true, 'no Reset button any more');
    const next = host.querySelector<HTMLElement>('a.nca-review-next');
    expect(next !== null && /Next to review/.test(next.textContent ?? ''), true, 'Next to review is a link');
    expect(Array.from(host.querySelectorAll('button')).some((b) => /Next/.test(b.textContent ?? '')), false,
      'no Next button any more');
    next!.click();
    expect(calls.next, 1);
    host.remove();
    const ro = renderPane(decided, {editable: false});
    expect(ro.host.querySelector('.fa-undo'), null, 'read-only: no Reset');
    ro.host.remove();
  });
  test('bulk: only units of the SAME category join (reasons are keyed by action AND category)', async () => {
    const notEst = unit({key: '1', category: 'not-estimated', allowedActions: ['acknowledge'], flagCodes: []});
    const noExp = unit({key: '2', category: 'no-exposure', allowedActions: ['acknowledge'], flagCodes: []});
    const notEst2 = unit({key: '3', category: 'not-estimated', allowedActions: ['acknowledge'], flagCodes: []});
    const mixed = renderPane(notEst, {selected: [noExp]});
    expect((mixed.host.textContent ?? '').includes('Decide for'), false,
      'no-exposure never joins a not-estimated batch');
    mixed.host.remove();
    const same = renderPane(notEst, {selected: [notEst2]});
    expect((same.host.textContent ?? '').includes('Decide for 2 selected profiles'), true, 'same category → bulk');
    same.host.remove();
  });
  test('the composite excluded-cell tooltip says it is the pooled composite that is withheld', async () => {
    const d = {action: 'exclude' as const, reasons: [{code: 'EXC_ERRATIC', text: 'Non-monotonic'}], user: 'a',
      at: '2026-09-30T12:00:00Z', fingerprint: 'f'};
    const c = excludedCellTipSpec(d, true);
    expect(c.title, 'excluded (λz)');
    expect(c.summary.includes('pooled composite'), true);
    expect((c.detail ?? '').includes('Non-monotonic'), true, 'the reasons ride along');
    expect(excludedCellTipSpec(d).summary.includes('every statistic'), true, 'the dense wording is unchanged');
  });
  test('nothing left to review: no "Nothing left" line and no Next link', async () => {
    const {host} = renderPane(unit({decision: {action: 'accept', reasons: [{code: 'X', text: 'x'}], user: 'a',
      at: '2026-09-30T12:00:00Z', fingerprint: 'f'}}), {pendingOrder: []});
    expect((host.textContent ?? '').includes('Nothing left'), false);
    expect(host.querySelector('.nca-review-next'), null);
    host.remove();
  });
});

category('Workflow stages: reason labels (QA round)', () => {
  test('every judgement reason has a short label distinct from its persisted full text', async () => {
    const pairs: [ReviewAction, 'flagged' | 'not-estimated' | 'no-exposure'][] = [
      ['accept', 'flagged'], ['exclude', 'flagged'], ['acknowledge', 'not-estimated'],
      ['acknowledge', 'no-exposure'],
    ];
    for (const [a, c] of pairs) {
      for (const r of judgementReasons(a, c)) {
        expect(r.label.length > 0 && r.label.length <= 28, true, `${r.code} label is short: "${r.label}"`);
        expect(r.text.length >= r.label.length, true, `${r.code} text is the full reason`);
      }
    }
  });
  test('a criterion reason is already short: its label IS its text', async () => {
    for (const r of criterionReasons(['ADJR2_LT', 'SPAN_LT'], T)) expect(r.label, r.text);
  });
});
