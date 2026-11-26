import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {labMarkText, perGroupText, planCallOf, planView} from '../plan-view';
import {LAB_SENSITIVE} from '../engine';
import {emptyMeta} from '../study-meta';
import {Arm, armName, armsOf} from '../arms-input';
import {PLAN_TAG, parsePlan} from '../plan-file';
import {openStage} from '../stages';
import {exampleStudy} from './example-study';

// the animals per group and the end of study are no longer form fields: set from their cards, they decide how many
// animals a study uses, so these tests click the cards and read the call the next plan takes
const arms = [{name: 'Cpd-A', therapy: 'small_molecule', metric: 'TGI', effect: 60}];

// the calendar hint counts the plans a user opens: these tests leave the tester's count as it was
const HINT = 'calendarHintViews';
let hintCount: string | null | undefined;

async function opened(inputs: {[param: string]: unknown}): Promise<DG.ViewBase> {
  hintCount = grok.userSettings.getValue('PmxStudio', HINT);
  const view = await planView({name: 'Card test', study: emptyMeta(),
    inputs: {line: 'MC38', endpoint: 1000, arms, power: 80, correction: 'holm', ...inputs}});
  grok.shell.addView(view);
  // the first plan: the cards carry its figures
  await awaitCheck(() => /\d/.test(card(view, 'Analyzable')?.textContent ?? ''), 'no plan in 60 s', 60000);
  return view;
}

const card = (view: DG.ViewBase, label: string) => [...view.root.querySelectorAll<HTMLElement>('.pmx-stat')]
  .find((c) => c.textContent!.includes(label));
const label = (view: DG.ViewBase, name: string) => card(view, name)!.textContent!.includes('· set');
const param = (view: DG.ViewBase, name: string) => planCallOf(view)!.inputs[name];

/** Clicks a card and sets `value` in the dialog it opens. */
async function setOn(view: DG.ViewBase, name: string, value: number): Promise<void> {
  card(view, name)!.click();
  await awaitCheck(() => !!document.querySelector('.u2-dialog input'), 'no dialog', 3000);
  const box = document.querySelector<HTMLInputElement>('.u2-dialog input')!;
  box.value = String(value);
  box.dispatchEvent(new Event('input', {bubbles: true}));
  box.dispatchEvent(new Event('change', {bubbles: true}));
  await delay(100);
  [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog-buttons button')]
    .find((b) => b.textContent!.trim().toLowerCase() === 'set')!.click();
}

/** Clicks a set card and picks the way back in the menu it opens. */
async function backOn(view: DG.ViewBase, name: string, item: string): Promise<void> {
  card(view, name)!.click();
  await awaitCheck(() => [...document.querySelectorAll('.u2-menu-item')].some((e) => e.textContent!.includes(item)),
    'no menu', 3000);
  [...document.querySelectorAll<HTMLElement>('.u2-menu-item')].find((e) => e.textContent!.includes(item))!.click();
}

/** Closes what a test left open (a dialog, a menu, when it failed half way) and the view; the hint count back. */
function closed(view: DG.ViewBase): void {
  document.querySelectorAll<HTMLButtonElement>('.u2-dialog-buttons button').forEach((b) => {
    if (b.textContent!.trim().toLowerCase() === 'cancel')
      b.click();
  });
  document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
  view.close();
  if (hintCount)
    grok.userSettings.add('PmxStudio', HINT, hintCount);
  else
    grok.userSettings.delete('PmxStudio', HINT);
  void grok.userSettings.flush();
}

category('Plan view', () => {
  test('the analyzable card says the animals per group, a range when the arms differ', async () => {
    expect(perGroupText([13, 13, 13, 13]), '13 per group');
    expect(perGroupText([14, 13, 12, 14]), '12–14 per group', '53 animals in 4 groups');
    expect(perGroupText([14, null]), null, 'a group no size reaches');
    expect(perGroupText([]), null);
  });

  test('the lab mark: on the analyzable card of Hepa1-6, PDX-NSCLC and PDX-PDAC only', async () => {
    expect([...LAB_SENSITIVE].sort().join(), 'Hepa1-6,PDX-NSCLC,PDX-PDAC');
    const none = ': it may need more animals.';
    expect(labMarkText(9, 6).endsWith('it would be 9 per group.'), true);
    expect(labMarkText(null, 6).endsWith('it would be over 40 per group.'), true);
    // the guard: a number only when above the card's own
    for (const [n, central] of [[6, 6], [5, 6], [null, null], [9, null]] as const)
      expect(labMarkText(n, central).endsWith(none), true, `${n} vs ${central}`);
    expect(labMarkText('failed', 6).endsWith(none), true);
    const shown = (view: DG.ViewBase) => {
      const mark = card(view, 'Analyzable')!.querySelector<HTMLElement>('.pmx-lab-mark');
      return !!mark && mark.style.display !== 'none';
    };
    for (const [line, want] of [['Hepa1-6', true], ['MC38', false]] as const) {
      const view = await opened({line});
      try {
        expect(shown(view), want, line);
      } finally {
        closed(view);
      }
    }
  }, {timeout: 240000});

  test('the animals per group and the end of study, set from their cards and back to the computed ones', async () => {
    const view = await opened({});
    try {
      expect(label(view, 'Analyzable'), false);
      await setOn(view, 'Analyzable', 12);
      await awaitCheck(() => label(view, 'Analyzable'), 'the card does not say set', 60000);
      expect(param(view, 'n'), 12);
      await backOn(view, 'Analyzable', 'Back to what the power needs');
      await awaitCheck(() => !label(view, 'Analyzable'), 'the card still says set', 60000);
      // empty, not 0: the target power sizes the groups again
      expect(param(view, 'n') == null, true, `n is ${param(view, 'n')}`);

      await setOn(view, 'End of study', 18);
      await awaitCheck(() => label(view, 'End of study'), 'the end card does not say set', 60000);
      expect(param(view, 'lastDay'), 18);
      expect(card(view, 'End of study')!.textContent!.includes('Day 18'), true);
      await backOn(view, 'End of study', 'Back to when the controls reach the endpoint');
      await awaitCheck(() => !label(view, 'End of study'), 'the end card still says set', 60000);
      expect(param(view, 'lastDay') == null, true, `lastDay is ${param(view, 'lastDay')}`);
    } finally {
      closed(view);
    }
  }, {timeout: 240000});

  test('an arm renamed in its card keeps what was typed in its details', async () => {
    hintCount = grok.userSettings.getValue('PmxStudio', HINT);
    // the arm's details: the route typed (its name gives none)
    const view = await planView({name: 'Rename test', study: {...emptyMeta(), groups: [{group: 'Cpd-A',
      compound: 'Cpd-A', dose: null, unit: '', route: 'i.v.', schedule: ''}]},
    inputs: {line: 'MC38', endpoint: 1000, arms, power: 80, correction: 'holm'}});
    grok.shell.addView(view);
    try {
      await awaitCheck(() => !!view.root.querySelector('.pmx-arm-title textarea'), 'no arm', 60000);
      const area = view.root.querySelector<HTMLTextAreaElement>('.pmx-arm-title textarea')!;
      area.value = 'Cpd-A 30 mg/kg QD';
      area.dispatchEvent(new Event('input', {bubbles: true}));
      area.dispatchEvent(new Event('change', {bubbles: true}));
      // the line under the name reads the details kept under the new name
      await awaitCheck(() => (view.root.querySelector('.pmx-arm')?.textContent ?? '')
        .includes('Cpd-A · 30 mg/kg · QD · i.v. (route from the group details)'), 'the details did not follow', 10000);
    } finally {
      closed(view);
    }
  }, {timeout: 120000});

  test('a saved plan with the animals and the end set opens with both cards saying so', async () => {
    const view = await opened({n: 6, lastDay: 18});
    try {
      await awaitCheck(() => label(view, 'Analyzable') && label(view, 'End of study'), 'not shown as set', 60000);
      expect(param(view, 'n'), 6);
      expect(param(view, 'lastDay'), 18);
    } finally {
      closed(view);
    }
  }, {timeout: 120000});

  test('Plan from a study\'s stage bar opens with the study\'s plan, and comes back to it', async () => {
    const text = (await exampleStudy()).getTag(PLAN_TAG)!;
    const want = (parsePlan(text).inputs.arms as Arm[]).map(armName).join();
    const before = new Set(grok.shell.views);
    const added = () => [...grok.shell.views].filter((v) => !before.has(v));
    await openStage('plan', text);
    const view = added()[0];
    try {
      expect(!!view, true, 'no Plan view');
      await awaitCheck(() => armsOf(param(view, 'arms') as DG.DataFrame, false).map(armName).join() === want,
        'the study\'s arms not loaded', 30000);
      await openStage('plan', text);
      expect(added().length, 1);
    } finally {
      view?.close();
    }
  }, {timeout: 60000});
});
