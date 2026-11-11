import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {planCallOf, planView} from '../plan-view';
import {emptyMeta} from '../study-meta';

// the animals per group and the end of study are no longer form fields: set from their cards, they decide how many
// animals a study uses, so these tests click the cards and read the call the next plan takes
const arms = [{name: 'Cpd-A', therapy: 'small_molecule', metric: 'TGI', effect: 60}];

async function opened(inputs: {[param: string]: unknown}): Promise<DG.ViewBase> {
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

category('Plan view', () => {
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
      view.close();
    }
  }, {timeout: 240000});

  test('a saved plan with the animals and the end set opens with both cards saying so', async () => {
    const view = await opened({n: 6, lastDay: 18});
    try {
      await awaitCheck(() => label(view, 'Analyzable') && label(view, 'End of study'), 'not shown as set', 60000);
      expect(param(view, 'n'), 6);
      expect(param(view, 'lastDay'), 18);
    } finally {
      view.close();
    }
  }, {timeout: 120000});
});
