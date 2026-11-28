import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {endpointProblem, labMarkText, perGroupText, planCallOf, planFileOf, planSaveOf, planView, plansMadeCount,
  reviewDays} from '../plan-view';
import {historyProblem} from '../history-dialog';
import {inWords} from '../ui';
import {OPENED_AS, linkRow, openStage, pathOfView} from '../stages';
import {cruxCallsOf} from '../crux';
import {PmxDay} from '../pmx-objects';
import {LAB_SENSITIVE} from '../engine';
import {emptyMeta} from '../study-meta';
import {Arm, armName, armsOf, convertedEffect, saidAs} from '../arms-input';
import {PLAN_TAG, framesOf, parsePlan, shownOf} from '../plan-file';
import {newStudyTable} from '../study-file';
import {exampleStudy} from './example-study';
import {linkedRow, memoryTables, open, save} from '../store';
import {tumorModels} from '../priors';
import {registerHistory} from '../history-models';

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
  // a plan the test changed asks first: its changes discarded, so no "Leave site?" stays armed for the next tests
  // (the platform's dialog, not u2's)
  [...document.querySelectorAll<HTMLElement>('.d4-dialog button, .d4-dialog .ui-btn')]
    .find((b) => b.textContent!.trim().toLowerCase() === 'discard')?.click();
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

  test('Save pressed while the plan is updated waits for it, says so, then saves', async () => {
    const view = await opened({});
    // (File › Save presses the view's Save and reads its words)
    const saveLabel = () => planSaveOf(view)!.querySelector<HTMLElement>('.u2-btn-label') ?? planSaveOf(view);
    const nameDialog = () => [...document.querySelectorAll('.u2-dialog')].some((d) => d.textContent!.includes('Name'));
    try {
      // (a size no other test plans, so not a plan the engine kept; Save pressed the moment the card's change starts
      // the plan, the next task)
      await setOn(view, 'Analyzable', 11);
      await delay(0);
      saveLabel()!.closest('button')!.click();
      expect(saveLabel()!.textContent, 'Saving when the plan is ready…');
      expect(nameDialog(), false, 'saved before the plan was ready');
      // the plan made: the save goes on (a plan never saved asks its name first)
      await awaitCheck(nameDialog, 'not saved once the plan was ready', 60000);
      expect(label(view, 'Analyzable'), true);
      expect(saveLabel()!.textContent, 'Save plan');
    } finally {
      closed(view);
    }
  }, {timeout: 120000});

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

  test('each arm’s swatch is its colour in the charts, after an arm is added and after one is renamed', async () => {
    const view = await opened({arms: [...arms, {name: 'Cpd-B', therapy: 'small_molecule', metric: 'TGI', effect: 50}]});
    // the charts' colour of a group: the result's swatches (colorsIn, the one source the charts are painted from)
    const charted = (name: string) => [...view.root.querySelectorAll<HTMLElement>('span')].find((s) =>
      s.children.length === 2 && s.children[1].textContent === name && !s.closest('.pmx-arm'))
      ?.querySelector<HTMLElement>('span')?.style.background;
    const cards = () => [...view.root.querySelectorAll<HTMLElement>('.pmx-arm')].map((c) => ({
      name: c.querySelector<HTMLTextAreaElement>('.pmx-arm-title textarea')!.value,
      color: c.querySelector<HTMLElement>('[data-pmx=arm-swatch]')!.style.background}));
    // (an arm added has no name yet, so no chart: its own colour, none of the others')
    const same = () => cards().every((c, i, all) => !!c.color && (c.name ? c.color === charted(c.name) :
      all.every((o, j) => j === i || o.color !== c.color)));
    try {
      expect(cards().length, 2);
      expect(same(), true, JSON.stringify(cards()));
      expect(cards()[0].color !== cards()[1].color, true, 'two arms, two colours');
      view.root.querySelector<HTMLElement>('[data-u2=pmx-add-arm]')!.click();
      await awaitCheck(() => cards().length === 3 && same(), `added: ${JSON.stringify(cards())}`, 60000);
      const area = view.root.querySelector<HTMLTextAreaElement>('.pmx-arm-title textarea')!;
      area.value = 'Cpd-C';
      area.dispatchEvent(new Event('input', {bubbles: true}));
      area.dispatchEvent(new Event('change', {bubbles: true}));
      await awaitCheck(() => cards()[0].name === 'Cpd-C' && !!charted('Cpd-C') && same(),
        `renamed: ${JSON.stringify(cards())}`, 60000);
    } finally {
      closed(view);
    }
  }, {timeout: 180000});

  test('a calendar day’s panel says it is measured exactly when its cell does, day 0 and the end too', async () => {
    const view = await opened({});
    const before = grok.shell.o;
    try {
      await awaitCheck(() => !!view.root.querySelector('.pmx-calendar-box [aria-label^="Day 0"]'), 'no calendar',
        60000);
      const cells = [...view.root.querySelectorAll<HTMLElement>('.pmx-calendar-box [aria-label^="Day "]')];
      expect(cells.length > 7, true);
      for (const cell of cells) {
        const day = Number(/^Day ([0-9]+)/.exec(cell.getAttribute('aria-label')!)![1]);
        const title = new RegExp(`^Day ${day}(?![0-9])`);
        const shown = () => grok.shell.o instanceof PmxDay ? grok.shell.o.details().textContent ?? '' : '';
        // (the platform takes no new current object for a moment after the last one: clicked again until it does)
        for (let i = 0; i < 10 && !title.test(shown()); i++) {
          cell.click();
          await delay(600);
        }
        expect(title.test(shown()), true, `day ${day}: no panel`);
        const text = shown();
        const marked = [...cell.children].some((c) => c.textContent === 'Measure');
        expect(/MeasureYes/.test(text), marked, `day ${day}: ${text}`);
      }
    } finally {
      grok.shell.o = before;
      closed(view);
    }
  }, {timeout: 240000});

  // what the cards say, in order: the numbers a plan shows
  test('a saved plan reopens on its saved tables, as computed; one saved without them is computed', async () => {
    const cards = (v: DG.ViewBase) => [...v.root.querySelectorAll('.pmx-stat')].map((c) => c.textContent).join(' | ');
    const view = await opened({});
    const views: DG.ViewBase[] = [view];
    try {
      // (saved once its schedules are compared, as a user saves it: the file keeps them)
      await awaitCheck(() => ![...view.root.querySelectorAll('*')].some((e) => e.children.length === 0 &&
        e.textContent === 'Comparing measurement schedules…'), 'schedules not compared', 120000);
      const fresh = cards(view);
      const file = parsePlan(JSON.stringify(await planFileOf(view, 'Reopen test')));
      expect(!!file.predicted.shown, true, 'no tables saved');
      // the same plan again: not designed again to show it (its charts' bands are drawn), the same numbers
      const designs = () => cruxCallsOf('stats.tumor_study_design');
      const before = designs();
      const again = await planView(file);
      views.push(again);
      grok.shell.addView(again);
      await awaitCheck(() => /[0-9]/.test(card(again, 'Analyzable')?.textContent ?? ''), 'not shown', 10000);
      expect(cards(again), fresh);
      expect(designs() - before, 0, 'the saved plan was designed again');
      // a file from before the tables were kept: computed as it opens, the same numbers
      const old = JSON.parse(JSON.stringify(file));
      delete old.predicted.shown;
      const computed = await planView(parsePlan(JSON.stringify(old)));
      views.push(computed);
      grok.shell.addView(computed);
      await awaitCheck(() => /[0-9]/.test(card(computed, 'Analyzable')?.textContent ?? ''), 'old file not shown',
        60000);
      expect(cards(computed), fresh);
      // tables computed by another engine, app or atlas: computed again, and no save while that runs (it would keep
      // the old tables as this engine's)
      const other = parsePlan(JSON.stringify(file));
      other.predicted.shown!.engine = 'another engine';
      const mid = plansMadeCount();
      const redone = await planView(other);
      views.push(redone);
      grok.shell.addView(redone);
      let refused = false;
      await planFileOf(redone, 'Too early').catch(() => refused = true);
      expect(refused, true, 'saved while computing');
      await awaitCheck(() => plansMadeCount() > mid && /[0-9]/.test(card(redone, 'Analyzable')?.textContent ?? ''),
        'not computed again', 60000);
      expect(cards(redone), fresh);
      // tables that lack what Plan reads first: computed
      const f = framesOf(file.predicted.shown)!;
      f.plan.columns.remove('study_end');
      const s = file.predicted.shown!;
      expect(framesOf(shownOf(f, s.asked, s.declared, s.engine)), null);
      // a study followed from it carries the plan without its tables
      const text = JSON.stringify(file);
      const tag = newStudyTable(file).getTag(PLAN_TAG) ?? '';
      expect(tag.includes('"shown"'), false);
      expect(JSON.stringify(parsePlan(tag).inputs), JSON.stringify(file.inputs));
      console.log(`PMX Studio: demo .pmxplan ${text.length} characters with its tables, ${tag.length} without`);
    } finally {
      views.reverse().forEach(closed);
    }
  }, {timeout: 240000});

  test('a plan of every model saves as a row and reads back, its controls reaching the endpoint or not', async () => {
    // (2000 mm³: slow models, Panc02 among them, never reach it in the study, and their end has no median; a model
    // of the user's own history, synthetic, too)
    const rows = [0, 3, 7, 10, 14, 17, 21].flatMap((day) => [0, 1, 2, 3].map((a) =>
      ({animal: `A${a}`, day, volume: 120 * Math.exp((0.08 + 0.005 * a) * day)})));
    const own = registerHistory({name: 'PDX-SAVE-TEST', rows, take: null, source: {sameLab: 'yes', sameHost: 'yes',
      passage: 'later', frozen: false, meansOnly: false, animalsPerGroup: null}});
    const t = memoryTables();
    const failed: string[] = [];
    for (const line of [...Object.keys(await tumorModels()), own]) {
      const view = await opened({line, endpoint: 2000});
      try {
        // (once the page is done with the plan: a save waits for it)
        let plan = null;
        for (let i = 0; i < 120 && !plan; i++)
          plan = await planFileOf(view, `${line} plan`).catch(() => delay(500).then(() => null));
        const back = await open('plan', (await save('plan', plan!, null, {}, t)).id, t);
        expect(parsePlan(JSON.stringify(back.doc)).predicted.line, plan!.predicted.line);
        if (line === own)
          expect(back.doc.history?.name, 'PDX-SAVE-TEST');
      } catch (e) {
        failed.push(`${line}: ${(e as Error)?.message ?? e}`);
      } finally {
        closed(view);
      }
    }
    expect(failed.join('; '), '');
  }, {timeout: 600000});

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

category('Plan checks (QA 0.1.3)', () => {
  test('the endpoint must be above the volume at randomization', async () => {
    expect(endpointProblem(500, 500) !== null, true, 'equal to V0');
    expect(endpointProblem(400, 500) !== null, true, 'below V0');
    expect(endpointProblem(1000, 100), null, 'above V0');
    expect(endpointProblem(null, 100), null, 'empty: the required check says so');
  });

  test('planned reviews: whole study days, once each', async () => {
    expect(JSON.stringify(reviewDays('14, 7')), JSON.stringify({days: [7, 14]}), 'sorted');
    expect(JSON.stringify(reviewDays('')), JSON.stringify({days: []}), 'none');
    expect(JSON.stringify(reviewDays(null)), JSON.stringify({days: []}), 'no text');
    for (const bad of ['7, 7', '-1', '10.5', 'abc', '0', '121'])
      expect('problem' in reviewDays(bad), true, bad);
  });

  test('TGI and T/C convert to the same treated volume, both ways', async () => {
    // V0 100, endpoint 1000: TGI 60 is a treated 460 mm³, T/C 46
    expect(convertedEffect(60, 'TC', 100, 1000), 46);
    expect(convertedEffect(46, 'TGI', 100, 1000), 60);
    expect(convertedEffect(95, 'TC', 100, 1000), 15);
    // a T/C under V0/C is a regression: TGI past 100%, up to just below complete regression (111.1% here: 111 is below)
    expect(convertedEffect(5, 'TGI', 100, 1000), 106);
    expect(convertedEffect(0, 'TGI', 100, 1000), 111, 'below complete regression');
    expect(convertedEffect(110, 'TC', 100, 1000), 1, 'TGI 110 is T/C 1 (no 5-95% clamp since 1527dbe6ec)');
    // said back and forth unchanged: the value entered, never a rounding of a rounding (52 drifted to 53)
    let arm: {metric: 'TGI' | 'TC', effect: number} = {metric: 'TGI', effect: 52};
    for (let k = 0; k < 3; k++) {
      const [tc, was] = saidAs(arm, 'TC', 100, 300);
      const [tgi] = saidAs({metric: 'TC', effect: tc}, 'TGI', 100, 300, was);
      expect(tgi, 52, `round ${k + 1}`);
      arm = {metric: 'TGI', effect: tgi};
    }
    // changed while said the other way: converted from the new value
    const [, was] = saidAs({metric: 'TGI', effect: 52}, 'TC', 100, 1000);
    expect(saidAs({metric: 'TC', effect: 40}, 'TGI', 100, 1000, was)[0], convertedEffect(40, 'TGI', 100, 1000));
  });

  test('an own model is checked in its dialog: growth, an atlas name, the take', async () => {
    const rows = [0, 3, 7, 10].map((day) => ({animal: 'A', day, volume: 100 * Math.exp(0.1 * day)}));
    expect(historyProblem(rows, 'PDX-117', ['MC38'], null), null, 'a good history');
    expect(historyProblem(rows, 'mc38', ['MC38'], null) !== null, true, 'an atlas model\'s name');
    expect(historyProblem(rows, 'PDX-117', [], 150) !== null, true, 'a take over 100%');
    expect(historyProblem(rows.slice(0, 2), 'PDX-117', [], null) !== null, true, 'no growth to read');
  });

  test('engine messages in words', async () => {
    expect(inWords('Arm A: stats.tumor_study_design: the endpoint volume must be above the volume at ' +
      'randomisation (> 0)'), 'Arm A: the endpoint volume must be above the volume at randomization');
  });

  test('a view showing a saved row names it in its URL, and lets go when it shows another', async () => {
    const root = document.createElement('div');
    const view = {root, path: ''} as unknown as DG.ViewBase;
    root.dataset.pmxStage = 'plan';
    linkRow(view, 'plan', 'abc');
    expect(view.path, '/plan?plan=abc', 'the app\'s own view: relative');
    expect(root.dataset[OPENED_AS], 'abc', 'the Overview row brings it up');
    root.dataset.pmxAdded = 'true';
    expect(pathOfView(view), '/apps/PmxStudio/plan?plan=abc', 'a view the app added: whole');
    linkRow(view, 'plan', null);
    expect(view.path, '/apps/PmxStudio/plan', 'a table not saved: no row');
    expect(root.dataset[OPENED_AS] === undefined, true, 'let go');
    expect(linkedRow('?plan=abc') === null || linkedRow('?plan=abc')!.id === 'abc', true, 'read back on reload');
  });
});
