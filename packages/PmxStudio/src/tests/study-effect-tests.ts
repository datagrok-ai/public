import * as DG from 'datagrok-api/dg';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Scope, signal} from '@datagrok-libraries/u2';

import {Arm, ArmsInput, NAME_EXAMPLES, armsFrame, armsOf} from '../arms-input';
import {parseArmName} from '../arm-name';
import {DrugPk, drugPkOf, proportionalOf} from '../exposure';
import {declarationOf} from '../frames';
import {armColors, groupColors} from '../colors';
import {PLAN_TAG, parsePlan} from '../plan-file';
import {Anchor, FromStudy, assumedText, PreviousStudy, anchorsIn, controlOf, effectFromStudy, effectText,
  linkedScenarios, nearestAnchor, noPkText, offerText, plainReason, silence, SIMPLIFIED, sizingEffect, undosedIn,
  undosedText}
  from '../study-effect';
import {emptyMeta} from '../study-meta';
import {isOwnStudy} from '../previous-studies';
import {exampleStudy} from './example-study';

const same = (a: unknown, b: unknown, msg?: string) => expect(JSON.stringify(a), JSON.stringify(b), msg);

/** An invented earlier study: a vehicle and two doses of Cpd-X, 4 animals each, measured on days 0, 7, 14 and 21;
 *  the treated tumors half the control's size on the last day (at 10 mg/kg). `ids`: the animals numbered 1..4 in
 *  every group, as xlsx studies number them. */
function earlier(line: string | null = 'MC38', ids: 'unique' | 'per group' = 'unique'): PreviousStudy {
  const rows: PreviousStudy['rows'] = [];
  for (const [group, f] of [['Vehicle', 1], ['Cpd-X 10 mg/kg i.p. QD', 0.5], ['Cpd-X 3 mg/kg i.p. QD', 0.8]] as const) {
    for (let a = 0; a < 4; a++) {
      for (const day of [0, 7, 14, 21]) {
        rows.push({animal: ids === 'unique' ? `${group}-${a}` : `${a + 1}`, group, day,
          volume: 100 * Math.exp(day / 10) * (day ? f : 1) + a});
      }
    }
  }
  // a fourth day-28 reading of the control only: not a day the effect can be read
  rows.push({animal: ids === 'unique' ? 'Vehicle-0' : '1', group: 'Vehicle', day: 28, volume: 1800});
  return {name: 'Earlier study', line, meta: {...emptyMeta(), start: '2026-01-05'}, rows};
}

/** NCA's PK of Cpd-X: 10 and 30 mg/kg, exposure in step with the dose. */
const cpdxPk: DrugPk = {doses: [10, 30], unit: 'mg/kg', auc: new Map([[10, '1,200 h*ng/mL']]), proportional: 'yes'};

const linked: FromStudy = {study: 'Earlier study', arm: 'Cpd-X 10 mg/kg i.p. QD', dosePrev: 10, doseNew: 20,
  unit: 'mg/kg', point: 65, range: [55, 80]};

category('Effect from an earlier study', () => {
  test('the groups that received the drug alone, read on the last day with 3 animals in each', async () => {
    const found = anchorsIn(earlier(), 'Cpd-X 20 mg/kg i.p. QD');
    same(found.map((a) => a.arm), ['Cpd-X 10 mg/kg i.p. QD', 'Cpd-X 3 mg/kg i.p. QD']);
    expect(found[0].control, 'Vehicle');
    expect(found[0].readDay, 21, 'day 28 has one control only');
    // 1 − ΔT/ΔC: (50 e^2.1 − 100) / (100 e^2.1 − 100) of the control's change from day 0
    expect(Math.round(found[0].tgi!), 57);
    // another drug, or a combination, has none
    expect(anchorsIn(earlier(), 'Cpd-Y 20 mg/kg').length, 0);
    expect(anchorsIn(earlier(), 'Cpd-X 20 mg/kg + anti-PD-1').length, 0);
    // animals counted as distinct animals of their group, also when numbered 1..n in every group
    expect(anchorsIn(earlier('MC38', 'per group'), 'Cpd-X 20 mg/kg')[0].readDay, 21);
  });

  test('the anchor is the dose nearest the arm\'s, whatever the order of the groups in the file', async () => {
    const found = anchorsIn(earlier(), 'Cpd-X 20 mg/kg');
    for (const order of [found, [...found].reverse()]) {
      expect(nearestAnchor(order, 'Cpd-X 20 mg/kg')!.arm, 'Cpd-X 10 mg/kg i.p. QD', '20: 2× from 10');
      expect(nearestAnchor(order, 'Cpd-X 4 mg/kg')!.arm, 'Cpd-X 3 mg/kg i.p. QD', '4: nearer 3 than 10');
    }
    expect(nearestAnchor([], 'Cpd-X 20 mg/kg'), null);
  });

  test('the note says the analysis\'s TGI when the study kept one, else the same TGI from its rows',
    async () => {
      const never = anchorsIn(earlier(), 'Cpd-X 20 mg/kg')[0];
      expect(offerText(never, 'MC38'), 'You have a study of Cpd-X at 10 mg/kg (TGI 57%, day 21).');
      const analyzed = {...earlier(), exposure: {effect: {day: 12, tgi: {'Cpd-X 10 mg/kg i.p. QD': 18.4}}}} as
        PreviousStudy;
      const a = anchorsIn(analyzed, 'Cpd-X 20 mg/kg')[0];
      expect([a.analyzed, a.readDay].join(), 'true,21', 'Crux still reads its own day');
      expect(offerText(a), 'You have a study of Cpd-X at 10 mg/kg in MC38 (TGI 18%, day 12).', 'no model planned yet');
      expect(offerText(a, 'MC38'), 'You have a study of Cpd-X at 10 mg/kg (TGI 18%, day 12).', 'same model: unsaid');
      expect(offerText(a, 'CT26').includes('10 mg/kg in MC38 ('), true, 'another model: said');
      const unrecorded = anchorsIn({...analyzed, line: null}, 'Cpd-X 20 mg/kg')[0];
      expect(offerText(unrecorded, 'CT26').includes('in a tumor model the study does not record'), true,
        'a study with no model: said, not taken for the plan\'s');
    });

  test('the control: the group its details say gets a vehicle, else the one named so; none: no anchor', async () => {
    const meta = {...emptyMeta(), groups: [{group: 'G1', compound: 'Vehicle', dose: null, unit: '', route: '',
      schedule: ''}, {group: 'G2', compound: 'Cpd-X', dose: 10, unit: 'mg/kg', route: 'i.p.', schedule: 'QD'}]};
    expect(controlOf(['G2', 'G1'], meta), 'G1');
    expect(controlOf(['Cpd-X 10 mg/kg', 'PBS'], null), 'PBS');
    expect(controlOf(['Cpd-X 10 mg/kg', 'Cpd-X 3 mg/kg'], null) === undefined, true);
    const study = earlier();
    study.rows = study.rows.filter((r) => r.group !== 'Vehicle');
    expect(anchorsIn(study, 'Cpd-X 20 mg/kg').length, 0, 'no control: silent, not the first group');
  });

  test('a study\'s own plan never offers that study: same ID, or same name before it has one', async () => {
    const s = {...earlier(), name: 'MC38 full study', meta: {...emptyMeta(), id: 'PMX-2026-014'}};
    expect(isOwnStudy(s, {id: 'PMX-2026-014', name: 'Another name'}), true, 'same ID');
    expect(isOwnStudy(s, {id: '', name: 'MC38 full study'}), true, 'same name, plan not saved yet');
    expect(isOwnStudy(s, {id: 'PMX-2026-015', name: 'Other'}), false, 'another study');
    expect(isOwnStudy({...s, meta: null}, {id: '', name: ''}), false, 'a new plan excludes nothing');
    expect(isOwnStudy(s, null), false);
  });

  test('a group keeps its colour whatever order a stage lists the groups in', async () => {
    const plan = armColors(['LD4172 20 mg/kg', 'anti-PD-1 100 ug', 'LD4172 + anti-PD-1']);
    const monitor = groupColors(['Vehicle', 'anti-PD-1 100 ug', 'LD4172 + anti-PD-1', 'LD4172 20 mg/kg'], 'Vehicle');
    for (const g of ['LD4172 20 mg/kg', 'anti-PD-1 100 ug', 'LD4172 + anti-PD-1'])
      expect(monitor[g], plan[g], g);
    expect(monitor['Vehicle'], plan['Control'], 'the control in the slate');
    same(Object.keys(armColors(['b', 'a'])), ['Control', 'b', 'a'], 'the listed order kept (legends, categories)');
  });

  test('the examples an empty arm name shows are names Plan reads whole', async () => {
    for (const ex of NAME_EXAMPLES) {
      const parts = parseArmName(ex);
      expect(parts.length > 0 && parts.every((p) => p.dose && p.route && p.schedule), true, ex);
    }
    expect(parseArmName(NAME_EXAMPLES[2]).length, 2, 'the combination');
  });

  test('arms added unnamed keep their examples; the table numbers them', async () => {
    const root = Scope.root('unnamed arms test');
    try {
      const input = root.build(() => new ArmsInput({}));
      input.add();
      input.add();
      await delay(50);
      const areas = [...input.root.querySelectorAll<HTMLTextAreaElement>('.pmx-arm textarea')];
      same(areas.map((a) => a.value), ['', ''], 'the first card, rebuilt by the second add, still empty');
      same(areas.map((a) => a.placeholder), NAME_EXAMPLES.slice(0, 2).map((e) => `e.g. ${e}`));
      same(armsOf(input.value).map((a) => a.name), ['Treatment 1', 'Treatment 2']);
    } finally {
      root.dispose();
    }
  });

  test('a group of the drug that says no dose: no anchor, but named so Plan can ask for it', async () => {
    const study = earlier();
    study.rows = study.rows.filter((r) => r.group !== 'Cpd-X 3 mg/kg i.p. QD')
      .map((r) => r.group === 'Cpd-X 10 mg/kg i.p. QD' ? {...r, group: 'Cpd-X'} : r);
    expect(anchorsIn(study, 'Cpd-X 20 mg/kg').length, 0);
    expect(undosedIn(study, 'Cpd-X 20 mg/kg'), 'Cpd-X');
    // its details give the dose: an anchor, nothing to ask
    study.meta = {...study.meta!, groups: [{group: 'Cpd-X', compound: 'Cpd-X', dose: 10, unit: 'mg/kg',
      route: 'i.p.', schedule: 'QD'}]};
    expect(anchorsIn(study, 'Cpd-X 20 mg/kg').length, 1);
    expect(undosedIn(study, 'Cpd-X 20 mg/kg'), null);
    // another drug, a dosed study, or no control: nothing to ask
    expect(undosedIn(earlier(), 'Cpd-X 20 mg/kg'), null);
    expect(undosedIn(study, 'Cpd-Y 20 mg/kg'), null);
    expect(undosedText('Cpd-X', 'Earlier study', 'Cpd-X').includes('does not say its dose'), true);
  });

  test('Plan stays silent without asking Crux: combination, dose, unit, model, day', async () => {
    const a: Anchor = anchorsIn(earlier(), 'Cpd-X 20 mg/kg')[0];
    expect(silence('Cpd-X 20 mg/kg i.p. QD', a, 'MC38', 21), null);
    expect(silence('Cpd-X 20 mg/kg', a, 'MC38', 22)!.includes('day 21'), true, 'a day past its last: silent');
    // any day the study measured, not only its last: the plan's own day is read there
    same(a.days, [0, 7, 14, 21]);
    expect(silence('Cpd-X 20 mg/kg', a, 'MC38', 14), null, 'a day it measured');
    expect(silence('Cpd-X 20 mg/kg', a, 'MC38', 13.6), null, 'the plan\'s day rounded');
    expect(silence('Cpd-X 20 mg/kg', a, 'MC38', 10)!.includes('days 7 and 14'), true, 'a day it did not measure');
    expect(silence('Cpd-X 20 mg/kg + anti-PD-1', a, 'MC38', 21)!.includes('combination'), true);
    expect(silence('Cpd-X', a, 'MC38', 21)!.includes('dose'), true);
    expect(silence('Cpd-X 20 mg', a, 'MC38', 21)!.includes('unit'), true);
    expect(silence('Cpd-X 20 mg/kg', a, 'CT26', 21)!.includes('MC38'), true);
    expect(silence('Cpd-X 20 mg/kg', anchorsIn(earlier(null), 'Cpd-X 20 mg/kg')[0], 'MC38', 21)!
      .includes('tumor model'), true);
    expect(silence('Cpd-X 20 mg/kg', a, 'MC38', 30)!.includes('day 21'), true);
    expect(silence('Cpd-X 20 mg/kg', a, 'MC38', null)!.includes('not known yet'), true);
    // another kind of host mouse
    const host = (l: string | null) => l === 'MC38' ? 'syngeneic' : 'humanized';
    expect(silence('Cpd-X 20 mg/kg', a, 'MC38', 21, host), null);
    const other = {...a, study: {...a.study, line: 'MC38-hu'}};
    expect(silence('Cpd-X 20 mg/kg', other, 'MC38', 21, host)!.includes('host'), true);
    // its days counted from the first dose: a day 0 in its details, no days before
    const noStart = {...a, study: {...a.study, meta: {...a.study.meta!, start: null}}};
    expect(silence('Cpd-X 20 mg/kg', noStart, 'MC38', 21)!.includes('day 0'), true);
    const before = {...a, study: {...a.study, rows: [...a.study.rows, {animal: 'x', group: 'Vehicle', day: -3,
      volume: 90}]}};
    expect(silence('Cpd-X 20 mg/kg', before, 'MC38', 21)!.includes('before day 0'), true);
  });

  test('the effect field: point and range, or why; T/C turned round', async () => {
    expect(effectText(linked, 'TGI'), '≈ 65%, between 55 and 80%, from your 10 mg/kg study');
    expect(effectText(linked, 'TC'), '≈ 35%, between 20 and 45%, from your 10 mg/kg study');
    expect(effectText({...linked, point: null, range: null}, 'TGI'), 'Estimating from your 10 mg/kg study…');
    expect(effectText({...linked, point: null, range: null, why: 'the route is not your study\'s'}, 'TGI'),
      'No basis to scale from your 10 mg/kg study: the route is not your study\'s');
    expect(assumedText('Cpd-X 20 mg/kg', linked), 'Assumes the same schedule and route as your 10 mg/kg study.');
    expect(assumedText('Cpd-X 20 mg/kg QD', linked), 'Assumes the same route as your 10 mg/kg study.');
    expect(assumedText('Cpd-X 20 mg/kg i.p. QD', linked), '');
    // each of Crux's eight reasons has its own sentence; an unknown one, the general
    const codes = ['schedule_differs', 'route_differs', 'ratio_out_of_range', 'pk_not_proportional',
      'dose_outside_pk_window', 'too_few_animals', 'anchor_effect_out_of_range', 'control_fit_failed'];
    expect(new Set(codes.map(plainReason)).size, 8);
    expect(codes.map(plainReason).includes('no basis to scale from your study'), false);
    expect(plainReason('ratio_out_of_range'), 'the new dose is more than 4 times, or less than a quarter of, ' +
      'your study\'s');
    expect(plainReason('something new'), 'no basis to scale from your study');
  });

  test('a linked arm is sized on the low end; the scenarios add its low end and point', async () => {
    expect(sizingEffect(linked, 'TGI'), 55);
    expect(sizingEffect(linked, 'TC'), 45);
    expect(sizingEffect({...linked, range: [2, 30]}, 'TGI'), 5, 'the field\'s floor');
    expect(sizingEffect({...linked, range: null, point: null, why: 'x'}, 'TGI'), null);
    expect(sizingEffect(undefined, 'TGI'), null);
    same(linkedScenarios([{fromStudy: linked}, {}]), [55, 65]);
    // the declaration Plan sizes: the low end, whatever effect the table carries
    const arms = armsFrame([{name: 'Cpd-X 20 mg/kg', therapy: 'small_molecule', metric: 'TGI', effect: 70,
      fromStudy: linked}, {name: 'Cpd-Z', therapy: 'small_molecule', metric: 'TGI', effect: 40}]);
    const d = declarationOf({line: 'MC38', endpoint: 1000, power: 80, correction: 'holm', arms});
    same(d.arms.map((a) => a.effect), [55, 40]);
  });

  test('the link survives the arms table, a saved plan and the cards', async () => {
    const arms: Arm[] = [{name: 'Cpd-X 20 mg/kg', therapy: 'small_molecule', metric: 'TGI', effect: 55,
      therapySource: 'user', fromStudy: linked}, {name: 'Cpd-Z', therapy: 'small_molecule', metric: 'TGI',
      effect: 40, therapySource: 'user'}];
    same(armsOf(armsFrame(arms)), arms, 'through the table');
    const measured: Arm = {...arms[0], fromStudy: {...linked, anchorTgi: 50}};
    same(armsOf(armsFrame([measured]))[0], measured, 'with the effect the study measured');
    // through a .pmxplan: the inputs' arms as a list, read back by parsePlan, then given to the form as a table
    const plan = parsePlan((await exampleStudy()).getTag(PLAN_TAG)!);
    const back = parsePlan(JSON.stringify({...plan, inputs: {...plan.inputs, arms}}));
    same(armsOf(armsFrame(back.inputs.arms)), arms, 'through a saved plan');
    // a damaged link is dropped, the arm kept
    expect(armsOf(armsFrame([{...arms[1], fromStudy: {study: 'x'} as any}]))[0].fromStudy === undefined, true);
    // through the cards: shown, sized on the low end, written back unchanged (no earlier study at hand)
    const root = Scope.root('linked arm test');
    try {
      const input = root.build(() => new ArmsInput({host: signal<string | null>(null),
        metric: signal<'TGI' | 'TC'>('TGI')}));
      document.body.append(input.root);
      input.value = armsFrame(back.inputs.arms);
      await delay(50);
      same(armsOf(input.value), arms, 'through the cards');
      const card = input.root.querySelector<HTMLElement>('.pmx-arm')!;
      // the effect field as always, holding the low end the arm is sized on; one line under it says where from
      const field = card.querySelector<HTMLInputElement>('input:not([type=checkbox])')!;
      expect(field.value, '55');
      const shown = card.querySelector<HTMLElement>('[data-pmx=from-study] span')!;
      expect(shown.textContent, 'Scaled from your 10 mg/kg study and its PK: between 55 and 80% (≈ 65%)');
      // every range says its model is the simplified one
      expect(shown.title.endsWith(SIMPLIFIED), true, shown.title);
      input.root.remove();
    } finally {
      root.dispose();
    }
  });

  test('the PK of a drug: its dose levels and the AUC at each, from the frames of NCA', async () => {
    const results = DG.DataFrame.fromObjects([{compound: 'Cpd-X', DOSE: 10}])!;
    results.setTag('nca.units', JSON.stringify({dose: 'mg/kg'}));
    const summary = DG.DataFrame.fromObjects([['Cpd-X', 30, 4100], ['Cpd-X', 10, 1200], ['Cpd-Y', 5, 300]]
      .map(([compound, DOSE, MEAN_AUC_LAST]) => ({compound, DOSE, MEAN_AUC_LAST, CV_PCT_AUC_LAST: 18})))!;
    summary.col('MEAN_AUC_LAST')!.setTag('units', 'h*ng/mL');
    summary.setTag('nca.keyColumns', JSON.stringify(['compound']));
    const pk = drugPkOf('cpdx', results, summary, null)!;
    same(pk.doses, [10, 30]);
    expect(pk.unit, 'mg/kg');
    expect(pk.auc.get(10), '1,200 (CV 18%) h*ng/mL');
    expect(pk.proportional, 'yes', '120 and 137 per mg/kg: within 25%');
    expect(drugPkOf('cpdz', results, summary, null), null);
  });

  test('the PK of a drug by the arm\'s route; profiles NCA left to review are said', async () => {
    const results = DG.DataFrame.fromObjects([{compound: 'Cpd-X', DOSE: 10}])!;
    results.setTag('nca.units', JSON.stringify({dose: 'mg/kg'}));
    const summary = DG.DataFrame.fromObjects([['IV-bolus', 10, 3000, 6, 0], ['IP', 10, 1200, 6, 2]]
      .map(([ROUTE, DOSE, MEAN_AUC_LAST, N, N_REVIEW]) => ({compound: 'Cpd-X', ROUTE, DOSE, MEAN_AUC_LAST,
        CV_PCT_AUC_LAST: 18, N, N_REVIEW})))!;
    summary.setTag('nca.keyColumns', JSON.stringify(['compound']));
    expect(drugPkOf('cpdx', results, summary, null, 'i.p.')!.auc.get(10),
      '1,200 (CV 18%); 2 of 6 profiles still to review in NCA', 'the i.p. row, not the first');
    expect(drugPkOf('cpdx', results, summary, null, 'i.v.')!.auc.get(10), '3,000 (CV 18%)', 'none to review');
    expect(drugPkOf('cpdx', results, summary, null, 'p.o.'), null, 'no PK by mouth');
  });

  test('exposure in step with the dose: AUC per dose within 25% for every pair', async () => {
    expect(proportionalOf(new Map([[10, 1200]])), 'unknown');
    expect(proportionalOf(new Map([[10, 1000], [30, 3700]])), 'yes');
    expect(proportionalOf(new Map([[10, 1000], [30, 4000]])), 'no');
    expect(proportionalOf(new Map([[3, 300], [10, 1000], [30, 2000]])), 'no', 'one pair apart is enough');
  });

  test('no PK of the drug: silent, saying what to add', async () => {
    const a = anchorsIn(earlier(), 'Cpd-X 20 mg/kg')[0];
    same(await effectFromStudy(a, 'Cpd-X 20 mg/kg', null), {point: null, range: null, why: noPkText('Cpd-X')});
    expect((await effectFromStudy(a, 'Cpd-X 20 mg/kg', {...cpdxPk, unit: 'mg'})).why, noPkText('Cpd-X'),
      'other unit');
  });

  test('within the PK: a point and a range from Crux, the same with animals numbered per group', async () => {
    const r = await effectFromStudy(anchorsIn(earlier(), 'Cpd-X 20 mg/kg')[0], 'Cpd-X 20 mg/kg', cpdxPk);
    if (r.why === 'this version of PMX Studio cannot scale from an earlier study yet')
      return 'skipped: the loaded Crux module has no stats.tumor_effect_from_study yet';
    expect(r.why === undefined, true, JSON.stringify(r));
    expect(Number.isFinite(r.point) && r.range !== null && r.range.every(Number.isFinite), true, JSON.stringify(r));
    const again = await effectFromStudy(anchorsIn(earlier('MC38', 'per group'), 'Cpd-X 20 mg/kg')[0],
      'Cpd-X 20 mg/kg', cpdxPk);
    same(again, r, 'ids 1..4 in every group are 12 animals, not 4');
    // read on the plan's day, one the study measured before its last: Crux reads the rows up to it
    const on14 = await effectFromStudy(anchorsIn(earlier(), 'Cpd-X 20 mg/kg')[0], 'Cpd-X 20 mg/kg', cpdxPk, 14);
    expect(on14.why === undefined && on14.range !== null, true, JSON.stringify(on14));
  }, {timeout: 60000});
});
