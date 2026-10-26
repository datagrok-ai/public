import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {OBSERVED, PREDICTED, exampleStudy} from '../analyze';
import {PLAN_TAG, parsePlan, powerAt} from '../plan-file';
import {noEffectNote} from '../plan-vs-observed';

const OUTPUTS = ['analysis', 'comparisons', 'animals', 'volumes', 'twins', 'effects', 'checks', 'versusPlan',
  'predicted', 'tolerability', 'weights'];

/** pmxAnalyzeStudy as the Analyze view calls it: the outputs by name. */
async function analyze(study: DG.DataFrame, params: {[k: string]: any} = {}): Promise<{[k: string]: DG.DataFrame}> {
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxAnalyzeStudy'})[0].prepare({study,
    animal: study.col('animal'), group: study.col('group'), day: study.col('day'), volume: study.col('volume'),
    line: 'MC38', endpoint: 1500, correction: 'dunnett', ...params});
  await fc.call();
  return Object.fromEntries(OUTPUTS.map((o) => [o, fc.outputs[o] as DG.DataFrame]));
}

category('Analyze', () => {
  test('tolerability: the example\'s body weight read per group, and nothing without the column', async () => {
    const study = await exampleStudy();
    const o = await analyze(study, {bodyWeight: study.col('body_weight')});
    const t = o.tolerability;
    // three groups, each as weighed and without the tumor's mass
    expect(t.rowCount, 6);
    const row = (group: string, basis: string) => {
      const i = Array.from({length: t.rowCount}, (_, k) => k)
        .find((k) => t.get('group', k) === group && t.get('basis', k) === basis)!;
      return (c: string) => t.get(c, i);
    };
    const [vehicle, strong] = ['Vehicle', 'Cpd TGI 70%'].map((g) => row(g, 'Measured'));
    // the stronger arm loses more weight than the vehicle, whose animals never lose 20%
    expect(strong('mean_nadir') < vehicle('mean_nadir'), true,
      `mean nadir ${strong('mean_nadir')} vs ${vehicle('mean_nadir')}`);
    expect(vehicle('grade_3'), 0);
    expect(vehicle('n'), o.analysis.get('n_control', 0));
    // the vehicle's tumors weigh the most: without them its weight change is lower on its last day
    const w = o.weights;
    const last = (basis: string) => {
      const ks = Array.from({length: w.rowCount}, (_, k) => k)
        .filter((k) => w.get('group', k) === 'Vehicle' && w.get('basis', k) === basis);
      return w.get('mean', ks.reduce((a, b) => w.get('day', b) > w.get('day', a) ? b : a));
    };
    expect(last('Tumor-corrected') < last('Measured'), true, `${last('Tumor-corrected')} vs ${last('Measured')}`);
    // the NCI limits: a mean nadir over 20% as weighed, or over 10% of the animals dead or removed for their weight
    for (const g of ['Vehicle', 'Cpd TGI 40%', 'Cpd TGI 70%'].map((x) => row(x, 'Measured'))) {
      expect(g('toxic'), g('mean_nadir') < -20 || (g('removed_death') + g('removed_weight')) / g('n') > 0.1,
        `nadir ${g('mean_nadir')}, removed ${g('removed_weight')} of ${g('n')}`);
    }
    expect(vehicle('toxic'), false);
    // a removal reason given: one vehicle animal found dead counts as a death
    const reason = DG.Column.fromStrings('removal', study.col('animal')!.toList().map((a, i) =>
      a === study.get('animal', 0) && i === 0 ? 'Found dead' : ''));
    study.columns.add(reason);
    const dead = (await analyze(study, {bodyWeight: study.col('body_weight'), removalReason: reason})).tolerability;
    const v = Array.from({length: dead.rowCount}, (_, k) => k).find((k) => dead.get('group', k) === 'Vehicle')!;
    expect(dead.get('removed_death', v), 1);
    study.columns.remove('removal');
    // no body weight column: empty tables, the rest of the analysis unchanged
    const bare = await analyze(study);
    expect(bare.tolerability.rowCount + bare.weights.rowCount, 0);
    expect(bare.comparisons.get('p_value', 1), o.comparisons.get('p_value', 1));
  }, {timeout: 300000});

  test('tolerability: the tumor\'s density, the line\'s, the one given or the study\'s weighed tumors\'', async () => {
    const study = await exampleStudy();
    const densities = async (params: {[k: string]: any}) => {
      const t = (await analyze(study, {bodyWeight: study.col('body_weight'), ...params})).tolerability;
      return Array.from({length: t.rowCount}, (_, i) => ({group: t.get('group', i), density: t.get('density', i),
        pairs: t.get('density_pairs', i)}));
    };
    // no tumor weights: MC38's typical density (the atlas), no pairs
    for (const g of await densities({}))
      expect(Math.abs(g.density - 0.749) < 1e-6 && g.pairs === 0, true, JSON.stringify(g));
    for (const g of await densities({density: 1.2}))
      expect(Math.abs(g.density - 1.2) < 1e-6, true, JSON.stringify(g));
    // the example's excised tumors, drawn at MC38's density: every vehicle tumor calibrates it, near the line's
    const weighed = await densities({tumorWeight: study.col('tumor_weight')});
    const vehicle = weighed.find((g) => g.group === 'Vehicle')!;
    expect(vehicle.pairs, study.col('group')!.toList().filter((g, i) => g === 'Vehicle' &&
      !study.col('tumor_weight')!.isNone(i)).length);
    for (const g of weighed)
      expect(g.density > 0.749 / 1.3 && g.density < 0.749 * 1.3, true, JSON.stringify(g));
  }, {timeout: 300000});

  test('a study Plan predicts: its end, its effects and its strong arm found', async () => {
    const study = await exampleStudy();
    const o = await analyze(study);
    expect(o.analysis.get('control', 0), 'Vehicle');
    expect(o.comparisons.rowCount, 2);
    // one lab of the atlas, with the animals Plan asks for: the evaluation day in the range a study run as planned
    // reaches the endpoint in, the 70% arm detected, every animal analyzed
    const plan = parsePlan(study.getTag(PLAN_TAG)!);
    const {end} = plan.predicted.checks;
    const day = o.analysis.get('evaluation_day', 0);
    expect(day >= end.q025 && (end.q975 === null || day <= end.q975), true, `day ${day} in ${end.q025}–${end.q975}`);
    expect(o.comparisons.get('activity', 1), 'Slows growth', `p ${o.comparisons.get('p_value', 1)}`);
    expect(o.animals.rowCount, plan.predicted.nControl! +
      plan.predicted.arms.reduce((t, a) => t + a.n!, 0));
  }, {timeout: 300000});

  test('the twin: the strong arm smaller than its twins, bands ordered, controls in their bands', async () => {
    const o = await analyze(await exampleStudy());
    const c = o.comparisons;
    // TGI 40% and 70% arms: the paired effect follows the declared one, more clear responders in the strong arm
    const [weak, strong] = [c.get('paired_effect', 0), c.get('paired_effect', 1)];
    expect(strong > weak && strong > 40, true, `paired effects ${weak} / ${strong}`);
    // compared as shares, since the arms have the animals each needs, not the same
    const share = (i: number) => c.get('responders', i) / c.get('twin_animals', i);
    expect(share(1) > share(0) && share(1) >= 0.5, true, `responders ${c.get('responders', 0)} of ` +
      `${c.get('twin_animals', 0)} / ${c.get('responders', 1)} of ${c.get('twin_animals', 1)}`);
    for (let i = 0; i < o.twins.rowCount; i++) {
      const [lo, mid, hi] = ['twin_q05', 'twin_median', 'twin_q95'].map((k) => o.twins.get(k, i));
      expect(lo < mid && mid < hi, true, `row ${i}: ${lo} ${mid} ${hi}`);
    }
    expect(o.effects.rowCount, c.get('twin_animals', 0) + c.get('twin_animals', 1));
    // the example's controls come from the line's own model: most of their measurements in their bands
    const coverage = o.analysis.get('control_coverage', 0);
    expect(coverage > 0.7, true, `control coverage ${coverage}`);
  }, {timeout: 300000});

  test('the corrections order the p-values: none ≤ Dunnett ≤ Bonferroni', async () => {
    const study = await exampleStudy();
    const p = async (correction: string) => (await analyze(study, {correction})).comparisons.col('p_value')!.toList();
    const [none, dunnett, bonferroni] = [await p('none'), await p('dunnett'), await p('bonferroni')];
    for (let i = 0; i < none.length; i++) {
      expect(none[i] <= dunnett[i] + 1e-12, true, `arm ${i}: ${none[i]} vs ${dunnett[i]}`);
      expect(dunnett[i] <= bonferroni[i] + 1e-12, true, `arm ${i}: ${dunnett[i]} vs ${bonferroni[i]}`);
    }
  }, {timeout: 300000});

  test('plan vs observed: the example against its own plan, and without one against the atlas', async () => {
    const study = await exampleStudy();
    // the plan the example carries reads back as a plan; anything else is refused
    const plan = parsePlan(study.getTag(PLAN_TAG)!);
    expect(plan.predicted.arms.length, 2);
    // the plan carries what a study run as planned shows: each arm's range holds its assumed TGI
    expect(plan.predicted.checks.arms.map((a) => a.name).join(), plan.predicted.arms.map((a) => a.name).join());
    for (const [k, a] of plan.predicted.checks.arms.entries()) {
      const assumed = plan.predicted.arms[k].effect;
      expect(a.tgi.q025 < assumed && assumed < (a.tgi.q975 ?? Infinity), true, `${a.name}: ${JSON.stringify(a.tgi)}`);
    }
    for (const junk of ['{}', 'not json', JSON.stringify({...plan, version: 2})]) {
      let refused = false;
      try {
        parsePlan(junk);
      } catch {
        refused = true;
      }
      expect(refused, true, `accepted ${junk.slice(0, 20)}`);
    }
    const o = await analyze(study);
    const a = (c: string) => o.analysis.get(c, 0);
    expect(a('plan'), plan.name);
    // one lab of the atlas, not picked by its outcome, so each verdict is what its figure and its range say
    const inside = (x: number, lo: number, hi: number | null) => x >= lo && (hi === null || x <= hi);
    expect(a('control_verdict') === 'As predicted', inside(a('control_slope'), a('control_slope_q025'),
      a('control_slope_q975') || null), `slope ${a('control_slope')} in ${a('control_slope_q025')}–` +
      `${a('control_slope_q975')}: ${a('control_verdict')}`);
    const g = o.versusPlan;
    expect(g.rowCount, 2);
    expect(g.col('assumed')!.toList().join(), '40,70');
    const c = o.comparisons;
    for (let i = 0; i < g.rowCount; i++) {
      const x = (k: string) => g.get(k, i);
      expect(x('lower') < x('observed') && x('observed') < x('upper'), true,
        `row ${i}: ${x('lower')} ${x('observed')} ${x('upper')}`);
      // the result is the test on the slopes, in both tabs
      const [lo, hi] = [c.get('slope_difference_lower', i), c.get('slope_difference_upper', i)];
      expect(x('activity'), hi < 0 ? 'Slows growth' : lo > 0 ? 'Grows faster' : 'No clear effect', `${lo}–${hi}`);
      expect(c.get('activity', i), x('activity'));
      // against the plan: inside its range, as planned only when the test tells the group from the control
      const held = inside(x('observed'), x('expected_q025'), x('expected_q975') || null);
      expect(x('verdict') === 'As planned' || x('verdict') === 'Inconclusive', held, `${x('verdict')}`);
      if (held)
        expect(x('verdict'), x('activity') === 'Slows growth' ? 'As planned' : 'Inconclusive');
      // each arm has the animals Plan asked for, and the power the plan gave it with them
      const arm = plan.predicted.arms[i];
      expect(x('n'), arm.n, `${arm.name}: ${x('n')} analysed, ${arm.n} planned`);
      expect(Math.abs(x('power') - powerAt(arm.power!, x('n'))!) < 1e-9 && x('power') >= 0.79, true,
        `${arm.name}: power ${x('power')}`);
    }
    expect(new Set(o.predicted.col('series')!.toList()).size, 2);

    // no plan: nothing assumed, the controls and the end against the atlas
    study.setTag(PLAN_TAG, '');
    const n = await analyze(study);
    // (a string column gives an empty string for a missing value)
    expect(!n.analysis.get('plan', 0), true);
    expect(n.versusPlan.col('assumed')!.toList().every((x) => x === null), true);
    expect(n.versusPlan.col('verdict')!.toList().every((x) => !x), true);
    // with no plan, still each group against the control
    expect(n.versusPlan.col('activity')!.toList().join(), o.versusPlan.col('activity')!.toList().join());
    expect([...new Set(n.predicted.col('series')!.toList())].sort().join(), [OBSERVED, PREDICTED].sort().join());
    expect(n.analysis.get('control_verdict', 0) !== null, true);
  }, {timeout: 300000});

  test('no clear effect: what the study could have detected, said with it', async () => {
    // with a plan, the power it gave the arm with these animals, and the animals it needed
    const row = (r: {[k: string]: unknown}) => (k: string) => r[k] ?? null;
    const planned = noEffectNote(row({group: 'A', activity: 'No clear effect', p_value: 0.31, n: 8, power: 0.55,
      n_needed: 26}), row({plan: 'P', power_target: 0.8}))!;
    expect(planned.includes('With 8 animals, a study like this one detects the planned effect in 55 of 100 labs; ' +
      'the plan needed 26 for 80 of 100.'), true, planned);
    expect(noEffectNote(row({group: 'A', activity: 'Slows growth', p_value: 0.01, n: 8}), row({})), null);
    // a plan saved before it kept the arms' power
    expect(noEffectNote(row({group: 'A', activity: 'No clear effect', p_value: 0.4, n: 8}), row({plan: 'P'}))!
      .includes('cannot tell this group from the control'), true);
    // no plan: a sham group, half of the example's controls, against the other half
    const study = await exampleStudy();
    study.setTag(PLAN_TAG, '');
    const sham = study.clone(DG.BitSet.create(study.rowCount, (i) => study.get('group', i) === 'Vehicle'));
    sham.setTag(PLAN_TAG, '');
    const ids = [...new Set(sham.col('animal')!.toList())];
    for (let i = 0; i < sham.rowCount; i++) {
      if (ids.indexOf(sham.get('animal', i)) % 2)
        sham.set('group', i, 'Sham');
    }
    const s = await analyze(sham, {control: 'Vehicle'});
    const x = (k: string) => s.versusPlan.get(k, 0);
    expect(x('activity'), 'No clear effect', `p ${x('p_value')}`);
    expect(x('mde') > 0 && x('mde') <= 100, true, `mde ${x('mde')}`);
    expect(noEffectNote(x, (k) => s.analysis.get(k, 0))!.includes(`a TGI of ${Math.round(x('mde'))}% or more`), true);
  }, {timeout: 300000});
});
