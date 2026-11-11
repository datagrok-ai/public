import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {armsFrame} from '../arms-input';
import {fixedValue} from '../ui';
import {ENDPOINT, HUMANE, LAB_WEEKDAYS, MIN_GROUP, StudyDeclaration, designStudy, finite} from '../engine';
import {declarationOf, tgisOf} from '../frames';
import {enginePriors} from '../priors';

/** pmxPlanStudy as the Plan view calls it: the outputs by name. */
async function plan(params: {[k: string]: any}): Promise<{[k: string]: DG.DataFrame}> {
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxPlanStudy'})[0].prepare({line: 'MC38', endpoint: 1000,
    power: 80, correction: 'holm', ...params});
  await fc.call();
  return Object.fromEntries(['plan', 'armResults', 'powerCurve', 'curves', 'duration', 'schedules', 'calendar',
    'scenarios'].map((o) => [o, fc.outputs[o] as DG.DataFrame]));
}

const arm = (name: string, effect = 60) => ({name, therapy: 'small_molecule' as const, metric: 'TGI' as const, effect});

category('Plan', () => {
  test('a value set from a card is checked as its field was: whole, in range; empty goes back', async () => {
    const v = (raw: unknown) => JSON.stringify(fixedValue(raw, 2, 100));
    expect(v(12), JSON.stringify({value: 12}));
    expect(v(null), JSON.stringify({value: null}), 'back to the computed value: empty, not 0');
    expect(v(''), JSON.stringify({value: null}));
    expect(v(undefined), JSON.stringify({value: null}));
    for (const bad of [7.5, 1, 1000, NaN, 0, -0, Infinity, -Infinity, 'abc'])
      expect('error' in fixedValue(bad, 2, 100), true, `${bad} accepted`);
  });

  test('one treated arm: Holm and no correction size it alike', async () => {
    const n = async (correction: string) => (await plan({arms: armsFrame([arm('Cpd-A')]), correction}))
      .armResults.get('n', 0);
    expect(await n('holm'), await n('none'));
  });

  test('TGI or T/C is set for the study: it sets every arm; left empty, each arm keeps its own', async () => {
    const arms = armsFrame([{name: 'Cpd-1', therapy: 'small_molecule', metric: 'TGI', effect: 40},
      {name: 'Cpd-2', therapy: 'small_molecule', metric: 'TC', effect: 40}]);
    const base = {line: 'MC38', endpoint: 1000, arms, power: 80, correction: 'holm'};
    expect(declarationOf({...base, effectAs: 'T/C'}).arms.map((a) => a.metric).join(), 'TC,TC');
    expect(declarationOf({...base, effectAs: 'TGI'}).arms.map((a) => a.metric).join(), 'TGI,TGI');
    expect(declarationOf(base).arms.map((a) => a.metric).join(), 'TGI,TC');
  });

  test('one arm: the study, its groups, the power curve and the volumes', async () => {
    const o = await plan({arms: armsFrame([arm('Cpd-1')])});
    expect(o.plan.rowCount, 1);
    expect(o.armResults.rowCount, 1);
    const n = o.armResults.get('n', 0);
    expect(n >= MIN_GROUP && n <= 40, true, `n ${n}`);
    // after randomization 10% are lost, before it 20%: each count covers its stage
    expect(o.armResults.get('randomise', 0) > n, true);
    expect(o.plan.get('implant_total', 0) > o.plan.get('randomise_total', 0), true);
    // the control as large as the arm; the totals are the groups' sums
    expect(o.plan.get('n_control', 0), n);
    expect(o.plan.get('n_total', 0), 2 * n);
    // the kill rate gives the declared effect: its median over the labs
    const realised = o.armResults.get('effect_realised', 0);
    expect(Math.abs(realised - 60) < 0.5, true, `${realised}`);
    const groups = new Set(o.curves.col('arm')!.toList());
    expect(groups.has('Control') && groups.has('Cpd-1'), true);
    expect(o.curves.col('median')!.meta.units, 'mm³');
  }, {timeout: 120000});

  test('three arms: each its own size, the control the largest, Holm at least one arm\'s', async () => {
    const one = (await plan({arms: armsFrame([arm('A', 50)])})).armResults.get('n', 0);
    const three = [arm('A', 50), arm('B', 60), arm('C', 70)];
    const corrected = await plan({arms: armsFrame(three)});
    const none = await plan({arms: armsFrame(three), correction: 'none'});
    const sizes = (o: {[k: string]: DG.DataFrame}) => [0, 1, 2].map((i) => o.armResults.get('n', i) as number);
    expect(sizes(corrected)[0] >= one, true, `${sizes(corrected)} vs ${one}`);
    // without correction the 50% arm needs what it needs alone; a stronger arm needs no more
    expect(sizes(none)[0], one);
    expect(sizes(none)[2] <= sizes(none)[0], true);
    expect(corrected.plan.get('n_control', 0), Math.max(...sizes(corrected)));
    const randomised = [0, 1, 2].map((i) => corrected.armResults.get('randomise', i) as number);
    expect(corrected.plan.get('randomise_total', 0),
      corrected.plan.get('randomise_control', 0) + randomised.reduce((a, b) => a + b));
    // every group at once keeps its animals less often than each one does
    const all = corrected.plan.get('confidence_study', 0);
    expect(all > 0.5 && all < 1, true, `every group at once: ${all}`);
  }, {timeout: 300000});

  test('a strong arm is held at the minimum group size, with a note', async () => {
    const o = await plan({arms: armsFrame([arm('Strong', 90)])});
    expect(o.armResults.get('n_needed', 0) < MIN_GROUP, true, `needs ${o.armResults.get('n_needed', 0)}`);
    expect(o.armResults.get('n', 0), MIN_GROUP);
    expect(String(o.plan.get('notes_groups', 0)).includes('minimum per group'), true);
  }, {timeout: 120000});

  test('a line the atlas does not pin down: sized for the power averaged over its fits, said when they differ',
    async () => {
      const marked = Object.entries((await enginePriors()).lines)
        .filter(([, l]) => (l.calibration?.identification?.solutions?.length ?? 0) > 1);
      for (const [line, l] of marked) {
        const o = await plan({line, arms: armsFrame([arm('Cpd-1')])});
        // each fit's own answer, the prior's first ('>40': none up to 40)
        const fits = String(o.armResults.get('n_fits', 0)).split(', ').map((x) => x === '>40' ? Infinity : Number(x));
        expect(fits.length, l.calibration!.identification!.solutions!.length, `${line}: ${fits}`);
        // (two more than the most any needs: a fit's simulated power can dip under the target past its own n, and the
        // other fits' each come from a block of the studies, noisier)
        const needed = finite(o.armResults.get('n_needed', 0)) ?? Infinity;
        expect(Math.min(...fits) - 1 <= needed && needed <= Math.max(...fits) + 2, true,
          `${line}: ${needed} of ${fits}`);
        expect(finite(o.armResults.get('n', 0)), Number.isFinite(needed) ? Math.max(needed, MIN_GROUP) : null);
        // said only when they are more than one animal apart (one is a block's noise)
        const sized = fits.map((x) => Math.max(x, MIN_GROUP));
        expect(String(o.plan.get('notes_groups', 0)).includes('do not pin down'),
          Math.max(...sized) - Math.min(...sized) > 1, `${line}: ${fits}`);
        // the advanced choice: the atlas's best fit alone, its own n and nothing said
        const best = await plan({line, arms: armsFrame([arm('Cpd-1')]), fits: 'Best fit only'});
        expect(String(best.armResults.get('n_fits', 0)), String(fits[0]), line);
        expect(finite(best.armResults.get('n_needed', 0)), Number.isFinite(fits[0]) ? fits[0] : null, line);
        expect(String(best.plan.get('notes_groups', 0)).includes('do not pin down'), false, line);
      }
      // a line with one fit: its own n, and nothing said
      const one = await plan({arms: armsFrame([arm('Cpd-1')])});
      expect(String(one.armResults.get('n_fits', 0)), String(one.armResults.get('n_needed', 0)));
      expect(String(one.plan.get('notes_groups', 0)).includes('do not pin down'), false);
    }, {timeout: 300000});

  test('several equally good fits: the power is their average, at the n found or set', async () => {
    const d: StudyDeclaration = {line: 'MC38', vRand: 100, endpoint: ENDPOINT, humane: HUMANE, lastDay: null,
      weekdays: LAB_WEEKDAYS, startWeekday: 0, power: 0.8, alpha: 0.05, correction: 'holm', n: null,
      arms: [arm('Cpd-1')], preLoss: 0.2, postLoss: 0.1, confidence: 0.9, lossRule: 'confident', schedules: false};
    const lp = (await enginePriors()).lines['MC38'];
    const own = {sd_within: {log_alpha: lp.sd_within.log_alpha, log_beta: lp.sd_within.log_beta},
      study_prior: {sd_log_alpha: lp.study_prior!.sd_log_alpha, sd_log_beta: lp.study_prior!.sd_log_beta},
      distance_mean6: 1};
    // a second fit with a study's animals twice as far apart: it needs more of them
    const wider = {...own, sd_within: {log_alpha: 2 * own.sd_within.log_alpha, log_beta: 2 * own.sd_within.log_beta}};
    const kept = lp.calibration;
    try {
      lp.calibration = {identification: {status: 'poorly_identified', rule: '', why: '',
        solutions: [{...own, best: true}, wider]}};
      const s = await designStudy(d);
      const [first, second] = s.arms[0].fits as number[];
      expect(s.arms[0].fits.length, 2);
      expect(second > first, true, `${s.arms[0].fits}`);
      // sized for the power averaged over the two: between what each needs alone
      const needed = s.arms[0].nAssurance!;
      expect(first <= needed && needed <= second, true, `${needed} of ${s.arms[0].fits}`);
      expect(s.arms[0].n, Math.max(needed, MIN_GROUP));
      expect(s.arms[0].assurance! >= d.power, true, `${s.arms[0].assurance}`);
      expect(s.nControl, s.arms[0].n);
      expect(s.notes.groups.some((x) => x.includes('do not pin down') &&
        x.includes(`need about ${Math.max(first, MIN_GROUP)} to ${second}`)), true, s.notes.groups.join(' | '));
      // a set number of animals: its power is the two fits' average too, below the prior's own fit's
      const [both, alone] = [await designStudy({...d, n: 8, mde: false}),
        await designStudy({...d, n: 8, mde: false, fits: false})];
      expect(both.arms[0].fits.length, 2);
      expect(alone.arms[0].fits.length, 1);
      expect(both.arms[0].assurance! < alone.arms[0].assurance!, true,
        `${both.arms[0].assurance} vs ${alone.arms[0].assurance}`);
      // the prior's own fit not marked: which one it carries is not known, and none counts
      lp.calibration = {identification: {status: 'poorly_identified', rule: '', why: '', solutions: [own, wider]}};
      const unmarked = await designStudy(d);
      expect(unmarked.arms[0].fits.length, 1);
      expect(unmarked.arms[0].nAssurance, first);
      // three identical fits: the two others share the studies (two blocks), and the average is the one fit's power
      lp.calibration = {identification: {status: 'poorly_identified', rule: '', why: '',
        solutions: [{...own, best: true}, own, own]}};
      const [three, single] = [await designStudy({...d, n: 8, mde: false}),
        await designStudy({...d, n: 8, mde: false, fits: false})];
      expect(three.arms[0].fits.length, 3);
      expect(Math.abs(three.arms[0].assurance! - single.arms[0].assurance!) < 0.05, true,
        `${three.arms[0].assurance} vs ${single.arms[0].assurance}`);
      const sized = await designStudy(d);
      const each = sized.arms[0].fits as number[];
      expect(Math.min(...each) - 1 <= sized.arms[0].nAssurance! && sized.arms[0].nAssurance! <= Math.max(...each) + 1,
        true, `${sized.arms[0].nAssurance} of ${each}`);
      expect(sized.notes.groups.some((x) => x.includes('do not pin down')), Math.max(...each) - Math.min(...each) > 1);
      // three different fits: the averaged power is the plain mean of each fit's own, so a fit weighed more than
      // another, or a block that is not its studies' mean, fails here. Not guarded: how many blocks, or the blocks
      // lost (each fit on every study): both give the same mean, only slower, which the timings show.
      // (0.05: a block of half Crux's 400 studies; the power varies across them by about 0.25, so its Monte Carlo error
      // is about 0.25 / sqrt(200) = 0.018, and 0.05 is under three of it)
      const mid = {...own,
        sd_within: {log_alpha: 1.5 * own.sd_within.log_alpha, log_beta: 1.5 * own.sd_within.log_beta}};
      lp.calibration = {identification: {status: 'poorly_identified', rule: '', why: '',
        solutions: [{...own, best: true}, wider, mid]}};
      const mixed = await designStudy({...d, n: 8, mde: false});
      expect(mixed.arms[0].fits.length, 3);
      const shipped = lp.sd_within;
      const byFit: number[] = [];
      try {
        for (const fit of [own, wider, mid]) {
          lp.sd_within = {...shipped, ...fit.sd_within};
          byFit.push((await designStudy({...d, n: 8, mde: false, fits: false})).arms[0].assurance!);
        }
      } finally {
        lp.sd_within = shipped;
      }
      const mean = byFit.reduce((t, x) => t + x, 0) / byFit.length;
      expect(Math.abs(mixed.arms[0].assurance! - mean) < 0.05, true, `${mixed.arms[0].assurance} vs ${byFit}`);
    } finally {
      lp.calibration = kept;
    }
  }, {timeout: 300000});

  test('no losses: as many to randomize and implant as analyzable; 10 at 20% loss need 15 (90% sure)', async () => {
    const none = await plan({arms: armsFrame([arm('Cpd-1')]), preLoss: 0, postLoss: 0});
    const n = none.armResults.get('n', 0);
    // the atlas's MC38 non-progressing share (0.19%) is the floor of the post-randomization loss
    expect(none.armResults.get('randomise', 0) - n <= 1, true);
    expect(none.plan.get('implant_total', 0), none.plan.get('randomise_total', 0));
    // P(Binom(15, 0.8) >= 10) = 0.94, P(Binom(14, 0.8) >= 10) = 0.87 (scipy): 15 is the fewest for 90%
    const fixed = await plan({arms: armsFrame([arm('Cpd-1')]), n: 10, preLoss: 0, postLoss: 20});
    expect(fixed.armResults.get('randomise', 0), 15);
  }, {timeout: 300000});

  test('a chosen end of study: the effects stay on the endpoint day, the study is evaluated then', async () => {
    const base = await plan({arms: armsFrame([arm('Cpd-1')])});
    const e = base.plan.get('endpoint_day', 0);
    // the study ends on the first measurement day from the endpoint day on (Mon, Wed, Fri: at most 3 days later)
    const end = base.plan.get('study_end', 0);
    expect(end >= e && end - e <= 3, true, `ends on day ${end}, the endpoint day is ${e}`);
    expect(base.duration.col('day')!.toList().includes(e), true);
    const later = await plan({arms: armsFrame([arm('Cpd-1')]), lastDay: e + 7});
    expect(later.plan.get('study_end', 0), e + 7);
    expect(later.plan.get('endpoint_day', 0), e);
    expect(later.armResults.get('kill_rate_per_day', 0), base.armResults.get('kill_rate_per_day', 0));
    // a later end lowers the assurance a little at every n (~0.008 at MC38's 60 % TGI, measured with 200 and 400
    // calibration labs): n may cross the power target one animal later, not more
    expect(later.armResults.get('n', 0) <= base.armResults.get('n', 0) + 1, true);
  }, {timeout: 300000});

  test('schedules: the plan\'s own is the design, fewer measurements need no fewer animals', async () => {
    const o = await plan({arms: armsFrame([arm('Cpd-1')]), measurementDays: 'Mon, Thu'});
    const s = o.schedules;
    const at = (k: number) => s.toJson().find((r: any) => r.per_week === k);
    expect(s.col('per_week')!.toList().join(), '1,2,3,5,7');
    expect(at(2).n_total, o.plan.get('n_total', 0));
    expect(at(2).measurements, o.plan.get('measurements', 0));
    // with v2 the labs' spread dominates: the schedules' assurances differ by little more than the simulation's error
    expect(at(1).n_total >= at(7).n_total && at(1).assurance < at(7).assurance + 0.02, true);
    // the calendar: from randomization to the end, every measurement day in it
    const days = o.calendar.col('day')!.toList();
    expect(days[0], 0);
    expect(days[days.length - 1], o.plan.get('study_end', 0));
    expect(String(o.plan.get('measurement_days', 0)).split(', ').every((t) => days.includes(Number(t))), true);
  }, {timeout: 120000});

  test('randomized on a Wednesday: the measurements stay on Mon / Wed / Fri', async () => {
    const o = await plan({arms: armsFrame([arm('Cpd-1')]), startWeekday: 'Wednesday'});
    expect(o.plan.get('start_weekday', 0), 2);
    expect(String(o.plan.get('measurement_days', 0)).startsWith('0, 2, 5, 7, 9, 12'), true,
      String(o.plan.get('measurement_days', 0)));
  }, {timeout: 120000});

  test('measured on Tuesdays only: the study ends on the first Tuesday from the endpoint day', async () => {
    const o = await plan({arms: armsFrame([arm('Cpd-1')]), measurementDays: 'tue'});
    const [end, e] = [o.plan.get('study_end', 0), o.plan.get('endpoint_day', 0)];
    expect(end >= e && end - e < 7 && end % 7 === 1, true, `end ${end}, endpoint day ${e}`);
    expect(o.plan.get('measurement_weekdays', 0), 'Tue');
  }, {timeout: 120000});

  test('Explore: a larger assumed TGI needs no more animals, whatever the form\'s n and metric', async () => {
    const declared = {n: 12, arms: armsFrame([{name: 'Cpd-1', therapy: 'small_molecule', metric: 'TC', effect: 40}])};
    const o = await plan({...declared, tgiScenarios: '80, 40, 60'});
    const rows = o.scenarios.toJson() as {[k: string]: number}[];
    const totals = rows.filter((r) => r.per_week === 3).map((r) => r.n_total);
    expect(rows.filter((r) => r.per_week === 3).map((r) => r.tgi).join(), '40,60,80');
    expect(totals[0] >= totals[1] && totals[1] >= totals[2] && totals[2] >= 2 * MIN_GROUP, true, `${totals}`);
    // the assumption changes the animals, not the controls' calendar
    expect(new Set(rows.map((r) => r.study_end)).size, 1);
    // and the plan itself stays what the form declares
    const alone = await plan(declared);
    expect(alone.scenarios.rowCount, 0);
    expect(o.plan.toCsv(), alone.plan.toCsv());
  }, {timeout: 300000});

  test('Explore: a TGI no group size reaches has no count, not the column\'s empty value', async () => {
    const s = (await plan({arms: armsFrame([arm('Cpd-1')]), tgiScenarios: '15, 80'})).scenarios;
    const at = (tgi: number, column: string) => {
      const row = Array.from({length: s.rowCount}, (_, i) => i)
        .find((i) => s.get('tgi', i) === tgi && s.get('per_week', i) === 3)!;
      return finite(s.get(column, row));
    };
    for (const column of ['n_total', 'randomise_total', 'implant_total']) {
      expect(at(15, column), null, column);
      expect((at(80, column) ?? 0) > 0, true, column);
    }
  }, {timeout: 120000});

  test('Explore: the assumed TGIs of a text', async () => {
    expect(tgisOf(null).length, 0);
    expect(tgisOf(' ').length, 0);
    expect(tgisOf('70%; 30, 30 50').join(), '30,50,70');
    for (const bad of ['0', '100', 'high', '10 20 30 40 50 60 70 80 90']) {
      let refused = false;
      try {
        tgisOf(bad);
      } catch {
        refused = true;
      }
      expect(refused, true, bad);
    }
  });

  test('without the schedules, the same plan and an empty schedules table', async () => {
    const all = await plan({arms: armsFrame([arm('Cpd-1')])});
    const lean = await plan({arms: armsFrame([arm('Cpd-1')]), compareSchedules: false});
    expect(lean.schedules.rowCount, 0);
    expect(all.schedules.rowCount > 0, true);
    for (const name of ['plan', 'armResults', 'powerCurve', 'duration', 'calendar'])
      expect(lean[name].toCsv(), all[name].toCsv(), name);
  }, {timeout: 120000});

  test('at a fixed n, without the detectable effect the same plan; calls in parallel as one by one', async () => {
    const one = {arms: armsFrame([arm('Cpd-1')]), n: 8, compareSchedules: false};
    const all = await plan(one);
    const lean = await plan({...one, detectableEffect: false});
    expect(typeof all.armResults.get('mde', 0), 'number');
    expect(lean.armResults.col('mde')!.isNone(0), true);
    // without it the design runs in shards (engine.ts): the same sums in another order, equal to rounding
    const rounded = (df: DG.DataFrame) => df.toCsv().replace(/\d+\.\d+/g, (x) => Number(x).toFixed(10));
    expect(rounded(lean.plan), rounded(all.plan));
    expect(rounded(lean.powerCurve), rounded(all.powerCurve));
    const [a, b] = await Promise.all([plan(one), plan(one)]);
    expect(a.armResults.toCsv(), all.armResults.toCsv());
    expect(b.armResults.toCsv(), all.armResults.toCsv());
  }, {timeout: 120000});

  test('the humane limit: declared, in the plan, and a note when the endpoint is not below it', async () => {
    const o = await plan({arms: armsFrame([arm('Cpd-1')]), humane: 1000});
    expect(o.plan.get('humane_limit', 0), 1000);
    expect(String(o.plan.get('notes_growth', 0)).includes('not below the humane limit'), true);
    // a limit over the guideline's 1.5 cm mean diameter is said; the default (1500 mm³) is under it
    const [high, usual] = await Promise.all([plan({arms: armsFrame([arm('Cpd-1')]), humane: 2000}),
      plan({arms: armsFrame([arm('Cpd-1')])})]);
    expect(String(high.plan.get('notes_growth', 0)).includes('Workman'), true);
    expect(usual.plan.get('humane_limit', 0), 1500);
    expect(String(usual.plan.get('notes_growth', 0)).includes('Workman'), false);
  }, {timeout: 120000});

  test('a fixed n gives each arm its power and smallest detectable effect', async () => {
    const free = await plan({arms: armsFrame([arm('Cpd-1')])});
    const n = free.armResults.get('n', 0);
    // two fewer: with one fewer the power is just under the target, and the effect searched (other simulated
    // animals) can land a fraction of a point on either side of the declared one
    const fewer = await plan({arms: armsFrame([arm('Cpd-1')]), n: n - 2});
    expect(fewer.armResults.get('n', 0), n - 2);
    expect(fewer.armResults.get('assurance', 0) < 0.8, true);
    expect(fewer.armResults.get('mde', 0) > 60, true, `mde ${fewer.armResults.get('mde', 0)}`);
  }, {timeout: 300000});
});
