import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {armsFrame} from '../arms-input';
import {MIN_GROUP} from '../engine';

/** pmxPlanStudy as the Plan view calls it: the outputs by name. */
async function plan(params: {[k: string]: any}): Promise<{[k: string]: DG.DataFrame}> {
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxPlanStudy'})[0].prepare({line: 'MC38', endpoint: 1500,
    power: 80, correction: 'dunnett', ...params});
  await fc.call();
  return Object.fromEntries(['plan', 'armResults', 'powerCurve', 'curves', 'duration', 'schedules', 'calendar']
    .map((o) => [o, fc.outputs[o] as DG.DataFrame]));
}

const arm = (name: string, effect = 60) => ({name, therapy: 'small_molecule' as const, metric: 'TGI' as const, effect});

category('Plan', () => {
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

  test('three arms: each its own size, the control the largest, Dunnett at least one arm\'s', async () => {
    const one = (await plan({arms: armsFrame([arm('A', 50)])})).armResults.get('n', 0);
    const three = [arm('A', 50), arm('B', 60), arm('C', 70)];
    const dunnett = await plan({arms: armsFrame(three)});
    const none = await plan({arms: armsFrame(three), correction: 'none'});
    const sizes = (o: {[k: string]: DG.DataFrame}) => [0, 1, 2].map((i) => o.armResults.get('n', i) as number);
    expect(sizes(dunnett)[0] >= one, true, `${sizes(dunnett)} vs ${one}`);
    // without correction the 50% arm needs what it needs alone; a stronger arm needs no more
    expect(sizes(none)[0], one);
    expect(sizes(none)[2] <= sizes(none)[0], true);
    expect(dunnett.plan.get('n_control', 0), Math.max(...sizes(dunnett)));
    const randomised = [0, 1, 2].map((i) => dunnett.armResults.get('randomise', i) as number);
    expect(dunnett.plan.get('randomise_total', 0),
      dunnett.plan.get('randomise_control', 0) + randomised.reduce((a, b) => a + b));
  }, {timeout: 300000});

  test('a strong arm is held at the minimum group size, with a note', async () => {
    const o = await plan({arms: armsFrame([arm('Strong', 90)])});
    expect(o.armResults.get('n_needed', 0) < MIN_GROUP, true, `needs ${o.armResults.get('n_needed', 0)}`);
    expect(o.armResults.get('n', 0), MIN_GROUP);
    expect(String(o.plan.get('notes_groups', 0)).includes('minimum per group'), true);
  }, {timeout: 120000});

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
    expect(base.plan.get('study_end', 0), e);
    expect(base.duration.col('day')!.toList().includes(e), true);
    const later = await plan({arms: armsFrame([arm('Cpd-1')]), lastDay: e + 7});
    expect(later.plan.get('study_end', 0), e + 7);
    expect(later.plan.get('endpoint_day', 0), e);
    expect(later.armResults.get('kill_rate_per_day', 0), base.armResults.get('kill_rate_per_day', 0));
    expect(later.armResults.get('n', 0) <= base.armResults.get('n', 0), true);
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

  test('without the schedules, the same plan and an empty schedules table', async () => {
    const all = await plan({arms: armsFrame([arm('Cpd-1')])});
    const lean = await plan({arms: armsFrame([arm('Cpd-1')]), compareSchedules: false});
    expect(lean.schedules.rowCount, 0);
    expect(all.schedules.rowCount > 0, true);
    for (const name of ['plan', 'armResults', 'powerCurve', 'duration', 'calendar'])
      expect(lean[name].toCsv(), all[name].toCsv(), name);
  }, {timeout: 120000});

  test('the humane limit: declared, in the plan, and a note when the endpoint is not below it', async () => {
    const o = await plan({arms: armsFrame([arm('Cpd-1')]), humane: 1500});
    expect(o.plan.get('humane_limit', 0), 1500);
    expect(String(o.plan.get('notes_growth', 0)).includes('not below the humane limit'), true);
  }, {timeout: 120000});

  test('a fixed n gives each arm its power and smallest detectable effect', async () => {
    const free = await plan({arms: armsFrame([arm('Cpd-1')])});
    const n = free.armResults.get('n', 0);
    const fewer = await plan({arms: armsFrame([arm('Cpd-1')]), n: n - 1});
    expect(fewer.armResults.get('n', 0), n - 1);
    expect(fewer.armResults.get('assurance', 0) < 0.8, true);
    expect(fewer.armResults.get('mde', 0) > 60, true, `mde ${fewer.armResults.get('mde', 0)}`);
  }, {timeout: 300000});
});
