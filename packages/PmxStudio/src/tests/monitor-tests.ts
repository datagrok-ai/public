import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {Scope} from '@datagrok-libraries/u2';

import {exampleStudy} from '../analyze';
import {allocateGroups, checkLevels, futureDays, outsideBand, parseLooks, promptsByChance, reasonsOf,
  weightReasons} from '../monitor';
import {PLAN_TAG, parsePlan} from '../plan-file';
import {LIMITS_TAG, LOOKS_TAG, WEIGHT_TAG, addLog, logOf, newStudyTable, parseStudy, studyOf, tableOf}
  from '../study-file';
import {panelTabs} from '../panels';

/** The example study up to `today`, as a follower has it. */
async function upTo(today: number): Promise<DG.DataFrame> {
  const whole = await exampleStudy();
  const df = whole.clone(DG.BitSet.create(whole.rowCount, (i) => whole.get('day', i) <= today));
  df.setTag(PLAN_TAG, whole.getTag(PLAN_TAG)!);
  return df;
}

async function monitor(study: DG.DataFrame, looks: string, bodyWeight: DG.Column | null = null):
Promise<{[k: string]: DG.DataFrame}> {
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxMonitorStudy'})[0].prepare({study,
    animal: study.col('animal'), group: study.col('group'), day: study.col('day'), volume: study.col('volume'), looks,
    line: 'MC38', endpoint: 1500, humane: 2000, correction: 'holm', bodyWeight, weightLimit: 20});
  await fc.call();
  return Object.fromEntries(['status', 'animals', 'days', 'decisions', 'tolerability', 'weights'].map((o) =>
    [o, fc.outputs[o] as DG.DataFrame]));
}

category('Monitor', () => {
  test('measure again by chance, in the words of a session', async () => {
    expect(promptsByChance(30, 0.037), 'about 1 a session');
    expect(promptsByChance(8, 0.037), 'about one every 3 sessions');
    // a line's calibrated levels go to Crux as they are; what they cost by chance is both tails of the 99% range
    const c = checkLevels({levels: {'90': {lo: 0.04, hi: 0.97}, '99': {lo: 0.002, hi: 0.998}}, papers: 12,
      points: 900, pooled: false, validated: {'99': {miss_lo: 0.004, miss_hi: 0.007, cover: 0.989}}, source: ''});
    expect(JSON.stringify(c.levels), '{"band_lo":0.04,"band_hi":0.97,"check_lo":0.002,"check_hi":0.998}');
    expect(Math.abs(c.chance! - 0.011) < 1e-12, true);
    // no calibration: Crux's own levels, and no rate to quote
    expect(JSON.stringify(checkLevels(undefined)), '{"levels":{},"chance":null}');
  });

  test('reviews and days to come', async () => {
    expectArray(parseLooks('10, 7;7'), [7, 10]);
    expectArray(parseLooks(null), []);
    let bad = false;
    try {
      parseLooks('day ten');
    } catch {
      bad = true;
    }
    expect(bad, true);
    expectArray(futureDays([0, 3, 6, 9]).slice(0, 2), [12, 15]);
  });

  test('a value entered is checked against its band, on its day only', async () => {
    const bands = new Map([['a', {check_lower: 100, check_upper: 200}],
      ['gone', {check_lower: null, check_upper: null}]]);
    const off = outsideBand([{animal: 'a', day: 12, volume: 300}, {animal: 'a', day: 12, volume: 150},
      {animal: 'a', day: 14, volume: 900}, {animal: 'a', day: 12, volume: NaN}, {animal: 'gone', day: 12, volume: 5},
      {animal: 'new', day: 12, volume: 5}], 12, bands);
    expect(off.length, 1);
    expect(off[0].volume, 300);
  });

  test('body weight: the limit, a margin before it, and an animal not weighed once', async () => {
    const r = (change: number | null, lastDay = 10, present = true) => weightReasons({change, lastDay}, 10, 7, 20,
      present);
    expect(r(-25)[0].startsWith('Over the 20%'), true);
    expect(r(-15)[0].startsWith('Losing weight'), true);
    expectArray(r(-14), []);
    expectArray(r(5), []);
    expectArray(r(null), []);
    expect(r(-1, 7)[0], 'Not weighed on day 10');
    // listed on the session it first misses, and not once it is gone from the study
    expectArray(r(-1, 4), []);
    expectArray(r(-1, 7, false), []);
  });

  test('body weight is optional, and when given it works', async () => {
    const study = await upTo(9);
    const without = await monitor(study, '');
    expect(without.tolerability.rowCount, 0);
    expect(without.animals.col('weight_change'), null);
    const withWeight = await monitor(study, '', study.col('body_weight'));
    expect(withWeight.tolerability.rowCount > 0, true);
    expect(withWeight.weights.rowCount > 0, true);
    expect(withWeight.animals.col('weight_change') !== null, true);
    expect(withWeight.animals.rowCount, without.animals.rowCount);
  }, {timeout: 120000});

  test('a saved study keeps its body weights and limits', async () => {
    const study = await upTo(9);
    study.setTag(WEIGHT_TAG, 'body_weight');
    const back = tableOf(parseStudy(JSON.stringify(studyOf(study, 1800, 15))));
    expect(back.getTag(WEIGHT_TAG), 'body_weight');
    expect(back.col('body_weight')!.stats.valueCount, study.col('body_weight')!.stats.valueCount);
    expect(JSON.parse(back.getTag(LIMITS_TAG)!).weightLimit, 15);
    expect(JSON.parse(back.getTag(LIMITS_TAG)!).humane, 1800);
    // a weight column the study does not use as its body weight travels as any other column: kept, not the weight
    const unused = tableOf(parseStudy(JSON.stringify(studyOf(await upTo(9), 2000))));
    expect(unused.col('body_weight') !== null, true);
    expect(unused.getTag(WEIGHT_TAG) ?? null, null);
    // and a table without one stays without
    const none = await upTo(9);
    none.columns.remove('body_weight');
    expect(tableOf(parseStudy(JSON.stringify(studyOf(none, 2000)))).col('body_weight'), null);
  });

  test('who needs attention and why', async () => {
    const r = (p: number, flag: number | null, lastDay = 10, v = 500) => reasonsOf({p_limit_next: p, last_flag: flag},
      10, 2000, lastDay, v, 7);
    // gone since before the previous session: listed once, not every session after
    expectArray(r(0.01, 0, 2), []);
    expectArray(r(0.01, 0), []);
    // a first value has no range to fall out of
    expectArray(r(0.01, null), []);
    expect(r(0.5, 0)[0].startsWith('May pass'), true);
    expect(r(0.01, -1)[0].startsWith('Last value unexpected'), true);
    expect(r(0.01, 1)[0].startsWith('Last value unexpected'), true);
    expect(r(0.01, 0, 7)[0], 'Not measured on day 10');
    expect(r(0.01, 0, 10, 2100)[0].startsWith('Over'), true);
  });

  test('a study followed: forecasts for every animal, the chance of success only on a review', async () => {
    const study = await upTo(9);
    const animals = new Set(study.col('animal')!.toList()).size;
    const onReview = await monitor(study, '9');
    expect(onReview.animals.rowCount, animals);
    expect(onReview.status.get('at_look', 0), true);
    expect(onReview.decisions.rowCount, 2);
    for (let i = 0; i < 2; i++) {
      const p = onReview.decisions.get('success', i);
      expect(p >= 0 && p <= 1, true, `chance ${p}`);
    }
    // the groups' growth reaches beyond today
    expect((onReview.days.col('day')!.toList() as number[]).some((d) => d > 9), true);
    // each animal's next value: the range to measure again outside lies around the expected one, and its last
    // value is flagged below (-1), inside (0) or above (1) the range it had
    let forecast = 0;
    for (let i = 0; i < onReview.animals.rowCount; i++) {
      if (onReview.animals.col('next_median')!.isNone(i))
        continue;
      forecast++;
      const a = (c: string): number => onReview.animals.get(c, i);
      expect(a('check_lower') < a('next_lower') && a('next_lower') < a('next_median') &&
        a('next_median') < a('next_upper') && a('next_upper') < a('check_upper'), true,
      `ranges of ${onReview.animals.get('animal', i)}`);
      expect([-1, 0, 1].includes(a('last_flag')), true, `flag ${a('last_flag')}`);
    }
    expect(forecast > 0, true);
    // what "measure again" costs by chance is quoted only where the line's ranges were validated
    const chance = onReview.status.col('check_chance');
    expect(!chance || chance.isNone(0) || (chance.get(0) > 0 && chance.get(0) < 0.2), true);
    const between = await monitor(study, '14');
    expect(between.status.get('at_look', 0), false);
    expect(between.decisions.rowCount, 0);
    expect(between.status.get('next_look', 0), 14);
  }, {timeout: 120000});

  test('groups alike in their tumors, the seed deciding who gets what', async () => {
    const animal = Array.from({length: 13}, (_, i) => `m${i}`);
    const volume = animal.map((_, i) => 80 + 3 * i);
    const groups = [{name: 'Vehicle', n: 6}, {name: 'A', n: 6}];
    const {allocation, balance} = await allocateGroups(animal, volume, groups, 3);
    expect(balance.rowCount, 2);
    expect(allocation.rowCount, 13);
    expect(allocation.col('group')!.toList().filter((g) => g).length, 12);
    const again = await allocateGroups(animal, volume, groups, 4);
    expect(allocation.col('group')!.toList().join() !== again.allocation.col('group')!.toList().join(), true);
  });

  test('a study file keeps the plan, the reviews, the measurements and the log', async () => {
    const whole = await exampleStudy();
    const plan = parsePlan(whole.getTag(PLAN_TAG)!);
    expect(newStudyTable(plan).rowCount, 0);
    const df = await upTo(5);
    df.name = 'S';
    df.setTag(LOOKS_TAG, '7, 10');
    addLog(df, 'Groups revealed');
    const back = tableOf(parseStudy(JSON.stringify(studyOf(df, 2000))));
    expect(back.rowCount, df.rowCount);
    expect(back.getTag(LOOKS_TAG), '7, 10');
    expect(logOf(back)[0].what, 'Groups revealed');
    expect(parsePlan(back.getTag(PLAN_TAG)!).name, plan.name);
    let bad = false;
    try {
      parseStudy('{"format": "other"}');
    } catch {
      bad = true;
    }
    expect(bad, true);
  });

  test('A panel on another tab opens that tab; the report shows every tab, then the open one alone again', async () => {
    const panel = (id: string) => ({id, title: id, hint: '', content: [document.createElement('div')]});
    const root = Scope.root('panels test');
    const t = root.build(() => panelTabs([panel('a'), panel('b')], [{name: 'One', layout: [{id: 'a', span: 12}]},
      {name: 'Two', layout: [{id: 'b', span: 12}]}], 'pmx-test', 'test'));
    const shown = () => [...t.root.querySelectorAll<HTMLElement>('.pmx-grid-host')].map((g) => !g.hidden).join();
    expect(shown(), 'true,false');
    t.show('b');
    expect(shown(), 'false,true');
    expect(await t.all(async () => shown()), 'true,true');
    expect(shown(), 'false,true');
    root.dispose();
  });
});
