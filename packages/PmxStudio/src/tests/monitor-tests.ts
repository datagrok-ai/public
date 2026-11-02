import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

import {exampleStudy} from '../analyze';
import {allocateGroups, futureDays, parseLooks, reasonsOf} from '../monitor';
import {PLAN_TAG, parsePlan} from '../plan-file';
import {LOG_TAG, LOOKS_TAG, addLog, logOf, newStudyTable, parseStudy, studyOf, tableOf} from '../study-file';

/** The example study up to `today`, as a follower has it. */
async function upTo(today: number): Promise<DG.DataFrame> {
  const whole = await exampleStudy();
  const df = whole.clone(DG.BitSet.create(whole.rowCount, (i) => whole.get('day', i) <= today));
  df.setTag(PLAN_TAG, whole.getTag(PLAN_TAG)!);
  return df;
}

async function monitor(study: DG.DataFrame, looks: string): Promise<{[k: string]: DG.DataFrame}> {
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxMonitorStudy'})[0].prepare({study,
    animal: study.col('animal'), group: study.col('group'), day: study.col('day'), volume: study.col('volume'), looks,
    line: 'MC38', endpoint: 1500, humane: 2000, correction: 'dunnett'});
  await fc.call();
  return Object.fromEntries(['status', 'animals', 'days', 'decisions'].map((o) => [o, fc.outputs[o] as DG.DataFrame]));
}

category('Monitor', () => {
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

  test('who needs attention and why', async () => {
    const r = (p: number, pit: number, lastDay = 10, v = 500) => reasonsOf({p_limit_next: p, last_pit: pit}, 10, 2000,
      lastDay, v, 7);
    // gone since before the previous session: listed once, not every session after
    expectArray(r(0.01, 0.5, 2), []);
    expectArray(r(0.01, 0.5), []);
    expect(r(0.5, 0.5)[0].startsWith('May pass'), true);
    expect(r(0.01, 0.001)[0].startsWith('Last value unexpected'), true);
    expect(r(0.01, 0.5, 7)[0], 'Not measured on day 10');
    expect(r(0.01, 0.5, 10, 2100)[0].startsWith('Over'), true);
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
    expect(back.getTag(LOG_TAG)!.length > 0, true);
    expect(parsePlan(back.getTag(PLAN_TAG)!).name, plan.name);
    let bad = false;
    try {
      parseStudy('{"format": "other"}');
    } catch {
      bad = true;
    }
    expect(bad, true);
  });
});
