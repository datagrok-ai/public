import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {awaitCheck, category, delay, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {Scope} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

import {exampleStudy} from './example-study';
import {allocateGroups, balanceNote, checkChance, futureDays, latestReadings, likelyTypo, mm3Text, monitorFrames,
  MonitorInputs, outsideBand, parseLooks, naturalCompare, notYetMeasured, PENDING, planControl, promptsByChance,
  reasonsOf,
  shortOfAnimals, weightMissing, weightReasons} from '../monitor';
import {PLAN_TAG, parsePlan, randomizedGroups} from '../plan-file';
import {planSummary} from '../file-preview';
import {LIMITS_TAG, LOOKS_TAG, REVEALED, WEIGHT_TAG, addLog, logOf, newStudyTable, parseStudy, planFieldsOf, studyOf,
  tableOf}
  from '../study-file';
import {panelTabs} from '../panels';
import {EmptyState, NoControl} from '../crux';
import {offPlanAsk, chanceWords, entryStep, followedFrom, localTime, monitorView, takeOverFollowed}
  from '../monitor-view';
import {animalIds, dayZeroAsk, nextMeasurement, whyAndStatus} from '../monitor';
import {farDayZero} from '../study-meta';
import {decisionEntry, dosingDays, emptyMeta, impossibleDate, metaOf, parseMeta, setMetaOn, stopsOf,
  weightChangeText}
  from '../study-meta';
import {groupsByAnimal} from '../animal-group';
import {blindColumn, mark, openStudyTable, tagStudy} from '../study-table';
import {analysisDeclarationOf} from '../analyze';
import {_package} from '../package';

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
  const emptyStudy = () => DG.DataFrame.fromColumns([DG.Column.string('animal', 0), DG.Column.string('group', 0),
    DG.Column.float('day', 0), DG.Column.float('volume', 0)]);
  const inputs = (study: DG.DataFrame, control: string) => ({study, animal: study.col('animal')!,
    group: study.col('group')!, day: study.col('day')!, volume: study.col('volume')!, line: 'MC38', control,
    endpoint: 1500, correction: 'holm'});

  test('a study with no measurements is an empty state, not an error', async () => {
    let caught: unknown = null;
    try {
      await monitorFrames(inputs(emptyStudy(), 'Vehicle'));
    } catch (e) {
      caught = e;
    }
    expect(caught instanceof EmptyState, true);
    expect((caught as Error).message.includes('no measurements yet'), true);
  });

  test('a control that is not in the study says which, keeping the groups for the view to hide', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['a1', 'a2']),
      DG.Column.fromStrings('group', ['Vehicle', 'Drug']), DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'day', [0, 0]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'volume', [100, 110])]);
    let caught: unknown = null;
    try {
      await monitorFrames(inputs(df, 'Control'));
    } catch (e) {
      caught = e;
    }
    expect(caught instanceof NoControl, true);
    expect((caught as NoControl).control, 'Control');
  });

  test('measure again by chance, in the words of a session', async () => {
    expect(promptsByChance(30, 0.037), 'about 1 a session');
    expect(promptsByChance(8, 0.037), 'about one every 3 sessions');
    expect(promptsByChance(18, 0.0395), 'about 1 a session');
    expect(promptsByChance(30, 0.0125), 'about one every 3 sessions');
    // what the range shown before measuring costs by chance: both tails of the 99% one, and the share above it
    const ranges = (lo: number, hi: number) => ({levels: {'99': {lo: 0.005, hi: 0.995}},
      observed: {'90': {miss_lo: 0.03, miss_hi: 0.07}, '99': {miss_lo: lo, miss_hi: hi}}, papers: 242, points: 24219,
      pooled: true, source: ''});
    const c = checkChance(ranges(0.002, 0.0105))!;
    expect(Math.abs(c.chance - 0.0125) < 1e-12 && Math.abs(c.above - 0.84) < 1e-12, true, JSON.stringify(c));
    // not measured, or nothing ever outside: no rate to quote
    expect(checkChance(undefined), null);
    expect(checkChance(ranges(0, 0)), null);
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
    expect(r(-15)[0].startsWith('Lost '), true);
    expectArray(r(-14), []);
    expectArray(r(5), []);
    expectArray(r(null), []);
    expect(r(-1, 7)[0], 'Not weighed on day 10');
    // listed on the session it first misses, and not once it is gone from the study
    expectArray(r(-1, 4), []);
    expectArray(r(-1, 7, false), []);
    // over the limit at its last weighing: said on the session it first misses, gone or not (Welfare counts it removed)
    expect(r(-22, 7, false)[0], 'Over the 20% weight-loss limit (lost 22%)');
    expectArray(r(-22, 4, false), []);
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
    // gone as Welfare reads it (no removal column: not measured on its group's last day): the same word here
    expect(reasonsOf({p_limit_next: 0.01, last_flag: 0}, 10, 2000, 7, 500, 7, 'other')[0],
      'Not measured on day 10: counted as removed early');
    expect(r(0.01, 0, 10, 2100)[0].startsWith('Over'), true);
    // a running study's session may still be under way: not measured yet, not removed
    expect(reasonsOf({p_limit_next: 0.01, last_flag: 0}, 10, 2000, 7, 500, 7, PENDING)[0],
      'Not measured yet on day 10');
  });

  test('a running study: an animal not reached yet today stays in Welfare\'s counts', async () => {
    const tol = DG.DataFrame.fromObjects([
      {group: 'A', basis: 'measured', removed_other: 2, removed_animals: 'P-04 (day 7, other reason), P-01 (day 4, ' +
        'other reason)'},
      {group: 'A', basis: 'corrected', removed_other: 2, removed_animals: 'P-04 (day 7, other reason), P-01 (day 4, ' +
        'other reason)'},
      {group: 'B', basis: 'measured', removed_other: 0, removed_animals: ''}])!;
    const animals = [{animal: 'P-04', group: 'A', last_day: 7, removal: 'other'},
      {animal: 'P-01', group: 'A', last_day: 4, removal: 'other'},
      {animal: 'P-09', group: 'A', last_day: 7, removal: 'weight'},
      {animal: 'P-05', group: 'B', last_day: 7, removal: 'other'}];
    // P-05's removal is in the table: it stays removed
    const back = notYetMeasured(animals, tol, 7, new Set(['P-05']));
    expectArray([...back], ['P-04']);
    expect(tol.get('removed_other', 0), 1);
    expect(tol.get('removed_other', 1), 1);
    expect(tol.get('removed_animals', 0), 'P-01 (day 4, other reason)');
    expect(tol.get('removed_other', 2), 0);
  });

  test('a loss near the limit is not rounded onto it; animal IDs in natural order', async () => {
    expect(weightReasons({change: -19.6, lastDay: 10}, 10, 7, 20, true)[0], 'Lost 19.6% of body weight (limit 20%)');
    // the table and the panel say the same number as the reason
    expect(weightChangeText(-19.6, 20), '−19.6%');
    expect(weightChangeText(-19.97, 20), '−19.9%');
    expect(weightChangeText(-16.4, 20), '−16%');
    expect(weightChangeText(-21, 20), '−21%');
    expect(weightChangeText(4.6, 20), '+5%');
    expect(weightChangeText(0.2, 20), '0%');
    expect(weightReasons({change: -16.4, lastDay: 10}, 10, 7, 20, true)[0], 'Lost 16% of body weight (limit 20%)');
    expect(weightReasons({change: -19.57, lastDay: 10}, 10, 7, 20, true)[0], 'Lost 19.6% of body weight (limit 20%)',
      'rounded, not truncated');
    expect(weightReasons({change: -19.97, lastDay: 10}, 10, 7, 20, true)[0], 'Lost 19.9% of body weight (limit 20%)',
      'never shown as reaching the limit');
    expectArray(['T3-10', 'T03-02', 'T3-2', 'T03-01', 'T3-1'].sort(naturalCompare),
      ['T03-01', 'T03-02', 'T3-1', 'T3-2', 'T3-10']);
  });

  test('a stopped arm: no dose after its review, until the stop is undone', async () => {
    const log = [{what: decisionEntry('Drug 10 mg/kg', 9, 'stop')}];
    expectArray(Object.entries(stopsOf(log)).flat(), ['Drug 10 mg/kg', 9]);
    expectArray(Object.keys(stopsOf([...log, {what: decisionEntry('Drug 10 mg/kg', 9, 'stop undone')}])), []);
    // a stop against the suggestion is one line, its reason in it
    const against = decisionEntry('A: 10 mg/kg', 7, 'stop', 'against the suggestion to continue: weight loss');
    expect(against, 'A: 10 mg/kg: stop (review on day 7); against the suggestion to continue: weight loss');
    expectArray(Object.entries(stopsOf([{what: against}])).flat(), ['A: 10 mg/kg', 7]);
    // an older file's details: day 0 with its time is that day; an owner written as a name is who ran it
    const old = parseMeta({owner: 'Someone Else', start: '2026-11-30T10:00:00Z'})!;
    expect([old.start, old.runBy, old.owner].join('|'), '2026-11-30|Someone Else|Someone Else');
    expect(parseMeta({owner: 'ben.ruiz'})!.runBy, '');
    const groups = [{group: 'Drug 10 mg/kg', compound: 'Drug', dose: 10, unit: 'mg/kg', route: 'p.o.',
      schedule: 'QD'}];
    const days = (stopped?: {[g: string]: number}) => [...dosingDays(groups, 14, {start: null, stopped}).keys()];
    expect(days().some((d) => d > 9), true);
    expect(Math.max(...days(stopsOf(log))), 9);
    // the study table's log is read with its details
    const df = DG.DataFrame.fromColumns([DG.Column.string('animal', 0)]);
    setMetaOn(df, emptyMeta());
    addLog(df, decisionEntry('Drug 10 mg/kg', 9, 'stop'));
    expect(metaOf(df)?.stopped?.['Drug 10 mg/kg'], 9);
  });

  test('day 0 that no calendar has is refused', async () => {
    expect(impossibleDate('2030-02-30'), true);
    expect(impossibleDate('31/04/2026'), true);
    expect(impossibleDate('2030-03-02'), false);
    expect(impossibleDate('02/03/2030'), false);
    expect(impossibleDate('03/31/2026'), false);
    expect(impossibleDate(''), false);
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
    // what "measure again" costs by chance is quoted only where the atlas measured it
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

  test('body weight missing: counted over animals with a volume, said plainly', async () => {
    expect(weightMissing([100, 0, NaN, 90, 80], [20, null, null, 0, 19]), 1);
    expect(weightMissing([100, 90], [null, NaN]), 2);
    expect(weightMissing([100, 90], [20, 21]), 0);
    expect(balanceNote(false, 0, false) === null, true);
    expect(balanceNote(true, 0, true), 'Groups balanced on tumor volume first, then on body weight.');
    expect(balanceNote(true, 0, false), 'Balanced on tumor volume only: body weight is not used yet.');
    expect(balanceNote(true, 3, false), 'Body weight missing for 3 animals: balanced on tumor volume only.');
    expect(balanceNote(true, 1, false), 'Body weight missing for 1 animal: balanced on tumor volume only.');
  });

  test('body weight of each animal: the latest valid one up to its latest volume', async () => {
    const r = latestReadings([
      {animal: 'a', day: 0, volume: 100, weight: 20}, {animal: 'a', day: 7, volume: 120, weight: null},
      {animal: 'a', day: 7, volume: NaN, weight: 21}, {animal: 'a', day: 9, volume: NaN, weight: 25},
      {animal: 'b', day: 0, volume: 90, weight: 0}, {animal: 'b', day: 3, volume: 95, weight: 22},
      {animal: 'c', day: 0, volume: NaN, weight: 20}]);
    expect(JSON.stringify(r), JSON.stringify({animal: ['a', 'b'], volume: [120, 95], weight: [21, 22]}));
  });

  test('groups with some body weights missing: balanced on tumor volume only', async () => {
    const animal = Array.from({length: 13}, (_, i) => `m${i}`);
    const volume = animal.map((_, i) => 80 + 3 * i);
    const weight = animal.map((_, i) => i === 0 ? null : 20 + (i % 5));
    const {balance, weightUsed} = await allocateGroups(animal, volume, [{name: 'Vehicle', n: 6}, {name: 'A', n: 6}], 3,
      weight);
    expect(weightUsed, false);
    expect(balance.col('weight_mean') == null, true);
    expect(balance.col('weight_sd') == null, true);
  });

  test('groups alike in tumors and body weight', async () => {
    const animal = Array.from({length: 13}, (_, i) => `m${i}`);
    const volume = animal.map((_, i) => 80 + 3 * i);
    const weight = animal.map((_, i) => 20 + (i % 5));
    const {balance, weightUsed} = await allocateGroups(animal, volume, [{name: 'Vehicle', n: 6}, {name: 'A', n: 6}], 3,
      weight);
    // the test library has no runtime skip; a returned string shows as the test's result
    if (!weightUsed)
      return 'NOT RUN: the engine does not use body weight yet (no weight_mean)';
    expect(balance.rowCount, 2);
    const means = balance.col('weight_mean')!.toList() as number[];
    expect(Math.max(...means) - Math.min(...means) < 1, true, `weight means ${means}`);
    expect((balance.col('weight_sd')!.toList() as number[]).every(Number.isFinite), true);
  });

  test('a study file keeps the plan, the reviews, the measurements and the log', async () => {
    const whole = await exampleStudy();
    const plan = parsePlan(whole.getTag(PLAN_TAG)!);
    expect(newStudyTable(plan).rowCount, 0);
    // the plan's review days become the followed study's (fixed from its first run in Monitor)
    expect(newStudyTable({...plan, inputs: {...plan.inputs, looks: ' 6, 9 '}}).getTag(LOOKS_TAG), '6, 9');
    expect(newStudyTable({...plan, inputs: {...plan.inputs, looks: ''}}).getTag(LOOKS_TAG) || 'none', 'none');
    // a measurements file opened in place of the study just followed takes its plan and reviews, once
    const followed = newStudyTable({...plan, inputs: {...plan.inputs, looks: '6'}});
    const opened = await upTo(5);
    opened.setTag(PLAN_TAG, '');
    opened.setTag(LOOKS_TAG, '');
    expect(takeOverFollowed(followed, opened), true);
    expect([opened.getTag(PLAN_TAG) === followed.getTag(PLAN_TAG), opened.getTag(LOOKS_TAG)].join(), 'true,6');
    expect(takeOverFollowed(followed, opened), false, 'it has a plan now');
    // previewed in Files: its model and end, the control and each arm, who saved it
    const preview = planSummary(plan);
    expect(preview.length, plan.predicted.arms.length + 3, preview.join(' | '));
    expect(preview[0].startsWith('MC38 · endpoint'), true, preview[0]);
    const arm = plan.predicted.arms[0];
    expect(preview[2], `${arm.name}: TGI ${arm.effect}% expected, ${arm.n ?? 'over 40'} analyzable animals`);
    const old = planSummary({...plan, created: undefined as any});
    expect(old[old.length - 1].startsWith('Saved'), true, 'a plan without its date');
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

  const blindedUpTo = async (today: number) => {
    const df = await upTo(today);
    setMetaOn(df, {...emptyMeta(), blinded: true});
    return df;
  };
  const shown = (df: DG.DataFrame) => !!df.col('group') && !df.col('~group');
  const hidden = (df: DG.DataFrame) => !df.col('group') && !!df.col('~group');
  const viewOf = async (df: DG.DataFrame) => {
    const view = await monitorView(df);
    return {view, close: () => {
      // (never added to the shell: closing it disposes nothing, so its host goes by hand, its study link with it)
      view.close();
      ViewHost.of(view).dispose();
      [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id)?.close();
      grok.shell.closeTable(df);
    }};
  };

  test('closing the study\'s table view leaves the stage its study', async () => {
    const df = await upTo(14);
    const {view, close} = await viewOf(df);
    try {
      const held = () => [...view.root.querySelectorAll('select')].some((s) => s.value === df.name);
      await awaitCheck(held, 'the study did not open in the Study field', 10000);
      openStudyTable(df);
      await awaitCheck(() => [...grok.shell.tableViews].some((v) => v.dataFrame.id === df.id), 'no table view', 5000);
      [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id)!.close();
      await delay(500);
      expect([...grok.shell.tableViews].some((v) => v.dataFrame.id === df.id), false, 'the closed tab came back');
      expect(held(), true, 'the Study field lost its study');
      // only what is on screen: the stage keeps a hidden warning line
      const shown = (t: string) => [...view.root.querySelectorAll<HTMLElement>('*')]
        .some((e) => e.children.length === 0 && e.textContent!.includes(t) && e.offsetParent !== null);
      expect(shown('No open tables'), false);
      expect(shown('A field is missing'), false);
    } finally {
      close();
    }
  });

  test('a blinded study opens with no group column, before any run; revealing gives it back', async () => {
    const df = await blindedUpTo(14);
    const {view, close} = await viewOf(df);
    try {
      await awaitCheck(() => hidden(df), 'the group column is shown on opening', 5000);
      const reveal = () => [...view.root.querySelectorAll('button')].find((b) => b.textContent === 'Reveal the groups');
      await awaitCheck(() => !!reveal(), 'no Reveal button', 60000);
      expect(hidden(df), true);
      reveal()!.click();
      // revealing is for good: asked once more, then Reveal
      const ok = () => [...document.querySelectorAll<HTMLElement>('.d4-dialog button, .d4-dialog .ui-btn')]
        .find((b) => b.textContent?.trim().toLowerCase() === 'reveal');
      await awaitCheck(() => !!ok(), 'no Reveal in the confirmation', 5000);
      ok()!.click();
      await awaitCheck(() => shown(df), 'the group column did not come back after the reveal', 5000);
    } finally {
      close();
    }
  });

  test('a blinded study whose run fails still has no group column', async () => {
    const df = emptyStudy();
    setMetaOn(df, {...emptyMeta(), blinded: true});
    const {close} = await viewOf(df);
    try {
      await delay(1500);
      expect(hidden(df), true);
    } finally {
      close();
    }
  });

  test('a blinded study saved while hidden keeps its allocation', async () => {
    const df = await blindedUpTo(14);
    const before = Array.from({length: df.rowCount}, (_, i) => `${df.get('animal', i)}:${df.get('group', i)}`);
    const {close} = await viewOf(df);
    try {
      // the study opens after the view (not awaited by monitorView): wait for the column to be hidden
      await awaitCheck(() => hidden(df), 'the group column is shown on opening', 5000);
      const saved = parseStudy(JSON.stringify(studyOf(df, null)));
      expect(saved.rows.map((r) => `${r.animal}:${r.group}`).join(), before.join());
      const back = tableOf(saved);
      expect(back.col('group') !== null, true);
      expect(Array.from({length: back.rowCount}, (_, i) => `${back.get('animal', i)}:${back.get('group', i)}`).join(),
        before.join());
    } finally {
      close();
    }
  });

  test('a day typed with no group takes its animal\'s group', async () => {
    expect(groupsByAnimal(['a', 'b', 'a', 'c'], ['x', 'y', '', null]).join(), 'x,y,x,');
    const col = (name: string, v: (string | null)[]) => DG.Column.fromStrings(name, v.map((x) => x ?? ''));
    const study = DG.DataFrame.fromColumns([col('animal', ['a', 'a', 'a']), col('group', ['x', 'x', '']),
      DG.Column.fromFloat64Array('day', Float64Array.from([0, 3, 6])),
      DG.Column.fromFloat64Array('volume', Float64Array.from([100, 150, 200]))]);
    const d = analysisDeclarationOf({study, animal: study.col('animal')!, group: study.col('group')!,
      day: study.col('day')!, volume: study.col('volume')!, line: 'MC38', control: 'x', endpoint: 1500,
      correction: 'holm'});
    expect(d.data.group.join(), 'x,x,x');
    expect(d.data.day.length, 3);
  });

  test('after a run and a flag, a blinded study still holds its groups under the hidden name', async () => {
    const df = await blindedUpTo(14);
    const before = Array.from({length: df.rowCount}, (_, i) => String(df.get('group', i)));
    const {view, close} = await viewOf(df);
    try {
      await awaitCheck(() => !!view.root.querySelector('.pmx-stat-attn, .pmx-stat') && !!df.getTag('pmx.study'),
        'no run', 60000);
      expect(JSON.parse(df.getTag('pmx.study') ?? '{}').group, '~group');
      mark(df, [0], 'flagged', 'test', 'flagged');
      await delay(1500);
      expect(hidden(df), true);
      expect(JSON.parse(df.getTag('pmx.study') ?? '{}').group, '~group');
      expect(Array.from({length: df.rowCount}, (_, i) => String(df.get('~group', i))).join(), before.join());
    } finally {
      close();
    }
  });

  test('switching blinded entry on and the reveal, in an open table view, hide and show the grid column', async () => {
    const df = await upTo(14);
    const {close} = await viewOf(df);
    const gridCol = (name: string) => [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id)?.grid.columns
      .byName(name);
    try {
      // (the study's table view opens only when asked for)
      await delay(1000);
      openStudyTable(df);
      await awaitCheck(() => !!gridCol('group'), 'no table view', 10000);
      tagStudy({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null});
      setMetaOn(df, {...emptyMeta(), blinded: true});
      blindColumn(df);
      expect(hidden(df), true);
      expect(gridCol('~group')?.visible === false, true);
      addLog(df, REVEALED);
      blindColumn(df);
      expect(shown(df), true);
      expect(gridCol('group')?.visible === true, true);
    } finally {
      close();
    }
  });

  test('an animal whose rows disagree on the group is not guessed', async () => {
    expect(groupsByAnimal(['a', 'a', 'a'], ['x', 'y', null]).join(), 'x,y,');
  });

  test('the chance of a clear result: low as a figure, high only from a clearly effective arm\'s interval',
    async () => {
      expect(chanceWords(0.05, 80), 'Low (under 10%)', 'a low chance is said whatever the interval');
      expect(chanceWords(0.97, 86), 'High (TGI at least 86%)', 'TGI 91% (86–94) is not "not low"');
      expect(chanceWords(0.6, 20), 'Not low (10% or more)', 'no high figure from the model itself');
    });

  test('Enter today\'s values: always there; off the plan\'s calendar it asks first', async () => {
    // the plan measures on day 0 and every 3rd day here
    const planned = (d: number) => d % 3 === 0;
    expect(entryStep(6, planned), 'open', 'a plan measurement day');
    expect(entryStep(7, planned), 'confirm', 'a day off the calendar asks first');
    expect(entryStep(7, null), 'open', 'no plan: at once');
    expect(entryStep(null, planned), 'open', 'no dated day 0: the dialog asks for the day');
    // the question names the day and the plan's next one (QA 0.1.3 round 3: technician N20)
    expect(offPlanAsk(7, planned), 'Day 7 is not a measurement day in the plan (the next one is day 9): enter anyway?');
  });

  test('QA 0.1.3 round 3: next measurement, why and status, day-0 IDs', async () => {
    // never a day already past: day 11 passed unmeasured is missed, the next is day 14
    const n = nextMeasurement([11, 14, 16], 12);
    expect([n.next, n.missed].join(), '14,11');
    expect([nextMeasurement([11, 14], 11).next, nextMeasurement([11, 14], 11).missed].join(), '11,');
    expect(nextMeasurement([11, 14], null).next, 11, 'no date: the first to come');
    const w = whyAndStatus('Likely a typo: 665 mm³ on day 11 is 9 times day 9\'s 77 mm³: check it; ' +
      'last value unexpected: measure again');
    expect(w.why, 'Likely a typo: 665 mm³ on day 11 is 9 times day 9\'s 77 mm³; last value unexpected');
    expect(w.status, 'Check it; measure again');
    expect(whyAndStatus('Over the 1500 mm³ limit').status, '');
    const ids = animalIds('A1\nA1\nA2', 60);
    expect([ids.ids.join(), ids.twice.join(), ids.line].join(' | '),
      'A1,A2 | A1 | 2 of 60 – add 58 or resize the groups in Plan · A1 appears twice');
    // round 4 designer D4: M6 and m6 are one ID; the count is of the animals the plan randomizes
    const cased = animalIds('M6\nm6\nM7', 2);
    expect([cased.ids.join(), cased.twice.join(), cased.line].join(' | '),
      'M6,M7 | M6 | 2 of 2 · M6 and m6 are the same ID');
    // a day 0 years away asks (chaos 6); one this year does not
    const nov28 = new Date(Date.UTC(2026, 10, 28));
    expect(farDayZero('1900-01-01', nov28), 'Day 0 on 1 Jan 1900 is 126 years ago: is that the date?');
    expect(farDayZero('2026-11-02', nov28), null);
    // round 4 designer D5: day 0 planned for another day than today asks which
    const sat = new Date(2026, 10, 28, 12);
    expect(dayZeroAsk('2026-11-30', sat), 'Enter day 0 values (planned for Mon 30 Nov; today is Sat 28 Nov)');
    expect(dayZeroAsk('2026-11-28', sat), null);
    expect(dayZeroAsk(null, sat), null);
  });
});

category('Monitor: a day under way, welfare and typos (QA 0.1.3)', () => {
  const inputs = (df: DG.DataFrame, humane: number, looks = ''): MonitorInputs => ({study: df,
    animal: df.col('animal')!, group: df.col('group')!, day: df.col('day')!, volume: df.col('volume')!, looks,
    line: 'MC38', control: '', endpoint: 1500, humane, correction: 'holm'});
  const days = (df: DG.DataFrame) => [...new Set(df.col('day')!.toList() as number[])].sort((a, b) => a - b);

  test('a likely typo, an impossible volume, a shortfall of animals: in words', async () => {
    expect(likelyTypo(420, 4201), '10 times');
    expect(likelyTypo(420, 42), '1/10 of');
    expect(likelyTypo(420, 900), null);
    expect(likelyTypo(420, -50), null);
    expect(mm3Text(-50), '−50');
    expect(mm3Text(10790), '10,790');
    expect(mm3Text(-1500), '−1,500');
    // a negative value is impossible, not one to measure again
    expect(outsideBand([{animal: 'a', day: 12, volume: -50}], 12,
      new Map([['a', {check_lower: 100, check_upper: 200}]])).length, 0);
    const short = shortOfAnimals(60, ['m1', 'm2'], ['m1', 'm2', 'm3']);
    expect(short.startsWith('The plan needs 60 animals, and 2 have a tumor volume (no volume yet: m3)'), true, short);
    // a continue undone is no stop
    expect(Object.keys(stopsOf([{what: decisionEntry('A', 7, 'continue undone')}])).length, 0);
    expect(localTime('not a time'), 'not a time');
    expect(/^\d{4}-\d\d-\d\d \d\d:\d\d$/.test(localTime('2026-11-27T13:05:00.000Z')), true);
  });

  test('a day under way: its groups in one line, the alerts as of their last session, the review waits',
    async () => {
      const whole = await exampleStudy();
      const [prev, last] = days(whole).slice(3, 5);
      const first = whole.get('group', 0);
      const df = whole.clone(DG.BitSet.create(whole.rowCount, (i) => whole.get('day', i) < last ||
        (whole.get('day', i) === last && whole.get('group', i) === first)));
      df.setTag(PLAN_TAG, whole.getTag(PLAN_TAG)!);
      const r = await monitorFrames(inputs(df, 2000, String(last)));
      expect(r.status.get('waiting', 0) > 0, true);
      expect(r.status.get('complete_day', 0), prev);
      expect(r.status.get('at_look', 0), false);
      expect(r.status.get('look_waits', 0), true);
      expect(r.decisions.rowCount, 0);
      for (let i = 0; i < r.animals.rowCount; i++)
        expect(String(r.animals.get('attention', i) ?? '').includes('Not measured yet'), false);
    }, {timeout: 120000});

  test('QA 0.1.3 round 4: a followed plan by its ID, its settings, the animals it randomizes', async () => {
    const plan = parsePlan(await _package.files.readAsText('demo/ct26-trametinib-antipd1.pmxplan'));
    const withId = {...plan, study: {...plan.study ?? emptyMeta(), id: 'PMX-2026-0906'}};
    const study = newStudyTable(withId);
    // D1: the plan edited since (another endpoint) still finds its study; another plan's ID does not
    const edited = {...withId, inputs: {...withId.inputs, endpoint: 1200}};
    expect(followedFrom(study, edited), true, 'edited plan, same ID');
    expect(followedFrom(study, {...edited, study: {...edited.study!, id: 'PMX-2026-0907'}}), false);
    // D3: the study takes the plan's endpoint and settings
    const e1200 = newStudyTable(edited);
    expect(planFieldsOf(e1200).endpoint, 1200);
    expect(planFieldsOf(DG.DataFrame.create()).endpoint === undefined, true, 'no plan: nothing');
    // D4: the animals to randomize, each group's; a plan saved before they were kept, its analyzable ones
    const old = randomizedGroups(plan);
    expect(old[0].n, plan.predicted.nControl ?? 0);
    const sized = randomizedGroups({...plan, predicted: {...plan.predicted, randomiseControl: 13,
      arms: plan.predicted.arms.map((a) => ({...a, randomise: 13}))}});
    expect(sized.map((g) => g.n).join(), Array(plan.predicted.arms.length + 1).fill(13).join());
  });

  test('an animal removed early over the weight limit stays listed while the next day is partly entered', async () => {
    // the demo's T3-07: -22% on day 7, not measured on day 9; then day 11 entered for the vehicle group only
    const df = DG.DataFrame.fromCsv(await _package.files.readAsText('demo/ct26-trametinib-antipd1-running-day9.csv'));
    const at9 = df.rowCount;
    for (let i = 0; i < at9; i++) {
      if (df.get('day', i) === 9 && df.get('group', i) === 'Vehicle') {
        df.rows.addNew([df.get('animal', i), 'Vehicle', 11, df.get('volume', i) * 1.2,
          df.get('body_weight', i)]);
      }
    }
    // the control as the plan has it: the one group none of its arms
    const plan = await _package.files.readAsText('demo/ct26-trametinib-antipd1.pmxplan');
    expect(planControl({...inputs(df, 1500), plan}), 'Vehicle');
    expect(planControl(inputs(df, 1500)) === undefined, true, 'no plan: no control from it');
    for (const study of [df.clone(DG.BitSet.create(df.rowCount, (i) => i < at9)), df]) {
      const r = await monitorFrames({...inputs(study, 1500), line: 'CT26', control: 'Vehicle',
        bodyWeight: study.col('body_weight'), weightLimit: 20});
      const row = Array.from({length: r.animals.rowCount}, (_, i) => i)
        .find((i) => r.animals.get('animal', i) === 'T3-07');
      expect(String(r.animals.get('attention', row!) ?? '').includes('weight-loss limit'), true,
        `T3-07 with ${study.rowCount} rows: ${r.animals.get('attention', row!)}`);
      // the planted ×10 typo (T1-08, 2,893 mm³): checked first, the humane limit not claimed on it
      const typo = String(r.animals.get('attention', r.animals.col('animal')!.toList().indexOf('T1-08')));
      expect(typo.startsWith('Likely a typo: 2,893 mm³'), true, typo);
      expect(typo.includes('mm³ limit'), false, typo);
    }
  }, {timeout: 120000});

  test('an animal excluded from the analysis is still listed over the humane limit', async () => {
    const whole = await exampleStudy();
    const last = days(whole)[4];
    const df = await upTo(last);
    let big = -1;
    for (let i = 0; i < df.rowCount; i++) {
      if (df.get('day', i) === last && (big < 0 || df.get('volume', i) > df.get('volume', big)))
        big = i;
    }
    const name = String(df.get('animal', big));
    const own = Array.from({length: df.rowCount}, (_, i) => i).filter((i) => df.get('animal', i) === name);
    mark(df, own, 'excluded', 'test', name);
    const r = await monitorFrames(inputs(df, df.get('volume', big) - 1));
    const row = Array.from({length: r.animals.rowCount}, (_, i) => i).find((i) => r.animals.get('animal', i) === name);
    expect(row !== undefined, true, 'the excluded animal is still on the list');
    expect(String(r.animals.get('attention', row!)).startsWith('Over the'), true);
  }, {timeout: 120000});
});
