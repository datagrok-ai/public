import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {AnalysisFrames, OBSERVED, PREDICTED, REANCHOR, controlsTo, doublingOf, endRange, endVerdict, intervalNote}
  from '../analyze';
import {HUMANE, LAB_WEEKDAYS, designStudy} from '../engine';
import {monitorEnd, planChecks, widenedEnd} from '../plan-checks';
import {bandFactors, humanReference} from '../priors';
import {drugKey} from '../arm-name';
import {exampleStudy} from './example-study';
import {PLAN_TAG, parsePlan, powerAt} from '../plan-file';
import {noEffectNote} from '../plan-vs-observed';
import {competingOf, projectedRows} from '../survival';
import {EXAMPLE, percentRange, plainError} from '../ui';
import {cellsOf} from '../wide-table';
import {readStudy} from '../study-reader';
import {_package} from '../package';
import {EngineError} from '../crux';
import {bandWords, combinationStep, independenceBand, notCombined, stepWords, synergyAlpha, synergyOf, synergyVerdict,
  synergyWords} from '../synergy';
import {arriveItems, receivesLine} from '../arrive';
import {ctdRows, ctdStudyOf} from '../ctd';
import {prepareItems} from '../prepare';
import {tumorModels} from '../priors';
import {ArmExposure, CurveRow, DrugPk, Exposure, HumanExposure, NcaByArm, Proportional, ReadoutGroup, ReferenceRow,
  aucFactor, concentrationFactor, concentrationMapping, curveSilence, dosesPerDay, drugMultiples, exposureAxis,
  exposureByArm, exposureByDrug, ncaByArm,
  exposureCurve, exposureHint, exposureJson, exposureRange, fromReference, humanMultiple, humansInto, medianDay,
  medianTgi, multipleText, ncaInstalled, parseExposure, patientsDose, patientsDoseOf, readout, referenceKept,
  referenceLine, referenceRow, renamedExposure, runNca, sourceText}
  from '../exposure';
import {MISSING, groupStatus, idMapCsv, modellingCsv, modellingDataset, noConcentrations} from '../modelling-dataset';
import {doseCount, readSchedule} from '../schedule';
import {EXCLUDED, EXPOSURE_TAG, FLAGGED, exposureOfTag, parseStudy, studyOf as studyFileOf, tableOf}
  from '../study-file';
import {armParts, detailsName, emptyMeta, ofGroups, parseMeta, renamedDetails, renamedGroups, withGroups}
  from '../study-meta';
import {groupDetails} from '../translation-report';
import type {StudyInputs} from '../explore';

const OUTPUTS = ['analysis', 'comparisons', 'animals', 'volumes', 'twins', 'effects', 'checks', 'versusPlan',
  'predicted', 'tolerability', 'weights'];

/** pmxAnalyzeStudy as the Analyze view calls it: the outputs by name. */
async function analyze(study: DG.DataFrame, params: {[k: string]: any} = {}): Promise<{[k: string]: DG.DataFrame}> {
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxAnalyzeStudy'})[0].prepare({study,
    animal: study.col('animal'), group: study.col('group'), day: study.col('day'), volume: study.col('volume'),
    line: 'MC38', endpoint: 1500, correction: 'holm', ...params});
  await fc.call();
  return Object.fromEntries(OUTPUTS.map((o) => [o, fc.outputs[o] as DG.DataFrame]));
}

category('Analyze', () => {
  test('the example (Yu 2024) reads the same TGI and p by its long file and as a column per day', async () => {
    const long = cellsOf(await _package.files.readAsText(EXAMPLE.file));
    // the long file as the app reads it, its body weight too; the same rows read as a column per day, no weight
    const days = [...new Set(long.slice(1).map((r) => r[2]))];
    const at = new Map(long.slice(1).map((r) => [`${r[0]}|${r[2]}`, r]));
    const wide = [['Animal', 'Group', ...days.map((d) => `Day ${d} (mm3)`)],
      ...[...new Set(long.slice(1).map((r) => r[0]))].map((a) =>
        [a, at.get(`${a}|${days[0]}`)![1], ...days.map((d) => at.get(`${a}|${d}`)?.[3] ?? '')])];
    const fromWide = DG.DataFrame.fromObjects(readStudy(wide).rows.map((r) => ({animal: r.animal, group: r.group,
      day: r.day, volume: r.volume})))!;
    const fromLong = DG.DataFrame.fromCsv(long.map((r) => r.join(',')).join('\n'));
    const result = async (study: DG.DataFrame, weight: boolean) => {
      const r = await analyze(study, {endpoint: 1000, control: 'Vehicle', withTwin: false,
        ...(weight ? {bodyWeight: study.col('body_weight')} : {})});
      const c = r.comparisons;
      return `day ${r.analysis.get('evaluation_day', 0)}; ` + Array.from({length: c.rowCount}, (_, i) =>
        `${c.get('arm', i)}: TGI ${Math.round(c.get('tgi', i))} p ${c.get('p_value', i).toFixed(3)}`).sort().join('; ');
    };
    const [byLong, byWide] = [await result(fromLong, true), await result(fromWide, false)];
    expect(byWide, byLong, 'the same study by either route');
    // QA W4 (build b0d6ecf4c0) noted anti-PD-1 p 0.22 off the screen: neither this code nor the engine nor the file has
    // changed since, and the test is deterministic, so it was a misreading (D702 c)
    expect(byLong, 'day 12; LD4172 20 mg/kg i.p. QD + anti-PD-1 100 ug i.p. Q3D: TGI 79 p 0.013; ' +
      'LD4172 20 mg/kg i.p. QD: TGI 18 p 0.291; anti-PD-1 100 ug i.p. Q3D: TGI 17 p 0.275', 'the example\'s numbers');
  }, {timeout: 120000});

  test('the predicted controls run to the evaluation day, between two of the plan\'s days', async () => {
    const plan = [0, 2, 4, 7].map((day) => ({day, q025: 50 * 2 ** (day / 3), median: 100 * 2 ** (day / 3),
      q975: 200 * 2 ** (day / 3)}));
    const band = controlsTo(plan, 6);
    expect(band.map((b) => b.day).join(), '0,2,4,6');
    expect(Math.abs(band[3].median - 400) < 1e-9, true, 'geometric between day 4 and day 7');
    expect(controlsTo(plan, 4).length, 3, 'a plan\'s own day: nothing added');
    expect(controlsTo(plan, 9).length, 4, 'past the last day: nothing made up');
  });

  test('without the twin, the same tests and effects, its tables empty', async () => {
    const study = await exampleStudy();
    const [all, lean] = [await analyze(study), await analyze(study, {withTwin: false})];
    expect(lean.analysis.get('twinned', 0), false);
    expect(all.analysis.get('twinned', 0), true);
    for (const name of ['twins', 'effects', 'checks']) {
      expect(lean[name].rowCount, 0, name);
      expect(all[name].rowCount > 0, true, name);
    }
    for (const name of ['animals', 'volumes', 'versusPlan', 'predicted'])
      expect(lean[name].toCsv(), all[name].toCsv(), name);
    for (const c of ['tgi', 'tc', 'p_value', 'slope_difference', 'activity'])
      expect(lean.comparisons.col(c)!.toList().join(), all.comparisons.col(c)!.toList().join(), c);
  }, {timeout: 120000});

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
    const [vehicle, strong] = ['Vehicle', 'Cpd-A 30 mg/kg'].map((g) => row(g, 'Measured'));
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
    for (const g of ['Vehicle', 'Cpd-A 10 mg/kg', 'Cpd-A 30 mg/kg'].map((x) => row(x, 'Measured'))) {
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

  test('doubling time from the growth rate; growth delay against the control', async () => {
    expect(Math.abs(doublingOf(Math.LN2 / 4)! - 4) < 1e-12, true);
    expect([doublingOf(0), doublingOf(-0.1), doublingOf(null)].join('|'), '||', 'no doubling unless it grows');
    const o = await analyze(await exampleStudy());
    for (let i = 0; i < o.animals.rowCount; i++) {
      const [s, d] = [o.animals.get('slope_per_day', i), o.animals.get('doubling_time', i)];
      if (s > 0)
        expect(Math.abs(d * s - Math.LN2) < 1e-9, true, `animal ${i}: ${s} ${d}`);
    }
    // the stronger arm reaches the endpoint later, or never (no delay then)
    const c = o.comparisons;
    const [weak, strong] = [c.get('growth_delay_days', 0), c.get('growth_delay_days', 1)];
    expect(c.col('growth_delay_days') !== null, true);
    expect(!Number.isFinite(strong) || !Number.isFinite(weak) || strong >= weak, true, `delays ${weak} / ${strong}`);
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

  test('the corrections order the p-values: none ≤ Holm ≤ Bonferroni; the interval is Bonferroni\'s', async () => {
    const study = await exampleStudy();
    const of = async (correction: string) => (await analyze(study, {correction})).comparisons;
    const [none, holm, bonferroni] = [await of('none'), await of('holm'), await of('bonferroni')];
    const p = (t: DG.DataFrame, i: number) => t.get('p_value', i) as number;
    for (let i = 0; i < none.rowCount; i++) {
      expect(p(none, i) <= p(holm, i) + 1e-12, true, `arm ${i}: ${p(none, i)} vs ${p(holm, i)}`);
      expect(p(holm, i) <= p(bonferroni, i) + 1e-12, true, `arm ${i}: ${p(holm, i)} vs ${p(bonferroni, i)}`);
      // 0 outside the interval at 1 - alpha / k says what Bonferroni's p-value says
      const out = holm.get('slope_difference_lower', i) > 0 || holm.get('slope_difference_upper', i) < 0;
      expect(out, p(bonferroni, i) < 0.05, `arm ${i}: interval and Bonferroni p ${p(bonferroni, i)}`);
      // and beside the test, what one variance for every group gives
      const common = holm.get('p_common_variance', i);
      expect(common >= 0 && common <= 1, true, `arm ${i}: p with a common variance ${common}`);
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
      // (to the column's single precision)
      expect(Math.abs(x('power') - powerAt(arm.power!, x('n'))!) < 1e-6 && x('power') >= 0.79, true,
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

  test('plan vs observed: a study randomized larger than planned, its controls predicted from its own start',
    async () => {
      const study = await exampleStudy();
      const plan = parsePlan(study.getTag(PLAN_TAG)!);
      const bandStart = (o: {[k: string]: DG.DataFrame}) => {
        const p = o.predicted;
        const i = Array.from({length: p.rowCount}, (_, k) => k)
          .find((k) => p.get('series', k) === PREDICTED && p.get('day', k) === 0)!;
        return p.get('median', i);
      };
      const notes = (o: {[k: string]: DG.DataFrame}) => String(o.analysis.get('plan_notes', 0));
      // as planned (the example's start near the plan's): the plan's own band, nothing said
      const asPlanned = await analyze(study);
      const planned = plan.predicted.checks.controls[0].median;
      expect(Math.abs(bandStart(asPlanned) - planned) < 1e-3, true, `${bandStart(asPlanned)} vs ${planned}`);
      expect(notes(asPlanned).includes('Randomized at'), false, notes(asPlanned));
      // every tumor 1.5 times larger: the band starts at the study's own start, the endpoint comes earlier
      const larger = study.clone();
      larger.setTag(PLAN_TAG, study.getTag(PLAN_TAG)!);
      larger.col('volume')!.init((i) => study.col('volume')!.get(i) * 1.5);
      const start = bandStart(asPlanned) * 1.5;
      const o = await analyze(larger);
      expect(Math.abs(Math.log(bandStart(o) / start)) < Math.log(REANCHOR), true, `band ${bandStart(o)} vs ${start}`);
      expect(o.analysis.get('end_median', 0) < plan.predicted.checks.end.median, true,
        `end ${o.analysis.get('end_median', 0)} vs the plan's ${plan.predicted.checks.end.median}`);
      expect(notes(o).includes('Randomized at'), true, notes(o));
    }, {timeout: 300000});

  test('the end range: Analyze and Monitor never call a day Plan printed earlier or later than predicted', async () => {
    const hull: string[] = [];
    for (const line of ['MC38', 'CT26', 'B16F10']) {
      const f = await bandFactors(line);
      const s = await designStudy({line, vRand: 100, endpoint: 1500, power: 0.8, alpha: 0.05, correction: 'holm',
        n: null, arms: [{name: 'x', therapy: 'small_molecule', metric: 'TGI', effect: 60}], preLoss: 0, postLoss: 0,
        lossRule: 'confident', confidence: 0.9, weekdays: LAB_WEEKDAYS, startWeekday: 0, humane: HUMANE,
        schedules: false, lastDay: null});
      const checks = await planChecks({line, vRand: 100, endpoint: 1500, humane: HUMANE, weekdays: LAB_WEEKDAYS,
        startWeekday: 0, endpointDay: s.endpointDay, nControl: 10, arms: []});
      const r = endRange(checks.end, f, s, false);
      const [q10, q90] = [Math.round(s.durationQ10), Math.round(s.durationQ90 ?? s.durationQ10)];
      let needed = 0;
      for (let day = q10; day <= q90; day++) {
        expect(endVerdict(day, day, r), 'As predicted',
          `${line} day ${day} in ${q10}-${q90} against ${JSON.stringify(r)}`);
        if (day < widenedEnd(checks.end, f).q025 || day > (widenedEnd(checks.end, f).q975 ?? Infinity))
          needed++;
      }
      hull.push(`${line}: ${needed} of ${q90 - q10 + 1} days needed the hull`);
      // Monitor shows a range with Plan's printed one inside it
      const m = monitorEnd({median: checks.end.median, q10: checks.end.q025, q90: checks.end.q975}, f, s);
      expect(m.q10 <= q10 && (m.q90 === null || m.q90 >= q90), true, `${line} Monitor ${JSON.stringify(m)}`);
    }
    console.log(hull.join('; '));
  }, {timeout: 300000});

  test('the end range: widened by the line factor, open stays open, no factor unchanged, reanchored without hull',
    async () => {
      const r = {q025: 8, median: 15, q975: 25 as number | null};
      const f = await bandFactors('MC38');
      const w = widenedEnd(r, f);
      // MC38's 1.21 in log days around the median, 15·(8/15)^1.21 = 7.01 to 15·(25/15)^1.21 = 27.8, rounded outward
      expect(JSON.stringify(w), JSON.stringify({q025: 7, median: 15, q975: 28}));
      const wide = widenedEnd({q025: 8.4, median: 15, q975: 25.4}, f);
      expect(Number.isInteger(wide.q025) && Number.isInteger(wide.q975!) && wide.q025 <= 8 && wide.q975! >= 25, true);
      expect(widenedEnd({...r, q975: null}, f).q975, null);
      expect(await bandFactors('4T1'), null);
      expect(widenedEnd(r, await bandFactors('4T1')), r);
      // a printed range wider than the widened one is kept; a reanchored study is widened only
      const p = {durationQ10: 5, durationQ90: 40 as number | null};
      expect(endRange(r, f, p, false).q025, 5);
      expect(endRange(r, f, p, true).q025 === w.q025 && endRange(r, f, p, true).q975 === w.q975, true);
      expect(endRange(r, null, p, false), r);
      // Monitor: the same factor and the same hull; no factor, no change
      const m = {median: 15, q10: 8, q90: 25 as number | null};
      expect(monitorEnd(m, f, null).q10 === w.q025 && monitorEnd(m, f, null).q90 === w.q975, true);
      expect(monitorEnd(m, f, null).q10 < m.q10 + 1e-9 && w.q975! > m.q90!, true, 'widened, not unchanged');
      expect(monitorEnd(m, f, p).q90, 40);
      // an open printed top stays open; a printed range that is not a number is not used; no plan, no hull
      expect(endRange(r, f, {durationQ10: 5, durationQ90: null}, false).q975, null);
      expect(monitorEnd(m, f, {durationQ10: 5, durationQ90: null}).q90, null);
      expect(JSON.stringify(endRange(r, f, {durationQ10: NaN, durationQ90: 40}, false)), JSON.stringify(w));
      expect(JSON.stringify(endRange(r, f, null, false)), JSON.stringify(w));
      // a median or q10 that is not a number: Monitor's range as it was
      expect(JSON.stringify(monitorEnd({median: NaN, q10: 8, q90: 25}, f, p)), JSON.stringify({q10: 8, q90: 25}));
      expect(monitorEnd({median: 15, q10: null as any, q90: 25}, f, p).q10, null);
      // every printed (rounded) day gets the verdict the printed range gives, with no plan and reanchored
      for (const re of [false, true]) {
        const e = endRange(r, f, null, re);
        for (let day = Math.round(e.q025); day <= Math.round(e.q975!); day++)
          expect(endVerdict(day, day, e), 'As predicted', `day ${day} reanchored ${re}`);
      }
      expect(JSON.stringify(monitorEnd(m, null, p)), JSON.stringify({q10: 8, q90: 25}));
    });

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
    // no plan: a sham group, half of the example's controls, against the other half: every other animal by its last
    // volume, so the halves are alike whatever example the atlas gives (every other by id is a chance split, and
    // one in twenty of those is significant)
    const study = await exampleStudy();
    study.setTag(PLAN_TAG, '');
    const sham = study.clone(DG.BitSet.create(study.rowCount, (i) => study.get('group', i) === 'Vehicle'));
    sham.setTag(PLAN_TAG, '');
    const last = new Map<string, {day: number, volume: number}>();
    for (let i = 0; i < sham.rowCount; i++) {
      const [a, day] = [sham.get('animal', i) as string, sham.get('day', i) as number];
      if (!(last.get(a)?.day! >= day))
        last.set(a, {day, volume: sham.get('volume', i)});
    }
    const ids = [...last.keys()].sort((a, b) => last.get(a)!.volume - last.get(b)!.volume);
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
  test('projected KM: each animal carried on at its last rate, to the endpoint or the horizon', async () => {
    const r = (animal: string, day: number, volume: number) => ({animal, group: 'G', day, volume});
    const rows = [r('A', 0, 100), r('A', 2, 200), r('A', 4, 400), r('B', 0, 100), r('B', 2, 80), r('B', 4, 60),
      r('C', 0, 900), r('C', 2, 2000), r('D', 0, 100)];
    const near = (x: number, y: number) => Math.abs(x - y) < 1e-6;
    const out = projectedRows(rows, 1600, 20);
    const a = out.find((x) => x.animal === 'A')!;
    expect(near(a.day, 8) && a.volume === 1600, true, 'A doubles every 2 days: 1600 on day 8');
    const b = out.find((x) => x.animal === 'B')!;
    expect(b.day === 20 && near(b.volume, 60), true, 'B shrinks: followed to the horizon as it is');
    expect(out.some((x) => x.animal === 'C' || x.animal === 'D'), false, 'past the endpoint, or one measurement');
    const short = projectedRows(rows, 1600, 6).find((x) => x.animal === 'A')!;
    expect(short.day === 6 && near(short.volume, 800), true, 'A by an earlier horizon: under the endpoint');
  });
  test('an example study for every tumor model, lines without a density of their own included', async () => {
    const arm = {name: 'Cpd 30 mg/kg', therapy: 'small_molecule' as const, metric: 'TGI' as const, effect: 70};
    const failed: string[] = [];
    await Promise.all(Object.keys(await tumorModels()).map((line) => exampleStudy({line, arms: [arm]})
      .then((df) => df.rowCount > 0 || failed.push(`${line}: no rows`),
        (e) => failed.push(`${line}: ${e?.message ?? e}`))));
    expect(failed.join('; '), '');
  }, {timeout: 180000});
  // Gompertz growth (0.25 − 0.03 ln V per day) less each animal's kill rate, 10 per group, days 0 to 21: A and B
  // 0.04 and 0.06, A+B `abKill`
  const combinationStudy = (abKill: number, perGroup = false) => {
    const kills: {[g: string]: number} = {'Control': 0, 'A': 0.04, 'B': 0.06, 'A+B': abKill};
    const rows: {[k: string]: any}[] = [];
    let k = 0;
    for (const [g, kill] of Object.entries(kills)) {
      for (let a = 0; a < 10; a++, k++) {
        const v0 = Math.log(100) + 0.2 * Math.sin(k);
        const kappa = kill * (1 + 0.3 * Math.sin(2 * k));
        for (let d = 0; d <= 21; d += 3) {
          const e = Math.exp(-0.03 * d);
          const lv = v0 * e + (0.25 / 0.03 - kappa / 0.03) * (1 - e) + 0.05 * Math.sin(5 * k + d);
          // a tumor under what calipers read recorded as 0, as studies do
          rows.push({animal: perGroup ? `${a + 1}` : `${g}-${a}`, group: g, day: d,
            volume: Math.exp(lv) < 3.4 ? 0 : Math.exp(lv)});
        }
      }
    }
    return {study: DG.DataFrame.fromObjects(rows)!, animal: 'animal', group: 'group', day: 'day', volume: 'volume',
      bodyWeight: null};
  };
  const combinationGroups = {control: 'Control', a: 'A', b: 'B', ab: 'A+B'};
  test('a combination on the ladder: they cancel, more than the two together', async () => {
    // A+B none (it cancels them), twice the two added up (no tumor near the floor) or far more
    const [study, groups] = [combinationStudy, combinationGroups];
    const cancel = await combinationStep(study(0), groups, 0.03);
    expect(cancel?.step, 'antagonism');
    expect(stepWords(cancel!.step), 'Worse than the better drug alone');
    const more = await combinationStep(study(0.2), groups, 0.03);
    expect(more?.step, 'beyond_independent_action', `floored ${more?.floored}`);
    // so strong its tumors go (0): more than the two together cannot be seen, so it reads better than each
    const gone = await combinationStep(study(0.6), groups, 0.03);
    expect(`${gone?.step} ${gone?.floored}`, 'better_than_each true');
  }, {timeout: 60000});
  test('Monitor\'s band of independence: kills less, as independent, kills more', async () => {
    // A+B none, as B (the better drug, about what independent action gives), far more
    const at = async (abKill: number) => (await independenceBand(combinationStudy(abKill), combinationGroups, 0.03))!;
    const none = await at(0);
    expect(`${none.verdict} ${none.day} ${none.days.length}`, 'less 21 7');
    expect(none.lower.every((x, j) => x <= none.mid[j] && none.mid[j] <= none.upper[j]), true);
    expect(bandWords('A+B', none), 'Day 21: A+B slows growth less than its two drugs acting independently would');
    expect((await at(0.06)).verdict, 'inside');
    expect((await at(0.2)).verdict, 'more');
    // animals numbered from 1 in each group, as studies do: the same
    expect((await independenceBand(combinationStudy(0.2, true), combinationGroups, 0.03))?.verdict, 'more');
    expect((await combinationStep(combinationStudy(0.2, true), combinationGroups, 0.03))?.step,
      'beyond_independent_action');
  }, {timeout: 60000});
  test('synergy, competing risks and ARRIVE on a made-up study', async () => {
    // four groups growing at known rates (ln V per day): A and B take 0.05 each off the control's 0.2, A+B 0.15
    const rates: {[g: string]: number} = {'Control': 0.2, 'A': 0.15, 'B': 0.15, 'A+B': 0.05};
    const rows: {[k: string]: any}[] = [];
    let k = 0;
    for (const [g, rate] of Object.entries(rates)) {
      for (let a = 0; a < 6; a++, k++) {
        const v0 = 100 * (1 + 0.1 * Math.sin(k));
        for (let d = 0; d <= 20; d += 2) {
          const v = v0 * Math.exp(rate * d + 0.02 * Math.sin(3 * k + d));
          // a treated animal removed (ulceration) on day 10, its last measurement
          const removed = g === 'A' && a < 2;
          if (removed && d > 10)
            break;
          rows.push({animal: `${g}-${a}`, group: g, day: d, volume: v, why: removed && d === 10 ? 'ulceration' : ''});
        }
      }
    }
    const df = DG.DataFrame.fromObjects(rows)!;
    const s = {study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null};
    const x = await synergyOf(s, {control: 'Control', a: 'A', b: 'B', ab: 'A+B'});
    expect(Math.abs(x.dif - 0.05) < 0.01, true, `dif ${x.dif}`);
    expect(x.p < 0.05 && synergyWords(x, 'A+B', 0.05).startsWith('Synergy'), true, `p ${x.p}`);
    // the dialog's verdict: more, about the same, less than the two drugs acting independently
    expect(synergyVerdict(x, 0.05), 'more');
    expect(synergyVerdict({...x, dif: 0.01, p: 0.3, pTwoSided: 0.6}, 0.05), 'about');
    expect(synergyVerdict({...x, dif: -0.05, p: 0.99, pTwoSided: 0.01}, 0.05), 'less');
    // at α 0.05 the plan's convention reads synergy when the one-sided p is under α / 2
    const at = (p: number) => synergyWords({...x, dif: 0.05, p, pTwoSided: 2 * p}, 'A+B', 0.05).startsWith('Synergy');
    expect(at(0.03), false, 'one-sided 0.03: not under 0.025');
    expect(at(0.01), true, 'one-sided 0.01');
    expect(synergyAlpha(null, 0.1), 0.1, 'no plan: the analysis α');
    expect(synergyAlpha('not a plan', 0.1), 0.1, 'a plan that cannot be read: the analysis α');
    expect(notCombined(null).length, 0, 'no plan: no arm unticked');
    expect(notCombined('not a plan').length, 0, 'a plan that cannot be read: no arm unticked');
    const c = await competingOf(s, 'why', 'Control', 1500);
    const a = c!.find((y) => y.group === 'A')!;
    expect(a.competing, 2);
    expect(a.shr < 1, true, `A against the control: SHR ${a.shr}`);
    expect(await competingOf(s, null, 'Control', 1500), null);
    // a group with no animal at the endpoint: not estimable (the fit gave 0.00, 0.00–0.00, p < 0.001)
    const still = DG.DataFrame.fromObjects(([['Control', 'C1', 2000, ''], ['Control', 'C2', 1800, ''],
      ['N', 'N1', 300, ''], ['N', 'N2', 200, 'ulceration']] as [string, string, number, string][])
      .flatMap(([group, animal, v, why]) => [{animal, group, day: 0, volume: 100, why: ''},
        {animal, group, day: 10, volume: v, why}]))!;
    const n = (await competingOf({...s, study: still}, 'why', 'Control', 1500))!.find((y) => y.group === 'N')!;
    expect([n.none, n.competing].join(), 'group,1');
    const items = arriveItems({control: 'Control', arms: ['A', 'B', 'A+B'], plan: null, exclusions: null,
      meta: null, randomized: false, endpoint: 1500, evaluationDay: 20, protocol: 'Slopes compared', alpha: 0.05,
      correction: 'holm', analyzed: []});
    expect(items.length, 10);
    expect(items.filter((i) => !i.done).map((i) => i.n).join(','), '2,4,5,8,9,10');
    // the methods in one line, the protocol behind it; what the app can add, linked to the study's details
    expect(items[6].says, 'Welch\'s t test on growth rates, Holm over 3 comparisons, α 0.05');
    expect(items[6].detail, 'Slopes compared');
    expect(items.filter((i) => i.fix === 'details').map((i) => i.n).join(','), '5,8,9');
    // the outcome on the day the study is evaluated
    expect(items[5].says, 'Tumor volume; endpoint 1500 mm³, evaluated on day 20');
    // not blinded is not met; blinded is
    const meta = {...emptyMeta(), animals: {species: 'mouse', strain: 'C57BL/6', sex: 'female', age: '8 weeks',
      weight: '20 g', source: 'Charles River'}};
    const facts = {control: 'Control', arms: ['A'], plan: 'MC38 plan', meta, randomized: true, endpoint: 1500,
      evaluationDay: 20, protocol: 'Slopes compared', analyzed: []};
    expect(arriveItems(facts)[4].status, 'todo');
    expect(arriveItems({...facts, meta: {...meta, blinded: true}})[4].status, 'met');
    // the animals: met once every field is filled, the line says what is and what is missing
    expect(arriveItems(facts)[7].says, 'mouse, C57BL/6, female, 8 weeks, 20 g, Charles River');
    const partly = arriveItems({...facts, meta: {...emptyMeta(), animals: {...emptyMeta().animals,
      strain: 'C57BL/6'}}})[7];
    expect(partly.status, 'todo');
    expect(partly.says, 'mouse, C57BL/6: add sex, age, weight, source');
    // files from before the animal fields: a mouse, the rest to add
    expect(parseMeta({id: 'PMX-2026-001'})!.animals.species, 'mouse');
    expect(parseMeta({id: 'PMX-2026-001'})!.animals.strain, '');
    // what the groups receive, in one line; the full list, routes too, on hover
    const groups = [{group: 'LD4172', compound: 'LD4172', dose: 20, unit: 'mg/kg', route: 'p.o.', schedule: 'QD'},
      {group: 'aPD1', compound: 'anti-PD-1', dose: 100, unit: 'µg', route: 'i.p.', schedule: 'Q3D'},
      {group: 'Combo', compound: 'LD4172', dose: 20, unit: 'mg/kg', route: 'p.o.', schedule: 'QD'},
      {group: 'Combo', compound: 'anti-PD-1', dose: 100, unit: 'µg', route: 'i.p.', schedule: 'Q3D'}];
    const given = arriveItems({...facts, meta: {...meta, groups}})[8];
    expect(given.says, 'LD4172 20 mg/kg p.o. QD · anti-PD-1 100 µg i.p. Q3D · both');
    expect(given.detail!.split('\n')[1], 'aPD1: anti-PD-1 100 µg i.p. Q3D');
    expect(receivesLine([groups[0], {...groups[1], dose: 50}, groups[2], groups[3]])!.says,
      'LD4172 20 mg/kg p.o. QD · anti-PD-1 50 µg i.p. Q3D · LD4172 20 mg/kg p.o. QD + anti-PD-1 100 µg i.p. Q3D');
    // from a plan: exclusions and results known once analyzed
    const planned = arriveItems({...facts, meta: null});
    expect(planned.filter((i) => !i.done).map((i) => i.n).join(','), '3,5,8,9,10');
    expect(planned.filter((i) => i.status === 'later').map((i) => i.n).join(','), '3,10');
    expect(planned[2].says, 'Once the study is analyzed');
  });
  test('PREPARE from a plan: what the plan settles, what the details add, the rest to do', async () => {
    const facts = {line: 'MC38', protocol: 'Slopes compared', humaneLimit: 2000, lastDay: 21, toImplant: 60,
      meta: null};
    const items = prepareItems(facts);
    expect(items.length, 15);
    expect([...new Set(items.map((i) => i.group))].length, 3);
    // the design and statistics are the study protocol, which the call-to-action copies
    expect(items.filter((i) => i.done).map((i) => i.n).join(','), '4');
    expect(items[3].fix, 'protocol');
    expect(items.filter((i) => i.fix === 'details').map((i) => i.n).join(','), '9,10');
    // what each group receives and the animals, once the details say them
    const meta = {...emptyMeta(), animals: {species: 'mouse', strain: 'C57BL/6', sex: 'female', age: '8 weeks',
      weight: '20 g', source: 'Charles River'},
    groups: [{group: 'A', compound: 'Cpd-A', dose: 10, unit: 'mg/kg', route: 'PO', schedule: 'QD'}]};
    const filled = prepareItems({...facts, meta});
    expect(filled.filter((i) => i.done).map((i) => i.n).join(','), '4,9,10');
    expect(filled[8].says, 'Cpd-A 10 mg/kg PO QD');
    expect(filled[9].says, 'MC38 · mouse, C57BL/6, female, 8 weeks, 20 g, Charles River');
    expect(items[9].says, 'MC38: add species, strain, sex, age, weight, source');
    // no protocol yet: nothing met
    expect(prepareItems({...facts, protocol: ''}).filter((i) => i.done).length, 0);
  });

  test('CTD table from an analysis: the control first, each group its effect and body weight', async () => {
    const r = await analyze(await exampleStudy()) as unknown as AnalysisFrames;
    const s = ctdStudyOf(r, null, 'example');
    expect(s.groups.length, r.comparisons.rowCount + 1);
    expect(s.groups[0].control && s.groups[0].name === r.analysis.get('control', 0), true);
    expect(s.groups[1].tgi === r.comparisons.get('tgi', 0) && s.groups[1].p === r.comparisons.get('p_value', 0),
      true);
    expect(s.endpointDay, r.analysis.get('evaluation_day', 0));
    // no details: the compound and dose read from the group's name, the route a field to fill
    expect(s.compounds.join(), 'Cpd-A');
    expect([s.groups[1].given, s.groups[1].route].join('|'), `${s.groups[1].name}|`);
    expect(ctdRows([s]).length, s.groups.length + 1);
  });

  test('exposure from NCA: matched by name, else by drug, dose and route; never guessed', async () => {
    // NCA's frames as it returns them, one row per profile (results) and per group (summary)
    const frames = (named: boolean) => {
      const groups = ['Cpd-A 10 mg/kg', 'Cpd-A 30 mg/kg', 'Combo'];
      const results = DG.DataFrame.fromObjects([
        {compound: 'Cpd-A', ROUTE: 'PO', DOSE: 10, VERDICT: 'Pass', REASONS: ''},
        {compound: 'Cpd-A', ROUTE: 'PO', DOSE: 30, VERDICT: 'Pass', REASONS: ''},
        {compound: 'Cpd-B', ROUTE: '', DOSE: null, VERDICT: 'Blocked',
          REASONS: 'sparse design with fewer than 3 animals, a composite needs 3 or more'},
      ].map((r, i) => named ? {...r, group: groups[i]} : r))!;
      results.setTag('nca.units', JSON.stringify({dose: 'mg/kg'}));
      const summary = DG.DataFrame.fromObjects([0, 1].map((i) => ({compound: 'Cpd-A',
        ...named ? {group: groups[i]} : {}, ROUTE: 'PO', DOSE: [10, 30][i], N: [6, 1][i],
        MEAN_AUC_LAST: [1200, 4100][i], CV_PCT_AUC_LAST: [18, null][i], MEAN_CMAX: [300, 900][i],
        CV_PCT_CMAX: [22, null][i]})))!;
      summary.col('MEAN_AUC_LAST')!.setTag('units', 'h*ng/mL');
      summary.setTag('nca.keyColumns', JSON.stringify(named ? ['compound', 'group'] : ['compound']));
      return {results, summary};
    };
    // by the group's own name
    const named = frames(true);
    const byName = exposureByArm(['Cpd-A 10 mg/kg', 'Combo', 'Missing'], named.results, named.summary, 'group');
    expect(byName.get('Cpd-A 10 mg/kg')!.exposure!.auc, 1200);
    expect(byName.get('Cpd-A 10 mg/kg')!.exposure!.aucUnit, 'h*ng/mL');
    expect(byName.get('Combo')!.why!, 'NCA could not compute exposure for Combo: too few animals sampled ' +
      '(a pooled profile needs at least 3)');
    expect(byName.get('Missing')!.why!.endsWith('the concentrations file has none for it'), true);
    // by drug, dose and route read in the arm's name; a combination, another unit or another route is no match
    const bare = frames(false);
    const byDrug = exposureByArm(['Cpd-A 30 mg/kg PO QD', 'Cpd-A 30 mg/kg i.p.', 'Cpd-A + Cpd-B', 'Cpd-A 10 mg'],
      bare.results, bare.summary, null);
    const x = byDrug.get('Cpd-A 30 mg/kg PO QD')!.exposure!;
    expect(x.auc === 4100 && x.aucCv === null && x.profiles === 1 && x.schedule === 'QD', true, JSON.stringify(x));
    expect(byDrug.get('Cpd-A 30 mg/kg i.p.')!.why!.endsWith('the concentrations file has none for it'), true);
    expect(byDrug.get('Cpd-A + Cpd-B')!.why!.includes('exact name only'), true);
    expect(byDrug.get('Cpd-A 10 mg')!.why!.includes('same unit'), true);
    // a dose NCA blocked has no dose of its own: by drug, the drug's reason is given, not "nothing found"
    const blockedByDrug = exposureByArm(['Cpd-B 5 mg/kg'], bare.results, bare.summary, null).get('Cpd-B 5 mg/kg')!;
    expect(blockedByDrug.why!.endsWith('too few animals sampled (a pooled profile needs at least 3)'), true,
      blockedByDrug.why!);
    // males and females given apart: not pooled here, said so
    const sexes = bare.summary.clone();
    sexes.columns.addNewString('SEX').init((i) => i === 0 ? 'M' : 'F');
    sexes.set('DOSE', 1, 10);
    expect(exposureByArm(['Cpd-A 10 mg/kg'], bare.results, sexes, null).get('Cpd-A 10 mg/kg')!.why!
      .endsWith('males and females apart'), true);
    // without the NCA package nothing is read; with it, two animals sampled once each give no exposure
    const sparse = DG.DataFrame.fromObjects(['m1', 'm2'].map((a, i) => ({subject: a, group: 'Cpd-A 10 mg/kg',
      time_h: [1, 4][i], conc_ng_ml: [120, 40][i], dose_mg_kg: 10, route: 'PO'})))!;
    const m = await concentrationMapping(sparse);
    if (!ncaInstalled()) {
      expect(m, null);
      return;
    }
    expect(m !== null, true, 'NCA reads the sparse table as concentrations');
    const run = await runNca(sparse, m);
    const said = exposureByArm(['Cpd-A 10 mg/kg'], run.results, run.summary, run.groupCol).get('Cpd-A 10 mg/kg')!;
    expect(said.why ?? '', 'NCA could not compute exposure for Cpd-A 10 mg/kg: too few animals sampled ' +
      '(a pooled profile needs at least 3)');
  });
});

category('Analyze: human exposure and modelling dataset', () => {
  test('doses a day from the schedule an arm\'s name writes', async () => {
    const cases: [string | null, number | null][] = [['QD', 1], ['BID', 2], ['qd x 14', 1], ['Q3D×4', 1 / 3],
      ['QW', 1 / 7], ['Q2W', 1 / 14], ['BIW', 2 / 7], ['QOD', 1 / 2], ['Q12H', 2], ['twice a day', 2],
      ['every 3 days', 1 / 3], ['once weekly', 1 / 7], ['x3', null], ['5 doses', null], [null, null]];
    for (const [s, n] of cases)
      expect(dosesPerDay(s), n, String(s));
  });

  test('exposure units: converted to mg·h/L and mg/L, or refused', async () => {
    const near = (a: number | null, b: number) => a !== null && Math.abs(a - b) <= 1e-12 * Math.max(1, b);
    expect(near(aucFactor('mg/L·h'), 1), true);
    expect(near(aucFactor('h*ng/mL'), 1e-3), true);
    expect(near(aucFactor('ng·h/mL'), 1e-3), true);
    expect(near(aucFactor('µg·h/mL'), 1), true);
    expect(near(aucFactor('ng*min/mL'), 1e-3 / 60), true);
    expect(aucFactor('nM·h'), null, 'molar: needs the molecular weight');
    expect(aucFactor('ng/mL'), null, 'no time');
    expect(near(concentrationFactor('ng/mL'), 1e-3), true);
    expect(near(concentrationFactor('µg/mL'), 1), true);
    expect(concentrationFactor('µM'), null);
  });

  test('the mouse\'s exposure as a multiple of the human\'s: daily, like with like, never across units', async () => {
    const x: Exposure = {auc: 10, aucCv: null, aucUnit: 'mg/L·h', cmax: 2, cmaxCv: null, cmaxUnit: 'mg/L',
      aucInf: 12, aucInfUnit: 'mg/L·h', dose: '10 mg/kg', route: 'PO', schedule: 'BID', profiles: 6, drug: 'Cpd-A',
      pooled: false};
    const near = (a: number | undefined, b: number) => a !== undefined && Math.abs(a - b) <= 1e-9 * Math.max(1, b);
    const human = (h: Partial<HumanExposure> = {}): Map<string, HumanExposure> => new Map([['cpda', {drug: 'Cpd-A',
      auc: 8000, cmax: 1000, unit: 'ng/mL', fuMouse: null, fuHuman: null, source: 'projected', ...h}]]);
    // AUCinf 12 mg·h/L × 2 doses a day = 24, against 8 mg·h/L: 3×, on total exposure
    const total = humanMultiple('Cpd-A 10 mg/kg BID', x, null, human())!;
    expect(near(total.ratio, 3), true, String(total.ratio));
    expect(total.tip!.includes('AUCinf') && total.tip!.includes('total exposure') &&
      total.tip!.includes('projected'), true, total.tip);
    // AUClast when NCA gives no AUCinf, said so
    const last = humanMultiple('Cpd-A 10 mg/kg QD', {...x, aucInf: null}, null, human())!;
    expect(near(last.ratio, 1.25), true);
    expect(last.tip!.includes('AUClast'), true);
    // unbound when both fractions are given
    const free = humanMultiple('Cpd-A 10 mg/kg QD', x, null, human({fuMouse: 0.1, fuHuman: 0.2, source: 'measured'}))!;
    expect(near(free.ratio, 0.75), true);
    expect(free.tip!.includes('unbound') && free.tip!.includes('measured'), true);
    // weekly: the AUC of one dose over 7 days
    expect(near(humanMultiple('Cpd-A 10 mg/kg QW', x, null, human())!.ratio, 1.5 / 7), true);
    // a combination, an unknown schedule, a unit NCA cannot convert: a reason, no ratio
    expect(humanMultiple('Cpd-A + Cpd-B', x, null, human())!.why!.includes('combination'), true);
    expect(humanMultiple('Cpd-A 10 mg/kg x3', x, null, human())!.why!.includes('schedule'), true);
    // continuous exposure: no multiple, said so (not an unread schedule)
    expect(humanMultiple('Cpd-A 10 mg/kg', {...x, schedule: 'drinking water'}, null, human())!.why!
      .includes('no multiple is given for continuous exposure'), true);
    expect(humanMultiple('Cpd-A 10 mg/kg QD', {...x, aucInfUnit: 'nM·h'}, null, human())!.why!.includes('unit'),
      true);
    // no human value for the drug, or the control: nothing to say
    expect(humanMultiple('Cpd-B 10 mg/kg QD', x, null, human()), null);
    expect(humanMultiple('Vehicle', undefined, null, human()), null);
    // no NCA group at the dose: scaled only when exposure rises in step with the dose
    // (proportionality is judged on the AUC used, the per-dose AUCs here; the drug's own flag is not read)
    const pk = (aucs: [number, number][], proportional: Proportional = 'yes'): DrugPk => ({doses: aucs.map(([d]) => d),
      unit: 'mg/kg', auc: new Map(), proportional,
      perDose: new Map(aucs.map(([d, a]) => [d, {auc: a, unit: 'mg/L·h', inf: true}]))});
    const linear = pk([[10, 10], [30, 30]]);
    expect(near(humanMultiple('Cpd-A 20 mg/kg QD', undefined, linear, human())!.ratio, 2.5), true);
    expect(near(humanMultiple('Cpd-A 30 mg/kg QD', undefined, pk([[10, 10], [30, 45]], 'no'), human())!.ratio,
      45 / 8), true, 'a dose with PK');
    expect(humanMultiple('Cpd-A 20 mg/kg QD', undefined, pk([[10, 10], [30, 45]], 'yes'), human())!.why!
      .includes('not scaled'), true, 'AUC 1.0 and 1.5 per mg/kg: not in step, whatever the flag says');
    // (one NCA dose: any other dose is outside it, so not extrapolated, which comes first)
    expect(humanMultiple('Cpd-A 20 mg/kg QD', undefined, pk([[10, 10]]), human())!.why!.includes('not extrapolated'),
      true);
    // within NCA's doses only: no extrapolation
    expect(humanMultiple('Cpd-A 60 mg/kg QD', undefined, linear, human())!.why!.includes('not extrapolated'), true);
    expect(humanMultiple('Cpd-A 5 mg/kg QD', undefined, linear, human())!.why!.includes('not extrapolated'), true);
    // a dose a hair off NCA's is NCA's
    expect(near(humanMultiple('Cpd-A 10.0000000001 mg/kg QD', undefined, linear, human())!.ratio, 1.25), true);
    // a dose of 0, an AUC of 0: a reason, no ratio
    expect(humanMultiple('Cpd-A 0 mg/kg QD', x, null, human())!.why!.includes('not above 0'), true);
    expect(humanMultiple('Cpd-A 10 mg/kg QD', {...x, aucInf: null, auc: 0}, null, human())!.why!
      .includes('not above 0'), true);
    expect(humanMultiple('Cpd-A 10 mg/kg QD', {...x, aucInf: null, auc: NaN}, null, human())!.why!
      .includes('not above 0'), true);
    // AUCinf from fewer animals than AUClast: said with the counts; the Cmax multiple says single-dose
    const fewer = humanMultiple('Cpd-A 10 mg/kg QD', {...x, nAucInf: 3, nAucLast: 6}, null, human())!;
    expect(fewer.tip!.includes('from 3 of the 6 animals') && fewer.tip!.includes('Single-dose Cmax'), true,
      fewer.tip);
    expect(multipleText(3), '3×');
    expect(multipleText(0.214), '0.21×');
    expect(multipleText(12.6), '13×');
    // outside the patients' range: CV 50 % puts the 95th percentile at exp(1.645·√ln 1.25) ≈ 2.18× the geometric
    // mean, which is the typed (arithmetic) mean / √1.25: ≈ 1.95× and ≈ 0.41× the typed mean
    expect(exposureRange(2.0, 50)?.dir, 'above');
    expect(exposureRange(2.0, 50)!.text.includes('the range seen in patients (95th percentile, CV 50%)'), true);
    expect(exposureRange(1.9, 50), null);
    expect(exposureRange(1.8, 50), null);
    expect(exposureRange(0.45, 50), null);
    expect(exposureRange(0.35, 50)?.dir, 'below');
    expect(exposureRange(0.35, 50)!.text.includes('5th percentile'), true);
    // without the CV, or with one typed as a fraction or out of range: 3-fold either way, said to be a convention
    expect(exposureRange(2.5, null), null);
    expect(exposureRange(2.5, 0.5), null);
    expect(exposureRange(2.5, 300), null);
    expect(exposureRange(3.5, null)!.text, 'more than 3-fold above the patients\' exposure (a convention, ' +
      'not patients\' ' +
      'spread)');
    expect(exposureRange(0.3, null)?.dir, 'below');
    // in the tooltip, qualified (projected, total) and never "unachievable"
    const high = humanMultiple('Cpd-A 10 mg/kg BID', x, null, human({auc: 2000, cv: 50}))!;
    expect(high.range, 'above');
    expect(high.tip!.includes('projected; total, not free') && !/unachievable|approved/i.test(high.tip!), true,
      high.tip);
    expect(humanMultiple('Cpd-A 10 mg/kg BID', x, null, human())!.range === undefined, true);
  });

  test('the reference table fills a drug it holds, says where from, and changes nothing for one it does not',
    async () => {
      const rows = (await humanReference()).drugs;
      // by a synonym, any case: Lynparza is olaparib (FDA label 12.3: 49.2 µg·h/mL (44%) at steady state with 300 mg
      // BID, window not stated: doubled, τ = 12 h assumed; 7.6 µg/mL, 82% bound)
      const ola = referenceRow(drugKey('LYNPARZA'), rows)!;
      expect(ola?.drug, 'olaparib');
      const h = fromReference('Lynparza', ola);
      expect([h.auc, h.cmax, h.unit, h.fuHuman, h.cv, h.source].join(), '98400,7600,ng/mL,0.18,44,measured');
      expect(ola.cycleDays, null);
      // the source as every view says it: the document, then the regimen the values are for
      expect(referenceLine(h).startsWith('Lynparza: FDA label 12.3 ('), true, referenceLine(h));
      expect(sourceText(h).endsWith(', 300 mg oral twice daily (tablets), continuous'), true, sourceText(h));
      expect(sourceText({...h, reference: undefined}), 'Your value');
      expect(sourceText({...h, reference: undefined, source: 'projected'}), 'Projected');
      // a value changed: still editable, now said to be edited
      expect(referenceKept(h, ola), true);
      expect(referenceKept({...h, auc: 50000}, ola), false);
      expect(referenceLine({...h, reference: {row: ola, edited: true}}).includes('continuous (edited)'), true);
      // not in the table: nothing filled, no caveat in the tooltip
      expect(referenceRow(drugKey('Cpd-Z'), rows) === undefined, true);
      const x: Exposure = {auc: 10, aucCv: null, aucUnit: 'mg/L·h', cmax: 2, cmaxCv: null, cmaxUnit: 'mg/L',
        aucInf: 12, aucInfUnit: 'mg/L·h', dose: '10 mg/kg', route: 'IV', schedule: 'QW', profiles: 6, drug: 'Cpd-A',
        pooled: false};
      const plain: HumanExposure = {drug: 'Cpd-A', auc: 4000, cmax: null, unit: 'ng/mL', fuMouse: null, fuHuman: null,
        source: 'measured'};
      const asIs = humanMultiple('Cpd-A 10 mg/kg QW', x, null, new Map([['cpda', plain]]))!;
      expect(!asIs.tip!.includes('Patients\' values') && !asIs.tip!.includes('over the cycle'), true, asIs.tip);
      // olaparib against a mouse dosed QD: daily against daily (AUCinf 12 mg·h/L against 98.4), not twice that
      const olaQd = humanMultiple('Lynparza 10 mg/kg QD', {...x, schedule: 'QD'}, null,
        new Map([[drugKey('Lynparza'), h]]))!;
      expect(Math.abs(olaQd.ratio! - 12 / 98.4) < 1e-9, true, String(olaQd.ratio));
      // trastuzumab (AUC per 21-day cycle) against a mouse dosed QW: the AUC ÷ 21 against the mouse's ÷ 7
      const tras = referenceRow(drugKey('Herceptin'), rows)!;
      expect(tras.cycleDays, 21);
      const trasQw = humanMultiple('Herceptin 10 mg/kg QW', x, null,
        new Map([[drugKey('Herceptin'), fromReference('Herceptin', tras)]]))!;
      expect(Math.abs(trasQw.ratio! - (12 / 7) / (tras.auc! / 1000 / 21)) < 1e-12, true, String(trasQw.ratio));
      expect(trasQw.tip!.includes('average daily exposure over the cycle') && trasQw.tip!.includes(
        'AUC per administration ÷ 21 days between administrations (on average)'), true, trasQw.tip);
      // an intermittent regimen also gives one dose against one administration (mouse 12 mg·h/L; patients' AUC as is)
      expect(trasQw.tip!.includes(`Per administration: the mouse's AUC per dose against the patients' AUC per ` +
        `administration is ${Number((12 / (tras.auc! / 1000)).toPrecision(3)).toLocaleString('en-US')}×`), true,
      trasQw.tip);
      expect(olaQd.tip!.includes('Per administration'), false, olaQd.tip);
      // what was measured, when it is more than the parent drug: intact oxaliplatin, not platinum; olaparib
      // says nothing
      const oxa = referenceRow(drugKey('Eloxatin'), rows)!;
      const oxaQw = humanMultiple('Eloxatin 10 mg/kg QW', x, null,
        new Map([[drugKey('Eloxatin'), fromReference('Eloxatin', oxa)]]))!;
      expect(oxaQw.tip!.includes('measured as oxaliplatin (intact'), true, oxaQw.tip);
      expect(olaQd.tip!.includes('measured as'), false, olaQd.tip);
      // a conjugate or a liposome is its own drug, never the antibody's or the free drug's row
      expect(referenceRow(drugKey('T-DXd'), rows)?.drug, 'trastuzumab-deruxtecan');
      expect(referenceRow(drugKey('trastuzumab emtansine'), rows)?.drug, 'trastuzumab-emtansine');
      expect(referenceRow(drugKey('liposomal doxorubicin'), rows) === undefined, true);
      // the same AUC rewritten in µg/mL keeps the cycle
      const inUg = {...fromReference('Herceptin', tras), unit: 'µg/mL', auc: tras.auc! / 1000};
      const trasUg = humanMultiple('Herceptin 10 mg/kg QW', x, null, new Map([[drugKey('Herceptin'), inUg]]))!;
      expect(Math.abs(trasUg.ratio! / trasQw.ratio! - 1) < 1e-9, true, String(trasUg.ratio));
      expect(referenceLine(inUg).includes('the unit was changed'), false);
      const row = (r: Partial<ReferenceRow>): ReferenceRow => ({drug: 'Cpd-A', synonyms: [], regimen: '', auc: 4000,
        aucWindow: 'AUC0-24 at steady state', cycleDays: null, aucCv: null, cmax: null, cmaxCv: null,
        fu: null, matrix: 'plasma', stat: 'arith_mean', tier: 'label', source: 'FDA label 12.3', ...r});
      const multiple = (r: Partial<ReferenceRow>) =>
        humanMultiple('Cpd-A 10 mg/kg QW', x, null, new Map([['cpda', fromReference('Cpd-A', row(r))]]))!;
      // whole blood, a median, a model estimate: said in the tooltip
      expect(multiple({matrix: 'whole_blood'}).tip!.includes('whole blood, not plasma'), true);
      expect(multiple({stat: 'median'}).tip!.includes('a median'), true);
      expect(multiple({stat: 'model'}).tip!.includes('model estimate'), true);
      // per 21-day infusion cycle: 4 mg·h/L ÷ 21 a day against the mouse's 12 ÷ 7, labelled
      const cycle = {cycleDays: 21, aucWindow: 'AUC(0-inf) after one 3-h infusion'};
      const infusion = multiple(cycle);
      expect(Math.abs(infusion.ratio! - (12 / 7) / (4 / 21)) < 1e-9, true, String(infusion.ratio));
      expect(infusion.tip!.includes('average daily exposure over the cycle') &&
        infusion.tip!.includes('one 3-h infusion'), true, infusion.tip);
      expect(Math.abs(multiple({}).ratio! - 3 / 7) < 1e-9, true);
      // an AUC edited by hand is a daily AUC0–24: no cycle
      const edited = {...fromReference('Cpd-A', row(cycle)), auc: 8000};
      const typed = humanMultiple('Cpd-A 10 mg/kg QW', x, null, new Map([['cpda', edited]]))!;
      expect(Math.abs(typed.ratio! - (12 / 7) / 8) < 1e-9, true, String(typed.ratio));
      expect(referenceLine({...edited, reference: {row: row(cycle), edited: true}})
        .includes('÷ 21 days between administrations (on average); your AUC is used as the daily AUC0–24'), true);
      // the unit changed but not the numbers of the table: said
      expect(referenceLine({...fromReference('Cpd-A', row({})), unit: 'µg/mL'}).includes('the unit was changed'), true);
    });

  test('the exposure is kept with the study: saved and opened again, both come back; an older study opens as before',
    async () => {
      const x: Exposure = {auc: 10, aucCv: 20, aucUnit: 'mg/L·h', cmax: 2, cmaxCv: null, cmaxUnit: 'mg/L',
        aucInf: 12, aucInfUnit: 'mg/L·h', dose: '10 mg/kg', route: 'PO', schedule: 'QD', profiles: 6, drug: 'Cpd-A',
        pooled: false};
      const pk: DrugPk = {doses: [10, 30], unit: 'mg/kg', auc: new Map([[10, '10'], [30, '30']]), proportional: 'yes',
        perDose: new Map([[10, {auc: 10, unit: 'mg/L·h', inf: true}], [30, {auc: 30, unit: 'mg/L·h',
          inf: true}]])};
      const row: ReferenceRow = {drug: 'Cpd-A', synonyms: [], regimen: '', auc: 8000, aucWindow: 'AUC0-24',
        cycleDays: null, aucCv: 40, cmax: null, cmaxCv: null, fu: null, matrix: 'plasma', stat: 'geo_mean',
        tier: 'label', source: 'FDA label 12.3'};
      const humans = new Map([['cpda', {...fromReference('Cpd-A', row), auc: 9000,
        reference: {row, edited: true}}]]);
      const nca: NcaByArm = {table: 'pk_conc', arms: new Map<string, ArmExposure>([['Cpd-A 10 mg/kg QD', {exposure: x}],
        ['Cpd-A 20 mg/kg QD', {why: 'no group at 20 mg/kg'}]]), pk: new Map([['Cpd-A 20 mg/kg QD', pk]])};
      const study = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['A1', 'A1']),
        DG.Column.fromStrings('group', ['Cpd-A 10 mg/kg QD', 'Cpd-A 10 mg/kg QD']),
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'day', [0, 3]), DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'volume',
          [100, 150])]);
      study.setTag(EXPOSURE_TAG, JSON.stringify(exposureJson(humans, nca)));
      const back = parseExposure(exposureOfTag(tableOf(parseStudy(JSON.stringify(studyFileOf(study, null))))));
      const h = back.humans.get('cpda')!;
      expect([h.auc, h.cv, h.source, h.reference?.edited, h.reference?.row.source].join(), '9000,40,measured,true,' +
        'FDA label 12.3');
      expect(back.nca?.table, 'pk_conc');
      expect(back.nca!.arms.get('Cpd-A 10 mg/kg QD')!.exposure!.auc, 10);
      expect(back.nca!.arms.get('Cpd-A 20 mg/kg QD')!.why, 'no group at 20 mg/kg');
      const pkBack = back.nca!.pk.get('Cpd-A 20 mg/kg QD')!;
      expect(pkBack.perDose!.get(30)!.auc, 30);
      // the multiple reads the same from what was kept as from what NCA gave
      expect(humanMultiple('Cpd-A 20 mg/kg QD', undefined, pkBack, back.humans)!.ratio,
        humanMultiple('Cpd-A 20 mg/kg QD', undefined, pk, humans)!.ratio);
      // an older study (no exposure kept) and a damaged tag: nothing, as before
      const old = studyFileOf(study, null);
      delete old.exposure;
      const plain = parseExposure(exposureOfTag(tableOf(parseStudy(JSON.stringify(old)))));
      expect(plain.humans.size === 0 && plain.nca === null, true);
      study.setTag(EXPOSURE_TAG, '{not json');
      expect(parseExposure(exposureOfTag(study)).nca, null);
      expect(parseExposure({humans: {a: {drug: 1}}, nca: {table: 't', arms: {g: {exposure: {auc: 'x'}}}}}).humans.size,
        0);
      // a group whose name no longer says the dose NCA read (doses swapped between names) is not restored; a renamed
      // group is not found by its new name
      const swapped = parseExposure({nca: {table: 't', pk: {}, arms: {'Cpd-A 30 mg/kg QD': {exposure: x},
        'Cpd-A 10 mg/kg QD': {exposure: x}}}}).nca!;
      expect(swapped.arms.has('Cpd-A 30 mg/kg QD'), false);
      expect(swapped.arms.has('Cpd-A 10 mg/kg QD'), true);
      expect(swapped.arms.has('Cpd-A 10mg/kg daily'), false);
      // a dose NCA read from a float32 column (0.3 as 0.30000001192092896) is the name's dose (Compare lost it)
      const float32 = parseExposure({nca: {table: 't', pk: {}, arms: {'Cpd-A 0.3 mg/kg QD':
        {exposure: {...x, dose: `${Math.fround(0.3)} mg/kg`}}}}}).nca!;
      expect(float32.arms.has('Cpd-A 0.3 mg/kg QD'), true);
      // a broken reference is dropped, the values kept
      const plainHuman: HumanExposure = {drug: 'Cpd-A', auc: 9000, cmax: null, unit: 'ng/mL', fuMouse: null,
        fuHuman: null, source: 'measured'};
      const broken = parseExposure({humans: {cpda: {...plainHuman, reference: {row: 'x', edited: false}}}}).humans;
      expect(broken.get('cpda')!.auc === 9000 && broken.get('cpda')!.reference === undefined, true);
      // the patients' exposure set into what a study kept (Compare, an open table): all else as it was
      const into = humansInto({...exposureJson(new Map(), nca), effect: {day: 21, tgi: {a: 40}}}, humans);
      expect([into.nca?.table, into.effect?.tgi.a, into.humans?.cpda.auc].join(), 'pk_conc,40,9000');
      expect(humansInto(into, new Map()).humans === undefined, true);
      // a group renamed takes its exposure and the analysis's TGI with it
      const renamed = renamedExposure({...exposureJson(humans, nca), effect: {day: 21,
        tgi: {'Cpd-A 10 mg/kg QD': 40}}}, new Map([['Cpd-A 10 mg/kg QD', 'Cpd-A 10 mg/kg QD×14']]));
      expect([renamed!.nca!.arms['Cpd-A 10 mg/kg QD×14'].exposure!.auc, renamed!.effect!.tgi['Cpd-A 10 mg/kg QD×14'],
        'Cpd-A 10 mg/kg QD' in renamed!.nca!.arms].join(), '10,40,false');
    });

  test('one sentence on exposure and effect, only when there is something to say', async () => {
    const g = (arm: string, tgi: number, auc: number | null, multiple: ReadoutGroup['multiple'] = null):
      ReadoutGroup => ({arm, tgi, auc, aucUnit: 'µg·h/mL', multiple});
    // with the patients' exposure
    const placed = readout([g('Cpd-A 10 mg/kg QD', 31, 4.1, {ratio: 0.46, range: 'below', cvKnown: true}),
      g('Cpd-A 30 mg/kg QD', 62, 12, {ratio: 1.4, cvKnown: true}),
      g('Cpd-A 100 mg/kg QD', 88, 41, {ratio: 4.6, range: 'above', cvKnown: true})]);
    expect(placed, 'At an exposure within the patients\' range (1.4×), 30 mg/kg gave a TGI of 62%; 100 mg/kg gave ' +
      '88% at 4.6× the patients\' exposure.');
    // without it: TGI against dose and AUC
    const dosed = readout([g('Cpd-A 10 mg/kg QD', 31, 4.1), g('Cpd-A 30 mg/kg QD', 62, 12),
      g('Cpd-A 100 mg/kg QD', 88, 41)]);
    expect(dosed, 'The groups gave a TGI of 31% at 10 and 88% at 100 mg/kg (AUC 4.1–41 µg·h/mL).');
    // neither: nothing
    expect(readout([g('Cpd-A 10 mg/kg QD', 31, null), g('Cpd-A 100 mg/kg QD', 88, null)]), null);
    // never a target, a clinical dose, a margin, a prediction
    for (const s of [placed!, dosed!])
      expect(/target|clinical dose|margin|predict|efficacious|did more|rose|fell/i.test(s), false, s);
  });

  test('Plan: the patients\' exposure as a mouse dose, from an earlier study that kept both', async () => {
    const x = (dose: number, aucInf: number): Exposure => ({auc: aucInf * 0.9, aucCv: null, aucUnit: 'mg/L·h',
      cmax: null, cmaxCv: null, cmaxUnit: 'mg/L', aucInf, aucInfUnit: 'mg/L·h', dose: `${dose} mg/kg`, route: 'PO',
      schedule: null, profiles: 6, drug: 'Cpd-A', pooled: false});
    const patients = (auc: number, h: Partial<HumanExposure> = {}) => new Map([['cpda', {drug: 'Cpd-A', auc,
      cmax: null, unit: 'ng/mL', fuMouse: null, fuHuman: null, source: 'measured' as const, ...h}]]);
    const nca = (doses: [number, number][], drug = 'Cpd-A'): NcaByArm => ({table: 'pk', pk: new Map(),
      arms: new Map(doses.map(([d, a]) => [`${drug} ${d} mg/kg QD`, {exposure: x(d, a)}]))});
    const two = nca([[10, 10], [30, 40]]);
    const at = (h: Map<string, HumanExposure>, n: NcaByArm, arm = 'Cpd-A 20 mg/kg QD') =>
      patientsDose(arm, exposureJson(h, n), 'CRC-017');
    // between 10 (10 mg·h/L) and 30 mg/kg (40), 20 mg·h/L: log-dose against log-AUC, 10 × 3^(ln 2 / ln 4) ≈ 17
    const between = at(patients(20000), two)!;
    // (the arm names no route: the study's single one, PO, is named)
    // said of the arm's own dose: 20 mg/kg gives 10 × 4^(ln 2 / ln 3) ≈ 24 mg·h/L a day, 1.2× the patients' 20
    expect(between.text, 'At 20 mg/kg, mice get about 1.2× the patients\' exposure; 17 mg/kg would match it');
    expect(between.tip.includes('Mouse doses with PK: 10, 30 mg/kg (QD PO in mice, from CRC-017\'s PK).'), true,
      between.tip);
    expect(between.tip.includes('AUCinf') && between.tip.includes('total exposure') &&
      !between.tip.includes('may underestimate'), true, between.tip);
    // one dose: proportional, said so; a projection is said to be one
    const one = at(patients(20000, {source: 'projected'}), nca([[10, 10]]))!;
    expect(one.text, 'At 20 mg/kg, mice get about 1× the projected patients\' exposure; 20 mg/kg would match it');
    expect(one.tip.includes('beyond them, its exposure taken in step with the dose'), true, one.tip);
    // beyond the doses with PK: no dose
    const above = at(patients(100000), two)!;
    expect(above.text, 'At 20 mg/kg, mice get about 0.24× the patients\' exposure');
    expect(at(patients(5000), two)!.text, 'At 20 mg/kg, mice get about 4.8× the patients\' exposure');
    // an arm with no dose: the patients' exposure placed on the doses with PK, as before
    expect(at(patients(20000), two, 'Cpd-A QD')!.text, 'Patients\' exposure ≈ 17 mg/kg QD PO in mice (between 10 ' +
      'and 30 mg/kg, from CRC-017\'s PK)');
    // a route-less arm against a study of two routes: nothing; naming one, its groups only
    const mixed: NcaByArm = {table: 'pk', pk: new Map(), arms: new Map([...two.arms, ['Cpd-A 10 mg/kg IV QD',
      {exposure: {...x(10, 30), route: 'IV-bolus'}}]])};
    expect(at(patients(20000), mixed), null);
    expect(at(patients(20000), mixed, 'Cpd-A 20 mg/kg PO QD')!.text.endsWith('17 mg/kg would match it'), true);
    // AUClast only: used, and said to read low
    const last: NcaByArm = {table: 'pk', pk: new Map(), arms: new Map([10, 30].map((d) => [`Cpd-A ${d} mg/kg QD`,
      {exposure: {...x(d, 1), aucInf: null, auc: d === 10 ? 10 : 40}}]))};
    const lastOnly = at(patients(20000), last)!;
    expect(lastOnly.text.endsWith('17 mg/kg would match it') && lastOnly.tip.includes('AUClast: may underestimate exposure'),
      true, lastOnly.tip);
    // an arm with no schedule: the study's, said so; a study with none either: nothing
    const studys = at(patients(20000), two, 'Cpd-A 20 mg/kg')!;
    expect(studys.text, 'At 20 mg/kg, mice get about 1.2× the patients\' exposure; 17 mg/kg would match it');
    expect(studys.tip.includes('PO in mice on CRC-017\'s QD schedule'), true, studys.tip);
    const unscheduled: NcaByArm = {table: 'pk', pk: new Map(), arms: new Map([10, 30].map((d) => [`Cpd-A ${d} mg/kg`,
      {exposure: x(d, d)}]))};
    expect(at(patients(20000), unscheduled, 'Cpd-A 20 mg/kg'), null);
    // the hint to add the schedule, only when a schedule would place the dose
    const earlier = (n: NcaByArm, h = patients(20000)) =>
      [{name: 'CRC-017', line: null, meta: null, rows: [], exposure: exposureJson(h, n)}];
    expect(patientsDoseOf('Cpd-A 20 mg/kg', earlier(unscheduled))!.text.startsWith('Add the schedule'),
      true);
    expect(patientsDoseOf('Cpd-A 20 mg/kg', earlier(two))!.text.startsWith('At 20 mg/kg'), true,
      'the study\'s schedule places it: no hint');
    expect(patientsDoseOf('Cpd-A', earlier(unscheduled)), null, 'no dose: a schedule places nothing');
    expect(patientsDoseOf('Cpd-B 20 mg/kg', earlier(unscheduled)), null, 'no exposure kept for the drug');
    // the study's group details first: groups named only "Low" and "High" placed by their details
    const coded: NcaByArm = {table: 'pk', pk: new Map(), arms: new Map([['Low', {exposure: x(10, 10)}],
      ['High', {exposure: x(30, 40)}]])};
    const details = {...emptyMeta(), groups: [10, 30].map((d, i) => ({group: i ? 'High' : 'Low', compound: 'Cpd-A',
      dose: d, unit: 'mg/kg', route: 'p.o.', schedule: 'QD'}))};
    expect(at(patients(20000), coded), null, 'without details the groups name no drug');
    expect(patientsDose('Cpd-A 20 mg/kg QD', exposureJson(patients(20000), coded), 'CRC-017', details)!.text
      .endsWith('17 mg/kg would match it'), true);
    // the doses a day follow the plan's rules (no dosing at weekends), never the earlier study's
    const noWeekends = {start: '2026-11-09', weekendsOff: true};
    const coded20 = (m: typeof details, rules: typeof noWeekends | null = null) =>
      patientsDose('Cpd-A 20 mg/kg QD', exposureJson(patients(20000), coded), 'CRC-017', m, rules)!.text;
    expect(coded20({...details, ...noWeekends}), coded20(details), 'the earlier study\'s weekends do not count');
    expect(coded20(details, noWeekends) !== coded20(details), true, 'the plan\'s weekends do');
    // unbound when both fractions are known
    expect(at(patients(20000, {fuMouse: 0.1, fuHuman: 0.1}), two)!.tip.includes('unbound'), true);
    // per 21-day cycle against a weekly schedule: 43 056 mg·h/L ÷ 21 a day against AUCs ÷ 7
    const tras: NcaByArm = {table: 'pk', pk: new Map(), arms: new Map([10, 30].map((d) => [`Herceptin ${d} mg/kg QW`,
      {exposure: {...x(d, d * 700), drug: 'Herceptin'}}]))};
    const row: ReferenceRow = {drug: 'trastuzumab', synonyms: ['Herceptin'], regimen: '', auc: 43056000,
      aucWindow: 'AUCss,0-21 days', cycleDays: 21, aucCv: null, cmax: null, cmaxCv: null, fu: null, matrix: 'plasma',
      stat: 'median', tier: 'label', source: 'FDA label 12.3'};
    const weekly = patientsDose('Herceptin 10 mg/kg QW', exposureJson(new Map([[drugKey('Herceptin'),
      fromReference('Herceptin', row)]]), tras), 'BC-2')!;
    expect(weekly.text.endsWith('21 mg/kg would match it'), true, weekly.text);
    // nothing kept, another drug: nothing
    expect(patientsDose('Cpd-A 20 mg/kg QD', undefined, 'old'), null);
    expect(at(patients(20000), two, 'Cpd-B 20 mg/kg QD'), null);
    for (const d of [between, one, above, weekly, lastOnly, studys]) {
      expect(/target|efficacious|recommended|clinical dose|margin|predict/i.test(`${d.text} ${d.tip}`), false,
        d.text + d.tip);
    }
  });

  test('the groups card asks for one missing thing at a time', async () => {
    expect(exposureHint({ncaInstalled: true, mousePk: false, patients: false, inTable: []})?.add, 'nca');
    expect(exposureHint({ncaInstalled: false, mousePk: false, patients: false, inTable: []}), null);
    const patients = exposureHint({ncaInstalled: true, mousePk: true, patients: false, inTable: ['olaparib']})!;
    expect(patients.add === 'patients' && patients.text.includes('olaparib is in PMX Studio\'s reference table'),
      true, patients.text);
    expect(exposureHint({ncaInstalled: true, mousePk: true, patients: true, placed: true, inTable: ['olaparib']}),
      null);
    // the patients' exposure is in but no group is placed: not asked again, said why
    const why = exposureHint({ncaInstalled: true, mousePk: true, patients: true, placed: false, inTable: [],
      why: 'Cpd-A 10 mg/kg QD: the human exposure\'s unit (nM) cannot be converted'})!;
    expect(why.add === null && why.text.includes('cannot be converted'), true, why.text);
  });

  test('modelling dataset: a dose per animal offered in mg/kg over the group\'s mean body weight', async () => {
    const rows = [['P1', 'anti-PD-1 100 ug', 0, 100, 19], ['P2', 'anti-PD-1 100 ug', 0, 100, 21],
      ['P1', 'anti-PD-1 100 ug', 3, 150, 25], ['C1', 'LD4172 20 + anti-PD-1 100 ug', 0, 90, 20]];
    const study = DG.DataFrame.fromObjects(rows.map(([animal, group, day, volume, bw]) =>
      ({animal, group, day, volume, bw})))!;
    const s: StudyInputs = {study, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: 'bw'};
    const pd1 = modellingDataset(s, null).lines.find((l) => l.label === 'anti-PD-1 100 ug')!;
    expect(pd1.status, 'no dose', 'still asked: only offered');
    expect(JSON.stringify(pd1.proposed), JSON.stringify({dose: 5, from: 'from 100 µg at 20.0 g mean weight'}),
      'the weights on the first dose day');
    const combo = modellingDataset(s, null).lines.filter((l) => l.group === 'LD4172 20 + anti-PD-1 100 ug');
    expect(combo.find((l) => l.label.endsWith('anti-PD-1'))?.proposed?.dose, 5, 'a combination\'s drug too');
    expect(combo.length === 2 && combo[0].proposed === undefined, true, 'a dose with no unit: nothing');
    const noBw = modellingDataset({...s, bodyWeight: null}, null).lines.find((l) => l.label === 'anti-PD-1 100 ug')!;
    expect(noBw.proposed === undefined, true, 'no body weights: as before');
  });

  test('modelling dataset: NONMEM rows, doses from the schedule, exclusions out, flags kept', async () => {
    const study = DG.DataFrame.fromObjects([
      ['V1', 'Vehicle', 0, 100, 20], ['V1', 'Vehicle', 3, 200, 21],
      ['A1', 'Cpd-A 10 mg/kg', 0, 110, 22], ['A1', 'Cpd-A 10 mg/kg', 3, 150, 21], ['A1', 'Cpd-A 10 mg/kg', 7, 160, 20],
      ['C1', 'Cpd-A + Cpd-B', 0, 90, 19], ['C1', 'Cpd-A + Cpd-B', 3, 95, 19],
    ].map(([animal, group, day, volume, bw]) => ({animal, group, day, volume, bw})))!;
    study.columns.addNewBool(EXCLUDED).init((i) => i === 4);
    study.columns.addNewBool(FLAGGED).init((i) => i === 3);
    const s: StudyInputs = {study, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: 'bw'};
    const meta = {...emptyMeta(), groups: [{group: 'Cpd-A 10 mg/kg', compound: 'Cpd-A', dose: 10, unit: 'mg/kg',
      route: 'p.o.', schedule: 'Q3D×4'}]};
    const d = modellingDataset(s, meta);
    expect(d.dosed, true);
    expect(d.ids.map((x) => `${x.id}:${x.animal}:${x.group}`).join(','), '1:V1:1,2:A1:2,3:C1:3');
    const a1 = d.rows.filter((r) => r.ID === 2);
    // day 7 excluded as a measurement, yet the animal was there: Q3D×4 doses on days 0, 3 and 6, none after day 7
    expect(a1.filter((r) => r.EVID === 0).map((r) => r.TIME).join(','), '0,3');
    expect(a1.filter((r) => r.EVID === 1).map((r) => `${r.TIME}:${r.AMT}:${r.CMT}`).join(','),
      '0:10:2,3:10:2,6:10:2');
    // at one time the observation first; the flagged measurement kept with FLAG 1
    expect(a1[0].EVID, 0);
    expect(a1.find((r) => r.TIME === 3 && r.EVID === 0)!.FLAG, 1);
    expect(a1.every((r) => r.DOSE === 10 && r.GROUP === 2 && r.COMBO === 0), true);
    // the control's DOSE is 0; a combination's is -99 with COMBO 1; neither is dosed
    expect(d.rows.filter((r) => r.ID === 1).every((r) => r.DOSE === 0 && r.EVID === 0), true);
    expect(d.rows.filter((r) => r.ID === 3).every((r) => r.DOSE === MISSING && r.COMBO === 1 && r.EVID === 0), true);
    expect(d.notes.some((n) => n.includes('Cpd-A + Cpd-B, Cpd-A') && n.includes('no dose')), true, d.notes.join(' | '));
    // the dialog's lines: the control, the dosed group, each drug of the combination without a dose
    expect(d.lines.map((l) => `${l.label}:${l.status}`).join('|'), 'Vehicle:control|Cpd-A 10 mg/kg:dosed|' +
      'Cpd-A + Cpd-B, Cpd-A:no dose|Cpd-A + Cpd-B, Cpd-B:no dose');
    // a dose typed in the dialog is the drug's dose; without a schedule the combination still gets no doses
    const typed = modellingDataset(s, meta, null, {dose: new Map([['Cpd-A + Cpd-B, Cpd-A', 5]])}).lines
      .find((l) => l.label === 'Cpd-A + Cpd-B, Cpd-A')!;
    expect(`${typed.dose}:${typed.status}`, '5:no schedule');
    // the dialog shows one row per group: the combination's status is its first missing thing over both drugs
    const combo = (x: typeof d) => groupStatus(x.lines.filter((l) => l.group === 'Cpd-A + Cpd-B'));
    expect(`${groupStatus(d.lines.slice(0, 1))}:${combo(d)}`, 'control:no dose');
    expect(combo(modellingDataset(s, meta, null, {dose: new Map([['Cpd-A + Cpd-B, Cpd-A', 5],
      ['Cpd-A + Cpd-B, Cpd-B', 5]])})), 'no schedule');
    const csv = modellingCsv(d).split('\n').filter((l) => !l.startsWith('@'));
    expect(csv[0], 'ID,TIME,DV,MDV,EVID,CMT,AMT,GROUP,DOSE,COMBO,BW,FLAG');
    expect(csv[1], '1,0,100,0,0,1,.,1,0,0,20,0');
    // a later first dose shifts the schedule: days 2 and 5, none past day 7
    const later = modellingDataset(s, meta, null, {firstDose: new Map([['Cpd-A 10 mg/kg', 2]])});
    expect(later.rows.filter((r) => r.ID === 2 && r.EVID === 1).map((r) => r.TIME).join(','), '2,5');
    // a first day the schedule writes wins over the dialog's: not added on top of it
    const written = {...meta, groups: [{...meta.groups[0], schedule: 'Q3D×4 from day 1'}]};
    expect(modellingDataset(s, written, null, {firstDose: new Map([['Cpd-A 10 mg/kg', 2]])}).rows
      .filter((r) => r.ID === 2 && r.EVID === 1).map((r) => r.TIME).join(','), '1,4,7');
    // a missing body weight is -99, never "." (NONMEM reads it as 0)
    const noBw = study.clone();
    noBw.set('bw', 0, null);
    expect(modellingDataset({...s, study: noBw}, meta).rows[0].BW, MISSING);
    expect(idMapCsv(d).split('\n')[3], '3,C1,3,Cpd-A + Cpd-B');
    // without a schedule: observations only, said so
    const bare = modellingDataset(s, null);
    expect(bare.dosed, false);
    expect(bare.rows.every((r) => r.EVID === 0), true);
    // NCA's concentrations of the study's own animals, in compartment 2, hours after the dose of the day given
    const conc = DG.DataFrame.fromObjects([{subject: 'A1', time_h: 2, conc: 1.5}, {subject: 'S9', time_h: 2,
      conc: 3}, {subject: 'A1', time_h: 24, conc: 0}])!;
    const mapping = {subject: 'subject', time: 'time_h', concentration: 'conc', units: {time: 'h',
      concentration: 'mg/L'}};
    const withPk = modellingDataset(s, meta, {table: conc, groupCol: null, mapping}, {pkDoseDay: 3});
    const pkRows = withPk.rows.filter((r) => r.CMT === 2 && r.EVID === 0);
    expect(pkRows.length, 1, 'the satellite animal S9 is not one of the study\'s; a 0 is left out');
    expect(pkRows[0].ID === 2 && Math.abs((pkRows[0].TIME as number) - (3 + 2 / 24)) < 1e-12 &&
      pkRows[0].DV === 1.5, true, JSON.stringify(pkRows[0]));
    expect(withPk.notes.some((n) => n.startsWith('1 concentration at or below 0')), true, withPk.notes.join(' | '));
    // none picked: the dialog says so, and names the study's own table when it is not open
    expect(noConcentrations(null).startsWith('No concentrations are included: pick a concentration table'), true);
    expect(noConcentrations('pk_conc').includes('pk_conc, the table the study\'s exposure was computed from, is not ' +
      'open'), true);
  });

  test('a renamed group, found by its animals, keeps its typed details; details of no group are not printed',
    async () => {
      const text = (m: Map<string, string>) => [...m].map(([a, b]) => `${a}>${b}`).join();
      const before = new Map([['Vehicle', ['V1', 'V2']], ['A', ['A1', 'A2']], ['B', ['B1', 'B2']]]);
      // two renamed at once: each pair by its animals (the order of the animals does not matter)
      expect(text(renamedGroups(before, new Map([['Vehicle', ['V1', 'V2']], ['A2', ['A2', 'A1']],
        ['B2', ['B1', 'B2']]]))), 'A>A2,B>B2');
      // one removed and another added with other animals: no rename
      expect(renamedGroups(before, new Map([['Vehicle', ['V1', 'V2']], ['A', ['A1', 'A2']],
        ['C', ['C1', 'C2']]])).size, 0);
      // an orphan row of details (Plan's "Control") moves nothing; the real rename A to A2 does
      const meta = withGroups(emptyMeta(), ['Control', 'A', 'Cpd-A 10 mg/kg QD']);
      const typed = {...meta, groups: meta.groups.map((g) => g.group === 'Cpd-A 10 mg/kg QD' ? {...g, route: 'p.o.'} :
        g)};
      const groups = ['Vehicle', 'A', 'Cpd-A 30 mg/kg QD'];
      const renames = renamedGroups(new Map([['Vehicle', ['V1']], ['A', ['A1']], ['Cpd-A 10 mg/kg QD', ['X1']]]),
        new Map([['Vehicle', ['V1']], ['A', ['A1']], ['Cpd-A 30 mg/kg QD', ['X1']]]));
      expect(text(renames), 'Cpd-A 10 mg/kg QD>Cpd-A 30 mg/kg QD');
      const after = withGroups(typed, groups, renames);
      // what was typed follows the renamed group, what its old name said does not
      const row = after.groups.find((g) => g.group === 'Cpd-A 30 mg/kg QD')!;
      expect([row.compound, row.dose, row.route, row.schedule].join(), 'Cpd-A,30,p.o.,QD');
      expect(after.groups.some((g) => g.group === 'Cpd-A 10 mg/kg QD'), false);
      // the orphan "Control" stays, not printed; with no renames given, nothing moves
      expect(ofGroups(after, groups).map((g) => g.group).join(), 'Vehicle,A,Cpd-A 30 mg/kg QD');
      expect(ofGroups(after).some((g) => g.group === 'Control'), true);
      expect(withGroups(typed, groups).groups.find((g) => g.group === 'Cpd-A 30 mg/kg QD')!.route, '');
      // Plan's explicit edit of an arm's name (old to new, one row): its typed route goes with it
      const planned = renamedDetails(typed, new Map([['Cpd-A 10 mg/kg QD', 'Cpd-A 30 mg/kg QD']]));
      expect(withGroups(planned, groups).groups.find((g) => g.group === 'Cpd-A 30 mg/kg QD')!.route, 'p.o.');
    });

  test('modelling dataset: animals named alike in two groups are two animals', async () => {
    const study = DG.DataFrame.fromObjects([['1', 'Vehicle', 0, 100], ['1', 'Cpd-A 10 mg/kg QD', 0, 110],
      ['2', 'Vehicle', 0, 90], ['1', 'Vehicle', 3, 150]]
      .map(([animal, group, day, volume]) => ({animal, group, day, volume})))!;
    const s: StudyInputs = {study, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null};
    const d = modellingDataset(s, null);
    expect(d.ids.map((x) => `${x.id}:${x.animal}:${x.groupName}`).join(','),
      '1:1:Vehicle,2:1:Cpd-A 10 mg/kg QD,3:2:Vehicle');
    expect(d.rows.filter((r) => r.ID === 1).length, 2);
    expect(d.rows.every((r) => r.BW === MISSING), true, 'no body weight column');
    // the concentrations join by group and animal: with NCA's group column, animal 1 of the treated group only
    const conc = DG.DataFrame.fromObjects([{subject: '1', grp: 'Cpd-A 10 mg/kg QD', time_h: 1, conc: 2}])!;
    const mapping = {subject: 'subject', time: 'time_h', concentration: 'conc', units: {time: 'h'}};
    const pk = modellingDataset(s, null, {table: conc, groupCol: 'grp', mapping}).rows
      .filter((r) => r.CMT === 2 && r.EVID === 0);
    expect(pk.length === 1 && pk[0].ID === 2, true, JSON.stringify(pk));
    // without it, a name two groups share is not placed
    expect(modellingDataset(s, null, {table: conc, groupCol: null, mapping}).rows
      .some((r) => r.CMT === 2 && r.EVID === 0), false);
  });

  test('modelling dataset: a day\'s doses spread over it, the last day\'s kept', async () => {
    const study = DG.DataFrame.fromObjects([['A1', 'Cpd-A 10 mg/kg BID', 0, 100], ['A1', 'Cpd-A 10 mg/kg BID', 1, 110],
      ['T1', 'Cpd-A 10 mg/kg TID', 0, 100]].map(([animal, group, day, volume]) => ({animal, group, day, volume})))!;
    const s: StudyInputs = {study, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null};
    const doses = (id: number) => modellingDataset(s, null).rows.filter((r) => r.ID === id && r.EVID === 1)
      .map((r) => r.TIME);
    expect(doses(1).join(','), '0,0.5,1,1.5', 'BID on the last measured day keeps its second dose');
    expect(doses(1).length, doseCount(readSchedule('BID'), 1));
    expect(doses(2).join(','), '0,0.3333,0.6667', 'TIME to 4 decimals');
  });

  test('exposure curve: between the doses tested only, a band of ±10; silent where it was not validated', async () => {
    // a vehicle growing by Gompertz (100 to a plateau of 3,000 mm³) and Cpd-X at the given kill rates κ (1/day), 4
    // animals a group (`few`: that dose 2), measured on days 0, 7, 14 and 21
    const arm = (d: number) => `Cpd-X ${d} mg/kg i.p. QD`;
    const study = (kappas: [number, number][], few: number | null = null): CurveRow[] => {
      const rows: CurveRow[] = [];
      const y = (t: number) => Math.log(3000) + (Math.log(100) - Math.log(3000)) * Math.exp(-0.1 * t);
      const g = (t: number) => (1 - Math.exp(-0.1 * t)) / 0.1;
      for (const [group, k, n] of [['Vehicle', 0, 4], ...kappas.map(([d, k]) => [arm(d), k, d === few ? 2 : 4])] as
        [string, number, number][]) {
        for (let a = 0; a < n; a++) {
          for (const day of [0, 7, 14, 21])
            rows.push({animal: `${group}-${a}`, group, day, volume: Math.exp(y(day) - k * g(day)) + a});
        }
      }
      return rows;
    };
    // κ in step with the dose: TGI about 21, 54 and 90% on day 21; NCA's AUC in step with the dose too
    const linear = (doses: number[], c = 0.00874) => doses.map((d) => [d, c * d] as [number, number]);
    const groupsOf = (doses: number[]) => doses.map((d) => ({arm: arm(d), auc: 100 * d}));
    const rows = study(linear([3, 10, 30]));
    const c = await exposureCurve(groupsOf([3, 10, 30]), rows, 'Vehicle', 21);
    if ('why' in c && c.why === 'this version of PMX Studio cannot draw it yet')
      return 'skipped: the loaded Crux module has no stats.tumor_effect_from_study yet';
    if ('why' in c)
      throw new Error(c.why);
    const xs = c.curve.map((p) => p.x);
    expect(xs[0] === 300 && xs[xs.length - 1] === 3000 && c.note === undefined, true,
      'from the lowest dose\'s AUC to the highest\'s, never beyond');
    expect(c.curve.every((p, i) => p.lo === Math.max(0, p.tgi - 10) && p.hi === Math.min(100, p.tgi + 10) &&
      (i === 0 || p.tgi > c.curve[i - 1].tgi)), true, JSON.stringify(c.curve));
    const mid = medianTgi(rows, arm(10), 'Vehicle', 21)!;
    expect(Math.abs(mid - 54) < 3, true, `1 − T/C of medians at 10 mg/kg: ${mid}`);
    // the day the medians are read on, which the chart's axis names
    expect(medianDay(rows, arm(10), 'Vehicle', 21), 21);
    // the axis: AUClast in step with the dose, else the dose
    expect(exposureAxis(groupsOf([3, 10, 30])).axis, 'auc');
    const notInStep = exposureAxis([300, 1000, 9000].map((auc, i) => ({arm: arm([3, 10, 30][i]), auc})));
    expect(notInStep.axis === 'dose' && notInStep.unit === 'mg/kg' && notInStep.x(10) === 10, true);
    expect(exposureAxis(groupsOf([3, 10, 30]).map((g, i) => ({...g, auc: i ? g.auc : null}))).axis, 'dose');
    // a group with a TGI under 5% counts as no kill, and the curve is still drawn
    const weak = await exposureCurve(groupsOf([3, 10, 30]), study([[3, 0.001], [10, 0.04], [30, 0.12]]), 'Vehicle',
      21);
    expect('curve' in weak, true, JSON.stringify(weak));
    // from the first dose over 95%: drawn below it when 3 doses remain, said so; else silent
    const four = await exposureCurve(groupsOf([3, 10, 30, 100]), study(linear([3, 10, 30, 100])), 'Vehicle', 21);
    expect('curve' in four && four.curve[four.curve.length - 1].x === 3000 && four.note!.startsWith('Drawn up to 30'),
      true, JSON.stringify(four));
    const stopped = await exposureCurve(groupsOf([3, 10, 30]), study(linear([3, 10, 30], 0.0133)), 'Vehicle', 21);
    expect('why' in stopped && stopped.why.startsWith('from 30 mg/kg the tumors stopped almost completely'), true,
      JSON.stringify(stopped));
    // doses too far apart for Crux's scaling (a kill rate more than 4 times another's)
    const apart = await exposureCurve(groupsOf([1, 10, 100]), study([[1, 0.001], [10, 0.03], [100, 0.25]]),
      'Vehicle', 21);
    expect(JSON.stringify(apart), JSON.stringify({why: 'the tested doses are too far apart to draw a curve between ' +
      'them'}));
    // fewer than 3 animals in a group
    const few = await exposureCurve(groupsOf([3, 10, 30]), study(linear([3, 10, 30]), 10), 'Vehicle', 21);
    expect(JSON.stringify(few), JSON.stringify({why: 'a group or the control has fewer than 3 animals measured on ' +
      'any one day'}));
    // silent from the names, saying why
    expect(curveSilence([arm(3), arm(10)]), 'fewer than 3 doses');
    expect(JSON.stringify(await exposureCurve(groupsOf([3, 10]), rows, 'Vehicle', 21)),
      JSON.stringify({why: 'fewer than 3 doses'}));
    expect(curveSilence([arm(3), arm(10), 'Cpd-X 30 mg/kg + Cpd-Y 5 mg/kg']), 'a group is a combination');
    expect(curveSilence(['anti-PD-1 3 mg/kg', 'anti-PD-1 10 mg/kg', 'anti-PD-1 30 mg/kg']),
      'the drug acts through the immune system');
    expect(curveSilence([arm(3), arm(10), 'Cpd-X 30 mg i.p. QD']), 'the doses are not all in one unit');
    expect(curveSilence([arm(3), arm(10), 'Cpd-X 30 mg/kg p.o. QD']),
      'the groups are not all dosed on the same schedule and route');
  }, {timeout: 120000});
});

category('Analyze: a combination inside translation', () => {
  const COMBO = 'Cpd-A 10 mg/kg Q3D + Cpd-B 5 mg/kg QD';

  test('a combination\'s details: one row per drug, read drug by drug', async () => {
    const meta = withGroups(emptyMeta(), ['Vehicle', COMBO]);
    const rows = meta.groups.filter((g) => g.group === COMBO);
    expect(rows.map((g) => `${g.compound}:${g.dose}:${g.schedule}`).join(','), 'Cpd-A:10:Q3D,Cpd-B:5:QD');
    // a route typed on the second drug's row: its part only, said so
    const typed = {...meta, groups: meta.groups.map((g) => g.group === COMBO && g.compound === 'Cpd-B' ?
      {...g, route: 'i.v.'} : g)};
    const parts = armParts(COMBO, typed);
    expect(parts.map((p) => `${p.drug}:${p.dose?.value}:${p.route}`).join(','), 'Cpd-A:10:null,Cpd-B:5:i.v.');
    expect(JSON.stringify(parts.from), '["route"]');
    expect(detailsName(COMBO, typed).includes(' + '), true);
    // renamed (Plan's edit of the name): the typed route goes with its drug's row, the rest from the new name
    const NEW = 'Cpd-A 20 mg/kg Q3D + Cpd-B 5 mg/kg QD';
    const after = withGroups(renamedDetails(typed, new Map([[COMBO, NEW]])), ['Vehicle', NEW]).groups
      .filter((g) => g.group === NEW);
    expect(after.map((g) => `${g.compound}:${g.dose}:${g.route}`).join(','), 'Cpd-A:20:,Cpd-B:5:i.v.');
    // an older study kept one empty row for the combination: the second drug's row is added
    const old = {...emptyMeta(), groups: [{group: COMBO, compound: '', dose: null, unit: '', route: '', schedule: ''}]};
    expect(withGroups(old, [COMBO]).groups.map((g) => g.compound).join(','), 'Cpd-A,Cpd-B');
    // ... and one whose row names the second drug gets the first one added, not the second twice
    const second = {...emptyMeta(), groups: [{...old.groups[0], compound: 'Cpd-B'}]};
    expect(withGroups(second, [COMBO]).groups.map((g) => g.compound).join(','), 'Cpd-B,Cpd-A');
    // a combination renamed to one drug: one row, one drug
    const MONO = 'Cpd-A 10 mg/kg Q3D';
    const mono = withGroups(renamedDetails(meta, new Map([[COMBO, MONO]])), ['Vehicle', MONO]);
    expect(mono.groups.filter((g) => g.group === MONO).length, 1);
    expect(armParts(MONO, mono).length, 1);
  });

  test('a combination\'s exposure and multiple, drug by drug; kept with the study', async () => {
    const results = DG.DataFrame.fromObjects(['Cpd-A', 'Cpd-B'].map((compound, i) => ({compound, group: COMBO,
      ROUTE: 'PO', DOSE: [10, 5][i], VERDICT: 'Pass', REASONS: ''})))!;
    results.setTag('nca.units', JSON.stringify({dose: 'mg/kg'}));
    const summary = DG.DataFrame.fromObjects(['Cpd-A', 'Cpd-B'].map((compound, i) => ({compound, group: COMBO,
      ROUTE: 'PO', DOSE: [10, 5][i], N: 6, MEAN_AUC_LAST: [16, 3][i], CV_PCT_AUC_LAST: 20, MEAN_CMAX: [4, 1][i],
      CV_PCT_CMAX: 25})))!;
    summary.col('MEAN_AUC_LAST')!.setTag('units', 'mg/L·h');
    summary.col('MEAN_CMAX')!.setTag('units', 'mg/L');
    summary.setTag('nca.keyColumns', JSON.stringify(['compound', 'group']));
    // the group alone is no match (two drugs); each drug is
    expect(exposureByArm([COMBO], results, summary, 'group').get(COMBO)!.why!.includes('given apart'), true);
    const each = exposureByDrug(COMBO, results, summary, 'group');
    expect(each.map((d) => `${d.part.drug}:${d.x.exposure?.auc}`).join(','), 'Cpd-A:16,Cpd-B:3');
    const nca = ncaByArm([COMBO], {results, summary, groupCol: 'group'}, 'pk_conc');
    expect(nca.arms.has(COMBO), false);
    expect(nca.drugs!.get(COMBO)!.length, 2);
    // Cpd-A: 16 mg·h/L every 3 days against 8 mg·h/L a day: 0.667×; Cpd-B has no patients' exposure: nothing said
    const humans = new Map<string, HumanExposure>([['cpda', {drug: 'Cpd-A', auc: 8000, cmax: null, unit: 'ng/mL',
      fuMouse: null, fuHuman: null, source: 'measured'}]]);
    const ms = drugMultiples(COMBO, nca.drugs!.get(COMBO)!, humans);
    expect(ms.map((d) => `${d.drug}:${d.m?.ratio === undefined ? null : Math.round(d.m.ratio * 1000) / 1000}`)
      .join(','), 'Cpd-A:0.667,Cpd-B:null');
    // the group's details: a block per drug, the drug named
    const det = groupDetails({group: COMBO, read: '', tgi: 60, day: 21, x: undefined, table: 'pk_conc', live: true,
      multiple: null, human: undefined, curve: null, drugs: ms.map((d) => ({drug: d.drug, multiple: d.m,
        human: humans.get(d.key), x: nca.drugs!.get(COMBO)!.find((x) => x.key === d.key)!.x}))});
    const heads = det.blocks.map((b) => b.heading);
    expect(['Mouse exposure, Cpd-A', 'Mouse exposure, Cpd-B', 'Against the patients, Cpd-A']
      .every((h) => heads.includes(h)), true, heads.join(' | '));
    // saved and opened again; renamed with its group
    const back = parseExposure(JSON.parse(JSON.stringify(exposureJson(humans, nca))));
    expect(back.nca!.drugs!.get(COMBO)!.map((d) => d.x.exposure?.auc).join(','), '16,3');
    const NEW = 'Cpd-A 10 mg/kg Q3D + Cpd-B 5 mg/kg QD p.o.';
    const moved = parseExposure(renamedExposure(exposureJson(humans, nca), new Map([[COMBO, NEW]])));
    expect(moved.nca!.drugs!.has(NEW), true);
    // a dose of one drug changed in the name since: that drug's exposure is left out
    const other = parseExposure(renamedExposure(exposureJson(humans, nca),
      new Map([[COMBO, 'Cpd-A 30 mg/kg Q3D + Cpd-B 5 mg/kg QD']])));
    expect(other.nca!.drugs!.get('Cpd-A 30 mg/kg Q3D + Cpd-B 5 mg/kg QD')!.map((d) => d.drug).join(), 'Cpd-B');
  });

  test('modelling dataset: each drug of a combination dosed into its own compartment', async () => {
    // radiation, never dosed in mg/kg, takes no compartment
    const groups = ['Vehicle', 'RT 2 Gy QD', 'Cpd-A 10 mg/kg Q3D', 'Cpd-B 5 mg/kg QD', COMBO];
    const study = DG.DataFrame.fromObjects(groups.flatMap((group, g) => [0, 3].map((day) =>
      ({animal: `M${g}`, group, day, volume: 100 + day, bw: 20}))))!;
    const s: StudyInputs = {study, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: 'bw'};
    const d = modellingDataset(s, null);
    expect(d.compartments.map((c) => `${c.drug}:${c.cmt}`).join(','), 'Cpd-A:2,Cpd-B:3');
    expect(d.dosedGroups.includes(COMBO), true);
    const doses = (id: number) => d.rows.filter((r) => r.ID === id && r.EVID === 1)
      .map((r) => `${r.CMT}@${r.TIME}:${r.AMT}`).join(',');
    expect(doses(2), '');
    expect(doses(3), '2@0:10,2@3:10');
    expect(doses(4), '3@0:5,3@1:5,3@2:5,3@3:5');
    expect(doses(5), '2@0:10,3@0:5,3@1:5,3@2:5,2@3:10,3@3:5');
    expect(d.rows.filter((r) => r.ID === 5).every((r) => r.DOSE === MISSING && r.COMBO === 1), true);
  });
});

category('Analyze: failures in plain words', () => {
  test('the engine\'s own text stays out; a sentence of ours is said', async () => {
    const said = 'The analysis could not be computed.';
    expect(plainError(new EngineError('panicked at src/stats.rs:42: index out of bounds'), said), said);
    expect(plainError(new Error('Pick the day column in the form.'), 'fallback'), 'Pick the day column in the form.');
    expect(plainError({weird: true}, 'fallback'), 'fallback');
  });
  test('an interval with a negative bound reads "−33 to 42%", never "(-33–42%)"', async () => {
    expect([percentRange(-33.2, 41.6), percentRange(12, 42), percentRange(-50, -10)].join('|'),
      '−33 to 42%|12–42%|−50 to −10%');
  });
});

category('Analyze: what the result and its numbers say together', () => {
  test('no dosing at weekends: the multiple is a week\'s average, said so', async () => {
    const x: Exposure = {auc: 10, aucCv: null, aucUnit: 'mg/L·h', cmax: null, cmaxCv: null, cmaxUnit: 'mg/L',
      aucInf: 12, aucInfUnit: 'mg/L·h', dose: '10 mg/kg', route: 'PO', schedule: 'QD', profiles: 6, drug: 'Cpd-A',
      pooled: false};
    const human = new Map<string, HumanExposure>([['cpda', {drug: 'Cpd-A', auc: 8000, cmax: null, unit: 'ng/mL',
      fuMouse: null, fuHuman: null, source: 'projected'}]]);
    const every = humanMultiple('Cpd-A 10 mg/kg QD', x, null, human)!;
    const weekdays = humanMultiple('Cpd-A 10 mg/kg QD', x, null, human,
      {...emptyMeta(), start: '2026-11-09', weekendsOff: true})!;
    expect(Math.abs(weekdays.ratio! / every.ratio! - 5 / 7) < 1e-9, true, `${weekdays.ratio} / ${every.ratio}`);
    expect(every.weekly === undefined, true);
    expect(weekdays.weekly!.includes('no dosing at weekends'), true, weekdays.weekly);
    const said = (weekly: boolean) => readout([{arm: 'Cpd-A 10 mg/kg QD', tgi: 60, auc: 12, aucUnit: 'mg/L·h',
      multiple: {ratio: weekdays.ratio!, cvKnown: false, weekly}}]) ?? '';
    expect(said(true).includes('averaged over the week'), true, said(true));
    expect(said(false).includes('averaged over the week'), false, said(false));
  });

  test('a result its T/C interval seems to contradict gets its note; one that agrees, none', async () => {
    const note = intervalNote([{activity: 'Slows growth', tcLower: 53, tcUpper: 113}], 21);
    expect(note, 'The result is decided on the growth curves over the whole study. The T/C interval describes ' +
      'day 21 alone, so it is wider.');
    expect(intervalNote([{activity: 'No clear effect', tcLower: 40, tcUpper: 80}], 21)!.endsWith('is not the test.'),
      true);
    expect(intervalNote([{activity: 'Slows growth', tcLower: 30, tcUpper: 70},
      {activity: 'No clear effect', tcLower: 60, tcUpper: 120}], 21) === null, true);
  });
});
