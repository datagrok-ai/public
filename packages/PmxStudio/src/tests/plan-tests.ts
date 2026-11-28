import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {anchorSource, armEffectProblem, armsFrame} from '../arms-input';
import {scheduleText} from '../schedule';
import {fixedValue} from '../ui';
import {ENDPOINT, HUMANE, LAB_WEEKDAYS, MIN_GROUP, StudyDeclaration, calibrated, combinationFailed, combinationN,
  combinationNote, combinationsOf, synergyClause, designCombination, designStudy, finite, LAB_FACTOR_HIGH,
  labSpreadGroup, endpointCandidates, TGI_ABOVE_100, tcOfTgi, tgiMax, tgiOfTc, tgiProblem,
  unreachedEffect} from '../engine';
import {cruxCalls} from '../crux';
import {declarationOf, tgisOf} from '../frames';
import {bandCalibration, enginePriors} from '../priors';
import {badValues} from '../plan-file';
import {combinationDrugs, demoStudies} from '../plan-view';
import {anchorOrder, findAnchor, previousOfTable, withSaved} from '../previous-studies';
import {FINISHED_TAG} from '../study-file';
import {DEMO_PLAN_ID} from '../hub-view';
import {parseArmName} from '../arm-name';

/** pmxPlanStudy as the Plan view calls it: the outputs by name. */
async function plan(params: {[k: string]: any}): Promise<{[k: string]: DG.DataFrame}> {
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxPlanStudy'})[0].prepare({line: 'MC38', endpoint: 1000,
    power: 80, correction: 'holm', ...params});
  await fc.call();
  return Object.fromEntries(['plan', 'armResults', 'powerCurve', 'curves', 'duration', 'schedules', 'calendar',
    'scenarios'].map((o) => [o, fc.outputs[o] as DG.DataFrame]));
}

const arm = (name: string, effect = 60) => ({name, therapy: 'small_molecule' as const, metric: 'TGI' as const, effect});

category('Plan files', () => {
  test('a plan file with values Plan cannot have saved is refused with what is wrong (QA 0.1.3 final #2, #11, #30)',
    async () => {
      const fileArm = (name: string, effect: unknown) => ({name, metric: 'TGI', effect});
      const file = (o: {arms?: unknown[], endpoint?: unknown, line?: string}) => ({inputs: {line: o.line ?? 'MC38',
        endpoint: o.endpoint ?? 1000, arms: o.arms ?? [fileArm('A', 60)]}, predicted: {line: 'MC38', endpoint: 1000,
        arms: [fileArm('A', 60)]}});
      expect(badValues(file({})), '');
      expect(badValues(file({endpoint: -5})), 'the endpoint volume is -5, not above 0');
      expect(badValues(file({arms: [fileArm('A', 200)]})), 'arm "A": expected TGI is 200%, it must be above 0 and ' +
        'below 111.1% (111.1% = complete regression with this start volume and endpoint)');
      // a regression (TGI past 100%) up to just below complete regression, 1000 / (1000 − 100) = 111.1%
      expect(badValues(file({arms: [fileArm('A', 105)]})), '');
      expect(badValues(file({arms: [fileArm('A', 111)]})), '');
      expect(badValues(file({arms: [fileArm('A', 111.2)]})).startsWith('arm "A": expected TGI is 111.2%'), true);
      expect(badValues(file({arms: [fileArm('A', 0)]})).includes('above 0'), true);
      // with the file's own volumes: from 300 to 600 mm³ the bound is 200%
      const own = {...file({arms: [fileArm('A', 150)]}), inputs: {line: 'MC38', endpoint: 600, vRand: 300,
        arms: [fileArm('A', 150)]}};
      own.predicted.arms = [fileArm('A', 150)];
      expect(badValues(own), '');
      // a T/C above 0, at most 100%
      expect(badValues(file({arms: [{name: 'A', metric: 'TC', effect: 3}]})), '');
      expect(badValues(file({arms: [{name: 'A', metric: 'TC', effect: 0}]})),
        'arm "A": expected T/C is 0%, it must be above 0 and at most 100%');
      expect(badValues(file({arms: [fileArm('A', 'NaN')]})).includes('TGI is NaN%'), true);
      expect(badValues(file({line: 'BOGUS'})), 'its model is BOGUS, but its predictions are for MC38');
      const many = badValues(file({arms: Array.from({length: 20}, (_, i) => fileArm(`T${i + 1}`, 60))}));
      expect(many.startsWith('20 treated arms'), true);
      expect(many.includes('"T16", "T17", "T18", "T19", "T20"'), true);
      // a dose below 0, in a group's details or an arm's name: refused with Plan's own words (QA 0.1.3 round 4 R4-3)
      const dosed = (dose: unknown) => ({...file({}), study: {groups: [{group: 'A', dose}]}});
      expect(badValues(dosed(-3)), 'group "A": a dose cannot be below 0 (it is -3)');
      expect(badValues(dosed('abc')).includes('a dose cannot be below 0'), true);
      expect(badValues(dosed(10)), '');
      expect(badValues(file({arms: [fileArm('Cpd -5 mg/kg', 60)]})).endsWith('a dose cannot be below 0'), true);
    });

  test('the protocol gives a combination\'s drugs their own regimen, each drug once (QA 0.1.3 final #5)',
    async () => {
      expect(combinationDrugs(parseArmName('cpd-a 30 mg/kg qd po + CPD-B 10 mg/kg ip q3dx4')).join(' | '),
        'cpd-a 30 mg/kg p.o. QD | CPD-B 10 mg/kg i.p. Q3D×4');
    });

  test('a running study open in the workspace is interim, never offered (QA 0.1.3 final #9, round 3 W3)',
    async () => {
      const df = DG.DataFrame.fromCsv('animal,group,day,volume\nV-01,Vehicle,0,100\nT-01,Cpd-A 5 mg/kg,0,100\n');
      expect(previousOfTable(df)?.interim, true);
      df.setTag(FINISHED_TAG, 'true');
      expect(previousOfTable(df)?.interim === undefined, true);
      const s = (name: string, line: string, interim?: true) =>
        ({name, line, meta: null, rows: [], ...interim ? {interim} : {}});
      expect(anchorOrder([s('run', 'CT26', true), s('other', 'MC38'), s('mine', 'CT26')], 'CT26')
        .map((x) => x.name).join(), 'mine,other');
    });

  test('an offer names its study by name, ID and day: two studies of one name told apart (QA round 3 W3)',
    async () => {
      const meta = {id: 'PMX-2026-0702'} as any;
      expect(anchorSource({study: {name: 'CT26 study', line: 'CT26', meta, rows: []}, tgiDay: 16, readDay: 14}),
        'CT26 study (PMX-2026-0702, finished, day 16)');
      expect(anchorSource({study: {name: 'CT26 study', line: 'CT26', meta: null, rows: []}, tgiDay: null,
        readDay: 14}), 'CT26 study (finished, day 14)');
      // round 4 designer wording 1: a name that says it is finished, not twice
      expect(anchorSource({study: {name: 'CT26 demo study, finished', line: 'CT26', meta: null, rows: []},
        tgiDay: 21, readDay: 21}), 'CT26 demo study (finished, day 21)');
    });

  test('a finished saved study is not hidden by an open copy still running (QA round 4 designer D6)', async () => {
    const s = (name: string, line: string, interim?: true) =>
      ({name, line, meta: null, rows: [], ...interim ? {interim} : {}});
    const all = withSaved([s('CT26 study', 'CT26', true)], [s('CT26 study', 'CT26'), s('MC38 study', 'MC38')]);
    expect(anchorOrder(all, 'CT26').map((x) => x.name).join(), 'CT26 study,MC38 study');
    expect(withSaved([s('CT26 study', 'CT26')], [s('CT26 study', 'CT26')]).length, 1, 'a finished open one stands');
  });

  test('a schedule is shown in one writing: codes in capitals, a count with × (QA round 3 U5)', async () => {
    expect(scheduleText('q3dx4'), 'Q3D×4');
    expect(scheduleText('Q3D x 4'), 'Q3D×4');
    expect(scheduleText('qd'), 'QD');
    expect(scheduleText('bid x14 from day 3'), 'BID×14 from day 3');
    expect(scheduleText('x3'), '×3');
    expect(scheduleText('twice weekly'), 'twice weekly');
    expect(scheduleText('every 3 days'), 'every 3 days');
    expect(parseArmName('Cpd-B 10 mg/kg ip q3dx4')[0].schedule, 'Q3D×4');
  });

  test('the demo plan scales its arms from the demo\'s own finished study only (QA 0.1.3 final #9)', async () => {
    const studies = await demoStudies();
    expect(studies.length, 1);
    expect(studies[0].interim === undefined, true);
    expect(studies[0].line, 'CT26');
    const a = await findAnchor('Trametinib 1 mg/kg p.o. QD', {id: DEMO_PLAN_ID, name: 'CT26 trametinib + anti-PD-1'},
      'CT26', demoStudies);
    expect(a?.study.name, studies[0].name);
    expect(a?.arm, 'Trametinib 1 mg/kg p.o. QD');
  });
});

category('Plan', () => {
  test('a lower endpoint is offered only where the field takes it: from 500 mm³, above the volume at randomization',
    async () => {
      expect(endpointCandidates(1000, 100).join(), '900,800,700,600,500');
      expect(endpointCandidates(600, 100).join(), '500');
      expect(endpointCandidates(1000, 650).join(), '900,800,700');
      expect(endpointCandidates(550, 100).length, 0, 'nothing below 500');
    });
  test('a value set from a card is checked as its field was: whole, in range; empty goes back', async () => {
    const v = (raw: unknown) => JSON.stringify(fixedValue(raw, 2, 100));
    expect(v(12), JSON.stringify({value: 12}));
    expect(v(null), JSON.stringify({value: null}), 'back to the computed value: empty, not 0');
    expect(v(''), JSON.stringify({value: null}));
    expect(v(undefined), JSON.stringify({value: null}));
    for (const bad of [7.5, 1, 1000, NaN, 0, -0, Infinity, -Infinity, 'abc'])
      expect('error' in fixedValue(bad, 2, 100), true, `${bad} accepted`);
  });

  test('the validated lines draw their bands as wide as the validation found them, others as the engine', async () => {
    const cal = await bandCalibration();
    const f = (line: string) => cal.lines[line] ?? null;
    expect([f('MC38')?.end80, f('CT26')?.end80, f('B16F10')?.end80].join(' '), '1.21 1.15 1.58');
    expect(f('4T1'), null, 'a line not validated has no factor');
    const curve = {day: 7, arm: 'Control', q05: 100, q25: 150, median: 200, q75: 260, q95: 400};
    const x = {endpointDay: 15, durationQ10: 9, durationQ90: 23 as number | null, curves: [curve]};
    // the end's 8-in-10 range, around day 15 in log days: day 9 to 23 becomes about 8.1 to 25.2 (MC38's 1.21)
    const mc38 = calibrated(x, f('MC38'));
    expect(mc38.durationQ10.toFixed(2), '8.08');
    expect(mc38.durationQ90!.toFixed(2), '25.16');
    expect([f('MC38')?.volume90, f('CT26')?.volume90, f('B16F10')?.volume90].join(' '), '1.25 1.21 1.28');
    // MC38's volume band of 100 to 400 mm³ around 200, 1.25 times as wide in log volume
    expect([mc38.curves[0].q05.toFixed(1), mc38.curves[0].q95.toFixed(1)].join(' '), '84.1 475.7');
    // a 90% volume band of 100 to 400 mm³ around 200, 1.2 times as wide in log volume; the quartiles untouched
    const v = calibrated(x, {end80: null, volume90: 1.2}).curves[0];
    expect([v.q05.toFixed(1), v.q95.toFixed(1), v.q25, v.median, v.q75].join(' '), '87.1 459.5 150 200 260');
    expect(calibrated(x, {end80: null, volume90: 1.2}).durationQ10, 9, 'the end untouched');
    // a range that never closes stays open; an end near day 0 stays above half a day; no factor: as it was
    expect(calibrated({...x, durationQ90: null}, f('B16F10')).durationQ90, null);
    expect(calibrated({...x, endpointDay: 2, durationQ10: 0.2}, f('B16F10')).durationQ10 >= 0.5 - 1e-9, true);
    expect(JSON.stringify(calibrated(x, null)), JSON.stringify(x));
    // only the drawn bands move: the animals and anything else the design carries pass through as they were
    const sized = calibrated({...x, nControl: 7, nTotal: 28}, f('B16F10'));
    expect(`${sized.nControl} ${sized.nTotal}`, '7 28');
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
    // a regression is a TGI past 100%, up to just below complete regression (111.1% from 100 to 1000 mm³)
    expect(tgisOf('50, 100, 105').join(), '50,100,105');
    expect(tgisOf('150', 600, 300).join(), '150', 'from 300 to 600 mm³ the bound is 200%');
    const refusal = (text: string, endpoint?: number, v0?: number) => {
      try {
        tgisOf(text, endpoint, v0);
        return '';
      } catch (e: any) {
        return String(e.message);
      }
    };
    for (const bad of ['0', '112', 'high', '10 20 30 40 50 60 70 80 90'])
      expect(refusal(bad) !== '', true, bad);
    expect(refusal('105', 500, 100), '', '125% from 100 to 500 mm³');
    expect(refusal('112'), 'Assumed TGIs: above 0 and below 111.1% (111.1% = complete regression with this start ' +
      'volume and endpoint), such as 30, 50, 70 (got "112")');
    expect(refusal('130', 500, 100).includes('below 125%'), true, refusal('130', 500, 100));
  });

  test('tgiMax: complete regression, endpoint / (endpoint − V0); each TGI said against it', async () => {
    expect(Math.abs(tgiMax(1000, 100) - 1000 / 9) < 1e-9, true);
    expect(tgiMax(600, 300), 200);
    expect(tgiMax(1000, 1000), 100, 'an endpoint not above V0: no regression bound');
    expect(tgiProblem(105, 1000, 100), null);
    expect(tgiProblem(111.2, 1000, 100), 'Expected TGI: above 0 and below 111.1% (111.1% = complete regression with ' +
      'this start volume and endpoint)');
    expect(tgiProblem(0, 1000, 100) !== null, true);
    expect(tgiProblem(null, 1000, 100) !== null, true);
    expect(TGI_ABOVE_100, 'Above 100% = the group\'s mean volume shrinks below its start (some animals may shrink ' +
      'earlier)');
    // the bound floored to one decimal: 111 is taken, so the message says 111.1, not 111 (QA 0.1.3 r4)
    expect(tgiProblem(111, 1000, 100), null);
    expect(tgiProblem(125, 500, 100)!.includes('below 125%'), true, 'a whole bound stays whole');
    // the arm card's message: TGI by the plan's volumes, T/C above 0 and below 100%
    expect(armEffectProblem(108, 'TGI', {v0: 100, endpoint: 1000}), null);
    expect(armEffectProblem(115, 'TGI', {v0: 100, endpoint: 1000}), tgiProblem(115, 1000, 100));
    expect(armEffectProblem(115, 'TGI', {v0: 100, endpoint: 500}), null, '125% from 100 to 500 mm³');
    expect(armEffectProblem(3, 'TC', {v0: 100, endpoint: 1000}), null);
    expect(armEffectProblem(0, 'TC', {v0: 100, endpoint: 1000}), 'Expected T/C: above 0 and below 100%');
    // Crux's refusal of an effect its kill rates cannot reach, in plain words
    expect(unreachedEffect('stats.tumor_study_design: with 100 % responders and onset day 3, the effect reaches at ' +
      'most TGI 104 % by day 14: lower the effect or change the archetype'), 'This effect needs a regression faster ' +
      'than the model allows (at most TGI 104% for this kind of treatment); try a lower TGI');
    expect(unreachedEffect('with 50 % responders and onset day 3, the effect reaches at most TGI 62 % by day 14: ' +
      'lower the effect or change the archetype')!.startsWith('This effect is stronger than the model allows'), true);
    expect(unreachedEffect('something else'), null);
  });

  test('a TGI of 110% plans (sent to Crux as its T/C), or says the model\'s "try a lower TGI"', async () => {
    const arms = armsFrame([arm('Cpd-A', 110)]);
    let message = '';
    let out: {[k: string]: DG.DataFrame} | null = null;
    try {
      out = await plan({arms});
    } catch (e: any) {
      message = String(e?.message ?? e);
    }
    // never Crux's own bound ("effect=110 is outside [0, 100]"): a plan, or why in words the user can act on
    expect(message === '' || message.includes('try a lower TGI'), true, message);
    if (out) {
      const r = out.armResults;
      expect(r.get('metric', 0), 'TGI');
      expect(r.get('effect', 0), 110);
      // the effect Crux gave said back as the TGI it was declared in, not a T/C
      expect(Math.abs(r.get('effect_realised', 0) - 110) < 5, true, String(r.get('effect_realised', 0)));
    }
  }, {timeout: 300000});

  test('TGI ≥ 100% and its T/C at the median lab: one formula, both ways', async () => {
    // V0 100, endpoint 1000: TGI 100 is T/C 10% (no growth), complete regression 111% is T/C 0
    expect(Math.abs(tcOfTgi(100, 1000, 100) - 10) < 1e-9, true);
    expect(Math.abs(tcOfTgi(1000 / 9, 1000, 100)) < 1e-9, true);
    expect(Math.abs(tcOfTgi(60, 1000, 100) - 46) < 1e-9, true, 'the same as convertedEffect');
    for (const tgi of [100, 105, 110, 30])
      expect(Math.abs(tgiOfTc(tcOfTgi(tgi, 1000, 100), 1000, 100) - tgi) < 1e-9, true, String(tgi));
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

  // MC38 on Plan's own model (A at 40 % TGI, B at 30 %, endpoint 1500, Mon / Wed / Fri). Crux's own checks pin the
  // numbers of the function; here its shape through this package's wasm: signs, order, nulls
  test('a combination against each drug alone, and what its synergy needs', async () => {
    const d: StudyDeclaration = {line: 'MC38', vRand: 100, endpoint: 1500, humane: 2000, lastDay: null,
      weekdays: LAB_WEEKDAYS, startWeekday: 0, power: 0.8, alpha: 0.05, correction: 'holm', n: null,
      arms: [arm('Alphanib', 40), arm('Betamab', 30), arm('Alphanib + Betamab', 59)], preLoss: 0.2, postLoss: 0.1,
      confidence: 0.9, lossRule: 'confident', schedules: false, mde: false, fits: false};
    const s = await designStudy(d);
    const [a, b] = s.arms;
    const added = await designCombination(d, a, b, 'TGI', 0, false);
    // asked again with the same inputs (a plan drawn again): no second engine call, the same row
    const calls = cruxCalls();
    expect(JSON.stringify(await designCombination(d, a, b, 'TGI', 0, false)), JSON.stringify(added));
    expect(cruxCalls() - calls, 0);
    // left to the caller (Plan, after the plan shows): the plan without the combination's sentence, the rest the same
    const quiet = await designStudy({...d, combinationNotes: false});
    expect(s.notes.groups.length - quiet.notes.groups.length, 1, s.notes.groups.join(' | '));
    expect(quiet.arms.map((x) => x.n).join(), s.arms.map((x) => x.n).join());
    expect(finite(added['effect_bliss']) !== null && finite(added['effect_ida']) !== null, true, 'the effects');
    expect(Number(added['effect_bliss']) > Number(added['effect_ida']), true, 'together, more than the best alone');
    const more = await designCombination(d, a, b, 'TGI', 0.5, false);
    expect(Number(more['effect_ab']) > Number(added['effect_ab']), true, 'a larger interaction, a larger effect');
    // more animals, more assurance; a larger interaction never needs more animals
    const powers = await Promise.all([10, 20, 30, 40].map((n) =>
      designCombination({...d, n}, a, b, 'TGI', 0, false).then((r) => Number(r['assurance_vs_best_at_n']))));
    expect(powers.every((p, k) => k === 0 || p >= powers[k - 1]), true, `assurance by n ${powers}`);
    const needs = await Promise.all([0.25, 0.5, 1].map((x) =>
      designCombination(d, a, b, 'TGI', x, false).then((r) => finite(r['n_vs_best']) ?? Infinity)));
    expect(needs.every((n, k) => k === 0 || n <= needs[k - 1]), true, `n_vs_best by interaction ${needs}`);
    // a larger group detects a smaller synergy
    const [at10, at36] = await Promise.all([10, 36].map((n) => designCombination({...d, n}, a, b, 'TGI', 0, true)));
    expect(Number(at36['mde_synergy']) < Number(at10['mde_synergy']), true,
      `mde ${at10['mde_synergy']} at 10, ${at36['mde_synergy']} at 36`);
    // the plan says it in its groups notes, for the arm named as the other two together: what its groups would give
    const said = s.notes.groups.find((x) => x.startsWith('Alphanib + Betamab: '));
    expect(!!said && said.includes('beats Alphanib and Betamab alone'), true, s.notes.groups.join(' | '));
  }, {timeout: 300000});

  // the sentences alone, on Crux's row as given (no engine): what the user reads in the groups notes
  test('the sentence of a combination', async () => {
    const K = ' (the group sizes stay as planned)';
    const note = (row: {[k: string]: unknown}) => combinationNote('A + B', row, 0.8, 'A and B');
    expect(note({n_vs_best: 12, n: 12, assurance_vs_best_at_n: 0.81}),
      'A + B would need 12 animals per group to show it beats A and B alone' + K);
    expect(note({n_vs_best: 18, n: 10, assurance_vs_best_at_n: 0.55}), 'A + B: showing it beats A and B alone ' +
      'takes 18 animals per group; with 10 per group it would show it in 55% of studies like this one' + K);
    // 0.796 against 0.8 is 79%, not "80%, short of the target"
    expect(note({n_vs_best: null, n: 40, assurance_vs_best_at_n: 0.796}), 'A + B: with 40 animals per group it ' +
      'would show it beats A and B alone in 79% of studies like this one, short of the target power' + K);
    expect(note({n_vs_best: null, n: 40, assurance_vs_best_at_n: 0.62}), 'A + B: with 40 animals per group it would ' +
      'show it beats A and B alone in 62% of studies like this one, short of the target power' + K);
    expect(note({n_vs_best: null, n: 50, assurance_vs_best_at_n: 0.85}), 'A + B: with 50 animals per group it would ' +
      'show it beats A and B alone in 85% of studies like this one' + K, 'above the search, at the target: not short');
    expect(note({n_vs_best: null, n: null, assurance_vs_best_at_n: null}),
      'A + B: showing it beats A and B alone takes more than 40 animals per group' + K);
    expect(combinationNote('A + B + C', {n_vs_best: 12, n: 12, assurance_vs_best_at_n: 0.81}, 0.8, 'A + B and C'),
      'A + B + C would need 12 animals per group to show it beats A + B and C alone' + K,
      'the arms named, not "each drug"');
    // the synergy: the effect the combination would need, in the plan's metric (never past 100% TGI)
    const mde = {mde_synergy: 1.36, synergy_testable: true};
    const needed = {effect_ab: 82.4, effect_bliss: 58.6};
    expect(synergyClause(mde, needed, 'TGI'), '; synergy would only show if it reached TGI 82%, against the 59% ' +
      'expected from the two drugs together');
    // T/C: lower is better, so the effect is a T/C the combination would have to fall to, below the expected one
    const tc = {effect_ab: 41.2, effect_bliss: 58.6};
    expect(synergyClause(mde, tc, 'TC'), '; synergy would only show if T/C fell to 41%, against the 59% expected ' +
      'from the two drugs together');
    expect(synergyClause(mde, {...tc, effect_ab: 0}, 'TC').includes('not tested for synergy'), true, 'T/C <= 0');
    expect(combinationNote('A + B', {n_vs_best: 12, n: 12, assurance_vs_best_at_n: 0.81}, 0.8, 'A alone and B alone'),
      'A + B would need 12 animals per group to show it beats A and B alone' + K, 'no "alone alone"');
    expect(synergyClause(mde, {...needed, effect_ab: 104}, 'TGI').includes('not tested for synergy'), true,
      'past 100% TGI: not testable');
    expect(synergyClause({...mde, synergy_testable: false}, needed, 'TGI').includes('not tested for synergy'), true);
    expect(synergyClause(mde, null, 'TGI').includes('not tested for synergy'), true, 'no row, no figure');
    expect(combinationFailed('A + B'), 'A + B: Plan could not work out how many animals it needs to show it beats ' +
      'the arms it is made of');
  });

  test('the combinations of a study, by the names of its arms', async () => {
    const d = (names: string[], no = -1) => ({arms: names.map((name, i) => ({...arm(name, 50),
      ...(i === no ? {noCombination: true} : {})}))}) as unknown as StudyDeclaration;
    const rates = [{killRate: 0.1}, {killRate: 0.1}, {killRate: 0.3}] as Parameters<typeof combinationsOf>[1];
    const found = combinationsOf(d(['Alphanib', 'Betamab', 'Alphanib + Betamab']), rates);
    expect(found.length, 1);
    expect(`${found[0].ab} ${found[0].a} ${found[0].b}`, '2 0 1');
    expect(Math.abs(found[0].interaction - 0.5) < 1e-9, true, `interaction ${found[0].interaction}`);
    expect(combinationsOf(d(['Alphanib', 'Betamab', 'Alphanib + Betamab'], 2), rates).length, 0, 'unticked');
    expect(combinationsOf(d(['Alphanib', 'Betamab', 'Gammastat']), rates).length, 0);
    expect(combinationsOf(d(['Alphanib', 'Betamab', 'Alphanib + Betamab']),
      [{killRate: 0}, {killRate: 0}, {killRate: 0}] as typeof rates).length, 0, 'no rates: no interaction');
  });
  test('combinations: doses, three drugs, controls, groups compared', async () => {
    const d = (names: string[]) => ({arms: names.map((name) => arm(name, 50))}) as unknown as StudyDeclaration;
    const rates = (n: number) => Array.from({length: n}, (_, i) => ({killRate: 0.1 + (i === n - 1 ? 0.1 : 0),
      n: 10 + i * 4})) as Parameters<typeof combinationsOf>[1];
    const pairs = (names: string[]) => combinationsOf(d(names), rates(names.length))
      .map((c) => `${c.ab}:${c.a}+${c.b}`).join(' ');
    // doses: the arm whose dose the combination gives; none when none does
    expect(pairs(['Cpd-X 5 mg/kg', 'Cpd-X 10 mg/kg', 'anti-PD-1', 'Cpd-X 10 mg/kg + anti-PD-1']), '3:1+2');
    expect(pairs(['Cpd-X 5 mg/kg', 'anti-PD-1', 'Cpd-X 10 mg/kg + anti-PD-1']), '', 'a dose no single arm has');
    expect(pairs(['Cpd-X', 'anti-PD-1', 'Cpd-X 10 mg/kg + anti-PD-1']), '2:0+1', 'a dose not written on one side');
    // three drugs: "A + B" is not a drug, the texts name the arms
    expect(pairs(['Cpd-A + Cpd-B', 'Cpd-C', 'Cpd-A + Cpd-B + Cpd-C']), '2:0+1');
    expect(combinationNote('Cpd-A + Cpd-B + Cpd-C', {n_vs_best: 9, n: 9, assurance_vs_best_at_n: 0.8}, 0.8,
      'Cpd-A + Cpd-B and Cpd-C').includes('beats Cpd-A + Cpd-B and Cpd-C alone'), true);
    // a control written into the name is not a drug
    expect(pairs(['Vehicle', 'saline', 'Vehicle + saline']), '');
    expect(pairs(['anti-PD-1', 'isotype', 'anti-PD-1 + isotype']), '');
    // the groups compared: the smallest of the three, not the control's
    const c = {ab: 2, a: 0, b: 1};
    const arms = rates(3).map((x, i) => ({...x, n: [14, 9, 20][i]})) as Parameters<typeof combinationsOf>[1];
    expect(combinationN({n: null} as StudyDeclaration, arms, c), 9);
    expect(combinationN({n: 12} as StudyDeclaration, arms, c), 12, 'the study\'s own');
    expect(combinationN({n: null} as StudyDeclaration, [{killRate: 1, n: null}, ...arms.slice(1)] as typeof arms, c) ===
      undefined, true, 'unsized: Crux own default');
  });

  test('the lab factor: the plan\'s own design, the central 0.90 its n exactly, 1.25 no fewer', async () => {
    const d: StudyDeclaration = {line: 'Hepa1-6', vRand: 100, endpoint: ENDPOINT, humane: HUMANE, lastDay: null,
      weekdays: LAB_WEEKDAYS, startWeekday: 0, power: 0.8, alpha: 0.05, correction: 'holm', n: null,
      arms: [arm('Cpd-1')], preLoss: 0.2, postLoss: 0.1, confidence: 0.9, lossRule: 'confident', schedules: false};
    expect((await enginePriors()).lines['Hepa1-6'].study_prior!.size_factor, 0.9, 'the central factor');
    const [central, high] = [await labSpreadGroup(d, 0.9), await labSpreadGroup(d, LAB_FACTOR_HIGH)];
    // the same path as the card's n (its fits averaged, its seeds): 0.90 / 0.90 is the plan itself
    expect(central, (await designStudy(d)).nControl, 'the central factor gives the plan\'s own n');
    expect(central !== null && high !== null, true, `no size: ${central}, ${high}`);
    expect(high! >= central!, true, `1.25 gives ${high}, 0.90 gives ${central}`);
  }, {timeout: 300000});
});
