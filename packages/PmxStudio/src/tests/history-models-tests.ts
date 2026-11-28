import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';

import {crux, rows as tableRows} from '../crux';
import {HUMANE, LAB_WEEKDAYS, designStudy} from '../engine';
import {SavedHistory, historyFit, historyKey, historyNotesOf, historyOf, registerHistory, savedHistory} from
  '../history-models';
import {controlRows, defaultName, historyDialog, takeOk} from '../history-dialog';
import {powerText} from '../plan-view';
import {HISTORY_LOG_BETA, historyNotes, otherLab} from '../model-history';
import {parsePlan} from '../plan-file';
import {bandFactors, historyModelName, hostClassOf, linePrior, modelName} from '../priors';
import {lackingImmunity, therapyFor} from '../arms-input';
import {HISTORY_TAG, registerHistoryOf, setHistoryOn} from '../study-file';

// synthetic histories only: no provider's data in tests
const DAYS = [0, 3, 7, 10, 14, 17, 21];
/** Five animals growing exponentially at about 0.09/day, each a little apart (deterministic). */
const ROWS = Array.from({length: 5}, (_, a) => DAYS.map((day) =>
  ({animal: `A${a + 1}`, day, volume: 120 * Math.exp((0.08 + 0.005 * a) * day) * (1 + 0.03 * Math.sin(a + day))})))
  .flat();
const OWN: SavedHistory = {name: 'PDX-TEST-1', rows: ROWS, take: null, source: {sameLab: 'yes', sameHost: 'yes',
  passage: 'later', frozen: false, meansOnly: false, animalsPerGroup: null}};

const design = (line: string) => designStudy({line, vRand: 100, endpoint: 1500, power: 0.8, alpha: 0.05,
  correction: 'holm', n: null, arms: [{name: 'x', therapy: 'small_molecule', metric: 'TGI', effect: 60}],
  preLoss: 0, postLoss: 0, lossRule: 'confident', confidence: 0.9, weekdays: LAB_WEEKDAYS, startWeekday: 0,
  humane: HUMANE, schedules: false, lastDay: null});

category('Tumor models: own history in the app', () => {
  test('the dialog keeps the control group only, and every row without groups', async () => {
    const rows = [{animal: 'a', group: 'Vehicle', day: 0, volume: 100}, {animal: 'b', group: 'Drug', day: 0,
      volume: 90}, {animal: 'a', group: 'Vehicle', day: 3, volume: null}];
    expect(controlRows(rows, 'Vehicle').length, 1);
    expect(controlRows(rows, null).length, 2, 'no group column: all measured rows');
  });

  test('the dialog opens with «Use this model», off until a history is read, and closes', async () => {
    historyDialog(() => {});
    const dialogs = () => [...document.querySelectorAll<HTMLElement>('.u2-dialog, .d4-dialog')]
      .filter((d) => d.textContent?.includes('Model not in the atlas'));
    const use = () => dialogs().flatMap((d) => [...d.querySelectorAll<HTMLButtonElement>('button')])
      .find((b) => b.textContent?.trim() === 'Use this model');
    await awaitCheck(() => !!use(), 'the dialog did not open with its Use this model button', 3000);
    expect(use()!.disabled, true, 'nothing read yet');
    const close = [...dialogs()[0].querySelectorAll<HTMLButtonElement>('button')]
      .find((b) => /^(cancel|close)$/i.test(b.textContent?.trim() ?? ''));
    expect(!!close, true, 'a way out');
    close!.click();
    await awaitCheck(() => dialogs().length === 0, 'the dialog did not close', 3000);
  });

  test('a history model is found by its key: its growth prior, its name, the engine\'s bands', async () => {
    const key = registerHistory(OWN);
    expect(key.startsWith('new:pdx:'), true, key);
    expect(registerHistory({...OWN}), key, 'the same history, the same key');
    expect(historyModelName(key), 'PDX-TEST-1');
    expect(savedHistory(key)?.name, 'PDX-TEST-1');
    const p = await linePrior(key);
    expect(p.log_beta, HISTORY_LOG_BETA);
    expect(p.papers, 0);
    expect(p.study_prior !== undefined && p.twin !== undefined, true, 'what the engine and the twin need');
    expect(await bandFactors(key), null, 'band factor 1: the engine\'s bands');
  });

  test('take: the pilot\'s is a loss before randomization, never after it (QA W2)', async () => {
    const ref = registerHistory(OWN);
    const pilot = registerHistory({...OWN, take: 50});
    const [r, q] = [await historyFit(ref)!, await historyFit(pilot)!];
    expect(q.prior.take.nonprogressing_map, r.prior.take.nonprogressing_map, 'the prior keeps the reference\'s');
    expect(q.prior.study_prior!.mean_log_alpha, r.prior.study_prior!.mean_log_alpha, 'the growth does not change');
    const [nr, nq] = [await historyNotesOf(ref, 21), await historyNotesOf(pilot, 21)];
    expect(nr.some((n) => n.includes('set in Advanced, as the history gives no take')), true, nr.join(' | '));
    expect(nq.some((n) => n.includes('from the pilot\'s take (50%') && n.includes('50% lost before randomization')),
      true, nq.join(' | '));
    // fewer tumors growing: more animals to implant, the same to randomize and analyse (design's losses are 0)
    const [dr, dq] = [await design(ref), await design(pilot)];
    expect(dq.preLoss, 0.5, 'the pilot\'s no-take, over the declared 0');
    expect(dq.postLoss, dr.postLoss, 'the loss after randomization stays');
    expect(dq.randomiseTotal, dr.randomiseTotal, 'the animals to randomize');
    expect((dq.implant ?? 0) > (dr.implant ?? 0), true, `implant ${dr.implant} / ${dq.implant}`);
    expect(JSON.stringify(dq.arms.map((a) => a.n)), JSON.stringify(dr.arms.map((a) => a.n)), 'n per group');
    expect(dq.nControl, dr.nControl, 'the control\'s n');
  }, {timeout: 300000});

  test('the notes the user decides on: another lab, and a plan past four weeks', async () => {
    const own = registerHistory(OWN);
    const other = registerHistory({...OWN, source: {...OWN.source, sameLab: 'no'}});
    const short = await historyNotesOf(own, 21);
    expect(short.some((n) => n.includes('another lab') || n.includes('assumed to keep its pace')), false,
      short.join(' | '));
    const long = await historyNotesOf(other, 35);
    expect(long.some((n) => n.includes('another lab')), true, long.join(' | '));
    expect(long.some((n) => n.includes('assumed to keep its pace')), true, long.join(' | '));
    expect((await historyNotesOf('MC38', 35)).length, 0, 'an atlas model: none');
  });

  test('the plan file keeps the history and reopens on its key; a damaged one is dropped', async () => {
    const key = registerHistory(OWN);
    const range = {q025: 1, median: 2, q975: 3};
    const plan = {format: 'pmx-plan', version: 1, name: 'p', created: '', author: '', history: OWN,
      inputs: {line: 'stale-key', arms: []}, predicted: {line: 'stale-key', endpoint: 1500, vRand: 100,
        studyEnd: 21, endpointDay: 21, durationQ10: 18, durationQ90: 25, nControl: 8, arms: [],
        checks: {end: range, controlSlope: range, controls: [{day: 0, median: 100}], arms: []}}};
    const p = parsePlan(JSON.stringify(plan));
    expect(p.history?.name, 'PDX-TEST-1');
    expect(p.predicted.line, key, 'the key the history gives');
    expect(p.inputs.line, key);
    expect(historyOf({...OWN, rows: []}), null, 'no rows: not a history');
    expect(historyOf({...OWN, take: 140}), null, 'a take above 100%');
    const damaged = parsePlan(JSON.stringify({...plan, history: {...OWN, source: {}}}));
    expect('history' in damaged, false);
  });

  test('a study table carries its model\'s history: reopened, the model is found again', async () => {
    const h = {...OWN, name: 'PDX-TEST-TAG'};
    const df = DG.DataFrame.fromObjects([{animal: 'a', group: 'Vehicle', day: 0, volume: 100}])!;
    df.setTag(HISTORY_TAG, JSON.stringify(h));
    registerHistoryOf(df);
    expect(historyModelName(registerHistory(h)), 'PDX-TEST-TAG');
    df.setTag(HISTORY_TAG, '{"broken": true}');
    registerHistoryOf(df); // a damaged tag: nothing registered, nothing thrown
  });

  test('a history model is shown by its name, from a file or table without «history» or the extension', async () => {
    expect(defaultName('PDX-SYN-01 history'), 'PDX-SYN-01');
    expect(defaultName('PDX-SYN-01_history.csv'), 'PDX-SYN-01');
    expect(defaultName('PDX-LU-117.xlsx'), 'PDX-LU-117');
    expect(defaultName('history.csv'), 'history.csv', 'nothing left: the name as it came');
    const h = {...OWN, name: 'PDX-TEST-NAME'};
    const key = registerHistory(h);
    expect(modelName(key), 'PDX-TEST-NAME');
    expect(modelName('MC38'), 'MC38');
  });

  test('an arm no size up to 40 powers shows the power at 40, never a 0% for none', async () => {
    const curve = DG.DataFrame.fromObjects([{arm: 'A', n: 39, implant: 40, assurance: 0.5},
      {arm: 'A', n: 40, implant: 41, assurance: 0.62}, {arm: 'B', n: 40, implant: 41, assurance: 0.9}])!;
    expect(powerText(null, DG.FLOAT_NULL, curve, 'A'), '62% with 40');
    expect(powerText(null, null, curve, 'C'), '—');
    expect(powerText(12, 0.81, curve, 'A'), '81%');
  });

  test('a history model is a PDX host: an unrecognised drug acts directly, not through immunity', async () => {
    const key = registerHistory({...OWN, name: 'PDX-TEST-HOST'});
    expect(hostClassOf(key, {}), 'immunodeficient');
    expect(hostClassOf('MC38', null), null, 'an atlas model not loaded: not known');
    expect(therapyFor('Drug A 50 mg/kg p.o. QD', hostClassOf(key, {})).therapy, 'small_molecule');
  });

  test('a drug acting through the immune system is flagged on an immunodeficient host, and only there', async () => {
    const pdx = hostClassOf(registerHistory({...OWN, name: 'PDX-TEST-IMMUNE'}), {});
    const pd1 = therapyFor('Anti-PD-1 10 mg/kg i.p. Q3D', pdx).therapy;
    expect(lackingImmunity('Anti-PD-1 10 mg/kg i.p. Q3D', pd1, pdx).join(), 'Anti-PD-1');
    expect(lackingImmunity('Anti-PD-1 10 mg/kg i.p. Q3D', pd1, 'immunocompetent').length, 0, 'a host that answers');
    expect(lackingImmunity('Anti-PD-1 10 mg/kg i.p. Q3D', pd1, null).length, 0, 'host not known yet: no warning');
    // a combination: only its immune drug; a drug acting on the tumor and a cell therapy are not flagged
    expect(lackingImmunity('Trametinib 1 mg/kg + anti-PD-1 10 mg/kg', 'immunotherapy', pdx).join(), 'anti-PD-1');
    expect(lackingImmunity('Trametinib 1 mg/kg', 'small_molecule', pdx).length, 0);
    expect(lackingImmunity('CAR-T cells 5e6 cells', 'cell_therapy', pdx).length, 0);
  });

  test('the pilot\'s take is a share of the tumors implanted: 1 to 100%, or none', async () => {
    expect([null, 1, 55.5, 100].every(takeOk), true);
    expect([0, 0.5, 100.1, 250, NaN].some(takeOk), false);
  });

  test('reopening a study never wipes the history its table carries (D694 A)', async () => {
    const h = {...OWN, name: 'PDX-TEST-KEEP'};
    const df = DG.DataFrame.fromObjects([{animal: 'a', group: 'Vehicle', day: 0, volume: 100}])!;
    df.setTag(HISTORY_TAG, JSON.stringify(h));
    // Monitor's remember() before any declaration has read the table: its own history is registered first
    setHistoryOn(df, historyKey(h));
    expect(historyOf(JSON.parse(df.getTag(HISTORY_TAG) || 'null'))?.name, 'PDX-TEST-KEEP');
    setHistoryOn(df, 'new:pdx:not-in-this-session');
    expect(historyOf(JSON.parse(df.getTag(HISTORY_TAG) || 'null'))?.name, 'PDX-TEST-KEEP', 'an unknown key: kept');
    setHistoryOn(df, 'MC38');
    expect(df.getTag(HISTORY_TAG) || '', '', 'an atlas model: taken off');
  });

  test('the same lab only with the same strain and site; not known counts as another lab (D694 B)', async () => {
    const lab = (sameHost: 'yes' | 'no' | 'unknown') => ({...OWN.source, sameHost});
    const wider = (sameHost: 'yes' | 'no' | 'unknown') =>
      historyNotes(lab(sameHost), 1, 0.7, 3).filter((n) => n.includes('ranges are wider'));
    expect(otherLab(lab('yes')), false);
    expect(wider('yes').length, 0, wider('yes').join(' | '));
    expect(otherLab(lab('unknown')), true);
    expect(wider('unknown')[0]?.startsWith('It is not known whether the history used the same mouse strain'), true,
      wider('unknown').join(' | '));
    expect(otherLab(lab('no')), true);
    expect(wider('no')[0]?.startsWith('The history used another mouse strain'), true, wider('no').join(' | '));
  });

  test('the history\'s growth through the app\'s engine (wasm) matches the contrast, case E-A (D694 D)', async () => {
    // analysis/out_pdx/history_validation/v_exp_engine.json, case E-A: the same table and parameters
    const table = {w: Float64Array.of(1), a: Float64Array.of(-2.2967795231752404),
      b: Float64Array.of(0.337259319390191)};
    const params = {log_beta: -14, log_v_ref: 4.678540596536774, sd_log_alpha: 0.2723895976827105, sd_log_beta: 0,
      sd_study_log_beta: 0.01, rho_study: 0, start_sd: 0.14579941798260548, nonprogressing: 0,
      me_c0: 0.054552708645625864, me_c1: 0.8588851445696877, noise_df: 2.851824039835491, sd_process: 0,
      sd_session: 0, v_rand: 150, n: 8, reps: 10000, last_day: 21, weekdays: 'mon,wed,fri', start_weekday: 0,
      humane: 1e9, study: 'prior'};
    // each study's mean ln V of its 8 controls on day 21, over 4 seeds
    const means: number[] = [];
    for (const seed of [20261124, 7, 11, 13]) {
      const sum = new Map<number, [number, number]>();
      for (const r of tableRows((await crux('stats.tumor_study_simulate', table, {...params, seed})).main as any)) {
        if (r.arm === 'control' && r.day === 21) {
          const s = sum.get(r.rep) ?? [0, 0];
          sum.set(r.rep, [s[0] + Math.log(r.volume), s[1] + 1]);
        }
      }
      for (const [t, n] of sum.values())
        means.push(t / n);
    }
    means.sort((x, y) => x - y);
    const q = (p: number) => {
      const i = p * (means.length - 1);
      const lo = Math.floor(i);
      return means[lo] + (means[Math.min(lo + 1, means.length - 1)] - means[lo]) * (i - lo);
    };
    const [p5, med, p95] = [q(0.05), q(0.5), q(0.95)];
    expect(means.length, 40000);
    expect(Math.abs(med - 7.1965) <= 0.02, true, `median ${med}`);
    expect(Math.abs(p5 - 6.2256) <= 0.05 && Math.abs(p95 - 8.9103) <= 0.05, true, `P5 ${p5}, P95 ${p95}`);
  }, {timeout: 600000});

  test('the twin works on a history model (Analyze, as Monitor\'s, through twinModel)', async () => {
    const key = registerHistory(OWN);
    // a synthetic study: controls at the history's pace, a treated arm at half of it
    const rows: {animal: string, group: string, day: number, volume: number}[] = [];
    for (const [group, rate] of [['Vehicle', 0.09], ['Drug', 0.045]] as const) {
      for (let a = 0; a < 6; a++) {
        for (const day of [0, 2, 4, 7, 9, 11, 14])
          rows.push({animal: `${group}-${a}`, group, day, volume: 110 * Math.exp((rate + 0.004 * (a - 2.5)) * day)});
      }
    }
    const study = DG.DataFrame.fromObjects(rows)!;
    // the history travels with the study, as Monitor puts it there: the function runs in the package's own bundle,
    // whose models this test's registration does not reach
    study.setTag(HISTORY_TAG, JSON.stringify(OWN));
    const fc =DG.Func.find({package: 'PmxStudio', name: 'pmxAnalyzeStudy'})[0].prepare({study,
      animal: study.col('animal'), group: study.col('group'), day: study.col('day'), volume: study.col('volume'),
      line: key, endpoint: 1500, correction: 'holm'});
    await fc.call();
    const analysis = fc.outputs['analysis'] as DG.DataFrame;
    const twins = fc.outputs['twins'] as DG.DataFrame;
    expect(analysis.get('twinned', 0), true);
    expect(twins.rowCount > 0, true, 'the twin gave its bands');
    for (let i = 0; i < twins.rowCount; i++) {
      const [lo, mid, hi] = ['twin_q05', 'twin_median', 'twin_q95'].map((k) => twins.get(k, i));
      expect(Number.isFinite(mid) && lo <= mid && mid <= hi, true, `row ${i}: ${lo} ${mid} ${hi}`);
    }
  }, {timeout: 300000});
});
