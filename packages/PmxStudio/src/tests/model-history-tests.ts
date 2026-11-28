import {category, expect, test} from '@datagrok-libraries/test/src/test';

import type {HistoryAnchor, LinePrior} from '../priors';
import {HistorySource, curveRate, historyNotes, historyPrior, historyStudy, otherLab,
  weighHistory} from '../model-history';

// synthetic answers only: no provider's data in tests
const OWN: HistorySource = {sameLab: 'yes', sameHost: 'yes', passage: 'later', frozen: false, meansOnly: false,
  animalsPerGroup: null};

category('Tumor models: own history', () => {
  test('the lab\'s shift applies unless the history is this lab\'s, same strain and site', async () => {
    expect(otherLab(OWN), false);
    expect(otherLab({...OWN, sameLab: 'no'}), true);
    expect(otherLab({...OWN, sameLab: 'unknown'}), true, 'not known counts as another lab');
    expect(otherLab({...OWN, sameHost: 'no'}), true, 'another strain or site');
    expect(otherLab({...OWN, sameHost: 'unknown'}), true, 'a strain or site not known counts as another lab (D694)');
  });

  test('an own pilot shows only where the growth comes from', async () => {
    const notes = historyNotes(OWN, 1, 0.738, 3);
    expect(notes.length, 1);
    expect(notes[0], 'Growth from this model’s own history (1 study), weighted 74% against PMX Studio’s PDX ' +
      'reference.');
  });

  test('a provider\'s sheet of means: every note, in order', async () => {
    const notes = historyNotes({sameLab: 'no', sameHost: 'no', passage: 'early', frozen: true, meansOnly: true,
      animalsPerGroup: null}, 10, 0.87, 3);
    expect(notes.length, 6);
    expect(notes[0].includes('(10 studies), weighted 87%'), true, notes[0]);
    expect(notes[1].startsWith('The history comes from another lab'), true);
    expect(notes[2].startsWith('The history used another mouse strain'), true);
    expect(notes[3].includes('early passage'), true);
    expect(notes[4].includes('frozen'), true);
    expect(notes[5].endsWith('3 per group assumed.'), true, notes[5]);
  });

  test('means with the animals stated assume nothing; an unknown lab says so', async () => {
    const notes = historyNotes({...OWN, sameLab: 'unknown', meansOnly: true, animalsPerGroup: 6}, 2, 0.5, 3);
    expect(notes.some((n) => n.includes('assumed')), false);
    expect(notes[1].startsWith('It is not known whether the history comes from this lab'), true);
  });

  // the reference's formula (history_anchor.py update) on this synthetic anchor: expected_history.py's values
  const A: HistoryAnchor = {mu: -2.5, tau_model: 0.45, tau_study: 0.25, sigma_animal: 0.27, window_volume_mm3: 357,
    tau_lab: 0.15, n_for_group_means: 3, provenance: {}};
  const close = (x: number, y: number, what: string) => expect(Math.abs(x - y) < 1e-9, true, `${what}: ${x} vs ${y}`);

  test('the history is weighed as the reference formula weighs it', async () => {
    const none = weighHistory(A, [], false);
    close(none.weight, 0, 'no history: weight');
    close(none.mean, -2.5, 'no history: the class mean');
    close(none.sd, 0.5147815070, 'no history: sd');
    for (const [what, studies, other, mean, sd, weight] of [
      ['own pilot, 5 animals', [{lnRate: -2, n: 5}], false, -2.1378496316, 0.3439899719, 0.7243007368],
      ['another lab', [{lnRate: -2, n: 5}], true, -2.1648238877, 0.3895557399, 0.6703522246],
      ['ten experiments of means, another lab', Array(10).fill({lnRate: -2, n: 3}), true, -2.0667151660,
        0.3346933555, 0.8665696679]] as [string, {lnRate: number, n: number}[], boolean, number, number, number][]) {
      const w = weighHistory(A, studies, other);
      expect(Math.abs(w.mean - mean) < 1e-8 && Math.abs(w.sd - sd) < 1e-8 && Math.abs(w.weight - weight) < 1e-8,
        true, `${what}: ${JSON.stringify(w)}`);
    }
  });

  test('a curve\'s growth rate: the slope of ln V, cut at a long gap; too short or flat gives none', async () => {
    const curve = (rate: number, days: number[]) => days.map((day) => ({day, volume: 100 * Math.exp(rate * day)}));
    close(curveRate(curve(0.08, [0, 3, 7, 10, 14, 17, 21]))!, 0.08, 'exponential');
    close(curveRate([...curve(0.08, [0, 3, 7, 10]), {day: 40, volume: 5}])!, 0.08, 'cut after a 30-day gap');
    expect(curveRate(curve(0.08, [0, 3, 7])), null, 'three points');
    expect(curveRate(curve(0, [0, 3, 7, 10])), null, 'no growth');
    expect(curveRate([...curve(0.08, [0, 3, 7]), {day: 10, volume: 0}]), null, 'a zero is not a point');
  });

  test('a study: the mean of its animals\' ln rates; group means take the stated or assumed animals', async () => {
    const rows = (animal: string, rate: number) => [0, 4, 8, 12, 16].map((day) =>
      ({animal, day, volume: 80 * Math.exp(rate * day)}));
    const s = historyStudy([...rows('a', 0.05), ...rows('b', 0.2), ...rows('c', 0)], false, null, 3)!;
    close(s.lnRate, (Math.log(0.05) + Math.log(0.2)) / 2, 'ln rate');
    expect(s.n, 2, 'the animal that did not grow is left out');
    expect(historyStudy(rows('exp 1', 0.1), true, null, 3)!.n, 3, 'means: assumed');
    expect(historyStudy(rows('exp 1', 0.1), true, 7, 3)!.n, 7, 'means: stated');
    expect(historyStudy(rows('a', 0), false, null, 3), null, 'nothing grew');
  });

  test('the growth prior: log alpha from the weighed rate, the rest the PDX lines\' mean', async () => {
    const line = (logVRef: number, beta: number, c0: number): LinePrior => ({papers: 9, animals: 90,
      log_alpha: {mixture: [], map_mean: 0, map_sd: 0, map_q025: 0, map_q975: 0, tau: 0.3, ess: 5},
      log_beta: beta, log_v_ref: logVRef, sd_within: {log_alpha: 0.1, log_beta: 0.1, log_v0: 0.5},
      study_prior: {mean_log_alpha: -2.2, sd_log_alpha: 0.5, mean_log_beta: beta, sd_log_beta: 1.4, rho: 0.8},
      noise: {c0, c1: 0.5, k: 1, df: 5, sd_session: 0.05},
      twin: {sd_log_alpha: 0.1, sd_log_beta: 0.1, sd_process: 0.02, first_factor: 1, start_df: 5, slow_share: 0.1,
        slow_shift: -0.5, skew: 0, validation: {loo_coverage_90: 0.9, measurements: 10, below: 1, above: 1},
        conformal: {levels: {}, papers: 1, points: 1, pooled: true, validated: false, source: 'atlas'}},
      start_sd: 0.2, sd_residual_log_v: 0.1,
      take: {nonprogressing_map: 0.05, nonprogressing: 0.05, animals: 90, definition: 'x'},
      facet_shifts_log_rate: [], map_settings: {}});
    const p = historyPrior([line(Math.log(100), -5, 0.1), line(Math.log(144), -4, 0.3)], A,
      {mean: -2.2, sd: 0.3}, 12);
    const sp = p.study_prior!;
    close(sp.mean_log_alpha, -2.2, 'log alpha: the weighed ln rate itself');
    close(sp.sd_log_alpha, 0.3, 'its sd');
    expect(JSON.stringify([sp.mean_log_beta, sp.sd_log_beta, sp.rho, p.log_beta]), '[-14,0.01,0,-14]', 'beta ~ 0');
    expect(JSON.stringify([p.sd_within.log_alpha, p.sd_within.log_beta]), '[0.27,0]', 'animals');
    close(p.noise!.c0, 0.2, 'noise: the lines\' mean');
    expect(p.noise!.sd_session, 0, 'no shift per session (the validated predictive had none)');
    close(p.start_sd, 0.2, 'start spread: the lines\'');
    expect(p.twin?.conformal === undefined, true, 'the lines\' validations are dropped');
    expect(p.animals, 12);
    expect(p.papers, 0);
    // beta ~ 0 makes the engine's closed form (Crux design.rs log_volume) the exponential over four weeks
    const [alpha, beta, u0, q] = [Math.exp(-2.2), Math.exp(-14), Math.log(150), p.log_v_ref];
    const top = q + alpha / beta;
    expect(Math.abs(top + (u0 - top) * Math.exp(-beta * 28) - (u0 + 28 * alpha)) < 1e-3, true, 'exponential');
  });

  test('a plan longer than four weeks says the growth is assumed sustained', async () => {
    expect(historyNotes(OWN, 1, 0.7, 3, 28).length, 1, '28 days: no note');
    const notes = historyNotes(OWN, 1, 0.7, 3, 35);
    expect(notes.length, 2);
    expect(notes[1].includes('up to about 4 weeks'), true, notes[1]);
  });
});
