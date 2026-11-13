/**
 * Occasion (study day / period) stratification of the derived analyses — topic
 * tk-study-day-stratification (F0–F6). Every expected number is read from the regenerated
 * PKNCA / R oracle CSVs in `files/occasion/` (rule 18: no transcribed literals); see
 * `files/REGEN.md` "Occasion fixtures" and
 * `docs/_internal/reference/pknca-validation/scripts/demo_study4_tk28d_dose_proportionality.R`.
 */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  isOccasionName, detectOccasionColumn, orderOccasionLevels, stripDayCode, resolveOccasion, splitsProfiles,
  computeDoseProportionality, computeBioavailability, isEntryStale, isDerivedBasisStale,
  computeRecipeFingerprint, computeRecipeHash, numAt, DP_PARAMETER_ID, VERDICT_DIFFERS_NOTE,
  appendAucTauColumn, AUC_TAU_COL, aggregateNcaResults, computeAccumulationRatios, computeSexFold,
} from '../adapter';
import type {ProfileRow, DoseProportionalityFit, DeltaBetaRow, AnalysisEntry} from '../adapter';
import {DERIVED_BASIS} from '../adapter/study-schema';
import {detectGroupByColumns} from '../pipeline/pipeline-mapping';
import {buildDemoStudy} from '../pipeline/demo-seed';
import {DEMO_STUDIES, DEMO_STUDIES_DIR} from '../adapter/demo-studies';
import {loadSampleByName} from '../views/nca-launch';

const DEMO_DEPS = {
  loadTable: (f: string) => loadSampleByName(`${DEMO_STUDIES_DIR}/${f}`),
  now: () => new Date().toISOString(),
};
import {summarizePipeline} from '../pipeline/pipeline-summary';

const DIR = 'System:AppData/Nca/occasion';

async function csv(name: string): Promise<DG.DataFrame> {
  return DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(`${DIR}/${name}`));
}

const relClose = (a: number, b: number, tol: number): boolean =>
  Math.abs(a - b) <= tol * Math.max(1e-12, Math.abs(b));

/** A string cell; '' when the column is absent (a single-route fixture carries no Route). */
const s = (df: DG.DataFrame, col: string, i: number): string => df.col(col)?.getString(i) ?? '';
const n = (df: DG.DataFrame, col: string, i: number): number => numAt(df.col(col)!, i);

/** The oracle per-profile frame as a "results" frame the adapter reads. */
async function oracleProfiles(name: string): Promise<DG.DataFrame> {
  return csv(name);
}

const DP_OPTS = {
  doseCol: 'Dose_mg_kg', sexCol: 'Sex', routeCol: 'Route', aucInfCol: 'aucinf.obs', aucLastCol: 'auclast',
  aucTauCol: 'auc0tau', cmaxCol: 'cmax', pctExtrapCol: 'aucpext.obs', tlastCol: 'tlast', subjectCol: 'SubjectID',
  occasionCol: 'StudyDay', tauH: 24,
} as const;

const ORACLE_BASIS: Record<string, string> = {AUCinf: DP_PARAMETER_ID.aucinf, AUC0tau: DP_PARAMETER_ID.auctau,
  AUClast: DP_PARAMETER_ID.auclast, Cmax: DP_PARAMETER_ID.cmax};
const fitId = (f: DoseProportionalityFit): string =>
  f.parameter === 'Cmax' ? DP_PARAMETER_ID.cmax : DP_PARAMETER_ID[f.basis!];

function animalsOf(df: DG.DataFrame): string[] {
  return Array.from({length: df.rowCount}, (_, i) => s(df, 'Animal', i));
}

// ── F1 — detection, animal, design ─────────────────────────────────────────────

category('Occasion: detection and animal (F1)', () => {
  test('one name list, both classes; a time-axis Day is never an occasion (AC-1.7, AC-1.8)', async () => {
    for (const nm of ['StudyDay', 'Study_Day', 'DAY', 'Day', 'VISITDY', 'VISITDAY', 'VISIT', 'NOMDY', 'Week',
      'Period', 'Occasion', 'OCC'])
      expect(isOccasionName(nm), true, nm);
    expect(isOccasionName('Dose_mg_kg'), false, 'dose');
    // Only profile-constant columns are candidates: a Day absent from the carried list is never picked.
    expect(detectOccasionColumn(['Sex', 'Dose_mg_kg']), null, 'no candidate');
    expect(detectOccasionColumn(['Period', 'StudyDay']), 'StudyDay', 'day-like name first');
    expect(detectOccasionColumn(['StudyDay', 'Week'], 'Week'), 'Week', 'the analyst pick wins');
  });

  test('levels: numeric order after a D/Day/W prefix, else first appearance', async () => {
    expect(orderOccasionLevels(['28', '1', '28']).levels.join(','), '1,28', 'numeric');
    expect(orderOccasionLevels(['D28', 'Day 1']).levels.join(','), 'Day 1,D28', 'prefixed');
    const r = orderOccasionLevels(['late', 'early']);
    expect(r.levels.join(',') + String(r.numeric), 'late,earlyfalse', 'first appearance');
  });

  test('rule (iii) day-coded id (AC-1.4)', async () => {
    expect(stripDayCode('T001_D28', '28'), 'T001', 'match');
    expect(stripDayCode('T001_D1', '2'), null, 'token != occasion');
    expect(stripDayCode('T001_DX', '1'), null, 'non-numeric token');
    expect(stripDayCode('R12 Day 1', 'Day 1'), 'R12', 'prefixed occasion');
  });

  test('AC-1.1: demo study 4 — StudyDay, rule (iii), repeat-dose; animal counts from the generator', async () => {
    const pp = await oracleProfiles('demo_study4_tk28d_nca_per_profile.csv');
    const rows: ProfileRow[] = Array.from({length: pp.rowCount}, (_, i) => ({
      subject: s(pp, 'SubjectID', i), occasion: s(pp, 'StudyDay', i), sex: s(pp, 'Sex', i),
      dose: n(pp, 'Dose_mg_kg', i), route: 'PO'}));
    const res = resolveOccasion('StudyDay', rows);
    expect(res.levels.join(','), '1,28', 'levels');
    expect(res.animal.rule, 'day-coded-id', 'rule (iii)');
    expect(res.animal.failing.length, 0, 'every animal verifies');
    expect(res.proposal.design, 'repeat-dose', 'proposed design');
    expect(res.stratify, true, 'DP / sex fold per occasion');
    expect(res.fPctPairing, 'across-occasions', 'F% waits for confirmation');
    // The oracle's own Animal column (rule iii in R) must equal ours, row for row.
    expect(res.animal.animals!.join(','), animalsOf(pp).join(','), 'animal per profile = oracle');
    const perSex = new Map<string, Set<string>>();
    res.animal.animals!.forEach((a, i) => perSex.set(rows[i].sex, (perSex.get(rows[i].sex) ?? new Set()).add(a)));
    // R2 F8: reconciled against the generator — the dosed animals per sex (the oracle drops vehicle).
    const doses = new Set(rows.map((r) => r.dose));
    expect(perSex.get('F')!.size, 3 * doses.size, 'F dosed animals = 3 per dose level');
    expect(perSex.get('M')!.size, 3 * doses.size, 'M dosed animals = 3 per dose level');
    const confirmed = resolveOccasion('StudyDay', rows, {occasionDesign: 'repeat-dose'});
    expect(confirmed.fPctPairing, 'same-occasion', 'confirmed → same-occasion F%');
  });

  test('AC-1.3 / AC-1.5: crossover, escalation, unknown, failing animals', async () => {
    const g = await csv('occasion_g_crossover_iv_po_per_profile.csv');
    const rowsOf = (df: DG.DataFrame, occ: string): ProfileRow[] => Array.from({length: df.rowCount}, (_, i) => ({
      subject: s(df, 'SubjectID', i), occasion: s(df, occ, i), sex: s(df, 'Sex', i),
      dose: n(df, 'Dose_mg_kg', i), route: s(df, 'Route', i)}));
    expect(resolveOccasion('Period', rowsOf(g, 'Period')).proposal.design, 'crossover', 'IV/PO crossover');
    const f = await csv('occasion_f_dose_escalation_per_profile.csv');
    const ef = resolveOccasion('StudyDay', rowsOf(f, 'StudyDay'));
    expect(ef.proposal.design + String(ef.stratify), 'escalationfalse', 'escalation pools DP');
    const i = await csv('occasion_i_crossover_po_po_per_profile.csv');
    const ui = resolveOccasion('Period', rowsOf(i, 'Period'));
    expect(ui.proposal.design + ui.fPctPairing, 'unknownacross-occasions', 'same-route PO/PO crossover stays unknown');
    const noAnimal: ProfileRow[] = [
      {subject: 'A', occasion: '1', sex: 'F', dose: 1, route: 'PO'},
      {subject: 'B', occasion: '28', sex: 'F', dose: 1, route: 'PO'}];
    expect(resolveOccasion('StudyDay', noAnimal).proposal.design, 'unknown', 'no animal key');
    const dup: ProfileRow[] = [
      {subject: 'X', occasion: '1', sex: 'F', dose: 1, route: 'PO', animal: 'A1'},
      {subject: 'Y', occasion: '1', sex: 'F', dose: 1, route: 'PO', animal: 'A1'},
      {subject: 'Z', occasion: '28', sex: 'F', dose: 1, route: 'PO', animal: 'A2'},
      {subject: 'W', occasion: '1', sex: 'F', dose: 1, route: 'PO', animal: 'A2'}];
    const rd = resolveOccasion('StudyDay', dup);
    expect(rd.animal.rule + rd.animal.failing.map((x) => x.animal).join(','), 'mappedA1', 'A1 fails, key stays');
  });

  test('AC-1.2: repeated subject id (rule ii) and NONMEM DAY + AnimalID (rule i)', async () => {
    const d = await csv('occasion_d_repeated_subject_per_profile.csv');
    const rows: ProfileRow[] = Array.from({length: d.rowCount}, (_, i) => ({
      subject: s(d, 'SubjectID', i), occasion: s(d, 'StudyDay', i), sex: s(d, 'Sex', i),
      dose: n(d, 'Dose_mg_kg', i), route: 'PO'}));
    expect(resolveOccasion('StudyDay', rows).animal.rule, 'repeated-subject', 'rule (ii)');
    const nm = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText('System:AppData/Nca/06_rat_tk28d_nonmem.csv'));
    expect(nm.col('DAY') !== null, true, 'NONMEM fixture has DAY');
    expect(detectOccasionColumn(['DAY', 'AnimalID']), 'DAY', 'DAY matched by the shared list');
  });

  test('AC-1.9: an occasion-named Day becomes a profile key only when it splits profiles', async () => {
    const k = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(`${DIR}/occasion_k_daily_profile.csv`));
    const keys = detectGroupByColumns(k, new Set(['SubjectID', 'Time_hr', 'Conc_ng_mL']),
      {subject: 'SubjectID', time: 'Time_hr'});
    expect(keys.includes('Day'), false, `per-row nominal Day must not split the profile (keys: ${keys})`);
    // splitsProfiles on its own: restarting 0..24 h blocks split; one cumulative series does not.
    expect(splitsProfiles(['a', 'a', 'a', 'a', 'a', 'a'], ['1', '1', '1', '2', '2', '2'], [0, 1, 2, 0, 1, 2]), true,
      'restarting blocks');
    expect(splitsProfiles(['a', 'a', 'a', 'a', 'a', 'a'], ['1', '1', '1', '2', '2', '2'], [0, 1, 2, 24, 25, 26]),
      false, 'cumulative time');
  });

  test('AC-1.6: derived choices never move the Setup fingerprint; they do move the hash', async () => {
    const base = {
      columns: {} as never, modes: {} as never, ncaRules: {} as never, units: {} as never,
      groupByCols: [], passthroughCols: [], acceptanceRSquared: 0.9,
    };
    const picked = {...base, derived: {occasionCol: 'StudyDay', tauH: 24}};
    expect(computeRecipeFingerprint(picked), computeRecipeFingerprint(base), 'fingerprint excludes derived');
    expect(computeRecipeHash(picked) !== computeRecipeHash(base), true, 'hash includes derived');
    expect(computeRecipeHash({...base, derived: undefined}), computeRecipeHash(base), 'absent derived = old hash');
  });
});

// ── F2 / F3a / F4 — dose proportionality per occasion ───────────────────────────

function oracleRow(dp: DG.DataFrame, sex: string, day: string, basis: string): number {
  for (let i = 0; i < dp.rowCount; i++)
    if (s(dp, 'Sex', i) === sex && s(dp, 'StudyDay', i) === day && s(dp, 'basis', i) === basis) return i;
  return -1;
}

async function checkDpAgainstOracle(stem: string, perProfile = `${stem}_per_profile.csv`): Promise<string[]> {
  const pp = await csv(perProfile);
  const dp = await csv(`${stem}_dp_per_day.csv`);
  const res = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp)});
  const checked: string[] = [];
  for (const f of res.fits) {
    const basis = Object.keys(ORACLE_BASIS).find((k) => ORACLE_BASIS[k] === fitId(f))!;
    const i = oracleRow(dp, f.sex, f.occasion, basis);
    expect(i >= 0, true, `oracle row for ${f.sex}/${f.occasion}/${basis}`);
    const tag = `${stem} ${f.sex} Day ${f.occasion} ${basis}`;
    expect(f.n, n(dp, 'n', i), `${tag}: N`);
    expect(f.nByDose, s(dp, 'n_by_dose', i), `${tag}: N per dose`);
    expect(f.df, n(dp, 'df', i), `${tag}: df`);
    expect(relClose(f.beta, n(dp, 'beta', i), 1e-6), true, `${tag}: β ${f.beta} vs ${n(dp, 'beta', i)}`);
    expect(relClose(f.slopeCI[0], n(dp, 'ci_low', i), 1e-6), true, `${tag}: CI low`);
    expect(relClose(f.slopeCI[1], n(dp, 'ci_high', i), 1e-6), true, `${tag}: CI high`);
    expect(f.carriesVerdict, s(dp, 'carries_verdict', i).toUpperCase() === 'TRUE', `${tag}: carries the verdict`);
    const oracleVerdict = s(dp, 'verdict', i);
    const ours = f.verdict === 'inconclusive-underpowered' ? 'inconclusive' : f.verdict;
    expect(ours, oracleVerdict, `${tag}: verdict class`);
    const differs = s(dp, 'verdict_differs', i).toUpperCase() === 'TRUE' && f.carriesVerdict;
    expect(f.notes.includes(VERDICT_DIFFERS_NOTE), differs, `${tag}: '${VERDICT_DIFFERS_NOTE}'`);
    checked.push(tag);
  }
  return checked;
}

function checkDeltaAgainstOracle(rows: ReadonlyArray<DeltaBetaRow>, db: DG.DataFrame, tag: string): void {
  for (let i = 0; i < db.rowCount; i++) {
    const r = rows.find((x) => x.sex === s(db, 'Sex', i))!;
    expect(r !== undefined && r.method === 'paired', true, `${tag} ${s(db, 'Sex', i)}: paired`);
    expect(r.nPaired, n(db, 'n_paired', i), `${tag}: paired N`);
    expect(relClose(r.deltaBeta, n(db, 'paired_delta_beta', i), 1e-6), true,
      `${tag}: Δβ ${r.deltaBeta} vs ${n(db, 'paired_delta_beta', i)}`);
    expect(relClose(r.ci[0], n(db, 'paired_ci_low', i), 1e-6), true, `${tag}: CI low`);
    expect(relClose(r.ci[1], n(db, 'paired_ci_high', i), 1e-6), true, `${tag}: CI high`);
    expect(r.df, n(db, 'paired_df', i), `${tag}: df`);
  }
}

category('Occasion: dose proportionality per occasion (F2, F3a, F4)', () => {
  test('AC-2.1(a) / AC-4.1: demo 4 per (Sex, Day, basis) on the oracle AUCs', async () => {
    const checked = await checkDpAgainstOracle('demo_study4_tk28d', 'demo_study4_tk28d_nca_per_profile.csv');
    // Day 1 (AUC∞ when the guard passes, AUC0–τ beside), Day 28 (AUC0–τ, AUC∞ beside), Cmax — × 2 sexes.
    expect(checked.length, 12, `fits checked: ${checked.join('; ')}`);
  });

  test('AC-2.1(b): demo 4 on the ENGINE\'s own AUCs — β within 0.002 and the same verdict class', async () => {
    const [b] = await buildDemoStudy(DEMO_DEPS, DEMO_STUDIES.find((x) => x.ref === 4)!, {review: true});
    const r = b.run.outcome.results;
    const dp = await csv('demo_study4_tk28d_dp_per_day.csv');
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    // AUC0–τ ≡ AUClast here: every dosed profile is sampled to exactly τ = 24 h (the oracle's tlast).
    for (let i = 0; i < pp.rowCount; i++) expect(n(pp, 'tlast', i), 24, `${s(pp, 'SubjectID', i)}: tlast = τ`);
    const rows: ProfileRow[] = Array.from({length: r.rowCount}, (_, i) => ({subject: s(r, 'SUBJECT', i),
      occasion: s(r, 'StudyDay', i), sex: s(r, 'Sex', i), dose: n(r, 'Dose_mg_kg', i), route: 'PO'}));
    const occ = resolveOccasion('StudyDay', rows);
    const res = computeDoseProportionality(r, {doseCol: 'Dose_mg_kg', sexCol: 'Sex', occasionCol: 'StudyDay',
      occasionLevels: occ.levels, aucTauCol: 'AUC_LAST', tauH: 24, animals: occ.animal.animals});
    const deviations: string[] = [];
    for (const f of res.fits) {
      const basis = Object.keys(ORACLE_BASIS).find((k) => ORACLE_BASIS[k] === fitId(f))!;
      const i = oracleRow(dp, f.sex, f.occasion, basis);
      const tag = `${f.sex} Day ${f.occasion} ${basis}`;
      const d = Math.abs(f.beta - n(dp, 'beta', i));
      deviations.push(`${tag} Δ=${d.toExponential(2)}`);
      expect(d <= 0.002, true, `${tag}: β ${f.beta} vs oracle ${n(dp, 'beta', i)}`);
      if (s(dp, 'knife_edge', i).toUpperCase() !== 'TRUE') {
        const ours = f.verdict === 'inconclusive-underpowered' ? 'inconclusive' : f.verdict;
        expect(ours, s(dp, 'verdict', i), `${tag}: verdict class`);
      }
    }
    return deviations.join('; ');
  }, {timeout: 300000});

  test('AC-3.1: demo 4 paired Δβ per sex = R lm(d ~ log(dose))', async () => {
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    const res = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp)});
    checkDeltaAgainstOracle(res.deltaBeta, await csv('demo_study4_tk28d_delta_beta.csv'), 'demo 4');
    // Without an animal key: two independent fits with the Welch–Satterthwaite df (R twofit_*).
    const two = computeDoseProportionality(pp, {...DP_OPTS});
    const db = await csv('demo_study4_tk28d_delta_beta.csv');
    for (let i = 0; i < db.rowCount; i++) {
      const r = two.deltaBeta.find((x) => x.sex === s(db, 'Sex', i))!;
      expect(r.method, 'two-fits', 'two fits');
      expect(relClose(r.df, n(db, 'twofit_df', i), 1e-6), true, `WS df ${r.df}`);
      expect(relClose(r.ci[1], n(db, 'twofit_ci_high', i), 1e-6), true, 'two-fit CI high');
    }
  });

  test('AC-3.3: attrition note fires only when the drop-out is concentrated', async () => {
    const variants = [['occasion_c_attrition', false], ['occasion_c2_attrition_topdose', true]] as const;
    for (const [stem, concentrated] of variants) {
      const pp = await csv(`${stem}_per_profile.csv`);
      const res = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp)});
      checkDeltaAgainstOracle(res.deltaBeta, await csv(`${stem}_delta_beta.csv`), stem);
      const m = res.deltaBeta.find((r) => r.sex === 'M')!;
      expect(m.notes.some((x) => x.startsWith('paired N')), true, `${stem}: paired N disclosed`);
      expect(m.notes.some((x) => x.includes('concentrated at')), concentrated, `${stem}: concentration note`);
    }
    await checkDpAgainstOracle('occasion_c_attrition');
  });

  test('AC-2.6: without an occasion the DP table keeps its single-basis layout', async () => {
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    const res = computeDoseProportionality(pp, {doseCol: 'Dose_mg_kg', sexCol: 'Sex', aucInfCol: 'aucinf.obs',
      cmaxCol: 'cmax'});
    expect(res.table.col('OCCASION'), null, 'no occasion column');
    expect(res.table.col('BASIS') !== null, true, 'legacy BASIS column');
    expect(res.deltaTable, null, 'no Δβ table');
    expect(res.fits.every((f) => f.carriesVerdict && f.occasion === ''), true, 'one verdict per stratum');
  });

  test('Δβ suppressed when no common basis exists (different sampling spans)', async () => {
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    // Pretend Day 28 was sampled to 48 h: no AUC0–τ, unequal spans → not comparable.
    const t = pp.col('tlast')!;
    const day = pp.col('StudyDay')!;
    for (let i = 0; i < pp.rowCount; i++) if (day.getString(i) === '28') t.set(i, 48);
    const res = computeDoseProportionality(pp, {...DP_OPTS, aucTauCol: undefined, animals: animalsOf(pp)});
    expect(res.deltaBeta.every((r) => r.method === 'suppressed'), true, 'suppressed');
    expect(res.deltaBeta[0].notes[0], 'not comparable: different sampling spans', 'label');
  });

  test('AUC0–τ per stratum: one profile not spanning τ → the whole stratum on labelled AUClast', async () => {
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    const tau = pp.col('auc0tau')!;
    const day = pp.col('StudyDay')!;
    let blanked = '';
    for (let i = 0; i < pp.rowCount && blanked === ''; i++) {
      if (day.getString(i) === '28' && s(pp, 'Sex', i) === 'F' && n(pp, 'Dose_mg_kg', i) > 0) {
        tau.set(i, null);
        blanked = s(pp, 'SubjectID', i);
      }
    }
    const res = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp)});
    const carried = (sex: string, occ: string): DoseProportionalityFit => res.fits.find((f) =>
      f.parameter === 'AUC' && f.sex === sex && f.occasion === occ && f.carriesVerdict)!;
    const f28 = carried('F', '28');
    expect(f28.basis, 'auclast', `F Day 28 (${blanked} not spanned) carries AUClast, not the spanned subset`);
    expect(f28.notes.some((x) => x.startsWith('AUClast: τ not spanned by 1 of ')), true, `note: ${f28.notes}`);
    expect(res.fits.some((f) => f.sex === 'F' && f.occasion === '28' && f.basis === 'auctau'), false,
      'no AUC0–τ fit on the partly spanned stratum');
    expect(carried('M', '28').basis, 'auctau', 'M Day 28 untouched');
  });

  test('sex homogeneity is withheld when the sexes carry different AUC bases (names the day)', async () => {
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    const base = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp)});
    expect(base.warnings.some((w) => w.code === 'SEX_BASIS_DIFFERS'), false, 'same bases: no warning');
    // One F Day-1 animal fails the AUC∞ guard → F Day 1 carries AUC0–τ, M Day 1 keeps AUC∞.
    const ext = pp.col('aucpext.obs')!;
    for (let i = 0; i < pp.rowCount; i++) {
      if (s(pp, 'StudyDay', i) === '1' && s(pp, 'Sex', i) === 'F' && n(pp, 'Dose_mg_kg', i) > 0) {
        ext.set(i, 30);
        break;
      }
    }
    const res = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp)});
    const w = res.warnings.filter((x) => x.code === 'SEX_BASIS_DIFFERS');
    expect(w.length, 1, 'one warning');
    expect(w[0].stratum.endsWith(', 1'), true, `stratum names the day: ${w[0].stratum}`);
    expect(res.sexHomogeneity.some((h) => h.occasion === '1'), false, 'no Day-1 sex claim');
  });

  test('AUC0–τ: a BLQ trough at τ (TLAST < τ, sampled to τ) takes the interval AUC; sampling stopped short → missing',
    async () => {
      // Rows: TLAST = τ (identity) · TLAST 12 sampled to 24 (BLQ trough) · TLAST 12 sampled to 12 · TLAST 30.
      const r = DG.DataFrame.fromColumns([
        DG.Column.fromFloat64Array('TLAST', new Float64Array([24, 12, 12, 30])),
        DG.Column.fromFloat64Array('AUC_LAST', new Float64Array([100, 80, 70, 120])),
      ]);
      const lastSample = [24, 24, 12, 30];
      const interval = [NaN, 81, 71, 110];
      appendAucTauColumn(r, 24, (i) => interval[i], (i) => lastSample[i]);
      const v = (i: number): number => numAt(r.col(AUC_TAU_COL)!, i);
      expect(v(0), 100, 'TLAST = τ: AUClast');
      expect(v(1), 81, 'BLQ trough at τ: interval AUC');
      expect(Number.isNaN(v(2)), true, 'sampling stopped before τ: missing');
      expect(v(3), 110, 'sampled past τ: interval AUC');
    });

  test('a confirmed crossover: no AUC0–τ, first-occasion rule on every period, Δβ labelled a period effect',
    async () => {
      const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
      const res = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp), singleDoseOccasions: true});
      expect(res.fits.some((f) => f.basis === 'auctau'), false, 'no AUC0–τ fit');
      for (const f of res.fits.filter((x) => x.parameter === 'AUC' && x.carriesVerdict))
        expect(f.basis === 'aucinf' || f.basis === 'auclast', true, `${f.sex} ${f.occasion}: ${f.basis}`);
      expect(res.deltaBeta.length > 0 && res.deltaBeta.every((r) =>
        r.notes.includes('period effect (single-dose periods), not accumulation')), true, 'Δβ label');
    });
});

// ── F2 — F%, sex comparison (cell summaries) ───────────────────────────────────

category('Occasion: F% and cell summaries (F2)', () => {
  test('AC-2.2: dose-normalised AUC∞ per day equals the oracle cell summary', async () => {
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    const cells = await csv('demo_study4_tk28d_cell_summary.csv');
    const res = computeDoseProportionality(pp, {...DP_OPTS, animals: animalsOf(pp)});
    const dn = res.dnTable;
    expect(dn.col('SUBJECT') !== null && dn.col('OCCASION') !== null && dn.col('ANIMAL') !== null, true,
      'DN table identifies each row (AC-2.5)');
    for (let i = 0; i < cells.rowCount; i++) {
      const sx = s(cells, 'Sex', i);
      const d = n(cells, 'Dose_mg_kg', i);
      const day = s(cells, 'StudyDay', i);
      const vals: number[] = [];
      for (let r = 0; r < dn.rowCount; r++) {
        if (s(dn, 'SEX', r) === sx && n(dn, 'DOSE', r) === d && s(dn, 'OCCASION', r) === day)
          vals.push(n(dn, 'DN_AUC', r));
      }
      expect(vals.length, n(cells, 'n_aucinf', i), `${sx} ${d} Day ${day}: N`);
      const mean = vals.reduce((a, b) => a + b, 0) / vals.length;
      expect(relClose(mean, n(cells, 'dn_aucinf', i), 1e-9), true, `${sx} ${d} Day ${day}: dn-AUC ${mean}`);
    }
  });

  test('AC-2.3 / AC-4.5: demo 4 engine run — sex comparison per day and declared-τ Rac equal the oracle',
    async () => {
      const [b] = await buildDemoStudy(DEMO_DEPS, DEMO_STUDIES.find((x) => x.ref === 4)!, {review: true});
      const r = b.run.outcome.results;
      const cells = await csv('demo_study4_tk28d_cell_summary.csv');
      // AC-2.3: the view's sex-comparison aggregate, [dose, sex, occasion] (rule 18: AUC within 0.1 %).
      const agg = aggregateNcaResults(r, ['Dose_mg_kg', 'Sex', 'StudyDay']);
      // AC-4.5: declared τ = 24 h; every profile is sampled to τ, so AUC0–τ ≡ AUClast.
      const acc = computeAccumulationRatios(r, {tau: 24, doseCol: 'Dose_mg_kg', studyDayCol: 'StudyDay',
        stratifyCols: ['Sex'], aucTauCol: 'AUC_LAST'});
      for (let i = 0; i < cells.rowCount; i++) {
        const sx = s(cells, 'Sex', i);
        const d = n(cells, 'Dose_mg_kg', i);
        const day = s(cells, 'StudyDay', i);
        const tag = `${sx} ${d} Day ${day}`;
        let a = -1;
        for (let k = 0; k < agg.rowCount; k++)
          if (s(agg, 'Sex', k) === sx && n(agg, 'Dose_mg_kg', k) === d && s(agg, 'StudyDay', k) === day) a = k;
        expect(a >= 0, true, `${tag}: aggregate row`);
        expect(n(agg, 'N', a), n(cells, 'n_auclast', i), `${tag}: N`);
        // The engine reports AUC in mg·h/L (canonical units); the oracle in ng·h/mL = 1000 × mg·h/L.
        const meanNg = n(agg, 'MEAN_AUC_LAST', a) * 1000;
        expect(relClose(meanNg, n(cells, 'mean_auclast', i), 1e-3), true,
          `${tag}: mean AUClast ${meanNg} vs ${n(cells, 'mean_auclast', i)} ng·h/mL`);
        expect(Math.abs(n(agg, 'CV_PCT_AUC_LAST', a) - n(cells, 'cv_auclast', i)) <= 0.1, true,
          `${tag}: %CV ${n(agg, 'CV_PCT_AUC_LAST', a)} vs ${n(cells, 'cv_auclast', i)}`);
        let c = -1;
        for (let k = 0; k < acc.table.rowCount; k++)
          if (s(acc.table, 'Sex', k) === sx && n(acc.table, 'DOSE', k) === d) c = k;
        expect(c >= 0, true, `${tag}: Rac row`);
        expect(relClose(n(acc.table, 'RAC_AUC', c), n(cells, 'rac_auc0tau', i), 1e-3), true,
          `${tag}: Rac ${n(acc.table, 'RAC_AUC', c)} vs ${n(cells, 'rac_auc0tau', i)}`);
      }
    }, {timeout: 300000});

  test('AC-2.4: a single-dose two-sex study gets a sex comparison (no study-day gate)', async () => {
    // Demo 5 (single-dose IV infusion, F and M): the view's sex-comparison path with no occasion.
    const [b] = await buildDemoStudy(DEMO_DEPS, DEMO_STUDIES.find((x) => x.ref === 5)!, {review: true});
    const r = b.run.outcome.results;
    const sexCol = ['Sex', 'SEX'].find((c) => r.col(c) !== null)!;
    const doseCol = ['Dose_mg_kg', 'DOSE'].find((c) => r.col(c) !== null)!;
    expect(new Set(Array.from({length: r.rowCount}, (_, i) => s(r, sexCol, i))).size, 2, 'two sexes');
    const fold = computeSexFold(aggregateNcaResults(r, [doseCol, sexCol]), {sexCol});
    expect(fold.table.rowCount > 0, true, 'sex comparison rows');
  }, {timeout: 300000});

  test('AC-4.2: τ unknown → the later-day verdict rides AUClast, labelled', async () => {
    const pp = await csv('demo_study4_tk28d_nca_per_profile.csv');
    const res = computeDoseProportionality(pp, {...DP_OPTS, aucTauCol: undefined, tauH: undefined,
      animals: animalsOf(pp)});
    for (const sx of ['F', 'M']) {
      const f = res.fits.find((x) => x.parameter === 'AUC' && x.sex === sx && x.occasion === '28' && x.carriesVerdict)!;
      expect(f.basis, 'auclast', `${sx} Day 28 basis`);
      expect(f.notes.includes('AUClast: τ unknown — AUC0–τ unavailable'), true, `${sx} Day 28 label: ${f.notes}`);
    }
    const ids = res.table.col('PARAMETER')!;
    expect(Array.from({length: res.table.rowCount}, (_, i) => ids.getString(i)).includes(DP_PARAMETER_ID.auclast),
      true, 'AUClast is the row identity in the verdict table');
  });

  test('P6: an EV cell with an IV reference only on another occasion takes it, with a warning', async () => {
    const e = await csv('occasion_e_iv_day1_only_per_profile.csv');
    const res = computeBioavailability(e, {doseCol: 'Dose_mg_kg', aucInfCol: 'aucinf.obs', routeCol: 'Route',
      stratifyCols: ['StudyDay', 'Sex'], occasionCol: 'StudyDay'});
    const w = res.warnings.filter((x) => x.code === 'CROSS_OCCASION_IV_REFERENCE');
    expect(w.length > 0, true, 'cross-occasion warning raised');
    expect(w[0].message.includes('StudyDay 28') && w[0].message.includes('StudyDay 1'), true, w[0].message);
    const t = res.table;
    for (let i = 0; i < t.rowCount; i++)
      expect(s(t, 'IV_OCCASION', i), '1', `row ${i}: the Day-1 IV reference`);
  });
});

// ── F5 — Pipeline F% pairing ────────────────────────────────────────────────────

category('Occasion: Pipeline F% (F5)', () => {
  test('AC-2.7: crossover pairs across periods; an unknown design pairs across with a note', async () => {
    const mk = (rows: [string, string, string, string, number, number][], occ = 'Period'): DG.DataFrame => {
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromStrings(occ, rows.map((r) => r[0])),
        DG.Column.fromStrings('SUBJECT', rows.map((r) => r[1])),
        DG.Column.fromStrings('SEX', rows.map((r) => r[2])),
        DG.Column.fromStrings('ROUTE', rows.map((r) => r[3])),
        DG.Column.fromFloat64Array('DOSE', Float64Array.from(rows.map((r) => r[4]))),
        DG.Column.fromFloat64Array('AUC_INF', Float64Array.from(rows.map((r) => r[5]))),
        DG.Column.fromStrings('VERDICT', rows.map(() => 'OK')),
      ]);
      df.setTag('nca.keyColumns', JSON.stringify([occ]));
      return df;
    };
    // IV/PO crossover: A and B get IV in one period, PO in the other.
    const cross = summarizePipeline(mk([['1', 'A', 'F', 'IV bolus', 1, 10], ['2', 'A', 'F', 'PO', 1, 5],
      ['1', 'B', 'F', 'PO', 1, 4], ['2', 'B', 'F', 'IV bolus', 1, 8]]));
    let poF = NaN;
    for (let i = 0; i < cross.rowCount; i++) if (s(cross, 'ROUTE', i) === 'PO') poF = n(cross, 'F_PCT', i);
    expect(Number.isFinite(poF), true, 'crossover F% estimated across periods');
    // Same-route PO/PO with a Period: design unknown → across-occasion pairing + a note.
    const unk = summarizePipeline(mk([['1', 'A', 'F', 'PO', 1, 5], ['2', 'A', 'F', 'PO', 1, 6],
      ['1', 'C', 'F', 'IV bolus', 1, 10]]));
    const notes = unk.getTag('nca.fPctNotes');
    expect(notes !== null && notes !== undefined && notes.includes('could not be classified'), true, `note: ${notes}`);
    // The note is also a visible column on the extravascular arms; the crossover summary has none.
    for (let i = 0; i < unk.rowCount; i++) {
      const ev = s(unk, 'ROUTE', i) === 'PO';
      expect(s(unk, 'F_PCT_NOTE', i).includes('could not be classified'), ev, `F_PCT_NOTE row ${i}`);
    }
    expect(cross.col('F_PCT_NOTE'), null, 'no note column on a classified crossover');

    // Repeat-dose (rule ii: the subject repeats across StudyDay, route and dose constant):
    // same-occasion pairing — Day 1 PO pairs Day 1 IV (5/10), Day 28 PO pairs Day 28 IV (8/20).
    const rep = summarizePipeline(mk([['1', 'A', 'F', 'PO', 1, 5], ['28', 'A', 'F', 'PO', 1, 8],
      ['1', 'C', 'F', 'IV bolus', 1, 10], ['28', 'C', 'F', 'IV bolus', 1, 20]], 'StudyDay'));
    const fOn = (day: string): number => {
      for (let i = 0; i < rep.rowCount; i++)
        if (s(rep, 'ROUTE', i) === 'PO' && s(rep, 'StudyDay', i) === day) return n(rep, 'F_PCT', i);
      return NaN;
    };
    expect(relClose(fOn('1'), 50, 1e-9), true, `Day 1 F% ${fOn('1')}`);
    expect(relClose(fOn('28'), 40, 1e-9), true, `Day 28 F% ${fOn('28')}`);
  });
});

// ── F6 — narrowed stale cue ─────────────────────────────────────────────────────

category('Occasion: derived-basis staleness (F6)', () => {
  const tier = (f: number | null): AnalysisEntry['summaryCache'] => ({
    analysisTier: {} as never, compoundTier: [{fPct: f === null ? null : {mean: f, cvPct: 1, n: 3}} as never],
  } as never);
  const base = {engineVersion: undefined, resultsBasis: undefined, hasCache: true, columns: {} as never,
    passthroughCols: []};

  test('AC-6.1: legacy multi-occasion %F is stale; single-occasion and fresh entries are not', async () => {
    expect(isDerivedBasisStale({...base, groupByCols: ['Dose_mg_kg', 'StudyDay'], summaryCache: tier(52)}), true,
      'legacy demo-4-like analysis with %F');
    expect(isDerivedBasisStale({...base, groupByCols: ['Dose_mg_kg'], summaryCache: tier(52)}), false,
      'no occasion named');
    expect(isDerivedBasisStale({...base, groupByCols: ['StudyDay'], summaryCache: tier(null)}), false, 'no %F');
    expect(isDerivedBasisStale({...base, groupByCols: ['StudyDay'], summaryCache: tier(52),
      derivedBasis: DERIVED_BASIS}), false, 'fresh stamp');
    // The ONE predicate every reader calls: isEntryStale carries the new term.
    expect(isEntryStale({engineVersion: 'x', resultsBasis: 'y', groupByCols: ['StudyDay'], columns: {} as never,
      passthroughCols: [], summaryCache: tier(52)} as never), true, 'isEntryStale reads it');
  });
});
