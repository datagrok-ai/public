/**
 * Hub demo studies (nca-hub-demo-studies F2 / F3 / F6): the planted verdicts on the shipped engine
 * through the seed's own build path (the authoritative seed screen, rule 16), the headless
 * workflow-state split, the seed protocol on a memory repository, and one live EMS run with the
 * reopen invariance (probe B1) and the View-only reading.
 */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  DEMO_STUDIES, DEMO_STUDIES_DIR, DEMO_THRESHOLDS, DEMO_STRUCTURE_COL, DEMO_F_ORDER,
  DEMO_SPARSE_VERDICTS,
} from '../adapter/demo-studies';
import type {DemoStudySpec} from '../adapter/demo-studies';
import {
  buildDemoStudy, seedDemoStudies, DemoDriftError, productionSeedDeps, loadDemoTable,
} from '../pipeline/demo-seed';
import type {BuiltDemoAnalysis, DemoSeedDeps} from '../pipeline/demo-seed';
import {appRunOf, appStateFrom} from '../pipeline/pipeline-save';
import {
  SPAN_RATIO_POOR, numAt, computeAccumulationRatios, computeBioavailability, buildSummaryCache,
  getNcaRepository, REVIEW_COL, REVIEW_EXCLUDED, NCA_STUDIES_CHANGED_EVENT, buildStudyGrid, studyCellExclusion,
  reportableResults, recipeOf, PASSTHROUGH_TAG, PASSTHROUGH_STRUCT,
} from '../adapter';
import type {AcceptanceThresholds, NcaRepository, ReviewUnit} from '../adapter';
import {loadSampleByName} from '../views/nca-launch';
import {openAnalysis} from '../views/nca-analysis-open';
import {RESULTS_GRID_TAG, currentRecipeFingerprint, liveWorkflowState} from '../views/nca-mapping-pane';
import type {NcaController} from '../controller/nca-controller';
import {MemoryNcaRepository} from './memory-nca-repository';
import {withRestrictedUser} from './restricted-user';

const T: AcceptanceThresholds = {
  lambdaAcceptanceRSquared: DEMO_THRESHOLDS.acceptanceRSquared, lambdaMinPoints: DEMO_THRESHOLDS.minPoints,
  extrapWarnPct: DEMO_THRESHOLDS.extrapWarnPct,
};

const DEPS = {
  loadTable: (f: string) => loadSampleByName(`${DEMO_STUDIES_DIR}/${f}`),
  now: () => new Date().toISOString(),
};

const study = (ref: number): DemoStudySpec => DEMO_STUDIES.find((s) => s.ref === ref)!;

/** Built once per (study, review) for the whole run. */
const memo = new Map<string, Promise<BuiltDemoAnalysis[]>>();
function built(ref: number, review = true): Promise<BuiltDemoAnalysis[]> {
  const k = `${ref}:${review}`;
  if (!memo.has(k)) memo.set(k, buildDemoStudy(DEPS, study(ref), {review}));
  return memo.get(k)!;
}

/** One unit of a run's current review snapshot. */
function units(b: BuiltDemoAnalysis): ReviewUnit[] {
  return b.run.controller.reviewSnapshot(T, {viable: b.run.viable, composites: b.run.composites}).units;
}
const unit = (b: BuiltDemoAnalysis, key: string): ReviewUnit => units(b).find((u) => u.key === key)!;
const codes = (u: ReviewUnit): string => [...u.flagCodes].sort().join('+');

/** A numeric results cell of one subject. */
function val(b: BuiltDemoAnalysis, subject: string, col: string): number {
  const res = b.run.outcome.results;
  const s = res.col('SUBJECT')!;
  for (let i = 0; i < res.rowCount; i++) if (s.getString(i) === subject) return numAt(res.col(col)!, i);
  return NaN;
}

const mean = (a: number[]): number => a.reduce((x, y) => x + y, 0) / a.length;
const sd = (a: number[]): number => {
  const m = mean(a);
  return Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / (a.length - 1));
};

// ─── F2: the planted verdicts on the shipped engine ─────────────────────────────

category('Demo studies: planted verdicts', () => {
  test('AC-2.1: every table builds; keys and annotations pass the app gates; thresholds pinned', async () => {
    expect(SPAN_RATIO_POOR, DEMO_THRESHOLDS.spanRatio, 'the span threshold the windows were designed for');
    for (const s of DEMO_STUDIES) {
      for (const b of await built(s.ref)) {
        const r = b.recipe;
        expect(r.ncaRules.lambdaZ.minRSquared, DEMO_THRESHOLDS.minRSquared, `${b.spec.label} minRSquared`);
        expect(r.ncaRules.lambdaZ.minPoints, DEMO_THRESHOLDS.minPoints, `${b.spec.label} minPoints`);
        expect(r.ncaRules.extrapWarnPct, DEMO_THRESHOLDS.extrapWarnPct, `${b.spec.label} extrap warn`);
        expect(r.acceptanceRSquared, DEMO_THRESHOLDS.acceptanceRSquared, `${b.spec.label} R_acc`);
        expect(r.groupByCols.join(','), s.groupBy.join(','), `${b.spec.label} Group by`);
        expect(r.passthroughCols.find((p) => p.name === DEMO_STRUCTURE_COL)?.isStructure, true,
          `${b.spec.label}: the structure annotation is a structure`);
        expect(b.input.df.name.endsWith('(simulated)'), true, `${b.spec.label}: the table says simulated`);
      }
    }
  }, {timeout: 300000});

  test('AC-1.3: the production load types SMILES as Molecule without Chem', async () => {
    const df = await loadDemoTable(productionSeedDeps(DEPS.loadTable), study(1));
    const col = df.col(DEMO_STRUCTURE_COL)!;
    expect(new Set(col.toList() as string[]).size, 8, 'one structure per compound');
    expect(col.semType, DG.SEMTYPE.MOLECULE, 'typed Molecule (grok test loads no Chem)');
  });

  test('AC-2.2 / AC-2.4: Study 1 — C1, C2, C3; nothing else flagged except N ≤ 3; F', async () => {
    const all = await built(1, false);
    const by = (cid: string): BuiltDemoAnalysis => all.find((b) => b.spec.label === cid)!;
    for (const k of ['R801', 'R802', 'R803', 'R804', 'R805', 'R806']) {
      const u = unit(by('DGK-0008'), k);
      expect(u.flagCodes.includes('SPAN_LT') && u.flagCodes.includes('EXTRAP_GT'), true, `C1 ${k}: ${codes(u)}`);
      const pct = val(by('DGK-0008'), k, 'PCT_EXTRAP');
      expect(pct > 20 && pct <= 35, true, `C1 ${k}: %AUCextrap ${pct} in (20, 35]`);
    }
    for (const k of ['R104', 'R105', 'R106']) expect(codes(unit(by('DGK-0001'), k)), 'NPTS_LE', `C2 ${k}`);
    expect(unit(by('DGK-0006'), 'R604').category, 'not-estimated', 'C3 double peak');
    expect(by('DGK-0006').run.outcome.results.col('STATUS')!.get(
      units(by('DGK-0006')).findIndex((u) => u.key === 'R604')), 'partial', 'C3 is engine partial (no window ≥ 0.85)');
    const planned = new Set(['R801', 'R802', 'R803', 'R804', 'R805', 'R806', 'R104', 'R105', 'R106', 'R604']);
    for (const b of all) {
      for (const u of units(b)) {
        if (planned.has(u.key)) continue;
        expect(u.category === 'accepted' || codes(u) === 'NPTS_LE', true,
          `unplanned ${b.spec.label} ${u.key}: ${u.category} ${codes(u)}`);
      }
    }
    // AC-2.5: F ≤ 100 % per compound; the §B order holds where adjacent compounds are distinguishable.
    const F = new Map<string, {f: number; se: number}>();
    for (const b of all) {
      const res = b.run.outcome.results;
      const fp = computeBioavailability(res, {doseCol: 'Dose_mg_kg', stratifyCols: ['Sex']}).table;
      const f = numAt(fp.col('F_PCT')!, 0);
      expect(f <= 100, true, `${b.spec.label}: F ${f} ≤ 100 %`);
      const route = res.col('ROUTE')!;
      const auc = (iv: boolean): number[] => Array.from({length: res.rowCount}, (_, i) => i)
        .filter((i) => route.getString(i).startsWith('IV') === iv)
        .map((i) => numAt(res.col('AUC_INF')!, i)).filter((x) => Number.isFinite(x));
      const po = auc(false); const iv = auc(true);
      // delta-method SE of the ratio of two independent group means (R2/N-5c)
      const rel = Math.sqrt((sd(po) / mean(po)) ** 2 / po.length + (sd(iv) / mean(iv)) ** 2 / iv.length);
      F.set(b.spec.label, {f, se: f * rel});
    }
    for (let i = 0; i + 1 < DEMO_F_ORDER.length; i++) {
      const a = F.get(DEMO_F_ORDER[i])!; const c = F.get(DEMO_F_ORDER[i + 1])!;
      if (Math.abs(a.f - c.f) > 2 * Math.sqrt(a.se ** 2 + c.se ** 2)) {
        expect(a.f > c.f, true,
          `F ${DEMO_F_ORDER[i]} (${a.f.toFixed(0)}) > ${DEMO_F_ORDER[i + 1]} (${c.f.toFixed(0)})`);
      }
    }
  }, {timeout: 300000});

  test('AC-2.2 / AC-2.5: Study 3 — C5 in [0.85, 0.90); F ≤ 100 % per dose', async () => {
    const [b] = await built(3, false);
    const u = unit(b, 'R316');
    const adj = val(b, 'R316', 'LAMBDAZ_ADJ_R2');
    expect(codes(u), 'ADJR2_LT', 'C5 flagged by adj. R² only');
    expect(adj >= 0.85 && adj < 0.90, true, `C5 adj. R² ${adj}`);
    for (const x of units(b)) {
      if (x.key !== 'R316')
        expect(x.category === 'accepted' || codes(x) === 'NPTS_LE', true, `unplanned ${x.key}: ${codes(x)}`);
    }
    const fp = computeBioavailability(b.run.outcome.results, {doseCol: 'Dose_mg_kg', stratifyCols: ['Sex']}).table;
    for (let i = 0; i < fp.rowCount; i++) expect(numAt(fp.col('F_PCT')!, i) <= 100, true, `F row ${i}`);
  }, {timeout: 300000});

  test('AC-2.2 / AC-2.6: Study 4 — C7, C8, controls, tlast = τ, Rac, t½', async () => {
    const [b] = await built(4, false);
    for (const k of ['T019_D1', 'T022_D1']) expect(codes(unit(b, k)), 'NPTS_LE', `C7 ${k}`);
    // C7: EXACTLY the two planted animals carry the delayed absorption at 100 mg/kg Day 1 — the latest
    // Tmax of the group (measured 4 h; the others 0.5–2 h). Other N ≤ 3 profiles there are the
    // ordinary 3-point fits the seed bulk-accepts (user decision).
    const d1 = ['T019', 'T020', 'T021', 'T022', 'T023', 'T024'].map((x) => `${x}_D1`);
    expect(d1.filter((k) => val(b, k, 'TMAX') >= 4).join(','), 'T019_D1,T022_D1',
      `C7: delayed absorption (Tmax ${d1.map((k) => `${k}=${val(b, k, 'TMAX')}`).join(' ')})`);
    const c8 = val(b, 'T016_D28', 'LAMBDAZ_ADJ_R2');
    expect(codes(unit(b, 'T016_D28')), 'ADJR2_LT', 'C8'); expect(c8 >= 0.85 && c8 < 0.90, true, `C8 ${c8}`);
    const c6 = val(b, 'T010_D28', 'LAMBDAZ_ADJ_R2');
    expect(codes(unit(b, 'T010_D28')), 'ADJR2_LT', 'C6 before'); expect(c6 >= 0.85, true, `C6 before ${c6}`);
    const res = b.run.outcome.results;
    const dose = res.col('Dose_mg_kg')!; const day = res.col('StudyDay')!; const sex = res.col('Sex')!;
    const subj = res.col('SUBJECT')!;
    const rows = Array.from({length: res.rowCount}, (_, i) => i);
    for (const i of rows) {
      const k = subj.getString(i);
      if (dose.getString(i) === '0') {
        expect(unit(b, k).category, 'no-exposure', `C4′ ${k}`);
        continue;
      }
      if (!['T019_D1', 'T022_D1', 'T016_D28', 'T010_D28'].includes(k)) {
        expect(unit(b, k).category === 'accepted' || codes(unit(b, k)) === 'NPTS_LE', true,
          `unplanned ${k}: ${codes(unit(b, k))}`);
      }
      expect(val(b, k, 'SPAN_RATIO') >= 2, true, `${k}: span ≥ 2`);
      expect(val(b, k, 'LAMBDAZ_TLAST'), 24, `${k}: measurable at 24 h (tlast = τ)`);
    }
    // Day-28 pre-dose quantifiable at every dosed level
    const obs = b.run.work;
    for (let i = 0; i < obs.rowCount; i++) {
      if (obs.col('StudyDay')!.getString(i) !== '28' || numAt(obs.col('Time_hr')!, i) !== 0 ||
        obs.col('Dose_mg_kg')!.getString(i) === '0') continue;
      expect(obs.col('BLQ')!.getString(i).toLowerCase(), 'false', `Day-28 trough row ${i} quantifiable`);
    }
    // Rac on AUCτ (Formula A) per sex × dose, as the TK panel computes it.
    const acc = computeAccumulationRatios(res, {tau: 24, doseCol: 'Dose_mg_kg', studyDayCol: 'StudyDay',
      stratifyCols: ['Sex'], tlastCol: 'LAMBDAZ_TLAST'});
    const t = acc.table;
    for (let i = 0; i < t.rowCount; i++) {
      const d = numAt(t.col('DOSE')!, i);
      const rac = numAt(t.col('RAC_AUC')!, i);
      const s = t.col('Sex')!.getString(i);
      if (d === 0) continue;
      if (d === 100) expect(rac >= 1.4 && rac <= 2.0, true, `Rac 100 ${s} = ${rac} in [1.4, 2.0]`);
      else expect(rac >= 0.85 && rac <= 1.15, true, `Rac ${d} ${s} = ${rac} in [0.85, 1.15] (linear identity)`);
    }
    const t28 = rows.filter((i) => dose.getString(i) === '100' && day.getString(i) === '28' && sex.getString(i) === 'M')
      .map((i) => numAt(res.col('HALF_LIFE')!, i));
    expect(mean(t28) >= 3.5 && mean(t28) <= 5.5, true, `Day-28 t½ 100 mg/kg (M) ${mean(t28)} in [3.5, 5.5]`);
  }, {timeout: 300000});

  test('AC-2.7 / AC-2.7a: Study 5 — infusion compute; C10 at 0.53 h; C9 auto-fit on 3 points', async () => {
    const [b] = await built(5, false);
    const res = b.run.outcome.results;
    for (let i = 0; i < res.rowCount; i++) {
      expect(res.col('ROUTE')!.getString(i), 'IV-infusion', 'computed as an infusion');
      expect(res.col('C0_EXTRAP')!.isNone(i), true, 'no C0 back-extrapolation');
    }
    expect(unit(b, 'R505').category, 'accepted', 'C10');
    const t = b.run.work.col('Time_hr')!; const s = b.run.work.col('SubjectID')!;
    expect(Array.from({length: t.length}, (_, i) => i).some((i) => s.getString(i) === 'R505' && numAt(t, i) === 0.53),
      true, 'C10 end-of-infusion sample at 0.53 h');
    // MRT with −T_inf/2 (0.25 h for the 0.5 h infusion)
    for (let i = 0; i < res.rowCount; i++) {
      const mrt = numAt(res.col('MRT')!, i);
      const ratio = numAt(res.col('AUMC_INF')!, i) / numAt(res.col('AUC_INF')!, i);
      if (Number.isFinite(mrt)) expect(Math.abs(mrt - (ratio - 0.25)) < 1e-6 * Math.max(1, mrt), true, `MRT row ${i}`);
    }
    expect(codes(unit(b, 'R510')), 'NPTS_LE', 'C9 auto-fit kept 3 points');
    expect(val(b, 'R510', 'LAMBDAZ_TFIRST'), 8, 'C9 auto window from 8 h');
  }, {timeout: 300000});

  test('AC-2.3 / AC-2.4 / C11: after the planted review — overrides, decisions, To review counts', async () => {
    const [b4] = await built(4);
    expect(unit(b4, 'T010_D28').category, 'accepted', 'C6 after force-out of 8 h');
    expect(val(b4, 'T010_D28', 'LAMBDAZ_NPOINTS') >= 4, true, 'C6 ≥ 4 fit points');
    const [b4raw] = await built(4, false);
    expect(val(b4, 'T010_D28', 'AUC_LAST'), val(b4raw, 'T010_D28', 'AUC_LAST'), 'C6 AUClast unchanged');
    expect(b4.pending, 3, 'Study 4: C7 × 2 + C8 to review');
    const [b5] = await built(5);
    expect(unit(b5, 'R510').category, 'accepted', 'C9 after force-in of 4 h, 6 h');
    expect(val(b5, 'R510', 'LAMBDAZ_NPOINTS'), 5, 'C9 5 points (4–24 h)');
    expect(val(b5, 'R510', 'LAMBDAZ_ADJ_R2') >= 0.90, true, 'C9 adj. R² ≥ 0.90');
    const [b3] = await built(3);
    expect(b3.run.outcome.results.col(REVIEW_COL)!.get(
      units(b3).findIndex((u) => u.key === 'R316')), REVIEW_EXCLUDED, 'C5 Excluded');
    for (const s of [1, 2, 3, 5])
      for (const b of await built(s)) expect(b.pending, 0, `${b.spec.label} fully reviewed`);

    for (const b of await built(2)) {
      const comp = units(b).filter((u) => u.kind === 'composite');
      expect(comp.length, 1, `${b.spec.label}: one composite`);
      expect(comp[0].decision?.action, DEMO_SPARSE_VERDICTS[b.spec.label].decision?.action, `C11 ${b.spec.label}`);
    }
  }, {timeout: 300000});
});

/** B0: a 0-dose group in a TK analysis, measured (the controls stay on this result; known gaps: two
 *  undisclosed blank rows, and SD = 0 at N = 0). The panels' own frames were read through the real app. */
category('Demo studies: zero-dose TK probe', () => {
  test('B0: the engine rows of the vehicle controls carry no non-finite value', async () => {
    const [b] = await built(4, false);
    const res = b.run.outcome.results;
    const dose = res.col('Dose_mg_kg')!;
    for (let i = 0; i < res.rowCount; i++) {
      if (dose.getString(i) !== '0') continue;
      for (const c of ['CL', 'VZ', 'CL_ABS', 'VZ_ABS', 'AUC_LAST', 'CMAX', 'HALF_LIFE']) {
        const col = res.col(c)!;
        expect(col.isNone(i) || Number.isFinite(col.get(i)), true, `${res.col('SUBJECT')!.getString(i)} ${c}`);
      }
      expect(res.col('DOSE_NORM')!.getString(i) !== 'NaN', true, 'DOSE_NORM is a basis word, never NaN');
    }
  }, {timeout: 300000});

  test('B0: the Accumulation table never prints a non-finite number for the 0-dose cell', async () => {
    const [b] = await built(4, false);
    const res = b.run.outcome.results;
    const acc = computeAccumulationRatios(res, {tau: 24, doseCol: 'Dose_mg_kg', studyDayCol: 'StudyDay',
      stratifyCols: ['Sex'], tlastCol: 'LAMBDAZ_TLAST'}).table;
    for (let i = 0; i < acc.rowCount; i++) {
      for (const c of acc.columns.toList()) {
        if (c.type !== DG.COLUMN_TYPE.FLOAT || c.isNone(i)) continue;
        const v = c.get(i) as number;
        expect(Number.isNaN(v) || Number.isFinite(v), true, `Rac ${c.name} row ${i}: ${v}`);
      }
    }
  }, {timeout: 300000});

  test('B0: the vehicle-control tiers of the study cache are null, never non-finite', async () => {
    const [b] = await built(4);
    const tiers = b.input.summaryCache!.compoundTier
      .filter((t) => t.groupValues.some((g) => g.key === 'Dose_mg_kg' && g.value === '0'));
    expect(tiers.length, 4, 'one tier per sex × day');
    for (const t of tiers) {
      for (const p of [t.cmax, t.aucInf, t.cl, t.vz, t.thalf]) {
        // n = 0 → the builder's "not estimated" mean (NaN in memory, null once stored as JSON), which
        // the study grid reads as NE — never a number and never ±Infinity.
        expect(p.mean === null || Number.isNaN(p.mean), true, `no mean (${p.mean})`);
        expect(p.n, 0, 'n = 0');
      }
    }
  }, {timeout: 300000});
});

// ─── F3: the headless split reproduces today's pipeline output (AC-3.1) ─────────

category('Demo studies: app run split', () => {
  test('AC-3.1: appRunOf + appStateFrom reproduce the pre-split cache, stage and confirmed recipe', async () => {
    // The pre-split `appStateOf` built the cache from the raw results before any snapshot; the split
    // writes the Review column first (CF-7). With no decision recorded the only difference is the
    // two Review columns in the per-tier key inventory — the same columns the app's own save carries.
    const strip = (c: unknown): string => JSON.stringify(c, (k, v) => (k === 'keyInventory' ?
      (v as Array<{key: string}>).filter((e) => e.key !== 'Review' && e.key !== '_lz_review') : v));
    for (const ref of [1, 2, 4, 5]) {
      for (const b of await built(ref, false)) {
        const r1 = appRunOf(b.input.df, b.mapping, b.rules, b.recipe);
        const r2 = appRunOf(b.input.df, b.mapping, b.rules, b.recipe);
        expect(r1.ok && r2.ok, true, `${b.spec.label}: run`);
        if (!r1.ok || !r2.ok) continue;
        const legacy = buildSummaryCache(r1.run.outcome.results, r1.run.outcome.perProfile, b.recipe, b.spec.label,
          r1.run.observations);
        const state = appStateFrom(r2.run, b.rules, b.recipe, b.spec.label, 'now');
        expect(strip(state.summaryCache), strip(legacy), `${b.spec.label}: identical cache`);
        expect(state.confirmed.fingerprint, b.input.workflow.confirmed!.fingerprint, `${b.spec.label}: fingerprint`);
        expect(state.stageCache.pending, b.pending, `${b.spec.label}: stage`);
      }
    }
  }, {timeout: 300000});
});

// ─── F6: the seed protocol (memory repository) ──────────────────────────────────

/** A memory repository with every demo study reachable, plus write-failure injection. */
function memDeps(repo: NcaRepository, over: Partial<DemoSeedDeps> = {}): DemoSeedDeps {
  return {
    repo, ...DEPS,
    ownerOf: async () => 'Jane Owner', shareGroupId: 'all-users', ...over,
  };
}
const col = (df: DG.DataFrame, name: string): string[] => df.col(name)!.toList() as string[];

/** A proxy of `repo` that answers `name` with `fn` (the other members pass through). */
function withMember<T extends object>(repo: T, name: string, fn: unknown): T {
  return new Proxy(repo, {get: (t, p, r) => {
    if (p === name) return fn;
    const v = Reflect.get(t, p, r) as unknown;
    return typeof v === 'function' ? (v as (...x: unknown[]) => unknown).bind(t) : v;
  }});
}

/** The analysis ids of one installed demo study. */
async function analysisIds(repo: NcaRepository, label: string): Promise<string[]> {
  return (await repo.listStudies()).find((x) => x.label === label)!.analyses.map((a) => a.id);
}

category('Demo studies: seed protocol', () => {
  test('any user installs the five studies, shared View to All users; a second call writes nothing', async () => {
    const repo = new MemoryNcaRepository(true);
    let changed = 0;
    const sub = grok.events.onCustomEvent(NCA_STUDIES_CHANGED_EVENT).subscribe(() => changed++);
    const shared: string[] = [];
    const spy = withMember(repo, 'share', async (id: string, g: string) => {
      shared.push(`${id}:${g}`);
      return repo.share(id, g);
    });
    const first = await seedDemoStudies(memDeps(spy));
    expect(col(first, 'Status').every((x) => x === 'installed'), true, col(first, 'Status').join('; '));
    const list = await repo.listStudies();
    expect(list.length, 5);
    expect(shared.length, 5, 'each study shared once');
    expect(shared.every((x) => x.endsWith(':all-users')), true, 'with All users');
    for (const s of DEMO_STUDIES) {
      const h = list.find((x) => x.label === s.label)!;
      expect(h.analyses.map((a) => a.label).sort().join(','), s.analyses.map((a) => a.label).sort().join(','));
      for (const e of (await repo.loadStudy(h.id)).analyses) {
        const full = (await repo.loadAnalysis(e.id)).entry;
        expect(full.confirmed !== null && full.stageCache !== null, true, `${e.label} confirmed with a stage`);
      }
    }
    const e4 = (await repo.loadAnalysis((await analysisIds(repo, study(4).label))[0])).entry;
    expect((e4.lzOverrides ?? []).length, 1, 'C6 override stored');
    expect(e4.stageCache!.pending, 3, 'Study 4 To review');
    const second = await seedDemoStudies(memDeps(repo));
    expect(col(second, 'Status').every((x) => x === 'up to date'), true, col(second, 'Status').join('; '));
    expect((await repo.listStudies()).length, 5, 'no duplicates');
    sub.unsubscribe();
    expect(changed > 0, true, 'the Hub and Browse tree are told (NCA_STUDIES_CHANGED_EVENT), no refresh needed');
  }, {timeout: 600000});

  test('another sci-comp version: the owner re-installs (planned first, new ids)', async () => {
    const repo = new MemoryNcaRepository();
    await seedDemoStudies(memDeps(repo, {specs: [study(5)]}));
    const before = (await repo.listStudies())[0].id;
    repo.bumpContent((await analysisIds(repo, study(5).label))[0], {engineVersion: '0.0.0-older'});
    const out = await seedDemoStudies(memDeps(repo, {specs: [study(5)]}));
    const st = col(out, 'Status');
    expect(st.length, 3, st.join('; '));
    expect(st[0].startsWith('planned: re-install') && st[0].includes('sci-comp'), true, st[0]);
    expect(st.slice(1).join(','), 'deleted,installed');
    const after = await repo.listStudies();
    expect(after.length, 1); expect(after[0].id !== before, true, 'a new id');
    expect(col(await seedDemoStudies(memDeps(repo, {specs: [study(5)]})), 'Status')[0], 'up to date');
  }, {timeout: 300000});

  test('newer bundled data: the owner re-installs', async () => {
    const repo = new MemoryNcaRepository();
    await seedDemoStudies(memDeps(repo, {specs: [study(3)]}));
    const id = (await analysisIds(repo, study(3).label))[0];
    const e = (await repo.loadAnalysis(id)).entry;
    repo.bumpContent(id, {sourceSignature: {...e.sourceSignature, hash: 'an-older-snapshot'}});
    const st = col(await seedDemoStudies(memDeps(repo, {specs: [study(3)]})), 'Status');
    expect(st[0].includes('older demo data'), true, st[0]);
    expect(st[st.length - 1], 'installed');
  }, {timeout: 300000});

  test('out of date, and the caller cannot delete it: told who can refresh it, nothing written', async () => {
    const repo = new MemoryNcaRepository();
    await seedDemoStudies(memDeps(repo, {specs: [study(5)]}));
    const before = (await repo.listStudies())[0].id;
    repo.bumpContent((await analysisIds(repo, study(5).label))[0], {engineVersion: '0.0.0-older'});
    const locked = withMember(repo, 'listStudies', async () =>
      (await repo.listStudies()).map((x) => ({...x, access: {...x.access, canDelete: false, canEdit: false}})));
    let writes = 0;
    const counted = withMember(locked, 'createStudy', async () => {writes++; throw new Error('no write expected');});
    const st = col(await seedDemoStudies(memDeps(counted, {specs: [study(5)]})), 'Status');
    expect(st.length, 1, st.join('; '));
    expect(st[0].startsWith('out of date') && st[0].includes('Jane Owner') && st[0].includes('administrator'),
      true, st[0]);
    expect(writes, 0, 'nothing written');
    expect((await repo.listStudies())[0].id, before, 'the study is untouched');
  }, {timeout: 300000});

  test('a same-label study with other analyses is refused and never deleted', async () => {
    const repo = new MemoryNcaRepository();
    const foreign = await repo.createStudy(study(3).label);
    const out = await seedDemoStudies(memDeps(repo, {specs: [study(3)]}));
    expect(col(out, 'Status')[0].startsWith('refused'), true, col(out, 'Status')[0]);
    const after = await repo.listStudies();
    expect(after.length, 1); expect(after[0].id, foreign.id, 'the user study is untouched');
  }, {timeout: 300000});

  test('AC-6.6: a failure mid-study leaves no partial study', async () => {
    const repo = new MemoryNcaRepository();
    let n = 0;
    const failing = withMember(repo, 'appendAnalysis', async (...a: Parameters<NcaRepository['appendAnalysis']>) => {
      if (++n === 2) throw new Error('injected');
      return repo.appendAnalysis(...a);
    });
    const out = await seedDemoStudies(memDeps(failing, {specs: [study(1)]}));
    expect(col(out, 'Status')[0].startsWith('not installed'), true, col(out, 'Status')[0]);
    expect((await repo.listStudies()).length, 0, 'the half-built study was deleted');
  }, {timeout: 300000});

  test('AC-3.2: drift — a planted decision the app would refuse throws', async () => {
    const s = {...study(3), analyses: [{...study(3).analyses[0],
      decisions: [{ref: 'X', units: ['R301'], action: 'exclude' as const, codes: [], other: 'x'}]}]};
    let err: unknown = null;
    try {await buildDemoStudy(DEPS, s);} catch (e) {err = e;}
    expect(err instanceof DemoDriftError, true, String(err));
  }, {timeout: 300000});
});

// ─── F6 live: EMS, reopen invariance (AC-3.3), View-only (AC-6.7) ────────────────

/** Renamed copies — a live run must never touch a real demo install on the server. */
const LIVE_SPECS: DemoStudySpec[] = DEMO_STUDIES.map((s) => ({...s, label: `zz-test ${s.label}`}));

/** Open one seeded analysis and assert the B1 contract (AC-3.3): the replay lands on the stored
 *  fingerprint, stage and pending count, drops no override, auto-resets no decision; the ✎ carries the
 *  seeding account and time; the source `Comment` sits at the overridden rows (AC-2.8); C5 is Excluded;
 *  and the stored cache equals the cache of the opened compute. A View-only open also has SAVE disabled. */
async function checkReopen(repo: NcaRepository, id: string, analyst: string, viewOnly: boolean): Promise<void> {
  const reader = viewOnly ? new (await import('../adapter')).EmsNcaRepository() : repo;
  const stored = (await reader.loadAnalysis(id)).entry;
  const label = `${stored.label} (${viewOnly ? 'View-only' : 'owner'})`;
  const view = await openAnalysis(id) as DG.TableView;
  try {
    const bag = view.temp as Record<string, unknown>;
    const t0 = Date.now();
    while (bag[RESULTS_GRID_TAG] == null && Date.now() - t0 < 60000) await DG.delay(250);
    await DG.delay(1500);
    expect(currentRecipeFingerprint(view), stored.confirmed!.fingerprint, `${label}: fingerprint`);
    const live = liveWorkflowState(view);
    expect(live.stage.stage, stored.stageCache!.stage, `${label}: stage`);
    expect(live.stage.pending, stored.stageCache!.pending, `${label}: pending`);
    const controller = bag['ncaController'] as NcaController;
    expect(controller.consumeReconciled(), 0, `${label}: 0 overrides pruned`);
    expect(controller.getLzOverrides().length, (stored.lzOverrides ?? []).length, `${label}: every override kept`);
    const current = (l: ReadonlyArray<{current: unknown}>): number => l.filter((x) => x.current !== null).length;
    expect(current(controller.getReviewLedgers()), current(stored.reviewDecisions ?? []),
      `${label}: 0 decisions auto-reset`);
    const src = view.dataFrame;
    const commentAt = (subject: string, t: number): string => {
      for (let i = 0; i < src.rowCount; i++) {
        if (src.col('SubjectID')!.getString(i) === subject && numAt(src.col('Time_hr')!, i) === t)
          return src.col('Comment')?.getString(i) ?? '';
      }
      return '';
    };
    const planted: Array<[string, number[], string]> =
      [['T010_D28', [8], 'Bioanalytical flag'], ['R510', [4, 6], 'Terminal phase']];
    for (const [k, times, text] of planted) {
      if (controller.profileIndex(k) === undefined) continue;
      expect(controller.overrideUserFor(k), analyst, `${label}: ✎ ${k} by the seeding account`);
      const at = controller.overrideTimeFor(k) ?? '';
      expect(at !== '' && at === (stored.lzOverrides ?? []).find((o) => o.profileKey === k)?.at, true,
        `${label}: ✎ ${k} time ${at}`);
      for (const t of times) expect(commentAt(k, t).includes(text), true, `${label}: Comment at ${k} ${t} h (AC-2.8)`);
    }
    const res = controller.results!;
    const subj = res.col('SUBJECT')!;
    for (let i = 0; i < res.rowCount; i++) {
      if (subj.getString(i) === 'R316')
        expect(res.col(REVIEW_COL)!.getString(i), REVIEW_EXCLUDED, `${label}: C5 Excluded`);
    }
    // Not vacuous: the cache of THIS open's compute, built as the app saves it, equals the stored one.
    const outcome = controller.lastOutcome as {perProfile: Parameters<typeof buildSummaryCache>[1]};
    const fresh = buildSummaryCache(reportableResults(res, controller.excludedKeys()).df, outcome.perProfile,
      {...recipeOf(stored), passthroughCols: stored.passthroughCols}, stored.label);
    const fields = ['n', 'cmax', 'aucInf', 'cl', 'vz', 'thalf', 'mrt', 'vss', 'fPct'] as const;
    const diffs: string[] = [];
    const tiers = stored.summaryCache?.compoundTier ?? [];
    expect(fresh.compoundTier.length, tiers.length, `${label}: tier count`);
    tiers.forEach((t, i) => {
      for (const f of fields) {
        const x = JSON.stringify(t[f]); const y = JSON.stringify(fresh.compoundTier[i]?.[f]);
        if (x !== y) diffs.push(`${t.compoundId.replace(//g, '/')}.${f}: stored ${x} vs open ${y}`);
      }
    });
    expect(diffs.slice(0, 3).join(' | '), '', `${label}: cached numbers equal the reopened compute`);
    if (viewOnly) {
      const save = view.root.ownerDocument.querySelector('[name="button-Save"]');
      expect(save?.classList.contains('d4-disabled'), true, `${label}: SAVE disabled`);
    }
  } finally {
    view.close();
  }
}

category('Demo studies: live seed', () => {
  test('install on EMS; a second call is up to date; reopen replays unchanged; View-only reads it', async () => {
    const repo = getNcaRepository();
    const me = await grok.dapi.users.current();
    const deps: DemoSeedDeps = {...memDeps(repo), shareGroupId: DG.Group.defaultGroupsIds['All users'],
      ownerOf: async () => me.friendlyName, specs: LIVE_SPECS};
    const ids: string[] = [];
    try {
      const out = await seedDemoStudies(deps);
      ids.push(...col(out, 'StudyId').filter((x) => x !== ''));
      expect(col(out, 'Status').every((x) => x === 'installed'), true, col(out, 'Status').join('; '));
      // The freshness check on EMS: the stored snapshots hash as the bundled tables, nothing re-installs.
      const again = await seedDemoStudies(deps);
      expect(col(again, 'Status').every((x) => x === 'up to date'), true, col(again, 'Status').join('; '));

      // AC-3.3 / AC-2.8 — every analysis reopens on its stored fingerprint, stage, overrides, decisions
      // and numbers, as the seeding account; then the same as a View-only user (AC-6.7).
      const analyses: string[] = [];
      for (const sid of ids) for (const x of (await repo.loadStudy(sid)).analyses) analyses.push(x.id);
      for (const id of analyses) await checkReopen(repo, id, me.friendlyName, false);
      const s3 = (await repo.listStudies()).find((x) => x.label === `zz-test ${study(3).label}`)!;
      const grid3 = buildStudyGrid(await repo.loadStudy(s3.id));
      let excl: string | null = null;
      for (const c of grid3.columns.names()) {
        for (let r = 0; r < grid3.rowCount; r++) {
          const x = studyCellExclusion(grid3, c, r);
          if (x !== null) excl = `n = ${x.n} of ${x.total}`;
        }
      }
      expect(excl, 'n = 2 of 3', 'C5: the study grid reads n = 2 of 3');

      const note = await withRestrictedUser('ncaview', async (u) => {
        await u.asUser(async () => {
          const r = new (await import('../adapter')).EmsNcaRepository();
          const seen = (await r.listStudies()).filter((x) => ids.includes(x.id));
          expect(seen.length, ids.length, 'every demo study is visible');
          for (const x of seen) {
            expect(x.access.canEdit, false, `${x.label}: View only`);
            const grid = buildStudyGrid(await r.loadStudy(x.id));
            // The overview gives a structure-tagged passthrough column the Molecule semType at render.
            expect(grid.columns.toList().some((c) => c.getTag(PASSTHROUGH_TAG) === PASSTHROUGH_STRUCT), true,
              `${x.label}: structures in the study grid`);
          }
        });
        for (const id of analyses) await u.asUser(() => checkReopen(repo, id, me.friendlyName, true));
        // Any user may call the function: on a current demo it writes nothing.
        await u.asUser(async () => {
          const r = new (await import('../adapter')).EmsNcaRepository();
          const st = col(await seedDemoStudies({...deps, repo: r}), 'Status');
          expect(st.every((x) => x === 'up to date'), true, `View-only call: ${st.join('; ')}`);
        });
        return 'ok';
      });
      return note ?? 'View-only part skipped: self-signup unavailable';
    } finally {
      for (const id of ids) await repo.deleteStudy(id).catch(() => {});
    }
  }, {timeout: 900000});
});
