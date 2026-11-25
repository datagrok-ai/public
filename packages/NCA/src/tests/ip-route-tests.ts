/**
 * nca-ip-route — intraperitoneal (i.p.) as a first-class route.
 *
 * F1 vocabulary + the derived code → label inverse (drift test over `nca.isKnownRouteCode`) + the
 * IP = PO numeric identity through the app's entry points; F2 `intraperitoneal infusion`; F3 the
 * "duration not used" note and its fingerprint exclusion; F6 review decisions kept across
 * `other` → `IP` (dense fingerprint alias, composite key alias, the effective ledger view).
 * F7 lives in `pipeline-save-tests.ts`, F8 in `adapter-event-driven-tests.ts`, the study-grid
 * IV/EV case in `study-grid-tests.ts`.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {
  normalizeRouteString, normalizeRouteColumn, IP_INFUSION_AS_SINGLE_DOSE, INFUSION_DURATION_NOT_USED,
  DERIVED_ROUTE_CODE_COL,
} from '../adapter/route-normalize';
import {ROUTE_STRING_TO_CODE, ROUTE_CODE_TO_STRING, routeStringToCode} from '../adapter/profile-inputs';
import type {ProfileScalars} from '../adapter/profile-inputs';
import type {ColumnContract} from '../adapter/column-contract';
import {computeBioavailability, isIntravascular, IP_F_MEANING, fPctMeaningFor} from '../adapter/bioavailability';
import {
  computeDecisionFingerprint, legacyOtherRouteFingerprint, NUMBER_FREE_NOTE_CODES,
} from '../adapter/fingerprint';
import type {DecisionFingerprintInput} from '../adapter/fingerprint';
import {
  legacyCompositeKey, compositeAccepted, compositeKeyAliases, sparseUnitKey,
} from '../adapter/review-units';
import {computeSparseForRecipe, computeSparseStrata, shortSparseNote, sparseCompositeUnits}
  from '../adapter/sparse-nca';
import type {SparseComputeMapping} from '../adapter/sparse-nca';
import {excludedCompositeRows} from '../adapter/reportable';
import {evaluateAcceptance} from '../adapter/acceptance';
import {decisionRowsOf, ledgersOf} from '../adapter/ems-nca-rows';
import {DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED} from '../adapter/defaults';
import {NcaController} from '../controller/nca-controller';
import {parseMapping} from '../pipeline/pipeline-mapping';
import {loadFixture, makeSyntheticObs, runPipelineOrThrow, singleSubjectContract} from './test-utils';

const TIMEOUT = 120000;
const T = {
  lambdaAcceptanceRSquared: DEFAULT_ACCEPTANCE_RSQUARED,
  lambdaMinPoints: DEFAULT_RULES.lambdaZ.minPoints,
  extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
};
const ACC = T;
/** An acceptance R² no real fit reaches: every dense profile is Flagged, so it needs review and
 *  can hold a decision (the F6 dense cases are about kept decisions, not about the verdict). */
const STRICT = {...T, lambdaAcceptanceRSquared: 0.999999999};

const AC_1_1_SPELLINGS = ['IP', 'i.p.', 'I.P.', 'ip', 'intraperitoneal', 'Intraperitoneally', 'intra-peritoneal',
  'INTRAPERITONEAL', 'I-PERITON', 'intraperitonial', 'IP bolus', 'i.p. injection', 'intraperitoneal bolus'];

function strCol(name: string, values: ReadonlyArray<string>): DG.Column {
  return DG.Column.fromStrings(name, [...values]);
}

/** Two-subject extravascular profiles (an absorption phase, then a clean decline). */
function evFrame(): DG.DataFrame {
  const t = [0, 0.5, 1, 2, 4, 8, 12, 24];
  const c1 = [0, 4.1, 6.0, 5.2, 3.1, 1.2, 0.48, 0.03];
  const c2 = [0, 3.0, 5.1, 4.8, 3.0, 1.3, 0.55, 0.04];
  return makeSyntheticObs({
    subjects: [...t.map(() => 'S1'), ...t.map(() => 'S2')],
    times: Float64Array.from([...t, ...t]), concs: Float64Array.from([...c1, ...c2]),
  });
}

function scalars(route: nca.RouteCode, infusionDuration: number | null = null): ProfileScalars {
  return {dose: 10, route, doseUnits: 'mg', concentrationUnits: 'mg/L', timeUnits: 'h', bodyWeight: null,
    bodyWeightUnit: 'kg', infusionDuration, lloq: 0};
}

/** Every numeric (float / int) column of two results frames, cell by cell. */
function expectSameNumbers(a: DG.DataFrame, b: DG.DataFrame, msg: string): void {
  let compared = 0;
  for (const col of a.columns) {
    if (col.type !== DG.COLUMN_TYPE.FLOAT && col.type !== DG.COLUMN_TYPE.INT) continue;
    const other = b.col(col.name);
    expect(other !== null, true, `${msg}: ${col.name} present`);
    for (let i = 0; i < a.rowCount; i++) {
      const x = col.isNone(i) ? 'none' : String(col.get(i));
      const y = other!.isNone(i) ? 'none' : String(other!.get(i));
      expect(x, y, `${msg}: ${col.name}[${i}]`);
      compared++;
    }
  }
  expect(compared > 20, true, `${msg}: the comparison is not vacuous (${compared} cells)`);
}

category('IP route: vocabulary (F1)', () => {
  test('AC-1.1 every listed i.p. spelling resolves to code 6; an all-IP column passes the pre-scan', async () => {
    for (const s of AC_1_1_SPELLINGS) {
      const r = normalizeRouteString(s, false);
      expect(r.status === 'mapped' && r.code === nca.ROUTE_IP, true, `'${s}' → IP`);
    }
    const res = normalizeRouteColumn(strCol('Route', AC_1_1_SPELLINGS),
      strCol('Subject', AC_1_1_SPELLINGS.map((_, i) => `A${i}`)));
    expect(res.errors.length, 0, 'no UNMAPPABLE_ROUTE_COLUMN on an all-i.p. column');
    expect([...res.codes].every((c) => c === nca.ROUTE_IP), true, 'every row → 6');
    expect([...res.warningsBySubject.values()].flat().some((w) => w.code === 'UNMAPPABLE_ROUTE'), false,
      'no UNMAPPABLE_ROUTE warning');
  });

  test('AC-1.2 neighbours stay unrecognised; `ip drip` still blocks', async () => {
    for (const s of ['intraportal', 'intrapleural', 'I-PLEURAL'])
      expect(normalizeRouteString(s, false).status, 'other', `'${s}' is not i.p.`);
    const res = normalizeRouteColumn(strCol('Route', ['intraportal', 'PO']), strCol('Subject', ['A', 'B']));
    expect((res.warningsBySubject.get('A') ?? []).some((w) => w.code === 'UNMAPPABLE_ROUTE'), true,
      'intraportal → UNMAPPABLE_ROUTE warning');
    expect(normalizeRouteString('ip drip', false).status, 'block', 'ip drip blocked by the IV-adjacency guard');
  });

  test('AC-1.5 drift: every known route code has a label that round-trips', async () => {
    let known = 0;
    for (let c = -1; c <= 32; c++) {
      if (!nca.isKnownRouteCode(c)) continue;
      known++;
      const label = ROUTE_CODE_TO_STRING[c];
      expect(label !== undefined, true, `code ${c} has a label`);
      expect(ROUTE_STRING_TO_CODE[label], c, `'${label}' round-trips to ${c}`);
    }
    expect(known >= 7, true, `codes 0..6 are known (${known})`);
    expect(ROUTE_CODE_TO_STRING[nca.ROUTE_IP], 'IP');
  });

  test('AC-1.7 the closed set offers IP; pipeline defaultRoute IP validates, `ip` throws', async () => {
    expect(routeStringToCode('IP'), nca.ROUTE_IP);
    const keys = Object.keys(ROUTE_STRING_TO_CODE);
    expect(keys.indexOf('IP') >= 0 && keys.indexOf('IP') < keys.indexOf('other'), true, 'IP listed before other');
    const base = {subject: 'Subject', time: 'Time', concentration: 'Conc', dose: 'Dose',
      units: {concentration: 'ng/mL', time: 'h', dose: 'mg/kg', bodyWeight: 'kg'}};
    expect(parseMapping(JSON.stringify({...base, defaultRoute: 'IP'})).defaultRoute, 'IP');
    let threw = false;
    try {parseMapping(JSON.stringify({...base, defaultRoute: 'ip'}));} catch {threw = true;}
    expect(threw, true, 'exact keys: `ip` is not a defaultRoute');
  });
});

category('IP route: IP = PO through the entry points (F1)', () => {
  test('AC-1.3 / AC-1.4 dense: identical numbers; ROUTE IP, CL_BASIS CL/F, no Vss / C0', async () => {
    const df = evFrame();
    const ip = runPipelineOrThrow(df, singleSubjectContract(df), scalars(nca.ROUTE_IP), DEFAULT_RULES).results;
    const po = runPipelineOrThrow(df, singleSubjectContract(df), scalars(nca.ROUTE_PO), DEFAULT_RULES).results;
    expectSameNumbers(ip, po, 'IP vs PO');
    for (let i = 0; i < ip.rowCount; i++) {
      expect(ip.col('ROUTE')!.getString(i), 'IP');
      expect(ip.col('CL_BASIS')!.getString(i), 'CL/F');
      expect(ip.col('VSS')!.isNone(i), true, 'Vss null on IP');
      expect(ip.col('C0_EXTRAP')!.isNone(i), true, 'C0 null on IP');
    }
  });

  test('AC-1.4 sparse: identical composites and fingerprints (route column IP vs PO)', async () => {
    const run = async (spelling: string): Promise<ReturnType<typeof computeSparseForRecipe>> => {
      const df = await loadFixture('04_mouse_sparse_destructive.csv');
      df.col('Route')!.init(() => spelling);
      const norm = normalizeRouteColumn(df.col('Route')!, df.col('MouseID')!);
      df.columns.add(DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, norm.codes));
      const m: SparseComputeMapping = {
        subjectName: 'MouseID', timeName: 'NominalTime_hr', concName: 'Conc_ng_mL', doseColName: 'Dose_mg_kg',
        bwColName: 'BodyWeight_g', units: {time: 'h', concentration: 'ng/mL', dose: 'mg/kg', bodyWeight: 'g'},
        route: {codesCol: DERIVED_ROUTE_CODE_COL, scalar: 'PO', rawCol: 'Route'},
        lloq: {colName: 'LLOQ_ng_mL', scalar: 0}, blqColName: 'BLQ', nominalTimeName: 'NominalTime_hr',
        animalIdName: 'MouseID', declaredTopology: null, acceptance: ACC, bootstrapIterations: 20,
      };
      return computeSparseForRecipe(df, m);
    };
    const ip = await run('i.p.');
    const po = await run('PO');
    expect(ip.error, null);
    expect(ip.result!.groups.length > 0 && ip.result!.groups.length === po.result!.groups.length, true);
    ip.result!.groups.forEach((g, i) => {
      const p = po.result!.groups[i];
      expect(g.route, 'IP');
      for (const k of ['aucLast', 'aucInf', 'cmax', 'lambdaZ', 'cl'] as const)
        expect(String(g[k]), String(p[k]), `${g.groupKey} ${k}`);
    });
    expect(sparseCompositeUnits(ip.result).map((u) => u.fingerprint).join(),
      sparseCompositeUnits(po.result).map((u) => u.fingerprint).join(), 'composite fingerprints hash route ""');
  }, {timeout: TIMEOUT});

  test('AC-1.6 / AC-1.8 mixed IV + IP + PO: %F has one IP and one PO stratum, IP is extravascular', async () => {
    expect(isIntravascular('IP'), false, 'the shared IV/EV owner reads IP as extravascular');
    const results = DG.DataFrame.fromColumns([
      strCol('ROUTE', ['IV-bolus', 'IV-bolus', 'IP', 'IP', 'PO', 'PO']),
      strCol('SEX', ['M', 'M', 'M', 'M', 'M', 'M']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([1, 1, 10, 10, 10, 10])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([10, 10, 80, 80, 30, 30])),
    ]);
    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['SEX']});
    const byRoute = new Map<string, number>();
    for (let i = 0; i < table.rowCount; i++)
      byRoute.set(table.col('ROUTE')!.getString(i), Number(table.col('F_PCT')!.get(i)));
    expect([...byRoute.keys()].sort().join(), 'IP,PO', 'one IP and one PO stratum');
    expect(Math.abs(byRoute.get('IP')! - 80) < 1e-9, true, `IP %F = 80 (got ${byRoute.get('IP')})`);
    expect(Math.abs(byRoute.get('PO')! - 30) < 1e-9, true, 'PO %F = 30');
  });

  test('AC-4.1 / AC-4.2 the %F meaning is on IP rows only, from one constant', async () => {
    expect(fPctMeaningFor(ROUTE_CODE_TO_STRING[nca.ROUTE_IP]), IP_F_MEANING, 'present on the IP row');
    for (const r of ['PO', 'SC', 'IM', 'other', 'IV-bolus', 'IV-infusion', ''])
      expect(fPctMeaningFor(r), null, `absent on route '${r}'`);
    expect(fPctMeaningFor(null), null, 'absent when the route is missing');
  });
});

category('IP route: intraperitoneal infusion and unused durations (F2, F3)', () => {
  test('AC-2.1 / AC-2.2 i.p. infusion → IP + warning; the guard is unchanged after the new set', async () => {
    for (const s of ['intraperitoneal infusion', 'IP infusion', 'i.p. infusion']) {
      const r = normalizeRouteString(s, false);
      expect(r.status === 'mapped' && r.code === nca.ROUTE_IP && r.ipInfusion, true, `'${s}' → IP + flag`);
    }
    const ivInf = normalizeRouteString('IV infusion', true);
    expect(ivInf.status === 'mapped' && ivInf.code === nca.ROUTE_IV_INFUSION, true, 'IV infusion unchanged');
    const inf = normalizeRouteString('infusion', true);
    expect(inf.status === 'mapped' && inf.code === nca.ROUTE_IV_INFUSION, true, 'bare infusion unchanged');
    const drip = normalizeRouteString('intravenous drip', true);
    expect(drip.status === 'mapped' && drip.code === nca.ROUTE_IV_INFUSION, true, 'intravenous drip unchanged');
    expect(normalizeRouteString('xyz infusion', false).status, 'block', 'unknown infusion still blocks');
    const res = normalizeRouteColumn(strCol('Route', ['intraperitoneal infusion', 'IP']),
      strCol('Subject', ['A', 'B']));
    expect((res.warningsBySubject.get('A') ?? []).some((w) => w.code === IP_INFUSION_AS_SINGLE_DOSE), true,
      'A warned');
    expect((res.warningsBySubject.get('B') ?? []).length, 0, 'plain IP: no warning');
  });

  test('AC-2.3 a sparse IP composite counts its i.p.-infusion members (k of n animals)', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    const animals = df.col('MouseID')!;
    const first = animals.getString(0);
    df.col('Route')!.init((i) => animals.getString(i) === first ? 'intraperitoneal infusion' : 'IP');
    const norm = normalizeRouteColumn(df.col('Route')!, animals);
    df.columns.add(DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, norm.codes));
    const strata = computeSparseStrata(df, {
      concCol: 'Conc_ng_mL', nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID', blqCol: 'BLQ',
      doseCol: 'Dose_mg_kg', bodyWeightCol: 'BodyWeight_g', lloq: {colName: 'LLOQ_ng_mL', scalar: 0},
      route: {codesCol: DERIVED_ROUTE_CODE_COL, scalar: 'PO', rawCol: 'Route'},
      units: {time: 'h', concentration: 'ng/mL', dose: 'mg/kg', bodyWeight: 'g'}, acceptance: ACC,
      bootstrapIterations: 20,
    });
    if (!strata.ok) throw new Error(strata.reason);
    const withWarn = strata.cells.filter((c) => c.group.warnings.some((w) => w.code === IP_INFUSION_AS_SINGLE_DOSE));
    expect(withWarn.length, 1, 'exactly the composite holding that animal');
    const msg = withWarn[0].group.warnings.find((w) => w.code === IP_INFUSION_AS_SINGLE_DOSE)!.message;
    const n = new Set(withWarn[0].rows.map((i) => animals.getString(i))).size;
    expect(msg.startsWith(`1 of ${n} animals`), true, `k of n wording: ${msg}`);
    // The NOTES cell reads words, not the raw code.
    expect(shortSparseNote(IP_INFUSION_AS_SINGLE_DOSE), 'i.p. infusion as single dose', 'short NOTES label');
    // Out of the composite fingerprint (SPARSE_NUMBER_FREE_CODES, R2 N-3).
    const without = computeSparseStrata(df, {
      concCol: 'Conc_ng_mL', nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID', blqCol: 'BLQ',
      doseCol: 'Dose_mg_kg', bodyWeightCol: 'BodyWeight_g', lloq: {colName: 'LLOQ_ng_mL', scalar: 0},
      route: {codesCol: DERIVED_ROUTE_CODE_COL, scalar: 'PO', rawCol: null},
      units: {time: 'h', concentration: 'ng/mL', dose: 'mg/kg', bodyWeight: 'g'}, acceptance: ACC,
      bootstrapIterations: 20,
    });
    if (!without.ok) throw new Error(without.reason);
    const fps = (s: typeof strata): string => sparseCompositeUnits({groups: s.cells.map((c) => c.group)} as never)
      .map((u) => u.fingerprint).join();
    expect(fps(strata), fps(without), 'the warning does not move a composite fingerprint');
  }, {timeout: TIMEOUT});

  test('AC-3.1 per-row duration on an EV row → INFUSION_DURATION_NOT_USED; F2 + F3 give one line', async () => {
    const res = normalizeRouteColumn(
      strCol('Route', ['PO', 'IP', 'IV-infusion', 'PO', 'PO', 'intraperitoneal infusion', 'topical']),
      strCol('Subject', ['A', 'B', 'C', 'D', 'E', 'F', 'G']),
      {infusionDurationCol: DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Dur', [2, 1, 1, 0, null, 1, 3])});
    const has = (s: string, code: string): boolean =>
      (res.warningsBySubject.get(s) ?? []).some((w) => w.code === code);
    expect(has('A', INFUSION_DURATION_NOT_USED), true, 'PO + duration');
    expect(has('B', INFUSION_DURATION_NOT_USED), true, 'IP + duration');
    expect(has('C', INFUSION_DURATION_NOT_USED), false, 'IV-infusion uses its duration');
    expect(has('D', INFUSION_DURATION_NOT_USED), false, 'duration 0');
    expect(has('E', INFUSION_DURATION_NOT_USED), false, 'duration null');
    expect(has('F', INFUSION_DURATION_NOT_USED), false, 'F2 + F3: only the i.p.-infusion line');
    expect(has('F', IP_INFUSION_AS_SINGLE_DOSE), true, 'F2 + F3: the i.p.-infusion line');
    expect(has('G', INFUSION_DURATION_NOT_USED), true, 'other (topical) + duration');
    // One subject with both `IP` and `intraperitoneal infusion` rows, each with a duration: still
    // one line, the i.p.-infusion warning (the per-subject pass removes the duration note).
    const mixed = normalizeRouteColumn(strCol('Route', ['IP', 'intraperitoneal infusion']),
      strCol('Subject', ['H', 'H']), {infusionDurationCol: DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Dur', [1, 1])});
    const h = (mixed.warningsBySubject.get('H') ?? []).map((w) => w.code);
    expect(h.includes(IP_INFUSION_AS_SINGLE_DOSE) && !h.includes(INFUSION_DURATION_NOT_USED), true,
      `mixed subject: ${h.join(', ')}`);
    const note = (res.warningsBySubject.get('A') ?? []).find((w) => w.code === INFUSION_DURATION_NOT_USED)!;
    expect(note.severity, 'info');
    // A positive SCALAR duration never raises the note (hidden unless the scalar route is IV-infusion).
    const scalar = normalizeRouteColumn(strCol('Route', ['PO']), strCol('Subject', ['A']),
      {infusionScalarPositive: true});
    expect((scalar.warningsBySubject.get('A') ?? []).some((w) => w.code === INFUSION_DURATION_NOT_USED), false,
      'scalar duration: no note');
  });

  test('AC-3.2 the note leaves the decision fingerprint byte-identical', async () => {
    expect(NUMBER_FREE_NOTE_CODES.has(INFUSION_DURATION_NOT_USED), true);
    const base: DecisionFingerprintInput = {points: [['1', '2']], dose: '10', route: 'PO', pointsUsedTimes: ['1'],
      lambdaZ: 0.1, adjR2: 0.99, aucInf: 50, pctExtrap: 5, verdict: 'Accepted', flagCodes: [], warnings: []};
    expect(computeDecisionFingerprint({...base, warnings: [{code: INFUSION_DURATION_NOT_USED, severity: 'info'}]}),
      computeDecisionFingerprint(base));
  });

  test('AC-3.3 extravascular numbers are identical with and without a duration', async () => {
    const df = evFrame();
    const a = runPipelineOrThrow(df, singleSubjectContract(df), scalars(nca.ROUTE_IP, 2), DEFAULT_RULES).results;
    const b = runPipelineOrThrow(df, singleSubjectContract(df), scalars(nca.ROUTE_IP, null), DEFAULT_RULES).results;
    expectSameNumbers(a, b, 'duration vs none');
  });
});

category('IP route: decisions survive other → IP (F6)', () => {
  const contract = (df: DG.DataFrame): ColumnContract => singleSubjectContract(df);
  /** The analysis as it was saved (i.p. read as `other` + UNMAPPABLE_ROUTE) and as it reopens (IP). */
  function pair(df = evFrame()): {legacy: NcaController; now: NcaController} {
    const warn = {code: 'UNMAPPABLE_ROUTE', severity: 'warning' as const,
      message: 'Route value \'i.p.\' is not a recognised route; treated as \'other\'.'};
    const legacy = new NcaController(df, contract(df), scalars(nca.ROUTE_OTHER), DEFAULT_RULES,
      {routeWarnings: new Map([['S1', [warn]], ['S2', [warn]]])});
    legacy.computeAll();
    const now = new NcaController(df, contract(df), scalars(nca.ROUTE_IP), DEFAULT_RULES);
    now.computeAll();
    return {legacy, now};
  }

  test('AC-6.1 a dense Exclude λz recorded under `other` is kept under IP (0 resets, still excluded)', async () => {
    const {legacy, now} = pair();
    const legacySnap = legacy.reviewSnapshot(STRICT, {viable: true});
    const s1 = legacySnap.units.find((u) => u.key === 'S1')!;
    // Decide on what the analyst saw then.
    const legacyFp = legacySnap.fingerprints.get('S1') ?? (() => {
      throw new Error(`S1 must need review for this test (category ${s1.category})`);
    })();
    legacy.decide('S1', 'exclude', [{code: 'EXC_ERRATIC', text: 'x'}], legacyFp);
    now.setReviewLedgers(legacy.getReviewLedgers());
    const snap = now.reviewSnapshot(STRICT, {viable: true});
    expect(snap.fingerprints.get('S1') !== legacyFp, true, 'the current fingerprint differs (route moved)');
    expect(now.reconcileDecisions(snap.accepted), 0, 'reset count 0');
    expect(now.decisionFor('S1')?.action, 'exclude', 'decision kept');
    expect(now.excludedKeys().has('S1'), true, 'Exclude λz still applied');
    expect(now.reviewSnapshot(STRICT, {viable: true}).units.find((u) => u.key === 'S1')!.decision !== null, true,
      'Review state: decided');
  });

  test('AC-6.2 a changed number still resets; AC-6.3 a PO profile gets no alias', async () => {
    const changed = evFrame();
    const conc = changed.col('conc')!;
    conc.set(5, (conc.get(5) as number) * 1.3); // S1, t = 8: the fit moves
    const legacyDf = evFrame();
    const {legacy} = pair(legacyDf);
    const fp = legacy.reviewSnapshot(STRICT, {viable: true}).fingerprints.get('S1')!;
    legacy.decide('S1', 'exclude', [{code: 'EXC_ERRATIC', text: 'x'}], fp);
    const now = new NcaController(changed, contract(changed), scalars(nca.ROUTE_IP), DEFAULT_RULES);
    now.computeAll();
    now.setReviewLedgers(legacy.getReviewLedgers());
    expect(now.reconcileDecisions(now.reviewSnapshot(STRICT, {viable: true}).accepted), 1, 'fit changed → reset');

    const po: DecisionFingerprintInput = {points: [], dose: '1', route: 'PO', pointsUsedTimes: [], lambdaZ: 1,
      adjR2: 1, aucInf: 1, pctExtrap: 1, verdict: 'Flagged', flagCodes: [], warnings: []};
    expect(legacyOtherRouteFingerprint(po), null, 'no alias for PO');
    expect(legacyOtherRouteFingerprint({...po, route: 'IP'}) !== null, true, 'alias for IP');
  });

  const flagged = evaluateAcceptance({spanRatio: 1.2, adjR2: 0.99, nPoints: 5, pctExtrap: 5}, T);
  const IP_KEY = 'ROUTE=IP|DOSE=10';
  const OTHER_KEY = 'ROUTE=other|DOSE=10';
  const IV_KEY = 'ROUTE=IV-bolus|DOSE=10';
  const composite = (groupKey: string, fingerprint: string) =>
    ({groupKey, acceptance: flagged, measurableTimepoints: 6, fingerprint});
  const legacyLedger = (fp: string) => ({unitKey: sparseUnitKey(OTHER_KEY), history: [], current: {
    action: 'exclude' as const, reasons: [{code: 'EXC_ERRATIC', text: 'x'}], user: 'u', at: '2026-11-01T00:00:00Z',
    fingerprint: fp}});
  function ctrl(): NcaController {
    const df = evFrame();
    const c = new NcaController(df, contract(df), scalars(nca.ROUTE_PO), DEFAULT_RULES);
    c.computeAll();
    return c;
  }

  test('pure: legacy key, aliases only when no current composite owns the legacy key', async () => {
    expect(legacyCompositeKey(sparseUnitKey(IP_KEY)), sparseUnitKey(OTHER_KEY));
    expect(legacyCompositeKey(sparseUnitKey('DOSE=10')), null, 'single-route keys do not change');
    expect(legacyCompositeKey('S1'), null);
    expect(compositeKeyAliases([sparseUnitKey(IP_KEY)]).get(sparseUnitKey(OTHER_KEY)), sparseUnitKey(IP_KEY));
    expect(compositeKeyAliases([sparseUnitKey(IP_KEY), sparseUnitKey(OTHER_KEY)]).size, 0, 'pooled: no alias');
    const acc = compositeAccepted([composite(IP_KEY, 'fpX')]);
    expect(acc.get(sparseUnitKey(OTHER_KEY))?.join(), 'fpX', 'legacy key accepts the new key\'s fingerprint');
  });

  test('AC-6.4 / AC-6.7 a composite Exclude λz under ROUTE=other masks the IP composite (view path)', async () => {
    const c = ctrl();
    c.setReviewLedgers([legacyLedger('fpX')]);
    const units = [composite(IP_KEY, 'fpX'), composite(IV_KEY, 'fpIV')];
    // The view's sparse reconcile (nca-mapping-pane): composite keys + the ONE builder.
    c.setCompositeKeys(units.map((u) => sparseUnitKey(u.groupKey)));
    expect(c.reconcileDecisions(compositeAccepted(units), (k) => k.startsWith('sparse:')), 0, 'reset count 0');
    expect(c.decisionFor(sparseUnitKey(IP_KEY))?.action, 'exclude', 'the IP composite shows the decision');
    expect(excludedCompositeRows([IP_KEY, IV_KEY], c.excludedKeys()).join(), 'true,false', 'masked');
    // The Confirm dry run path: reviewSnapshot with the composites → its own `accepted`.
    const snap = c.reviewSnapshot(T, {viable: true, composites: units});
    expect(snap.units.find((u) => u.key === sparseUnitKey(IP_KEY))!.decision?.action, 'exclude', 'classified decided');
    expect(c.reconcileDecisions(snap.accepted), 0, 'dry run predicts 0 resets');
    expect(c.decisionFor(sparseUnitKey(IP_KEY))?.action, 'exclude', 'still kept');
  });

  test('AC-6.4 / AC-6.8 a remaining `other` composite owns the legacy key; its changed members reset it', async () => {
    const c = ctrl();
    c.setReviewLedgers([legacyLedger('fpPooled')]);
    const units = [composite(IP_KEY, 'fpIP'), composite(OTHER_KEY, 'fpTopicalOnly'), composite(IV_KEY, 'fpIV')];
    c.setCompositeKeys(units.map((u) => sparseUnitKey(u.groupKey)));
    expect(c.decisionFor(sparseUnitKey(IP_KEY)), null, 'the IP composite has no decision');
    expect(c.reconcileDecisions(compositeAccepted(units), (k) => k.startsWith('sparse:')), 1,
      'the pooled `other` decision resets (fewer members), counted');
    expect(c.decisionFor(sparseUnitKey(OTHER_KEY)), null);
  });

  test('AC-6.5 / AC-6.6 a new decision on the IP composite leaves one current record; rows round-trip', async () => {
    const c = ctrl();
    c.setReviewLedgers([legacyLedger('fpX')]);
    const units = [composite(IP_KEY, 'fpX'), composite(IV_KEY, 'fpIV')];
    c.setCompositeKeys(units.map((u) => sparseUnitKey(u.groupKey)));
    c.decide(sparseUnitKey(IP_KEY), 'accept', [{code: 'ACC_VISUAL_OK', text: 'ok'}], 'fpX');
    const stored = c.getReviewLedgers();
    const currents = stored.filter((l) => l.current !== null);
    expect(currents.map((l) => l.unitKey).join(), sparseUnitKey(IP_KEY), 'one current record, on the new key');
    const legacy = stored.find((l) => l.unitKey === sparseUnitKey(OTHER_KEY))!;
    expect(legacy.history[0].supersededBy, 'replaced', 'the legacy record superseded `replaced`');
    expect(c.reconcileDecisions(compositeAccepted(units), (k) => k.startsWith('sparse:')), 0, 'reset notice 0');
    // AC-6.6: the append-only rows round-trip (no immutable-core change).
    const back = ledgersOf(decisionRowsOf('a', 's', stored));
    expect(JSON.stringify(back.map((l) => [l.unitKey, l.current?.fingerprint ?? null, l.history.length]).sort()),
      JSON.stringify(stored.map((l) => [l.unitKey, l.current?.fingerprint ?? null, l.history.length]).sort()));
  });

  test('a Reset on the IP composite resets the legacy decision it shows', async () => {
    const c = ctrl();
    c.setReviewLedgers([legacyLedger('fpX')]);
    c.setCompositeKeys([sparseUnitKey(IP_KEY), sparseUnitKey(IV_KEY)]);
    c.resetDecision(sparseUnitKey(IP_KEY));
    expect(c.decisionFor(sparseUnitKey(IP_KEY)), null);
    expect(c.getReviewLedgers().find((l) => l.unitKey === sparseUnitKey(OTHER_KEY))!.history[0].supersededBy, 'reset');
  });
});
