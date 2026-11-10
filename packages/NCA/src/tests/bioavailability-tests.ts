/**
 * Adapter — bioavailability (F%) category (per-row-route Slice D, AD-6).
 *
 * @see files/rat_multiroute_bioavailability.csv (PKNCA reference)
 *
 * AC-D1 is pinned to the **stored reference columns** (rule 16 / R2-N3): the
 * fixture CSV carries per-stratum `dn_AUC`, `iv_dn_AUC` and `F_pct`. We build a
 * per-subject results frame from those stored `dn_AUC` / `iv_dn_AUC` numbers
 * (one subject per non-IV stratum, one IV reference subject per sex) and assert
 * `computeBioavailability` reproduces the stored `F_pct` — so the ratio math is
 * validated against the real reference output, not a re-derivation.
 *
 * Also covers AC-D2 (no-IV stratum → NaN + warning), AC-D3 (F% > 100 uncapped),
 * arithmetic-vs-geometric `denominatorStat`, and AC-D4 (`aggregateNcaResults`
 * unmodified — no regression to per-compound aggregation).
 *
 * **Footgun gate:** `expect(cond === expected, true, msg)` form throughout.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat, before} from '@datagrok-libraries/test/src/test';

import {computeBioavailability} from '../adapter/bioavailability';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {unmappedGridColumns} from '../views/results-view';

const FIXTURE = 'System:AppData/Nca/rat_multiroute_bioavailability.csv';

/** Reference stratum rows parsed from the bundled CSV. */
interface RefRow {
  route: string; sex: string; dose: number;
  dnAuc: number; ivDnAuc: number; fPct: number;
}

async function loadReferenceRows(): Promise<RefRow[]> {
  const csv = await grok.dapi.files.readAsText(FIXTURE);
  const df = DG.DataFrame.fromCsv(csv);
  const rows: RefRow[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    rows.push({
      route: String(df.col('Route')!.get(i)),
      sex: String(df.col('Sex')!.get(i)),
      dose: df.col('Dose_mg_kg')!.get(i) as number,
      dnAuc: df.col('dn_AUC')!.get(i) as number,
      ivDnAuc: df.col('iv_dn_AUC')!.get(i) as number,
      fPct: df.col('F_pct')!.get(i) as number,
    });
  }
  return rows;
}

/**
 * Build a per-subject NCA results frame from the reference stratum rows:
 *  - one subject per non-IV stratum with `AUC_INF = dn_AUC × Dose` (so its
 *    dose-normalized AUCinf equals the stored stratum `dn_AUC`),
 *  - one IV reference subject per sex with `AUC_INF = iv_dn_AUC × 1` at IV
 *    dose 1 (so the IV stratum mean equals the stored `iv_dn_AUC`).
 */
function resultsFromReference(rows: RefRow[]): DG.DataFrame {
  const route: string[] = [];
  const sex: string[] = [];
  const dose: number[] = [];
  const aucInf: number[] = [];

  const ivBySex = new Map<string, number>();
  for (const r of rows) ivBySex.set(r.sex, r.ivDnAuc);
  for (const [s, ivDn] of ivBySex) {
    route.push('IV-bolus'); sex.push(s); dose.push(1); aucInf.push(ivDn * 1);
  }
  for (const r of rows) {
    route.push(r.route); sex.push(r.sex); dose.push(r.dose); aucInf.push(r.dnAuc * r.dose);
  }

  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('ROUTE', route),
    DG.Column.fromStrings('SEX', sex),
    DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from(dose)),
    DG.Column.fromFloat64Array('AUC_INF', Float64Array.from(aucInf)),
  ]);
}

category('Adapter: bioavailability (F%)', () => {
  let refRows: RefRow[] = [];
  before(async () => {
    refRows = await loadReferenceRows();
  });

  // BUG-12 drift guard: every column of the REAL computed table must resolve in
  // BIOAVAIL_FIXED, or configureBioavailabilityGrid leaves it with a raw header
  // and DG's generic stats tooltip. Only the caller's stratify columns are
  // exempt (user-named, self-describing).
  test('BUG-12: every F% table column has a presentation entry', async () => {
    const results = resultsFromReference(refRows);
    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['SEX']});
    const unmapped = unmappedGridColumns(table, 'bioavailability', ['SEX']);
    expect(unmapped.length, 0,
      `F% columns with no BIOAVAIL_FIXED entry: [${unmapped.join(', ')}]`);
  });

  test('F-S5.2: N and N (IV ref.) count each side of every %F ratio', async () => {
    const results = resultsFromReference(refRows);
    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['SEX']});
    const n = table.col('N')!;
    const ivN = table.col('IV_N')!;
    expect(n !== null && ivN !== null, true, 'both N columns present');
    for (let i = 0; i < table.rowCount; i++) {
      expect(Number(n.get(i)) >= 1, true, `row ${i}: an EV cell has ≥ 1 contributor`);
      const f = Number(table.col('F_PCT')!.get(i));
      if (Number.isFinite(f)) expect(Number(ivN.get(i)) >= 1, true, `row ${i}: a finite %F has an IV reference`);
    }
  });

  test('AC-D1: reproduces stored stratum F_pct (sex-stratified, arithmetic)', async () => {
    expect(refRows.length > 0, true, 'reference rows loaded');
    const results = resultsFromReference(refRows);
    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['SEX']});

    // Index the output by route|sex|dose for lookup.
    const idx = new Map<string, number>();
    for (let i = 0; i < table.rowCount; i++) {
      const k = `${table.col('ROUTE')!.get(i)}|${table.col('SEX')!.get(i)}|${table.col('DOSE')!.get(i)}`;
      idx.set(k, i);
    }
    for (const r of refRows) {
      const k = `${r.route}|${r.sex}|${r.dose}`;
      const row = idx.get(k);
      expect(row !== undefined, true, `output row for ${k} missing`);
      const got = table.col('F_PCT')!.get(row!) as number;
      // Pure arithmetic from the same stored numbers → tight tolerance.
      expectFloat(got, r.fPct, Math.abs(r.fPct) * 1e-6, `F_pct for ${k}`);
      // Provenance: stored IV denominator reproduced (R2-N3).
      expectFloat(table.col('IV_DN_AUC')!.get(row!) as number, r.ivDnAuc,
        Math.abs(r.ivDnAuc) * 1e-6, `iv_dn_AUC for ${k}`);
    }
    // Spot-check the headline strata explicitly.
    const f5 = idx.get('PO|F|5')!;
    expectFloat(table.col('F_PCT')!.get(f5) as number, 43.0677728553402, 0.01, 'PO/F/5 = 43.07');
    const m5 = idx.get('PO|M|5')!;
    expectFloat(table.col('F_PCT')!.get(m5) as number, 31.9451449763666, 0.01, 'PO/M/5 = 31.95');
  });

  test('AC-D3: F% > 100 reported uncapped (SC/F/25 = 117.5)', async () => {
    const results = resultsFromReference(refRows);
    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['SEX']});
    let found = false;
    for (let i = 0; i < table.rowCount; i++) {
      if (table.col('ROUTE')!.get(i) === 'SC' && table.col('SEX')!.get(i) === 'F' &&
          table.col('DOSE')!.get(i) === 25) {
        found = true;
        expectFloat(table.col('F_PCT')!.get(i) as number, 117.505634161016, 0.01, 'SC/F/25 uncapped > 100');
      }
    }
    expect(found, true, 'SC/F/25 stratum present');
  });

  test('AC-D2: stratum with no IV arm → F% NaN + NO_IV_REFERENCE_FOR_STRATUM', async () => {
    // PO/M with no IV-M subject present.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ROUTE', ['IV-bolus', 'PO']),
      DG.Column.fromStrings('SEX', ['F', 'M']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([1, 5])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([2.0, 4.0])),
    ]);
    const {table, warnings} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['SEX']});
    expect(warnings.some((w) => w.code === 'NO_IV_REFERENCE_FOR_STRATUM'), true, 'no-IV warning raised');
    // The PO/M row F% must be NaN (no pooled fallback).
    let nanSeen = false;
    for (let i = 0; i < table.rowCount; i++) {
      if (table.col('SEX')!.get(i) === 'M')
        nanSeen = table.col('F_PCT')!.isNone(i);
    }
    expect(nanSeen, true, 'PO/M F% is NaN with no IV-M reference');
  });

  test('denominatorStat: geometric mean differs from arithmetic on skewed IV arm', async () => {
    // Two IV-F subjects dn = {1, 4}: arith mean = 2.5, geo mean = 2.0.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ROUTE', ['IV-bolus', 'IV-bolus', 'PO']),
      DG.Column.fromStrings('SEX', ['F', 'F', 'F']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([1, 1, 1])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([1, 4, 2])),
    ]);
    const base = {doseCol: 'Dose_mg_kg', stratifyCols: ['SEX']} as const;
    const arith = computeBioavailability(results, {...base, denominatorStat: 'mean'});
    const geo = computeBioavailability(results, {...base, denominatorStat: 'gmean'});
    // PO dn = 2 → F% arith = 2/2.5*100 = 80; geo = 2/2*100 = 100.
    expectFloat(arith.table.col('F_PCT')!.get(0) as number, 80, 0.001, 'arithmetic denominator');
    expectFloat(geo.table.col('F_PCT')!.get(0) as number, 100, 0.001, 'geometric denominator');
    expect(arith.denominatorStat, 'mean', 'arith stat surfaced');
    expect(geo.denominatorStat, 'gmean', 'geo stat surfaced');
  });

  test('AC-D4: aggregateNcaResults unmodified — single-key per-group mean still works', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('COMPOUND', ['X', 'X', 'Y']),
      DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 3.0, 5.0])),
    ]);
    const agg = aggregateNcaResults(df, 'COMPOUND', {params: ['CL']});
    expect(agg.col('MEAN_CL') !== null, true, 'MEAN_CL present');
    // X mean = 2.0, Y mean = 5.0 — order is groupBy-native, so look both up.
    const means = new Map<string, number>();
    for (let i = 0; i < agg.rowCount; i++)
      means.set(String(agg.col('COMPOUND')!.get(i)), agg.col('MEAN_CL')!.get(i) as number);
    expectFloat(means.get('X')!, 2.0, 1e-9, 'X mean unchanged');
    expectFloat(means.get('Y')!, 5.0, 1e-9, 'Y mean unchanged');
  });
});
