/**
 * Adapter — sparse / destructive-sampling NCA (UC-04 / FR-301..306).
 *
 * @see src/adapter/sparse-nca.ts (implementation under test)
 *
 * Two layers:
 *  - **Reference parity (AC-B4 / AC-D1, rule 18):** the real mouse-sparse fixture
 *    (`04_mouse_sparse_destructive.csv`) end-to-end through `computeSparseNca`,
 *    asserted against the PKNCA golden values in `mouse_sparse-reference.json`
 *    (AUC ≤ 0.1 %, SE/df ≤ 1 %).
 *  - **Behavioural (AC-B1..B3):** topology cross-check, linear-only SE, and the
 *    group-mean-BW CL/Vss + `SPARSE_COMPOSITE_BW` surface on synthetic frames.
 */

import {expectNoRawNan} from './test-utils';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat, before} from '@datagrok-libraries/test/src/test';

import {computeSparseNca} from '../adapter/sparse-nca';
import type {SparseNcaResult} from '../adapter/sparse-nca';
import {unmappedGridColumns} from '../views/results-view';

const TIMEOUT = 60000;
const FIXTURE_PATH_PREFIX = 'System:AppData/Nca';

interface DestructiveGolden {
  dose_mg_kg: number; sparse_auclast: number; sparse_auc_se: number; sparse_auc_df: number;
  composite: {nominal_time: number[]; n: number[]; n_blq: number[]; mean: number[]; sd: number[]};
}
interface CompositeNcaGolden {
  dose_mg_kg: number; auclast: number; aucinf_obs: number; lambda_z: number; half_life: number;
  cl_obs: number; vz_obs: number; aucpext_obs: number; lambda_z_n_points: number; adj_r_squared: number;
  mrt_obs: number;
}
interface SparseGolden {destructive: DestructiveGolden[]; composite_nca: {groups: CompositeNcaGolden[]};}

/**
 * Load the mouse-sparse CSV and build the `computeSparseNca` input frame:
 * Conc (mg/L = ng/mL ÷ 1000), a BOOL BLQ column from the "True"/"False" strings,
 * and BodyWeight_kg (= g ÷ 1000). Mirrors the units the fixture golden uses.
 */
async function loadMouseFrame(): Promise<DG.DataFrame> {
  const csv = await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/04_mouse_sparse_destructive.csv`);
  const src = DG.DataFrame.fromCsv(csv);
  const n = src.rowCount;
  const concNg = src.col('Conc_ng_mL')!;
  const blqStr = src.col('BLQ')!;
  const bwG = src.col('BodyWeight_g')!;
  const concMgL = new Float64Array(n);
  const blqBool: boolean[] = [];
  const bwKg = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    concMgL[i] = (concNg.isNone(i) ? 0 : Number(concNg.get(i))) / 1000;
    blqBool.push(blqStr.getString(i).trim().toLowerCase() === 'true');
    bwKg[i] = (bwG.isNone(i) ? NaN : Number(bwG.get(i))) / 1000;
  }
  src.columns.add(DG.Column.fromFloat64Array('Conc', concMgL));
  src.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'BLQflag', blqBool));
  src.columns.add(DG.Column.fromFloat64Array('BodyWeight_kg', bwKg));
  return src;
}

function runMouse(df: DG.DataFrame): SparseNcaResult {
  return computeSparseNca(df, {
    concCol: 'Conc',
    nominalTimeCol: 'NominalTime_hr',
    animalIdCol: 'MouseID',
    blqCol: 'BLQflag',
    doseCol: 'Dose_mg_kg',
    bodyWeightCol: 'BodyWeight_kg',
    samplingTypeCol: 'SamplingType',
    groupCol: 'Dose_mg_kg',
    route: 'PO',
    concUnit: 'mg/L',
    lloq: 0.001,
    bootstrapIterations: 200,
  });
}

category('Adapter: sparse NCA (UC-04)', () => {
  let golden: SparseGolden;
  let mouse: DG.DataFrame;

  before(async () => {
    const txt = await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/mouse_sparse-reference.json`);
    golden = JSON.parse(txt) as SparseGolden;
    mouse = await loadMouseFrame();
  });

  // BUG-12 drift guard: both emitted tables must resolve fully in their
  // presentation maps, or configureSparseNcaGrid / configureSparseCompositeGrid
  // leaves a column with a raw header and DG's generic stats tooltip. NOTHING is
  // exempt here: both tables are assembled from fixed `cols` arrays, and the
  // caller's `groupCol` is used only to PARTITION rows (sparse-nca.ts:288-315),
  // never emitted as a passthrough column — so any unmapped name is a real gap.
  test('BUG-12: every sparse / composite column has a presentation entry', async () => {
    const res = runMouse(mouse);
    // AC-13 (F4 Arm 1): this producer's REAL output must carry no raw IEEE NaN.
    // Asserted here rather than in a duplicate fixture — the frame already exists.
    expectNoRawNan(res.parameterTable, 'sparse-nca.parameterTable');
    expectNoRawNan(res.compositeTable, 'sparse-nca.compositeTable');
    const p = unmappedGridColumns(res.parameterTable, 'sparseNca');
    expect(p.length, 0, `sparse columns with no SPARSE_NCA_FIXED entry: [${p.join(', ')}]`);
    const c = unmappedGridColumns(res.compositeTable, 'sparseComposite');
    expect(c.length, 0,
      `composite columns with no SPARSE_COMPOSITE_FIXED entry: [${c.join(', ')}]`);
  }, {timeout: TIMEOUT});

  // --- AC-B4 / AC-D1: reference parity against PKNCA golden ---------------
  test('AC-D1: destructive AUClast / SE / df match PKNCA within NFR-04/05', async () => {
    const res = runMouse(mouse);
    for (const g of golden.destructive) {
      const grp = res.groups.find((x) => Math.abs(x.doseMgKg - g.dose_mg_kg) < 1e-9);
      expect(grp !== undefined, true, `group ${g.dose_mg_kg} present`);
      expect(grp!.topology, 'destructive', `dose ${g.dose_mg_kg} topology`);
      expectFloat(grp!.aucLast, g.sparse_auclast, Math.abs(g.sparse_auclast) * 0.001,
        `dose ${g.dose_mg_kg}: AUClast out of 0.1%`);
      expectFloat(grp!.aucSe, g.sparse_auc_se, Math.abs(g.sparse_auc_se) * 0.01,
        `dose ${g.dose_mg_kg}: AUC SE out of 1%`);
      expectFloat(grp!.aucDf, g.sparse_auc_df, Math.abs(g.sparse_auc_df) * 0.01,
        `dose ${g.dose_mg_kg}: AUC df out of 1%`);
    }
  }, {timeout: TIMEOUT});

  // --- IMPL-GAP-SPARSE-AUC-METHOD-MIX: one (linear) AUC method per composite row ------------------
  // PKNCA `pk.nca` on the composite mean with auc.method = "linear" (golden `composite_nca`).
  // Tolerances: AUC 0.1 %, λz 0.5 %, derived 1 % (NFR-04 / NFR-05, rule 18).
  test('AUC method mix: composite AUCinf / λz / t½ / CL / Vz / %extrap match PKNCA (linear)', async () => {
    const res = runMouse(mouse);
    for (const g of golden.composite_nca.groups) {
      const grp = res.groups.find((x) => Math.abs(x.doseMgKg - g.dose_mg_kg) < 1e-9);
      expect(grp !== undefined, true, `group ${g.dose_mg_kg} present`);
      const d = `dose ${g.dose_mg_kg}`;
      expectFloat(grp!.aucLast, g.auclast, g.auclast * 0.001, `${d}: AUClast out of 0.1%`);
      expectFloat(grp!.aucInf, g.aucinf_obs, g.aucinf_obs * 0.001, `${d}: AUCinf out of 0.1%`);
      expectFloat(grp!.lambdaZ, g.lambda_z, g.lambda_z * 0.005, `${d}: λz out of 0.5%`);
      expect(grp!.lambdaZNPoints, g.lambda_z_n_points, `${d}: λz points`);
      expectFloat(grp!.halfLife, g.half_life, g.half_life * 0.01, `${d}: t½ out of 1%`);
      expectFloat(grp!.pctExtrap, g.aucpext_obs, g.aucpext_obs * 0.01, `${d}: %AUCextrap out of 1%`);
      const clMlMinKg = g.cl_obs * 1000 / 60; // L/h/kg → mL/min/kg (the per-kg CL unit)
      expectFloat(grp!.cl, clMlMinKg, clMlMinKg * 0.01, `${d}: CL out of 1%`);
      expectFloat(grp!.vz, g.vz_obs, g.vz_obs * 0.01, `${d}: Vz out of 1%`);
      expectFloat(grp!.mrt, g.mrt_obs, g.mrt_obs * 0.01, `${d}: MRT (linear AUMC) out of 1%`);
      expectFloat(grp!.lambdaZAdjR2, g.adj_r_squared, Math.abs(g.adj_r_squared) * 0.005, `${d}: adj R² out of 0.5%`);
      // A fragile fit ships flagged, never hidden (thesis #4): adj R² 0.099 on the 10 mg/kg composite.
      if (g.adj_r_squared < 0.9) {
        expect(grp!.acceptance.verdict, 'Flagged', `${d}: fragile λz fit is Flagged`);
        expect(grp!.acceptance.flagCodes.includes('ADJR2_LT'), true, `${d}: adj R² criterion violated`);
      }
    }
  }, {timeout: TIMEOUT});

  // REGRESSION_PIN for the old bug (the PKNCA test above is the ground truth). Scoped to this
  // extravascular fixture: on an IV-bolus composite the dense run back-extrapolates C0 while
  // `sparseAuc` anchors at (0, 0), so %AUCextrap is not relative to the shown AUClast there.
  test('AUC method mix: AUCinf ≥ AUClast and %AUCextrap rests on the shown AUClast', async () => {
    const res = runMouse(mouse);
    for (const g of res.groups) {
      const d = `dose ${g.doseMgKg}`;
      expect(g.aucInf >= g.aucLast, true, `${d}: AUCinf ${g.aucInf} < AUClast ${g.aucLast}`);
      const pct = (g.aucInf - g.aucLast) / g.aucInf * 100;
      expectFloat(g.pctExtrap, pct, 1e-9, `${d}: %AUCextrap is relative to the shown AUClast`);
      // Direct: the dense run's AUClast (AUCinf − Clast/λz, Clast = last positive composite mean)
      // is the shown closed-form AUClast.
      const tps = g.composite.timepoints;
      const clast = [...tps].reverse().find((t) => t.n > 0 && t.nBlq / t.n <= 0.5 && t.mean > 0)!.mean;
      expectFloat(g.aucInf - clast / g.lambdaZ, g.aucLast, g.aucLast * 1e-9, `${d}: dense AUClast = shown AUClast`);
    }
  }, {timeout: TIMEOUT});

  // REGRESSION_PIN of the disclosed IV-bolus residual (TODO IMPL-GAP-SPARSE-IV-BOLUS-ANCHOR, user
  // decision 2026-11-09: disclose now, decide the anchor in its own cycle). Measured on this composite
  // as IV bolus: closed-form AUClast 3.821 from (0, 0), dense AUClast 3.983 from the back-extrapolated C0.
  test('AUC method mix: an IV-bolus composite discloses its two AUClast anchors; PO does not', async () => {
    const po = runMouse(mouse);
    for (const g of po.groups)
      expect(g.warnings.some((w) => w.code === 'SPARSE_IV_BOLUS_C0_ANCHOR'), false, `PO dose ${g.doseMgKg}: no note`);
    const iv = computeSparseNca(mouse, {
      concCol: 'Conc', nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID', blqCol: 'BLQflag',
      doseCol: 'Dose_mg_kg', bodyWeightCol: 'BodyWeight_kg', groupCol: 'Dose_mg_kg',
      route: 'IV-bolus', concUnit: 'mg/L', lloq: 0.001, bootstrapIterations: 200,
    });
    for (const g of iv.groups) {
      const d = `IV dose ${g.doseMgKg}`;
      const note = g.warnings.find((w) => w.code === 'SPARSE_IV_BOLUS_C0_ANCHOR');
      expect(note !== undefined, true, `${d}: the anchor note is present`);
      expect(note!.message.includes(g.aucLast.toPrecision(4)), true, `${d}: the note names the shown AUClast`);
      expect(g.aucInf >= g.aucLast, true, `${d}: AUCinf ${g.aucInf} < AUClast ${g.aucLast}`);
    }
    expect(iv.parameterTable.col('NOTES')!.get(0).includes('IV bolus: AUClast from (0, 0)'), true,
      'the short note reaches the Notes column');
  }, {timeout: TIMEOUT});

  test('AC-D1: composite means/SD/n/%BLQ match the oracle', async () => {
    const res = runMouse(mouse);
    for (const g of golden.destructive) {
      const grp = res.groups.find((x) => Math.abs(x.doseMgKg - g.dose_mg_kg) < 1e-9)!;
      const tps = grp.composite.timepoints;
      expect(tps.length, g.composite.nominal_time.length, `dose ${g.dose_mg_kg} timepoint count`);
      tps.forEach((tp, i) => {
        expect(tp.n, g.composite.n[i], `dose ${g.dose_mg_kg} t${tp.time} n`);
        expect(tp.nBlq, g.composite.n_blq[i], `dose ${g.dose_mg_kg} t${tp.time} nBlq`);
        expectFloat(tp.mean, g.composite.mean[i], 1e-6, `dose ${g.dose_mg_kg} t${tp.time} mean`);
        expectFloat(tp.sd, g.composite.sd[i], 1e-6, `dose ${g.dose_mg_kg} t${tp.time} sd`);
      });
    }
  }, {timeout: TIMEOUT});

  // --- AC-B3: CL/Vss ship; the group-mean BW forms the absolute companions (SPARSE_COMPOSITE_BW) ---
  test('AC-B3: CL/Vss ship; CL abs. on group-mean BW with the SPARSE_COMPOSITE_BW note', async () => {
    const res = runMouse(mouse);
    const g10 = res.groups.find((x) => x.doseMgKg === 10)!;
    expect(Number.isFinite(g10.meanBodyWeightKg!), true, 'group-mean BW present');
    expect(Number.isFinite(g10.cl) && g10.cl > 0, true, 'CL ships (finite, positive)');
    expect(g10.warnings.some((w) => w.code === 'SPARSE_COMPOSITE_BW'), true,
      'SPARSE_COMPOSITE_BW provenance note present');
  }, {timeout: TIMEOUT});

  // --- nca-cl-bw-normalised-units D-d / AD-5 (amending uc-04 SF-1) -----------------------------
  // The oracle carries no CL, so the per-kg composite is asserted as the identity D-d defines on the
  // engine's OWN composite AUCinf: per kg = dose per kg / AUCinf (× 1000/60), and CL_ABS = dose per kg ×
  // group-mean BW / AUCinf — exactly the absolute value the pre-change group-mean-BW run reported.
  test('AC-4.1: /kg composite → CL = D_kg / AUCinf per kg; CL_ABS = D_kg × mean BW / AUCinf', async () => {
    const res = runMouse(mouse);
    for (const g of res.groups) {
      expect(Number.isFinite(g.aucInf), true, `dose ${g.doseMgKg}: AUCinf estimated`);
      const perKgLhKg = g.doseMgKg / g.aucInf;
      expectFloat(g.cl, perKgLhKg * 1000 / 60, perKgLhKg * 1000 / 60 * 1e-9, `dose ${g.doseMgKg}: CL per kg`);
      const abs = g.doseMgKg * g.meanBodyWeightKg! / g.aucInf;
      expectFloat(g.clAbs, abs, abs * 1e-12, `dose ${g.doseMgKg}: CL abs. = the group-mean-BW absolute`);
      expect(g.engineDoseNorm, 'per_kg', 'the engine ran on the dose per kg');
      expect(g.warnings.some((w) => w.code === 'SPARSE_PER_KG_WEIGHT_BASED'), true, 'assumption stated');
    }
    const t = res.parameterTable;
    expect(t.col('CL')!.getTag('units'), 'mL/min/kg', 'AC-4.4: per-kg unit tag');
    expect(t.col('CL_ABS')!.getTag('units'), 'L/h', 'companion unit tag');
    expect(t.col('BW_KG')!.getTag('units'), 'kg', 'group-mean BW unit tag');
  }, {timeout: TIMEOUT});

  test('AC-4.1: without a body weight the per-kg composite is still formed; CL_ABS is not', async () => {
    const noBw = computeSparseNca(mouse, {
      concCol: 'Conc', nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID', blqCol: 'BLQflag',
      doseCol: 'Dose_mg_kg', groupCol: 'Dose_mg_kg', route: 'PO', concUnit: 'mg/L', lloq: 0.001,
      bootstrapIterations: 50,
    });
    const withBw = runMouse(mouse);
    for (const g of noBw.groups) {
      const twin = withBw.groups.find((x) => x.doseMgKg === g.doseMgKg)!;
      expectFloat(g.cl, twin.cl, Math.abs(twin.cl) * 1e-12, `dose ${g.doseMgKg}: per kg needs no weight`);
      expect(Number.isNaN(g.clAbs), true, `dose ${g.doseMgKg}: no weight → no CL abs.`);
      expect(g.warnings.some((w) => w.code === 'SPARSE_NO_DOSE_NO_CL'), false, 'a dose is enough for CL');
    }
  }, {timeout: TIMEOUT});

  test('AC-4.3: the bootstrap CL interval is per kg — the point estimate sits inside it', async () => {
    const res = runMouse(mouse);
    for (const g of res.groups) {
      const b = g.bootstrapCl;
      if (b.suppressed || !Number.isFinite(b.ciLow)) continue;
      expect(b.ciLow <= g.cl && g.cl <= b.ciHigh, true,
        `dose ${g.doseMgKg}: CL ${g.cl} in [${b.ciLow}, ${b.ciHigh}] (mL/min/kg)`);
    }
  }, {timeout: TIMEOUT});

  test('AC-4.2: absolute-mg dose column → per kg = abs / mean BW; no BW → PER_KG_NOT_FORMED', async () => {
    const df = mouse.clone();
    const d = df.col('Dose_mg_kg')!; const bw = df.col('BodyWeight_kg')!;
    df.columns.add(DG.Column.fromFloat64Array('Dose_mg', Float64Array.from({length: df.rowCount},
      (_, i) => Number(d.get(i)) * (bw.get(i) as number))));
    const base = {concCol: 'Conc', nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID', blqCol: 'BLQflag',
      doseMgCol: 'Dose_mg', groupCol: 'Dose_mg_kg', route: 'PO', concUnit: 'mg/L' as const, lloq: 0.001,
      bootstrapIterations: 50};
    const withBw = computeSparseNca(df, {...base, bodyWeightCol: 'BodyWeight_kg'});
    const perKgRun = runMouse(mouse);
    for (const g of withBw.groups) {
      expect(g.engineDoseNorm, 'absolute', 'the engine ran on the mg dose');
      expectFloat(g.cl, g.clAbs / g.meanBodyWeightKg! * 1000 / 60, Math.abs(g.cl) * 1e-6,
        `dose ${g.doseLevel.value}: per kg = abs / mean BW`);
      const twin = perKgRun.groups.find((x) => x.doseMgKg === Number(g.doseLevel.value))!;
      // Weight-based mg doses: mean(D·BWᵢ)/mean(BWᵢ) = D, so the two per-kg values agree.
      expectFloat(g.cl, twin.cl, Math.abs(twin.cl) * 1e-9, `dose ${g.doseLevel.value}: same per kg as the /kg run`);
      expect(g.warnings.some((w) => w.code === 'COMPOSITE_NOT_WEIGHT_BASED'), false, 'mg ∝ BW here');
    }
    const noBw = computeSparseNca(df, base);
    for (const g of noBw.groups) {
      expect(Number.isNaN(g.cl), true, 'absolute mg, no BW → per kg not formed');
      expect(Number.isFinite(g.clAbs), true, 'CL abs. reported');
      expect(g.warnings.some((w) => w.code === 'PER_KG_NOT_FORMED'), true, 'the reason is stated');
    }
  }, {timeout: TIMEOUT});

  // --- AC-B2: closed-form SE is the linear-trapezoidal Bailer (destructive) -
  test('AC-B2: destructive design uses the Bailer closed-form CI method', async () => {
    const res = runMouse(mouse);
    const g10 = res.groups.find((x) => x.doseMgKg === 10)!;
    expect(g10.aucCiMethod, 'Bailer', 'destructive → Bailer closed-form');
    expect(g10.aucCi[0] < g10.aucLast && g10.aucCi[1] > g10.aucLast, true, 'CI brackets AUClast');
  }, {timeout: TIMEOUT});

  // --- AC-B1: declared-destructive but paired animals → topology mismatch --
  test('AC-B1: declared destructive with paired animals warns SPARSE_TOPOLOGY_MISMATCH', async () => {
    // 3 timepoints, 4 animals, each sampled twice (batch), labelled "Destructive".
    const animal = ['A', 'B', 'A', 'C', 'B', 'C'];
    const time = [1, 1, 2, 2, 4, 4];
    const conc = [10, 12, 6, 8, 2, 3];
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Animal', animal),
      DG.Column.fromFloat64Array('NomTime', Float64Array.from(time)),
      DG.Column.fromFloat64Array('Conc', Float64Array.from(conc)),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SType', animal.map(() => 'Destructive')),
      DG.Column.fromFloat64Array('Dose', Float64Array.from(animal.map(() => 10))),
    ]);
    const res = computeSparseNca(df, {
      concCol: 'Conc', nominalTimeCol: 'NomTime', animalIdCol: 'Animal',
      samplingTypeCol: 'SType', doseCol: 'Dose', groupCol: 'Dose', route: 'PO',
      bootstrapIterations: 50,
    });
    const g = res.groups[0];
    expect(g.topology, 'batch', 'derived topology is batch (paired animals)');
    expect(g.warnings.some((w) => w.code === 'SPARSE_TOPOLOGY_MISMATCH'), true,
      'SPARSE_TOPOLOGY_MISMATCH raised against the "Destructive" label');
  }, {timeout: TIMEOUT});

  // --- UC-04 mapping UI: scalar declaredTopology override (AC-A2) ----------
  // The Sampling-type form choice declares the topology WITHOUT a SamplingType
  // column; it routes through the `declaredTopology` option to the same AD-2
  // cross-check. Paired animals + a 'destructive' declaration ⇒ mismatch, proving
  // the scalar path alone (no samplingTypeCol) surfaces SPARSE_TOPOLOGY_MISMATCH.
  test('AC-A2: scalar declaredTopology=destructive on paired animals warns mismatch', async () => {
    const animal = ['A', 'B', 'A', 'C', 'B', 'C'];
    const time = [1, 1, 2, 2, 4, 4];
    const conc = [10, 12, 6, 8, 2, 3];
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Animal', animal),
      DG.Column.fromFloat64Array('NomTime', Float64Array.from(time)),
      DG.Column.fromFloat64Array('Conc', Float64Array.from(conc)),
      DG.Column.fromFloat64Array('Dose', Float64Array.from(animal.map(() => 10))),
    ]);
    const res = computeSparseNca(df, {
      concCol: 'Conc', nominalTimeCol: 'NomTime', animalIdCol: 'Animal',
      declaredTopology: 'destructive', // no samplingTypeCol — scalar path only
      doseCol: 'Dose', groupCol: 'Dose', route: 'PO', bootstrapIterations: 50,
    });
    const g = res.groups[0];
    expect(g.topology, 'batch', 'derived topology is batch (paired animals)');
    expect(g.warnings.some((w) => w.code === 'SPARSE_TOPOLOGY_MISMATCH'), true,
      'SPARSE_TOPOLOGY_MISMATCH raised against the scalar "destructive" declaration');
  }, {timeout: TIMEOUT});

  // Zero-regression: a serial declaration consistent with paired data raises NO
  // mismatch (the override is cross-checked, not blindly trusted).
  test('AC-A2: scalar declaredTopology=batch consistent with paired data → no mismatch', async () => {
    const animal = ['A', 'B', 'A', 'C', 'B', 'C'];
    const time = [1, 1, 2, 2, 4, 4];
    const conc = [10, 12, 6, 8, 2, 3];
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Animal', animal),
      DG.Column.fromFloat64Array('NomTime', Float64Array.from(time)),
      DG.Column.fromFloat64Array('Conc', Float64Array.from(conc)),
      DG.Column.fromFloat64Array('Dose', Float64Array.from(animal.map(() => 10))),
    ]);
    const res = computeSparseNca(df, {
      concCol: 'Conc', nominalTimeCol: 'NomTime', animalIdCol: 'Animal',
      declaredTopology: 'batch', doseCol: 'Dose', groupCol: 'Dose', route: 'PO',
      bootstrapIterations: 50,
    });
    expect(res.groups[0].warnings.some((w) => w.code === 'SPARSE_TOPOLOGY_MISMATCH'), false,
      'consistent batch declaration → no mismatch warning');
  }, {timeout: TIMEOUT});

  // --- Honest bootstrap framing (SF-2) ------------------------------------
  test('exploratory bootstrap AUC CI is annotated, never the headline', async () => {
    const res = runMouse(mouse);
    const g10 = res.groups.find((x) => x.doseMgKg === 10)!;
    // n=3 design is above the gate → bootstrap runs and is flagged exploratory.
    const exploratory = res.warnings.some((w) =>
      w.group === g10.groupKey && w.code === 'SPARSE_BOOTSTRAP_AUC_EXPLORATORY');
    const suppressed = res.warnings.some((w) =>
      w.group === g10.groupKey && w.code === 'SPARSE_BOOTSTRAP_N_TOO_SMALL');
    expect(exploratory || suppressed, true, 'bootstrap AUC is either exploratory-flagged or suppressed');
    // The closed-form CI is always present regardless of the bootstrap.
    expect(Number.isFinite(g10.aucCi[0]) && Number.isFinite(g10.aucCi[1]), true,
      'closed-form CI is the headline, always present');
  }, {timeout: TIMEOUT});
});
