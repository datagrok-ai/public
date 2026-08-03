/**
 * Sparse NCA validation (UC-04 / FR-301..306, AC-U1..U4).
 *
 * Validates `sparseAuc` / `buildCompositeProfile` / `summarizeBootstrap` against
 * the rule-18 oracle `__tests__/fixtures/05_mouse_sparse.json`:
 *   - DESTRUCTIVE: PKNCA 0.12.1 `pk.calc.sparse_auclast` (auc/se/df)
 *   - BATCH:       hand-derived Holder A1+A3 + Nedelman-Jia df (PKNCA df=NA)
 *
 * Tolerances: AUC ≤ 0.1 % (NFR-04); SE, df ≤ 1 % (NFR-05).
 */

import {readFileSync} from 'fs';
import {join} from 'path';
import {sparseAuc, buildCompositeProfile} from '../sparse';
import type {SparseInput} from '../sparse';
import {summarizeBootstrap} from '../bootstrap';

const DATASETS_DIR = join(__dirname, '..', '..', '__tests__', 'datasets');
const FIXTURES_DIR = join(__dirname, '..', '..', '__tests__', 'fixtures');

interface CompositeBlock {
  nominal_time: number[]; n: number[]; n_blq: number[];
  mean: number[]; sd: number[]; cv_pct: Array<number | string>;
}
interface DestructiveBlock {
  dose_mg_kg: number; sparse_auclast: number; sparse_auc_se: number;
  sparse_auc_df: number; composite: CompositeBlock;
}
interface BatchBlock {
  times: number[];
  samples: Array<{time: number; animal: string[]; conc: number[]}>;
  sparse_auclast: number; sparse_auc_se: number; sparse_auc_df: number;
}
interface SparseFixture {
  destructive: DestructiveBlock[]; batch: BatchBlock;
}

const fixture: SparseFixture = JSON.parse(
  readFileSync(join(FIXTURES_DIR, '05_mouse_sparse.json'), 'utf8'));

// --- Parse the mouse dataset into per-dose SparseInput (conc in mg/L). --------
interface MouseRow {dose: number; mouseId: number; nominalTime: number; concMgL: number; blq: number;}

function loadMouseRows(): MouseRow[] {
  const csv = readFileSync(join(DATASETS_DIR, '04_mouse_sparse_destructive.csv'), 'utf8');
  const lines = csv.trim().split(/\r?\n/);
  const rows: MouseRow[] = [];
  for (let i = 1; i < lines.length; i++) {
    const f = lines[i].split(',');
    rows.push({
      mouseId: Number(f[1].replace(/\D/g, '')),
      dose: Number(f[7]),
      nominalTime: Number(f[8]),
      concMgL: Number(f[10]) / 1000,
      blq: f[11] === 'True' ? 1 : 0,
    });
  }
  return rows;
}

function sparseInputForDose(rows: MouseRow[], dose: number): SparseInput {
  const sub = rows.filter((r) => r.dose === dose);
  return {
    nominalTime: Float64Array.from(sub.map((r) => r.nominalTime)),
    conc: Float64Array.from(sub.map((r) => r.concMgL)),
    blqMask: Uint8Array.from(sub.map((r) => r.blq)),
    lloq: 1.0 / 1000,
    animalId: Int32Array.from(sub.map((r) => r.mouseId)),
  };
}

const rel = (a: number, b: number): number => Math.abs(a - b) / Math.abs(b);

describe('sparseAuc — destructive (AC-U1, PKNCA oracle)', () => {
  const rows = loadMouseRows();
  for (const block of fixture.destructive) {
    it(`dose ${block.dose_mg_kg}: auc/se/df match PKNCA within NFR-04/05`, () => {
      const res = sparseAuc(sparseInputForDose(rows, block.dose_mg_kg), {blqRule: 'set-zero'});
      expect(res.topology).toBe('destructive');
      expect(rel(res.auc, block.sparse_auclast) < 1e-3).toBe(true); // 0.1% AUC
      expect(rel(res.se, block.sparse_auc_se) < 1e-2).toBe(true); // 1% SE
      expect(rel(res.df, block.sparse_auc_df) < 1e-2).toBe(true); // 1% df
    });
  }

  it('matches PKNCA to floating-point round-off (not just within tolerance)', () => {
    const block = fixture.destructive[0];
    const res = sparseAuc(sparseInputForDose(rows, block.dose_mg_kg), {blqRule: 'set-zero'});
    expect(rel(res.auc, block.sparse_auclast) < 1e-9).toBe(true);
    expect(rel(res.se, block.sparse_auc_se) < 1e-9).toBe(true);
    expect(rel(res.df, block.sparse_auc_df) < 1e-9).toBe(true);
  });
});

describe('buildCompositeProfile — composite means/SD/%CV (AC-U3)', () => {
  const rows = loadMouseRows();
  for (const block of fixture.destructive) {
    it(`dose ${block.dose_mg_kg}: per-timepoint mean/SD/n/%BLQ match the oracle`, () => {
      const comp = buildCompositeProfile(sparseInputForDose(rows, block.dose_mg_kg), 'set-zero');
      const c = block.composite;
      expect(comp.timepoints.length).toBe(c.nominal_time.length);
      comp.timepoints.forEach((tp, i) => {
        expect(tp.time).toBeCloseTo(c.nominal_time[i], 10);
        expect(tp.n).toBe(c.n[i]);
        expect(tp.nBlq).toBe(c.n_blq[i]);
        expect(Math.abs(tp.mean - c.mean[i]) < 1e-9).toBe(true);
        expect(Math.abs(tp.sd - c.sd[i]) < 1e-9).toBe(true);
      });
    });
  }
});

describe('sparseAuc — batch Holder A1+A3 + Nedelman-Jia df (AC-U2)', () => {
  function batchInput(): SparseInput {
    const b = fixture.batch;
    const nominalTime: number[] = []; const conc: number[] = []; const animalId: number[] = [];
    const codeOf = new Map<string, number>();
    for (const s of b.samples) {
      s.animal.forEach((a, k) => {
        if (!codeOf.has(a)) codeOf.set(a, codeOf.size + 1);
        nominalTime.push(s.time); conc.push(s.conc[k]); animalId.push(codeOf.get(a)!);
      });
    }
    return {
      nominalTime: Float64Array.from(nominalTime),
      conc: Float64Array.from(conc),
      blqMask: new Uint8Array(nominalTime.length),
      lloq: 0,
      animalId: Int32Array.from(animalId),
    };
  }

  it('reproduces the hand-derived batch auc/se/df oracle', () => {
    const res = sparseAuc(batchInput(), {blqRule: 'set-zero'});
    expect(res.topology).toBe('batch');
    expect(rel(res.auc, fixture.batch.sparse_auclast) < 1e-9).toBe(true);
    expect(rel(res.se, fixture.batch.sparse_auc_se) < 1e-9).toBe(true);
    expect(rel(res.df, fixture.batch.sparse_auc_df) < 1e-9).toBe(true);
  });

  it('Holder reduces to Bailer when r_ij = 0 (destructive ids on the same profile)', () => {
    // Same concentrations/times, but every observation a distinct animal → r_ij=0.
    const b = batchInput();
    const destructive: SparseInput = {...b, animalId: Int32Array.from(b.nominalTime.map((_, i) => i + 1))};
    const withOverlap = sparseAuc(b, {blqRule: 'set-zero'});
    const noOverlap = sparseAuc(destructive, {blqRule: 'set-zero'});
    expect(noOverlap.topology).toBe('destructive');
    // AUC identical (same means); SE differs (no covariance term).
    expect(rel(noOverlap.auc, withOverlap.auc) < 1e-12).toBe(true);
    expect(noOverlap.se).not.toBeCloseTo(withOverlap.se, 6);
  });
});

describe('sparseAuc — scientific guards', () => {
  const rows = loadMouseRows();

  it('AD-2: declared topology mismatch warns', () => {
    const res = sparseAuc(sparseInputForDose(rows, 10), {declaredTopology: 'batch'});
    expect(res.warnings.some((w) => w.code === 'SPARSE_TOPOLOGY_MISMATCH')).toBe(true);
  });

  it('AD-11: absent animal-ID assumes destructive + warns', () => {
    const inp = sparseInputForDose(rows, 10);
    const res = sparseAuc({...inp, animalId: null});
    expect(res.topology).toBe('destructive');
    expect(res.warnings.some((w) => w.code === 'SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE')).toBe(true);
  });

  it('AD-7: a singleton timepoint gets a modeled variance + SPARSE_VARIANCE_MODELED', () => {
    // 3 timepoints: two with n=3 (real variance), one with n=1 (modeled).
    const input: SparseInput = {
      nominalTime: Float64Array.from([1, 1, 1, 2, 2, 2, 4]),
      conc: Float64Array.from([10, 12, 11, 5, 6, 5.5, 2]),
      blqMask: new Uint8Array(7),
      lloq: 0,
      animalId: Int32Array.from([1, 2, 3, 4, 5, 6, 7]),
    };
    const withBorrow = sparseAuc(input, {varianceBorrow: true});
    const noBorrow = sparseAuc(input, {varianceBorrow: false});
    expect(withBorrow.warnings.some((w) => w.code === 'SPARSE_VARIANCE_MODELED')).toBe(true);
    // Borrowed variance widens the SE vs silently zeroing the singleton's variance.
    expect(withBorrow.se).toBeGreaterThan(noBorrow.se);
  });

  it('AD-3: a steep wide-gap terminal phase flags SPARSE_TERMINAL_OVEREST', () => {
    const res = sparseAuc(sparseInputForDose(rows, 10), {estimatedHalfLife: 1.5});
    expect(res.warnings.some((w) => w.code === 'SPARSE_TERMINAL_OVEREST')).toBe(true);
  });
});

describe('summarizeBootstrap (AC-U4)', () => {
  const rows = loadMouseRows();

  it('is reproducible at a fixed master seed', () => {
    const inp = sparseInputForDose(rows, 10);
    const stat = (r: SparseInput): number => sparseAuc(r).auc;
    const a = summarizeBootstrap(inp, stat, {iterations: 300, masterSeed: 7});
    const b = summarizeBootstrap(inp, stat, {iterations: 300, masterSeed: 7});
    expect(a.ci[0]).toBe(b.ci[0]);
    expect(a.ci[1]).toBe(b.ci[1]);
    expect(a.median).toBe(b.median);
  });

  it('median ≈ point estimate and CI brackets it (n=3 design, not suppressed)', () => {
    const inp = sparseInputForDose(rows, 10);
    const stat = (r: SparseInput): number => sparseAuc(r).auc;
    const s = summarizeBootstrap(inp, stat, {iterations: 800, masterSeed: 42});
    expect(s.suppressed).toBe(false);
    expect(rel(s.median, s.estimate) < 0.25).toBe(true);
    expect(s.ci[0]).toBeLessThanOrEqual(s.estimate);
    expect(s.ci[1]).toBeGreaterThanOrEqual(s.estimate);
  });

  it('suppresses a 5-timepoint × 2-animal destructive design unconditionally (NB-2)', () => {
    // 5 timepoints, 2 animals each → h = 3^5 = 243 ≤ 360.
    const nominalTime: number[] = []; const conc: number[] = []; const animalId: number[] = [];
    let id = 1;
    [1, 2, 4, 8, 12].forEach((t, ti) => {
      for (let k = 0; k < 2; k++) {nominalTime.push(t); conc.push(10 / (ti + 1) + k); animalId.push(id++);}
    });
    const inp: SparseInput = {
      nominalTime: Float64Array.from(nominalTime),
      conc: Float64Array.from(conc),
      blqMask: new Uint8Array(nominalTime.length),
      lloq: 0,
      animalId: Int32Array.from(animalId),
    };
    const s = summarizeBootstrap(inp, (r) => sparseAuc(r).auc, {minNPerTimepoint: 2});
    expect(s.suppressed).toBe(true);
    expect(Number.isNaN(s.ci[0])).toBe(true);
  });
});
