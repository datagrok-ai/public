/**
 * Test-only loaders for the sparse fixtures: the PO mouse dataset (04 / 05) and the
 * simulated IV-bolus dataset (07, `regen-sparse-fixture.R` block 4). Not a test file
 * (no `.test.ts`), shared by the sparse suites.
 */

import {readFileSync} from 'fs';
import {join} from 'path';
import type {SparseInput} from '../sparse';

export const DATASETS_DIR = join(__dirname, '..', '..', '__tests__', 'datasets');
export const FIXTURES_DIR = join(__dirname, '..', '..', '__tests__', 'fixtures');

/** One CSV row of the 04 / 07 layout (the two datasets share their columns). */
interface Row {mouseId: number; dose: number; nominalTime: number; concMgL: number; blq: number}

function loadRows(file: string): Row[] {
  const lines = readFileSync(join(DATASETS_DIR, file), 'utf8').trim().split(/\r?\n/);
  return lines.slice(1).map((line) => {
    const f = line.split(',');
    return {
      mouseId: Number(f[1].replace(/\D/g, '')),
      dose: Number(f[7]),
      nominalTime: Number(f[8]),
      concMgL: Number(f[10]) / 1000,
      blq: f[11] === 'True' ? 1 : 0,
    };
  });
}

function inputFor(rows: Row[], dose: number, lloq: number): SparseInput {
  const sub = rows.filter((r) => r.dose === dose);
  return {
    nominalTime: Float64Array.from(sub.map((r) => r.nominalTime)),
    conc: Float64Array.from(sub.map((r) => r.concMgL)),
    blqMask: Uint8Array.from(sub.map((r) => r.blq)),
    lloq,
    animalId: Int32Array.from(sub.map((r) => r.mouseId)),
  };
}

/** The PO mouse dataset, one input per dose group (mg/L, LLOQ 1 ng/mL). */
export function mouseInputs(): Array<{dose: number; input: SparseInput}> {
  const rows = loadRows('04_mouse_sparse_destructive.csv');
  return [...new Set(rows.map((r) => r.dose))].sort((a, b) => a - b)
    .map((dose) => ({dose, input: inputFor(rows, dose, 1 / 1000)}));
}

/** One composite block of `fixtures/07_iv_bolus_sparse.json`. */
export interface IvComposite {
  readonly dose_mg_kg: number;
  readonly times: number[];
  readonly path: 'logslope' | 'c1';
  readonly c0: number;
  readonly c0_method: string;
  readonly sparse_auclast_c0: number;
  readonly sparse_auclast: number;
  readonly sparse_auc_se: number;
  readonly sparse_auc_df: number;
  readonly sparse_auc_se_const_c0: number;
  readonly delta_se: number;
  readonly delta_df: number;
}

export const ivFixture: {composites: IvComposite[]} = JSON.parse(
  readFileSync(join(FIXTURES_DIR, '07_iv_bolus_sparse.json'), 'utf8'));

/** The IV-bolus dataset, one input per composite (mg/L, LLOQ 0.001 ng/mL). */
export function ivInputs(): Array<{block: IvComposite; input: SparseInput}> {
  const rows = loadRows('07_iv_bolus_sparse.csv');
  return ivFixture.composites.map((block) => ({block, input: inputFor(rows, block.dose_mg_kg, 1e-6)}));
}
