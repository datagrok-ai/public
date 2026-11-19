/** A temporary method to check molecule features end to end; delete it with its tests (core TODO "Forge: delete the
 * check method"). */
import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {isRecord} from '../preparation/preparation-options';
import {isReadableNumber} from '../training/default-features';

export const CHECK_METHOD = 'Forge Check (molecules)';
const MOLECULE = 'molecule';
const RIDGE = 1e-9;
const ZERO_PIVOT = 1e-12;

interface CheckModel { features: string[]; coefficients: number[]; intercept: number }

/** A numerical [target], at least one Molecule column in [df] and every other column numerical. */
export function isCheckApplicable(df: DG.DataFrame, target: DG.Column): boolean {
  const columns = df.columns.toList();
  const isMolecule = (c: DG.Column) => c.semType?.toLowerCase() === MOLECULE;
  return target.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME) && columns.some(isMolecule) &&
    columns.every((c) => isMolecule(c) || isReadableNumber(c));
}

/** Ordinary least squares with an intercept on the numerical columns of [df] (molecules ignored); without one, a
 * model that predicts 0. */
export function checkTrain(df: DG.DataFrame, target: DG.Column): Uint8Array {
  const columns = df.columns.toList().filter(isReadableNumber);
  const [intercept, ...coefficients] = columns.length > 0 ? leastSquares(columns, target, df.rowCount) : [0];
  const model: CheckModel = {features: columns.map((c) => c.name), coefficients, intercept};
  return new TextEncoder().encode(JSON.stringify(model));
}

/** The predictions of the model [blob] ({@link checkTrain}) for the rows of [df], in a column `prediction`. */
export function checkApply(df: DG.DataFrame, blob: unknown): DG.DataFrame {
  const model = checkModelOf(blob);
  const columns = model.features.map((name) => df.getCol(name));
  const values = new Float64Array(df.rowCount);
  for (let row = 0; row < df.rowCount; row++) {
    let value = model.intercept;
    for (let k = 0; k < columns.length; k++)
      value += model.coefficients[k] * columns[k].getNumber(row);
    values[row] = value;
  }
  return DG.DataFrame.fromColumns([DG.Column.fromFloat64Array('prediction', values)]);
}

/** The intercept, then the coefficients of [columns], of the ridge-stabilized normal equations for [target]. */
function leastSquares(columns: DG.Column[], target: DG.Column, rowCount: number): number[] {
  const n = columns.length + 1;
  const a = Array.from({length: n}, () => new Array<number>(n).fill(0));
  const b = new Array<number>(n).fill(0);
  const x = new Array<number>(n);
  for (let row = 0; row < rowCount; row++) {
    x[0] = 1;
    for (let k = 0; k < columns.length; k++)
      x[k + 1] = columns[k].getNumber(row);
    const y = target.getNumber(row);
    for (let i = 0; i < n; i++) {
      b[i] += x[i] * y;
      for (let j = 0; j < n; j++)
        a[i][j] += x[i] * x[j];
    }
  }
  for (let i = 1; i < n; i++)
    a[i][i] += RIDGE;
  return solve(a, b);
}

function checkModelOf(blob: unknown): CheckModel {
  let parsed: unknown;
  try {
    parsed = blob instanceof Uint8Array ? JSON.parse(new TextDecoder().decode(blob)) : undefined;
  } catch {
    parsed = undefined;
  }
  const isNumbers = (v: unknown): v is number[] => Array.isArray(v) && v.every((n) => typeof n === 'number');
  if (!isRecord(parsed) || !Array.isArray(parsed.features) || !parsed.features.every((f) => typeof f === 'string') ||
    !isNumbers(parsed.coefficients) || typeof parsed.intercept !== 'number' ||
    parsed.coefficients.length !== parsed.features.length)
    throw new ForgeError('The model file of Forge Check is damaged.');
  return {features: parsed.features, coefficients: parsed.coefficients, intercept: parsed.intercept};
}

/** The solution of [a] x = [b] by Gaussian elimination with partial pivoting; a zero pivot gives its unknown 0. */
function solve(a: number[][], b: number[]): number[] {
  const n = b.length;
  const m = a.map((row, i) => [...row, b[i]]);
  const isFree = new Array<boolean>(n).fill(false);
  for (let col = 0; col < n; col++) {
    let pivot = col;
    for (let r = col + 1; r < n; r++) {
      if (Math.abs(m[r][col]) > Math.abs(m[pivot][col]))
        pivot = r;
    }
    if (Math.abs(m[pivot][col]) < ZERO_PIVOT) {
      isFree[col] = true;
      continue;
    }
    [m[col], m[pivot]] = [m[pivot], m[col]];
    for (let r = col + 1; r < n; r++) {
      const factor = m[r][col] / m[col][col];
      for (let j = col; j <= n; j++)
        m[r][j] -= factor * m[col][j];
    }
  }
  const x = new Array<number>(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    if (isFree[i])
      continue;
    let sum = m[i][n];
    for (let j = i + 1; j < n; j++)
      sum -= m[i][j] * x[j];
    x[i] = sum / m[i][i];
  }
  return x;
}
