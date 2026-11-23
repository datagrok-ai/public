import * as DG from 'datagrok-api/dg';

export type ValidationScheme = 'kfold' | 'holdout';

/** How a training measures Validation: k-fold cross-validation on [folds] folds, or one hold out of [testFraction] of
 * the rows; [stratified] by class (a regression target is never stratified); [shuffle] by the split seed, else in row
 * order. */
export interface ValidationSettings {
  scheme: ValidationScheme;
  folds: number;
  testFraction: number;
  stratified: boolean;
  shuffle: boolean;
}

export const FOLDS = 5;
export const FOLDS_RANGE = [2, 20];
export const TEST_FRACTION = 0.2;
export const TEST_FRACTION_RANGE = [0.05, 0.5];
export const DEFAULT_SEED = 42;
export const DEFAULT_VALIDATION: ValidationSettings = {scheme: 'kfold', folds: FOLDS, testFraction: TEST_FRACTION,
  stratified: true, shuffle: true};

/** The folds a training predicts: k for k-fold, the one held out part for a hold out. */
export function validationFolds(validation: ValidationSettings): number {
  return validation.scheme === 'kfold' ? validation.folds : 1;
}

/** The rows a training needs: two per fold, or two held out rows. */
export function minimumRows(validation: ValidationSettings): number {
  return validation.scheme === 'kfold' ? 2 * validation.folds : Math.ceil(2 / validation.testFraction);
}

/** The rows a hold out of [testFraction] takes of [rowCount]: at least one, and one row left to fit. */
export function heldOutCount(rowCount: number, testFraction: number): number {
  return Math.min(rowCount - 1, Math.max(1, Math.round(testFraction * rowCount)));
}

/** Every row's fold for [validation] ({@link kFold}, {@link holdOutFold}): the folds 0..k-1, or 0 for the held out rows
 * and -1 for the others. */
export function splitRows(rowCount: number, validation: ValidationSettings, seed: number, stratify?: DG.Column):
  Int32Array {
  const {scheme, folds, testFraction, shuffle} = validation;
  return scheme === 'kfold' ? kFold(rowCount, folds, seed, stratify, shuffle) :
    holdOutFold(rowCount, testFraction, seed, stratify, shuffle);
}

/** Every row's fold: the rows shuffled by [seed] and dealt round-robin, so fold sizes differ by at most one; with
 * [stratify] (a classification target), each class is dealt in turn, so every fold holds each class's share (a class's
 * counts differ by at most one across the folds). Without [shuffle], the rows keep their order and [seed] is unused:
 * contiguous folds, or each class dealt round-robin in row order. */
export function kFold(rowCount: number, folds: number, seed: number, stratify?: DG.Column, shuffle = true):
  Int32Array {
  const order = rowsOf(rowCount, seed, shuffle);
  const fold = new Int32Array(rowCount);
  if (stratify === undefined) {
    for (let j = 0; j < rowCount; j++)
      fold[order[j]] = shuffle ? j % folds : Math.floor(j * folds / rowCount);
    return fold;
  }
  let k = 0;
  for (const rows of classGroups(order, stratify)) {
    for (const r of rows)
      fold[r] = k++ % folds;
  }
  return fold;
}

/** 0 for the {@link heldOutCount} rows held out, -1 for the others: the last ones of the rows shuffled by [seed] (in
 * row order without [shuffle]); with [stratify], each class's share of them from its own rows. */
export function holdOutFold(rowCount: number, testFraction: number, seed: number, stratify?: DG.Column,
  shuffle = true): Int32Array {
  const order = rowsOf(rowCount, seed, shuffle);
  const groups = stratify === undefined ? [[...order]] : classGroups(order, stratify);
  const quotas = groupQuotas(groups.map((g) => g.length), heldOutCount(rowCount, testFraction));
  const fold = new Int32Array(rowCount).fill(-1);
  for (const [i, group] of groups.entries()) {
    for (const r of group.slice(group.length - quotas[i]))
      fold[r] = 0;
  }
  return fold;
}

/** Each group's share of [total] by its size (the largest remainders round up, so the shares sum to [total]). */
export function groupQuotas(sizes: number[], total: number): number[] {
  const rows = sizes.reduce((sum, s) => sum + s, 0);
  const exact = sizes.map((s) => s * total / rows);
  const quotas = exact.map(Math.floor);
  const byRemainder = sizes.map((_, i) => i).sort((a, b) => (exact[b] - quotas[b]) - (exact[a] - quotas[a]) || a - b);
  const short = total - quotas.reduce((sum, q) => sum + q, 0);
  for (let k = 0; k < short; k++)
    quotas[byRemainder[k]]++;
  return quotas;
}

function rowsOf(rowCount: number, seed: number, shuffle: boolean): Int32Array {
  return shuffle ? shuffledRows(rowCount, seed) : Int32Array.from({length: rowCount}, (_, i) => i);
}

/** The rows 0..[rowCount] - 1 shuffled by [seed] (Fisher-Yates on {@link mulberry32}). */
export function shuffledRows(rowCount: number, seed: number): Int32Array {
  const random = mulberry32(seed);
  const order = new Int32Array(rowCount);
  for (let i = 0; i < rowCount; i++)
    order[i] = i;
  for (let i = rowCount - 1; i > 0; i--) {
    const j = Math.floor(random() * (i + 1));
    [order[i], order[j]] = [order[j], order[i]];
  }
  return order;
}

/** The rows of [order] per class of [target], each class's rows in that order, the classes by their first row. */
export function classGroups(order: Int32Array, target: DG.Column): number[][] {
  const groups = new Map<string, number[]>();
  for (const r of order) {
    const key = String(target.get(r));
    const rows = groups.get(key);
    if (rows === undefined)
      groups.set(key, [r]);
    else
      rows.push(r);
  }
  return [...groups.values()];
}

export function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
