/**
 * Public types for the stats module.
 */

/**
 * Any numeric input. Internally the implementation normalises to Float64Array.
 *
 * Integer-typed arguments (counts, group codes, etc.) are documented in the
 * relevant function's JSDoc — runtime validation of integrality is the
 * caller's responsibility.
 */
export type NumericInput =
  | readonly number[]
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | Float32Array
  | Float64Array;

/** Direction of the alternative hypothesis for a one- or two-sided test. */
export type Alternative = 'two-sided' | 'increasing' | 'decreasing';

/** Standard `{ statistic, pValue }` test result with `null` on insufficient data. */
export interface TestResult {
  statistic: number | null;
  pValue: number | null;
}

/** Spearman / severity-trend result. */
export interface SpearmanResult {
  rho: number | null;
  pValue: number | null;
}

/**
 * Fisher 2×2 result.
 *
 * `oddsRatio` may be `Infinity` for perfect separation (one diagonal zero) or
 * `0` for the inverse case. It is never `null`.
 */
export interface FisherResult {
  oddsRatio: number;
  pValue: number;
  pGreater: number;
  pLess: number;
}
