import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';

/** @type {*} A dictionary of basic binary operations. */
const _operations: {[name: string]: Function} = {
  '+': (a: number, b: number) => (a + b),
  '-': (a: number, b: number) => (a - b),
  '*': (a: number, b: number) => (a * b),
  '/': (a: number, b: number) => (a / b),
  'min': (a: number, b: number) => (a < b ? a : b),
  'max': (a: number, b: number) => (a > b ? a : b),
};

/**
 * Returns the indices that would sort an array.
 *
 * @param {Float32Array} values Array to sort.
 * @return {Int32Array} Array of indices that sort values along the first axis.
 */
function _argsort(values: Float32Array): Int32Array {
  const array = Array.from(values);
  return Int32Array.from(ArrayUtils.argSort(array));
}

/**
 * Take elements from an array.
 *
 * @param {Int32Array} order The indices of the values to extract.
 * @param {Float32Array} values The source array.
 * @return {Float32Array} The returned array has the same type as values.
 */
function _take(order: Int32Array, values: Float32Array): Float32Array {
  // TODO: Implement a general function for TypedArray.
  return Float32Array.from(values).map((_, i) => values[order[i]]);
}

/**
 * Assign elements of an array following the order given (floating-point version).
 *
 * @param {Float32Array} values The source array.
 * @param {Int32Array} order The order given.
 * @return {Float32Array} The returned array has the same type as values.
 */
function _give(values: Float32Array, order: Int32Array): Float32Array {
  const v = Float32Array.from(values);

  for (let i = 0; i < order.length; ++i) {
    v[order[i]] = values[i];
  }
  return v;
}

/**
 * Assign elements of an array following the order given (boolean version).
 *
 * @param {Array<boolean>} values The source array.
 * @param {Int32Array} order The order given.
 * @return {Array<boolean>} The returned array has the same type as values.
 */
function _giveb(values: Array<boolean>, order: Int32Array): Array<boolean> {
  const v = Array.from(values);

  for (let i = 0; i < order.length; ++i) {
    v[order[i]] = values[i];
  }
  return v;
}

/**
 * No frills empirical cdf used in fdrcorrection.
 *
 * @param {Float32Array} x The source array to take a dimension from.
 * @return {Float32Array} Empirical cdf.
 */
function _ecdf(x: Float32Array): Float32Array {
  const nobs = x.length;
  return Float32Array.from(x).map((_, i) => (i+1)/nobs);
}

/**
 * cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))
 *
 * @param {number} n The number given.
 * @return {number} cm value.
 */
function _cm(n: number): number {
  let sum = 0;
  for (let i = 0; i < n; ++i) {
    sum += 1 / (i+1);
  }
  return sum;
}

/**
 * Basic operation under a vector and a scalar.
 *
 * @param {Float32Array} values The source vector.
 * @param {number} scale The scalar.
 * @param {string} [op='*'] The operation to perform.
 * @return {Float32Array} New vector as a result of the operation.
 */
function _factor(values: Float32Array, scale: number, op = '*'): Float32Array {
  return Float32Array.from(values).map((v, _) => _operations[op](v, scale));
}

/**
 * Basic operation under two vectors.
 *
 * @param {Float32Array} values The first vector.
 * @param {Float32Array} scale The second vector.
 * @param {string} [op='*'] The operation to perform.
 * @return {Float32Array} New vector as a result of the operation.
 */
function _vfactor(values: Float32Array, scale: Float32Array, op = '*'): Float32Array {
  return Float32Array.from(values).map((v, i) => _operations[op](v, scale[i]));
}

/**
 * Accumulate the result of applying the min operator to all elements.
 *
 * @param {Float32Array} values The array to act on.
 * @return {Float32Array} The accumulated values.
 */
function _minimumAccumulate(values: Float32Array): Float32Array {
  const nItems = values.length;
  const r = Float32Array.from(values);

  for (let i = 0; i < nItems; ++i) {
    r[i] = values.slice(0, i+1).reduce((a, b, _, __) => (_operations['min'](a, b)));
  }
  return r;
}

/**
 * pvalue correction for false discovery rate
 *
 * @export
 * @param {Float32Array} pvals Set of p-values of the individual tests.
 * @param {number} [alpha=0.05] Family-wise error rate. Defaults to 0.05.
 * @param {string} [method='n'] {'i', 'indep', 'p', 'poscorr', 'n', 'negcorr'}, optional
 * Which method to use for FDR correction.
 * ``{'i', 'indep', 'p', 'poscorr'}`` all refer to ``fdr_bh``
 * (Benjamini/Hochberg for independent or positively
 * correlated tests). ``{'n', 'negcorr'}`` both refer to ``fdr_by``
 * (Benjamini/Yekutieli for general or negatively correlated tests).
 * Defaults to 'n'.
 * @param {boolean} [isSorted=false] If False (default), the p_values will be sorted, but the corrected
 * pvalues are in the original order. If True, then it assumed that the
 * pvalues are already sorted in ascending order.
 * @return {[Array<boolean>, Float32Array]} rejected : ndarray, bool
 * True if a hypothesis is rejected, False if not
 * pvalue-corrected : ndarray
 * pvalues adjusted for multiple hypothesis testing to limit FDR
 * @see
 * If there is prior information on the fraction of true hypothesis, then alpha
 * should be set to ``alpha * m/m_0`` where m is the number of tests,
 * given by the p-values, and m_0 is an estimate of the true hypothesis.
 * (see Benjamini, Krieger and Yekuteli)
 *
 * The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
 * of false hypotheses will be available (soon).
 *
 * Both methods exposed via this function (Benjamini/Hochberg, Benjamini/Yekutieli)
 * are also available in the function ``multipletests``, as ``method="fdr_bh"`` and
 * ``method="fdr_by"``, respectively.
 */
export function fdrcorrection(
  pvals: Float32Array,
  alpha: number = 0.05,
  method: string ='n',
  isSorted: boolean = false): [Array<boolean>, Float32Array]
// eslint-disable-next-line brace-style
{
  const nItems = pvals.length;
  let pvalsSorted: Float32Array;
  let pvalsSortind: Int32Array;
  let cm = 0;

  if (!isSorted) {
    pvalsSortind = _argsort(pvals);
    pvalsSorted = _take(pvalsSortind, pvals);
  } else {
    pvalsSortind = new Int32Array(nItems).fill(0).map((_, i) => (i));
    pvalsSorted = pvals; // alias
  }

  let ecdffactor = _ecdf(pvalsSorted);

  if (['i', 'indep', 'p', 'poscorr'].includes(method)) {
    ;
  } else if (['n', 'negcorr'].includes(method)) {
    cm = _cm(nItems);
    ecdffactor = _factor(ecdffactor, cm, '/');
  } else {
    throw new Error('only indep and negcorr implemented');
  }

  const reject: boolean[] = new Array(nItems).fill(false);
  let rejectmax = -1;

  for (let i = 0; i < nItems; ++i) {
    if (pvalsSorted[i] <= ecdffactor[i]*alpha) {
      rejectmax = i;
    }
  }

  if (rejectmax >= 0) {
    for (let i = 0; i < rejectmax; ++i) {
      reject[i] = true;
    }
  }

  let pvalsCorrected = _vfactor(pvalsSorted, ecdffactor, '/');
  pvalsCorrected = _minimumAccumulate(pvalsCorrected.reverse()).reverse();

  for (let i = 0; i < nItems; ++i) {
    if (pvalsCorrected[i] > 1) {
      pvalsCorrected[i] = 1;
    }
  }

  if (!isSorted) {
    const pvalsCorrected_ = _give(pvalsCorrected, pvalsSortind);
    const reject_ = _giveb(reject, pvalsSortind);
    return [reject_, pvalsCorrected_];
  }
  return [reject, pvalsCorrected];
}
