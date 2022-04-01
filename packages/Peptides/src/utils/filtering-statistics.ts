import * as DG from 'datagrok-api/dg';

import {tTest} from '@datagrok-libraries/statistics/src/tests';

/** Column statistics helper. */
export class FilteringStatistics {
  private data?: Float32Array;
  private stats: Stats;

  /**
   * Creates an instance of FilteringStatistics.
   * @param {Float32Array} [data] Numeric values to consider.
   */
  constructor(data?: Float32Array) {
    this.data = data;
    this.stats = {
      count: 0,
      pValue: 1.,
      meanDifference: 0.,
    };
  }

  /**
   * Sets values to make statistical analysis.
   * @param {Float32Array} data Those values.
   */
  setData(data: Float32Array) {
    this.data = data;
  }

  /**
   * Sets bit mask to split population into two groups.
   * @param {DG.BitSet} mask The mask to perform splitting.
   */
  setMask(mask: DG.BitSet) {
    if (this.data) {
      const selected = this.data.filter((_, i) => mask.get(i));
      const rest = this.data.filter((_, i) => !mask.get(i));
      this.stats = this.calcStats(selected, rest);
    }
  }

  /**
   * Calculates simple statistics on two samples.
   * @param {Float32Array} selected First sample.
   * @param {Float32Array} rest Second sample.
   * @return {Stats} Statistics.
   */
  calcStats(selected: Float32Array, rest: Float32Array): Stats {
    const testResult = tTest(selected, rest);
    const currentMeanDiff = testResult['Mean difference']!;
    return {
      count: selected.length,
      pValue: testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'],
      meanDifference: currentMeanDiff,
    };
  }

  /** Returns calculated statistics. */
  get result(): Stats {
    return this.stats;
  }
}

export interface Stats {
  count: number,
  pValue: number,
  meanDifference: number,
}
