import * as DG from 'datagrok-api/dg';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {RawData} from './types';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {ClusterType} from '../viewers/logo-summary';

export type StatsItem = {
  count: number,
  pValue: number | null,
  meanDifference: number,
  ratio: number,
  mask: BitArray,
  mean: number,
};

export type PositionStats = { [monomer: string]: StatsItem } & { general: SummaryStats };
export type MonomerPositionStats = { [position: string]: PositionStats } & { general: SummaryStats };
export type MonomerPositionStatsCache = {[colName: string]: {filter: Int32Array, stats: MonomerPositionStats, selection: Int32Array}};
export type ClusterStats = { [cluster: string]: StatsItem };
export type ClusterTypeStats = { [clusterType in ClusterType]: ClusterStats };
export type MasksInfo = {
  [monomerOrCluster: string]: { [positionOrClusterType: string | ClusterType]: { mask: BitArray } }
};

export type SummaryStats = {
  minCount: number, maxCount: number,
  minMeanDifference: number, maxMeanDifference: number,
  minPValue: number, maxPValue: number,
  minRatio: number, maxRatio: number,
  minMean: number, maxMean: number,
};

export const getAggregatedColName = (aggF: string, colName: string): string => `${aggF}(${colName})`;
export type AggregationColumns = { [col: string]: DG.AggregationType };

export function getStats(data: RawData | number[], bitArray: BitArray): StatsItem {
  if (data.length !== bitArray.length && data.some((v, i) => i >= bitArray.length ? v !== 0 : false))
    throw new Error('PeptidesError: Data and bit array have different lengths');
  if (bitArray.falseCount() === 0 || bitArray.trueCount() === 0)
    throw new Error('PeptidesError: One of the samples is empty');

  const selected = new Float32Array(bitArray.trueCount());
  const rest = new Float32Array(bitArray.falseCount());

  let selectedIndex = 0;
  let restIndex = 0;
  for (let i = 0; i < bitArray.length; ++i) {
    if (bitArray.getBit(i))
      selected[selectedIndex++] = data[i];
    else
      rest[restIndex++] = data[i];
  }

  const selectedMean = selected.reduce((a, b) => a + b, 0) / selected.length;
  if (selected.length === 1 || rest.length === 1) {
    const restMean = rest.reduce((a, b) => a + b, 0) / rest.length;
    return {
      count: selected.length,
      pValue: null,
      mean: selectedMean,
      meanDifference: selectedMean - restMean,
      ratio: selected.length / (bitArray.length),
      mask: bitArray,
    };
  }

  const testResult = tTest(selected, rest);
  const currentMeanDiff = testResult['Mean difference']!;
  return {
    count: selected.length,
    pValue: testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'],
    mean: selectedMean,
    meanDifference: currentMeanDiff,
    ratio: selected.length / (bitArray.length),
    mask: bitArray,
  };
}

export function getAggregatedValue(col: DG.Column<number>, agg: DG.AggregationType, mask?: DG.BitSet): number {
  const stat = DG.Stats.fromColumn(col, mask);
  if (!(agg in stat))
    throw new Error(`Aggregation type ${agg} is not supported`);
  //@ts-ignore: this is a hack to avoid using switch to access the getters
  return stat[agg] as number;
}

export function getAggregatedColumnValues(df: DG.DataFrame, columns: AggregationColumns,
  options: {
                                            filterDf?: boolean,
                                            mask?: DG.BitSet,
                                            fractionDigits?: number
                                          } = {}): StringDictionary {
  options.filterDf ??= false;
  options.fractionDigits ??= 3;

  const filteredDf = options.filterDf && df.filter.anyFalse ? df.clone(df.filter) : df;

  const colResults: StringDictionary = {};
  const aggrColsEntries = Object.entries(columns);
  aggrColsEntries.forEach((it) => { it[0] = it[0].substring(0, it[0].indexOf("%"))});
  for (const [colName, aggFn] of aggrColsEntries) {
    const newColName = getAggregatedColName(aggFn, colName);
    const value = getAggregatedValue(filteredDf.getCol(colName), aggFn, options.mask);
    colResults[newColName] = value.toFixed(options.fractionDigits);
  }
  return colResults;
}
