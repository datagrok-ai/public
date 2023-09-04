import * as DG from 'datagrok-api/dg';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {RawData} from './types';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export type Stats = {
  count: number,
  pValue: number | null,
  meanDifference: number,
  ratio: number,
};

export function getStats(data: RawData | number[], bitArray: BitArray): Stats {
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

  if (selected.length === 1 || rest.length === 1) {
    const selectedMean = selected.reduce((a, b) => a + b, 0) / selected.length;
    const restMean = rest.reduce((a, b) => a + b, 0) / rest.length;
    return {
      count: selected.length,
      pValue: null,
      meanDifference: selectedMean - restMean,
      ratio: selected.length / (bitArray.length),
    };
  }

  const testResult = tTest(selected, rest);
  const currentMeanDiff = testResult['Mean difference']!;
  return {
    count: selected.length,
    pValue: testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'],
    meanDifference: currentMeanDiff,
    ratio: selected.length / (bitArray.length),
  };
}

export function getAggregatedValue(col: DG.Column<number>, agg: DG.AggregationType, mask?: DG.BitSet): number {
  const stat = DG.Stats.fromColumn(col, mask);
  if (!(agg in stat))
    throw new Error(`Aggregation type ${agg} is not supported`);
  //@ts-ignore: this is a hack to avoid using switch to access the getters
  return stat[agg] as number;
}
