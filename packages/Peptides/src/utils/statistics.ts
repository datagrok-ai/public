import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {RawData} from './types';

export type Stats = {
  count: number,
  pValue: number,
  meanDifference: number,
  ratio: number,
};

export type MaskInfo = {
  trueCount: number,
  falseCount: number,
  mask: boolean[] | Int32Array,
};

export function getStats(data: RawData | number[], maskInfo: MaskInfo): Stats {
  const selected = new Float32Array(maskInfo.trueCount);
  const rest = new Float32Array(maskInfo.falseCount);

  let selectedIndex = 0;
  let restIndex = 0;
  for (let i = 0; i < data.length; ++i) {
    if (maskInfo.mask[i])
      selected[selectedIndex++] = data[i];
    else
      rest[restIndex++] = data[i];
  }

  const testResult = tTest(selected, rest);
  const currentMeanDiff = testResult['Mean difference']!;
  return {
    count: selected.length,
    pValue: testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'] || 0,
    meanDifference: currentMeanDiff || 0,
    ratio: selected.length / data.length,
  };
}
