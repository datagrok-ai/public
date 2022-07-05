import * as DG from 'datagrok-api/dg';

import {tTest} from '@datagrok-libraries/statistics/src/tests';

export type Stats = {
  count: number,
  pValue: number,
  meanDifference: number,
  ratio: number,
};

export function getStats(data: Float32Array | number[], mask: DG.BitSet): Stats {
  const selected = new Float32Array(mask.trueCount);
  const rest = new Float32Array(mask.falseCount);
  let selectedIndex = 0;
  let restIndex = 0;
  data.forEach((v, i) => mask.get(i) ? selected[selectedIndex++] = v : rest[restIndex++] = v);

  const testResult = tTest(selected, rest);
  const currentMeanDiff = testResult['Mean difference']!;
  return {
    count: selected.length,
    pValue: testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'],
    meanDifference: currentMeanDiff,
    ratio: selected.length / data.length,
  };
}
