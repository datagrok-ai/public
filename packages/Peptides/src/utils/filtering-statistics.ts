import * as DG from 'datagrok-api/dg';

import {tTest} from '@datagrok-libraries/statistics/src/tests';

export type Stats = {
  count: number,
  pValue: number,
  meanDifference: number,
  ratio: number,
};

type StatsData = Float32Array | Float64Array | Int32Array | Uint32Array | number[];

export function getStats(data: StatsData, mask: DG.BitSet): Stats {
  const selected = new Float32Array(mask.trueCount);
  const rest = new Float32Array(mask.falseCount);
  let selectedIndex = 0;
  let restIndex = 0;
  data.forEach((v, i) => mask.get(i) ? selected[selectedIndex++] = v : rest[restIndex++] = v);

  const testResult = tTest(selected, rest);
  const currentMeanDiff = testResult['Mean difference']!;
  const realCount = selected.length || data.length;
  return {
    count: realCount,
    pValue: testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'],
    meanDifference: currentMeanDiff,
    ratio: realCount / data.length,
  };
}
