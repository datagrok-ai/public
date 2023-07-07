/* eslint-disable no-tabs */
export type BoxPlotStatistics = {
  q1: number,
	q2: number,
	q3: number,
	lowerAdjacentValue: number,
	upperAdjacentValue: number,
};


export function calculateBoxPlotStatistics(values: number[]): BoxPlotStatistics {
  values.sort();

  const q1 = values[Math.trunc(values.length / 4)];
  const q2 = values.length % 2 === 0 ?
    (values[Math.trunc(values.length / 2)] + values[Math.trunc(values.length / 2) - 1]) / 2 :
    values[Math.trunc(values.length / 2)];
  const q3 = values[Math.trunc(values.length * 3 / 4)];

  const hSpread = q3 - q1;
  const step = hSpread * 1.5;

  let upperIdx = arrayBinarySearch(values, q3 + step);
  if (upperIdx < 0)
    upperIdx = ~upperIdx;
  let upperAdjacentValue = values[Math.max(0, upperIdx - 1)];
  if (upperAdjacentValue < q3)
    upperAdjacentValue = q3;

  let lowerIdx = arrayBinarySearch(values, q1 - step);
  if (lowerIdx < 0)
    lowerIdx = ~lowerIdx;
  let lowerAdjacentValue = values[lowerIdx];
  if (lowerAdjacentValue > q1)
    lowerAdjacentValue = q1;

  return {
    q1: q1,
    q2: q2,
    q3: q3,
    lowerAdjacentValue: lowerAdjacentValue,
    upperAdjacentValue: upperAdjacentValue,
  };
}

function arrayBinarySearch(items: number[], value: number, left: number = 0, right?: number): number {
  if (items.length === 0)
    return -1; //~0;

  if (right === undefined)
    right = items.length - 1;

  if (value < items[left])
    return -1; //~left;
  if (value > items[right])
    return -(right + 1); //~(right+1);

  while (right! - left > 1) {
    const mid = Math.trunc((right + left) / 2);
    if (value === items[mid])
      return mid;
    if (value < items[mid])
      right = mid;
    else
      left = mid;
  }

  if (items[left] === value)
    return left;
  if (items[right] === 0)
    return right;

  // no match - returning inverted index
  if (value < items[left])
    return -(left + 1); //~left;
  return -(right + 1); //~right;
}
