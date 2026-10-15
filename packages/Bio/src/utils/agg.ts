import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type AggValueListType = (number | null)[] | Float32Array | Int32Array;
export type AggFunc = (valueList: AggValueListType) => number | null;

/** Whether a list value is present (missing: null, NaN or the column null sentinels). */
const present = (x: number | null): x is number =>
  x != null && !Number.isNaN(x) && x !== DG.FLOAT_NULL && x !== DG.INT_NULL;

/** Aggregations computed here: building a platform column for each position and monomer
 * (WebLogo makes thousands of them on long sequences) took seconds on 100,000 rows. */
const DIRECT: {[agg: string]: AggFunc} = {
  [DG.AGG.AVG]: (valueList) => {
    let sum = 0;
    let n = 0;
    for (let i = 0; i < valueList.length; i++) {
      const x = valueList[i];
      if (present(x)) {
        sum += x;
        n++;
      }
    }
    return n > 0 ? sum / n : null;
  },
  [DG.AGG.SUM]: (valueList) => {
    let sum = 0;
    for (let i = 0; i < valueList.length; i++) {
      const x = valueList[i];
      if (present(x)) sum += x;
    }
    return sum;
  },
  [DG.AGG.MIN]: (valueList) => {
    let min: number | null = null;
    for (let i = 0; i < valueList.length; i++) {
      const x = valueList[i];
      if (present(x) && (min === null || x < min)) min = x;
    }
    return min;
  },
  [DG.AGG.MAX]: (valueList) => {
    let max: number | null = null;
    for (let i = 0; i < valueList.length; i++) {
      const x = valueList[i];
      if (present(x) && (max === null || x > max)) max = x;
    }
    return max;
  },
};

export function getAgg(agg: DG.AggregationType): AggFunc {
  if (agg in DIRECT) return DIRECT[agg];

  function buildCol(valueList: AggValueListType): DG.Column<number> {
    let resCol: DG.Column<number>;
    const resColName = `agg`;
    if (valueList instanceof Float32Array)
      resCol = DG.Column.fromFloat32Array(resColName, valueList as Float32Array);
    else if (valueList instanceof Int32Array)
      resCol = DG.Column.fromInt32Array(resColName, valueList as Int32Array);
    else
      resCol = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, resColName, valueList as (number | null)[]);

    return resCol;
  }

  return (valueList: AggValueListType): number | null => {
    const aggCol = buildCol(valueList);
    const res = aggCol.aggregate(agg);
    return res;
  };
}
