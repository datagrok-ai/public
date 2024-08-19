import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type AggValueListType = (number | null)[] | Float32Array | Int32Array;
export type AggFunc = (valueList: AggValueListType) => number | null;

export function getAgg(agg: DG.AggregationType): AggFunc {
  let res: AggFunc;

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
