import * as DG from 'datagrok-api/dg';


export namespace ts {
  /** A type guard function.
   * See https://stackoverflow.com/questions/64616994/typescript-type-narrowing-not-working-for-in-when-key-is-stored-in-a-variable*/
  export function hasProp<T>(obj: T, key: PropertyKey): key is keyof T {
    return key in obj;
  }
}

export namespace data {
  export function mapToRange(x: number, min1: number, max1: number, min2: number, max2: number): number {
    const range1 = max1 - min1;
    const range2 = max2 - min2;
    return (((x - min1) * range2) / range1) + min2;
  }

  export function aggToStat(dataframe: DG.DataFrame, columnName: string,
    aggregation: DG.AggregationType): number | null {
    const colStatsCall = 'dataframe.getCol(columnName).stats.';
    const stats = {
      avg: colStatsCall + 'avg',
      count: colStatsCall + 'totalCount',
      kurt: colStatsCall + 'kurt',
      max: colStatsCall + 'max',
      med: colStatsCall + 'med',
      min: colStatsCall + 'min',
      nulls: colStatsCall + 'missingValueCount',
      q1: colStatsCall + 'q1',
      q2: colStatsCall + 'q2',
      q3: colStatsCall + 'q3',
      skew: colStatsCall + 'skew',
      stdev: colStatsCall + 'stdev',
      sum: colStatsCall + 'sum',
      unique: colStatsCall + 'uniqueCount',
      values: colStatsCall + 'valuesCount',
      variance: colStatsCall + 'variance',
      '#selected': 'dataframe.selection.trueCount',
      first: 'dataframe.getCol(columnName).get(0)',
    };
    return ts.hasProp(stats, aggregation) ? eval(stats[aggregation]) : null;
  }
}
