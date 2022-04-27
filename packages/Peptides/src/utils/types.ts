import * as DG from 'datagrok-api/dg';

export type SubstitutionsTooltipData = {[aar: string]: number[][][]};
export type DataFrameDict = {[key: string]: DG.DataFrame};

export namespace BarChart {
  export type BarPart = {colName : string, aaName : string};
  export type BarStatsObject = {name: string, count: number, selectedCount: number};
}
