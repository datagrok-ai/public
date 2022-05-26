import * as DG from 'datagrok-api/dg';

export type SubstitutionCases = {[aar: string]: number[][][]};
export type SubstitutionTooltips = { [aar: string]: {}[][]; };
export type DataFrameDict = {[key: string]: DG.DataFrame};

export namespace BarChart {
  export type BarPart = {colName : string, aaName : string};
  export type BarStatsObject = {name: string, count: number, selectedCount: number};
}

type UTypedArray = Uint8Array | Uint16Array | Uint32Array;
//AAR: (Position: (index: indexList))
export type SubstitutionsInfo = Map<string, Map<string, Map<number, number[] | UTypedArray>>>;
