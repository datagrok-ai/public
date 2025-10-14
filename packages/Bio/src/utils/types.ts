import * as DG from 'datagrok-api/dg';
import {pepseaMethods} from './pepsea';

export type DataFrameDict = { [key: string]: DG.DataFrame };

export namespace BarChart {
  export type BarPart = { colName: string, aaName: string };
  export type BarStatsObject = { name: string, count: number, selectedCount: number };
}

export type UTypedArray = Uint8Array | Uint16Array | Uint32Array;
//AAR: (Position: (index: indexList))
export type SubstitutionsInfo = Map<string, Map<string, Map<number, number[] | UTypedArray>>>;
export type SelectionObject = { [position: string]: string[] };

export type MultipleSequenceAlignmentUIOptions = {
  col?: DG.Column<string> | null, clustersCol?: DG.Column | null,
  pepsea?: { method?: typeof pepseaMethods[number], gapOpen?: number, gapExtend?: number },
  kalign?: { gapOpen?: number, gapExtend?: number, terminalGap?: number }
};
