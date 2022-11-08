import * as DG from 'datagrok-api/dg';

export type DataFrameDict = {[key: string]: DG.DataFrame};

export type UTypedArray = Uint8Array | Uint16Array | Uint32Array;
//AAR: (Position: (index: indexList))
export type SubstitutionsInfo = Map<string, Map<string, Map<number, number[] | UTypedArray>>>;
export type PositionToAARList = {[postiton: string]: string[]};

export type MonomerColStats = {[monomer: string]: {count: number, selected: number}};
export type MonomerDfStats = {[position: string]: MonomerColStats};

export type ScalingMethods = 'none' | 'lg' | '-lg';
export type PeptidesSettings =
  {scaling?: ScalingMethods, isBidirectional?: boolean, maxMutations?: number, minActivityDelta?: number};

export type DrawOptions = {
  fontStyle?: string,
  upperLetterHeight?: number,
  upperLetterAscent?: number,
  bounds?: DG.Rect,
  textAlign?: CanvasTextAlign,
  textBaseline?: CanvasTextBaseline,
  marginVertical?: number,
  marginHorizontal?: number,
};

export type StatsInfo = {
  monomerCol: DG.Column<string>,
  countCol: DG.Column<number>,
  orderedIndexes: Int32Array,
}