import * as DG from 'datagrok-api/dg';

export type DataFrameDict = {[key: string]: DG.DataFrame};

export type UTypedArray = Uint8Array | Uint16Array | Uint32Array;
//AAR: (Position: (index: indexList))
export type SubstitutionsInfo = Map<string, Map<string, Map<number, number[] | UTypedArray>>>;
export type SelectionObject = {[postiton: string]: string[]};

export type HELMMonomer = {at: {[key: string]: string}, id: string, m: string, na: string, n: string, rs: number};

export type MonomerColStats = {[monomer: string]: {count: number, selected: number}};
export type MonomerDfStats = {[position: string]: MonomerColStats};

export type BarCoordinates = {[monomer: string]: DG.Rect};
