import * as DG from 'datagrok-api/dg';
import {SCALING_METHODS} from './constants';
import {ClusterType} from '../model';

export type DataFrameDict = {[key: string]: DG.DataFrame};

export type RawData = Int32Array | Uint32Array | Float32Array | Float64Array;
export type UTypedArray = Uint8Array | Uint16Array | Uint32Array;
//AAR: (Position: (index: indexList))
export type MutationCliffs = Map<string, Map<string, Map<number, number[] | UTypedArray>>>;
export type Selection = {[positionOrClusterType: string | ClusterType]: string[]};
export type SelectionItem = {positionOrClusterType: string | ClusterType, monomerOrCluster: string};
export type SelectionStats = {[positionOrClusterType: string | ClusterType]: { [monomerOrCluster: string]: number}};

export type PeptidesSettings = {
  sequenceColumnName?: string,
  activityColumnName?: string,
  clustersColumnName?: string,
  targetColumnName?: string,
  scaling?: SCALING_METHODS,
  isBidirectional?: boolean,
  showMonomerPosition?: boolean,
  showMostPotentResidues?: boolean,
  showLogoSummaryTable?: boolean,
  showDendrogram?: boolean,
  maxMutations?: number,
  minActivityDelta?: number,
  columns?: {[col: string]: DG.AggregationType},
};

export type DrawOptions = {
  symbolStyle?: string,
  upperLetterHeight?: number,
  upperLetterAscent?: number,
  bounds?: DG.Rect,
  textAlign?: CanvasTextAlign,
  textBaseline?: CanvasTextBaseline,
  marginVertical?: number,
  marginHorizontal?: number,
  headerStyle?: string,
  textHeight?: number,
};

export type StatsInfo = {
  monomerCol: DG.Column<string>,
  countCol: DG.Column<number>,
  orderedIndexes: Int32Array,
}

export type RawColumn = {name: string, rawData: RawData, cat?: string[]};

export type SelectionOptions = {shiftPressed: boolean, ctrlPressed: boolean};
