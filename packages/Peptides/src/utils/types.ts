import * as DG from 'datagrok-api/dg';
import {ClusterType} from '../viewers/logo-summary';
import {SCALING_METHODS} from './constants';
import {AggregationColumns, StatsItem} from './statistics';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

export type RawData = Int32Array | Uint32Array | Float32Array | Float64Array;
type MONOMER = string;
type POSITION = string;
type INDEX = number;
type INDEXES = number[] | UTypedArray;
export type UTypedArray = Uint8Array | Uint16Array | Uint32Array;
//Monomer: (Position: (index: indexList))
export type MutationCliffs = Map<MONOMER, Map<POSITION, Map<INDEX, INDEXES>>>;
//Monomer: (Position: (Stats))
export type MutationCliffStats = {
  stats: Map<MONOMER, Map<POSITION, StatsItem>>;
  minDiff: number; maxDiff: number; maxCount: number; minCount: number;
};
export type Selection = { [positionOrClusterType: string | ClusterType]: string[] };
export type SelectionItem = { positionOrClusterType: string | ClusterType, monomerOrCluster: string };
export type SelectionStats = { [positionOrClusterType: string | ClusterType]: { [monomerOrCluster: string]: number } };

export interface PeptidesSettings {
  sequenceColumnName: string,
  activityColumnName: string,
  activityScaling: SCALING_METHODS,
  showMonomerPosition?: boolean,
  showMostPotentResidues?: boolean,
  showLogoSummaryTable?: boolean,
  showDendrogram?: boolean,
  showClusterMaxActivity?: boolean,
  showSequenceSpace?: boolean,
  columns?: AggregationColumns,
  sequenceSpaceParams: SequenceSpaceParams,
  mclSettings: MCLSettings,
}

export class MCLSettings {
  maxIterations: number = 5;
  threshold: number = 80;
  distanceF: MmDistanceFunctionsNames = MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH;
  gapOpen: number = 1;
  gapExtend: number = 0.6;
  fingerprintType: string = 'Morgan';
}

export class SequenceSpaceParams {
  distanceF: MmDistanceFunctionsNames = MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH;
  gapOpen: number = 1;
  gapExtend: number = 0.6;
  clusterEmbeddings: boolean = true;
  epsilon: number = 0.01;
  minPts: number = 4;
  fingerprintType: string = 'Morgan';
  constructor(clusterEmbeddings?: boolean) {
    this.clusterEmbeddings = !!clusterEmbeddings;
  }
}

export type PartialPeptidesSettings = Partial<PeptidesSettings>;

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
  selectionWidth?: number,
};

export type RawColumn = { name: string, rawData: RawData, cat?: string[] };

export type SelectionOptions = { shiftPressed: boolean, ctrlPressed: boolean, notify?: boolean };

export type CachedWebLogoTooltip = { bar: string, tooltip: HTMLDivElement | null };
