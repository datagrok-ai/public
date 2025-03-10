import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {BiostructureData, BiostructureDataJson} from '../pdb/types';
import {IViewer} from './viewer';
import {PdbResDataFrameType} from '../pdb/pdb-helper';
import {TAGS as pdbTAGS} from '../pdb/index';

export enum RepresentationType {
  Cartoon = 'cartoon',
  Backbone = 'backbone',
  BallAndStick = 'ball+stick',
  Licorice = 'licorice',
  Hyperball = 'hyperball',
  Surface = 'surface'
}

export enum RegionStateOptionsType {
  FULL = 'full',
  COLLAPSED = 'collapsed',
  HIDDEN = 'hidden'
}

export enum SimpleRegionStateOptionsType {
  FULL = 'full',
  HIDDEN = 'hidden'
}

export enum PluginLayoutControlsDisplayType {
  OUTSIDE = 'outside',
  PORTRAIT = 'portrait',
  LANDSCAPE = 'landscape',
  REACTIVE = 'reactive'
}

export const BiostructurePropsDefault = new class {
  // -- Data --
  dataJson: string = BiostructureDataJson.empty;
  pdb: string | null = null;
  pdbTag: string | null = pdbTAGS.PDB;

  biostructureIdColumnName: string | null = null;
  /** DG.Func nqName */
  biostructureDataProvider: string | null = null;

  ligandColumnName: string | null = null;
  // pdbProvider: string = 'rcsb';
  // emdbProvider: string = 'rcsb';

  // -- Style --
  representation: RepresentationType = RepresentationType.Cartoon;

  // -- Layout --
  layoutIsExpanded: boolean = false;
  layoutShowControls: boolean = false;
  layoutRegionStateLeft: RegionStateOptionsType = RegionStateOptionsType.FULL;
  layoutRegionStateTop: SimpleRegionStateOptionsType = SimpleRegionStateOptionsType.FULL;
  layoutRegionStateRight: SimpleRegionStateOptionsType = SimpleRegionStateOptionsType.FULL;
  layoutRegionStateBottom: SimpleRegionStateOptionsType = SimpleRegionStateOptionsType.FULL;
  layoutControlsDisplay: PluginLayoutControlsDisplayType = PluginLayoutControlsDisplayType.OUTSIDE;

  showImportControls: boolean = false;
  showSessionControls: boolean = false;
  showStructureSourceControls: boolean = true;
  // showMeasurementsControls: true,
  // showStrucmotifSubmitControls: true,
  showSuperpositionControls: boolean = true;

  // showQuickStylesControls: false,
  // showStructureComponentControls: true,
  // showVolumeStreamingControls: true,
  // showAssemblySymmetryControls: true,
  showValidationReportControls: boolean = true;

  showMembraneOrientationPreset: boolean = false;

  // showNakbColorTheme: false,
  /**
   * Needed when running outside of sierra. If set to true, the strucmotif UI will use an absolute URL to sierra-prod.
   * Otherwise, the link will be relative on the current host.
   */
  detachedFromSierra: boolean = false;

  layoutShowRemoteState: boolean = false;
  layoutShowSequence: boolean = false;
  layoutShowLog: boolean = false;
  layoutShowLeftPanel: boolean = false;

  collapseLeftPanel: boolean = true;
  collapseRightPanel: boolean = true;
  viewportShowExpand: boolean = false;
  viewportShowControls: boolean = false;

  viewportShowSelectionMode: boolean = true;
  volumeStreamingServer: string = 'https://maps.rcsb.org/';
  backgroundColor: number = 0xffffff; // white
  showWelcomeToast: boolean = false;

  modelUrlProviders: string = '';
  showExportControls: boolean = false;
}();

export type BiostructureProps = typeof BiostructurePropsDefault;

export interface IBiostructureViewer extends IViewer {
  get dataEff(): BiostructureData | null;
  setOptions(options: Partial<BiostructureProps>): void;
}


export interface ISaguaroViewer extends IViewer {
  setData(df: PdbResDataFrameType | DG.DataFrame): void;
}

declare module 'datagrok-api/dg' {
  interface DataFramePlotHelper {
    // eslint-disable-next-line max-len
    fromType(viewerType: 'Biostructure', options: Partial<BiostructureProps>): Promise<DG.Viewer<BiostructureProps> & IBiostructureViewer>;
  }
}

export async function viewBiostructure(content: string, format: string = 'pdb', name?: string): Promise<void> {
  const packageName: string = 'BiostructureViewer';
  const funcName: string = 'viewBiostructure';
  await grok.functions.call(`${packageName}:${funcName}`, {content, format, name});
}
