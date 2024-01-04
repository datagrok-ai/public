import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';
import {MjolnirPointerEvent} from 'mjolnir.js';
import {PickingInfo} from '@deck.gl/core/typed';

import {TAGS as pdbTAGS} from '../pdb/index';
import {IViewer} from './viewer';
import {BiostructureDataJson} from '../pdb/types';


export type NodeStyleType = { [propName: string]: any };
export type StylesType = { [nodeName: string]: NodeStyleType };

export enum RepresentationType {
  Cartoon = 'cartoon',
  Backbone = 'backbone',
  BallAndStick = 'ball+stick',
  Licorice = 'licorice',
  Hyperball = 'hyperball',
  Surface = 'surface'
}

export const NglPropsDefault = new class {
  // -- Data --
  dataJson: string = BiostructureDataJson.empty;
  pdb: string | null = null;
  pdbTag: string | null = pdbTAGS.PDB;
  ligandColumnName: string | null = null;

  // -- Style --
  representation: RepresentationType = RepresentationType.Cartoon;

  // -- Behaviour --

  showSelectedRowsLigands: boolean = false;
  showCurrentRowLigand: boolean = true;
  showMouseOverRowLigand: boolean = true;
}();


export type NglProps = typeof NglPropsDefault;

export interface INglViewer extends IViewer {
  setOptions(options: Partial<NglProps>): void;

  get onAfterBuildView(): Observable<void>;
}

declare module 'datagrok-api/dg' {
  interface DataFramePlotHelper {
    fromType(viewerType: 'NGL', options: Partial<NglProps>): Promise<DG.Viewer<NglProps> & INglViewer>;
  }
}
