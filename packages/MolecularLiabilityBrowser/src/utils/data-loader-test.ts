import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {
  CdrMapType,
  DataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType, ObsPtmType,
  PtmMapType
} from './data-loader';
import {ArgumentOutOfRangeError} from 'rxjs';
import {VdRegion} from '@datagrok-libraries/bio/src/vd-regions';

export class DataLoaderTest extends DataLoader {

  private _schemes: string[];
  private _cdrs: string[];
  private _antigens: DG.DataFrame;
  private _vids: string[];
  private _vidsObsPtm: string[];
  private _filterProperties: FilterPropertiesType;
  private _mutcodes: MutcodesDataType;

  private _predictedPtmMap: PtmMapType;
  private _predictedCdrMap: CdrMapType;
  private _observedPtmMap: PtmMapType;
  private _observedCdrMap: CdrMapType;

  get schemes(): string[] { return this._schemes; }

  get cdrs(): string[] { return this._cdrs; }

  get antigens(): DG.DataFrame { return this._antigens; }

  get vids(): string[] { return this._vids; }

  get vidsObsPtm(): string[] { return this._vidsObsPtm; }

  get filterProperties(): FilterPropertiesType { return this._filterProperties; }

  get mutcodes(): MutcodesDataType { return this._mutcodes; }

  get predictedPtmMap(): PtmMapType { return this._predictedPtmMap; }

  get predictedCdrMap(): CdrMapType { return this._predictedCdrMap; }

  get observedPtmMap(): PtmMapType { return this._observedPtmMap; }

  get observedCdrMap(): CdrMapType { return this._observedCdrMap; }


  async init(startInit: number): Promise<void> {
    return Promise.resolve(undefined);

    await Promise.all([
      //load numbering schemes
      () => {
        this._schemes = ['imgt', 'aho'];
      },
      //load cdr definition list
      () => {
        this._cdrs = ['chothia', 'aroop'];
      },
      //load antigen list
      () => {

      },

    ]);
  }

  getLayoutBySchemeCdr(scheme: string, cdr: string): Promise<VdRegion[]> {
    return Promise.resolve([]);
  }

  getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    return Promise.resolve(undefined);
  }

  getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    return Promise.resolve(undefined);
  }

  getAnarci(scheme: string, chain: string, antigen: string): Promise<DG.DataFrame> {
    return Promise.resolve(undefined);
  }

  loadMlbDf(): Promise<DG.DataFrame> {
    return Promise.resolve(undefined);
  }

  loadTreeDf(): Promise<DG.DataFrame> {
    return Promise.resolve(undefined);
  }

  // -- PTM --

  getObservedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    throw new Error('Not implemented');
    throw new ArgumentOutOfRangeError();
  }

  getPredictedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    throw new Error('Not implemented');
    throw new ArgumentOutOfRangeError();
  }

  // -- 3D --

  load3D(vid: string): Promise<[JsonType, string, NumsType, ObsPtmType]> {
    throw new Error('Not implemented');
  }
}