import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

// import {encode as encodeToBase64, decode as decodeFromBase64} from 'uint8-to-base64';

import {
  catchToLog,
  CdrMapType,
  DataLoader, FilesForDataLoader, QueriesForDataLoader, DataQueryDict,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType,
  ObsPtmType,
  PtmMapType
} from './data-loader';
import {MlbDatabase} from './mlb-database';
import {channel} from 'diagnostics_channel';
import {VdRegion} from '@datagrok-libraries/bio/src/vd-regions';
import {packageName} from '../package';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

export class DataLoaderDb extends DataLoader {
  static _files = class {
    static tree = 'tree.csv';
  };

  private _startInit: number;

  private readonly _cache: MlbDatabase;
  private readonly _dlFiles: FilesForDataLoader;
  private readonly _dlQueries: QueriesForDataLoader;

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

  private fromStartInit(): string {
    return ((Date.now() - this._startInit) / 1000).toString();
  }

  public constructor(mlbQueries: DataQueryDict, serverListVersionDf?: DG.DataFrame) {
    super();

    this._cache = new MlbDatabase(serverListVersionDf);
    this._dlQueries = new QueriesForDataLoader(this._cache, mlbQueries);
    this._dlFiles = new FilesForDataLoader(this._cache);
  }

  async init(startInit: number) {
    this._startInit = startInit;
    // Here we should load files from src/externalData
    // But if we will use require(), commit will fail
    console.debug(`MLB: DataLoaderDb.init(), ${this.fromStartInit()} s`);

    // Checking files is disabled because it takes too long
    // await this.check_files(this._files);
    console.debug(`MLB: DataLoaderDb.init() check_files, ${this.fromStartInit()} s`);

    console.debug('MLB: DataLoaderDb.init() start, ' + `${this.fromStartInit()} s`);

    await Promise.all([
      this._dlQueries.listSchemes().then((value) => { this._schemes = value; }),
      this._dlQueries.listCdrs().then((value) => { this._cdrs = value; }),
      this._dlQueries.listAntigens().then((value) => { this._antigens = value; }),
      this._dlQueries.getVids().then((value) => { this._vids = value; }),
      this._dlQueries.getObservedPtmVids().then((value) => { this._vidsObsPtm = value; }),
      this._dlFiles.getFilterProperties().then((value) => { this._filterProperties = value; }),
      this._dlFiles.getMutcodes().then((value) => { this._mutcodes = value; }),
      this._dlFiles.getPredictedPtmMap().then((value) => { this._predictedPtmMap = value; }),
      this._dlFiles.getPredictedCdrMap().then((value) => { this._predictedCdrMap = value; }),
      this._dlFiles.getObservedPtmMap().then((value) => { this._observedPtmMap = value; }),
      this._dlFiles.getObservedCdrMap().then((value) => { this._observedCdrMap = value; }),
    ]);

    console.debug('MLB: DataLoaderDb.init() end, ' + `${this.fromStartInit()} s`);
  }

  async getLayoutBySchemeCdr(scheme: string, cdr: string): Promise<VdRegion[]> {
    return this._dlQueries.getLayoutBySchemeCdr(scheme, cdr);
  }

  async getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    return this._dlQueries.getMlbByAntigen(antigen);
  }

  async getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    return this._dlQueries.getTreeByAntigen(antigen);
  }

  async getAnarci(scheme: string, chain: string, antigen: string): Promise<DG.DataFrame> {
    return this._dlQueries.getAnarci(scheme, chain, antigen);
  }

  async loadMlbDf(): Promise<DG.DataFrame> {
    return this._dlQueries.loadMlbDf();
  }

  /** deprecated */
  async loadTreeDf(): Promise<DG.DataFrame> {
    const dfTxt: string = await grok.dapi.files
      .readAsText(`System:Data/${packageName}/${DataLoaderDb._files.tree}`);
    return DG.DataFrame.fromCsv(dfTxt);
  }

  // -- PTM --

  async getPredictedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    return this._dlQueries.getPredictedPtmByAntigen(antigen);
  }

  async getObservedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    return this._dlQueries.getObservedPtmByAntigen(antigen);
  }

  // -- 3D --

  async load3D(vid: string): Promise<[JsonType, string, NumsType, ObsPtmType]> {
    try {
      // const df: DG.DataFrame = await grok.functions.call(`${this._pName}:get3D`, {vid: vid});
      const df: DG.DataFrame = await this._dlQueries.load3D(vid);

      const jsonTxt = df.get('json', 0);
      const pdbTxt = df.get('pdb', 0);
      const realNumsTxt = df.get('real_nums', 0);
      const obsPtmTxt = df.get('obs_ptm', 0);

      const jsonValue: JsonType = jsonTxt ? JSON.parse(jsonTxt) : null;
      const pdbValue: string = pdbTxt;
      const realNumsValue: NumsType = realNumsTxt ? JSON.parse(realNumsTxt) : null;
      const obsPtmValue: ObsPtmType = obsPtmTxt ? JSON.parse(obsPtmTxt)['ptm_observed'] : null;

      return [jsonValue, pdbValue, realNumsValue, obsPtmValue];
    } catch (err: any) {
      const errStr = errorToConsole(err);
      console.error(`MLB: query get3D('${vid}' error: ${errStr}`);
      throw err;
    }
  }
}
