import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {encode as encodeToBase64, decode as decodeFromBase64} from 'uint8-to-base64';

import {
  catchToLog,
  CdrMapType,
  DataLoader, DataQueryDict, FilesForDataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType,
  ObsPtmType,
  PtmMapType, QueriesForDataLoader,
} from './data-loader';
import {IAntigen, ICdr, IScheme, IVid, IVidObsPtm, MlbDatabase} from './mlb-database';
import {VdRegion} from '@datagrok-libraries/bio/src/vd-regions';
import {packageName} from '../package';

export class DataLoaderFiles extends DataLoader {
  private _startInit: number;

  private readonly _cache: MlbDatabase;
  private readonly _dlFiles: FilesForDataLoader;
  private readonly _dlQueries: QueriesForDataLoader;

  static _files = class {
    static mlb = 'mlb.csv';
    static tree: 'tree.csv';
    static example = 'example.json';
    static examplePDB = 'examplePDB.json';
    static exampleOptm = 'exampleOptm.json';
    static realNums = 'exampleNums.json';
  };

  // private _files: { [name: string]: string } = {
  //   filterProps: 'properties.json',
  //   mutcodes: 'mutcodes.json',
  //   predictedPtmMap: 'ptm_map.json',
  //   predictedCdrMap: 'cdr_map.json',
  //   observedPtmMap: 'obs_ptm_map.json',
  //   observedCdrMap: 'obs_cdr_map.json',
  //   predictedPtm: 'ptm_in_cdr.d42',

  // };

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

  private _realNums: any;


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

  public constructor(mlbQueries: DataQueryDict, serverListVersionDf: DG.DataFrame) {
    super();

    this._cache = new MlbDatabase(serverListVersionDf);
    this._dlQueries = new QueriesForDataLoader(this._cache, mlbQueries);
    this._dlFiles = new FilesForDataLoader(this._cache);
  }

  async init(startInit: number): Promise<void> {
    this._startInit = startInit;

    console.debug(`MLB: DataLoaderFiles.init(), ${this.fromStartInit()} s`);
    // Check files disabled while too much time-consuming
    // await this.check_files(this._files);
    console.debug('MLB: DataLoaderFiles.init() check_files, ' +
      `${this.fromStartInit()} s`);

    await Promise.all([
      this._dlQueries.listSchemes().then((value) => { this._schemes = value;}),
      this._dlQueries.listCdrs().then((value) => {this._cdrs = value; }),
      this._dlQueries.listAntigens().then((value) => { this._antigens = value; }),

      //load available Vids
      Promise.resolve<string[]>(['VR000000008', 'VR000000043', 'VR000000044'])
        .then((value) => { this._vids = value;}),

      //load observed PTM data
      Promise.resolve<string[]>(['VR000000044'])
        .then((value) => this._vidsObsPtm = value),

      this._dlFiles.getFilterProperties().then((value) => { this._filterProperties = value; }),
      this._dlFiles.getMutcodes().then((value) => { this._mutcodes = value; }),
      this._dlFiles.getPredictedPtmMap().then((value) => { this._predictedPtmMap = value; }),
      this._dlFiles.getPredictedCdrMap().then((value) => { this._predictedCdrMap = value; }),
      this._dlFiles.getObservedPtmMap().then((value) => { this._observedPtmMap = value; }),
      this._dlFiles.getObservedCdrMap().then((value) => { this._observedCdrMap = value; }),
    ]);

    console.debug('MLB: DataLoaderFiles.init() preload_data, ' + `${this.fromStartInit()} s`);
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
    const dfTxt: string = await grok.dapi.files
      .readAsText(`System:AppData/${packageName}/${DataLoaderFiles._files.mlb}`);
    const df: DG.DataFrame = DG.DataFrame.fromCsv(dfTxt);

    // 'ngl' column have been removed from query 2022-04
    df.columns.remove('ngl');

    // Convert 'antigen ncbi id' to string as it is a column with lists of ids
    // TODO: Convert column type does not work as expected
    /* Examples: text -> loaded -> converted
       "3592,51561,3593" -> 3592515613593 -> "3592515613593.00"
       "952" -> 952 -> "952.00"
       "3553,3590,335,4045,105805883,102116898,3543,27132,100423062,389812,10892,9080,5284,3904,8740"
         -> Infinity -> "Infinity"
     */
    // df.changeColumnType('antigen_ncbi_id', DG.COLUMN_TYPE.STRING);

    return df;
  }

  /** deprecated */
  async loadTreeDf(): Promise<DG.DataFrame> {
    const dfTxt: string = await grok.dapi.files
      .readAsText(`System:AppData/${packageName}/${DataLoaderFiles._files.tree}`);
    return DG.DataFrame.fromCsv(dfTxt);
  }

  private async loadFileJson(path: string): Promise<Object> {
    const jsonTxt = await grok.dapi.files.readAsText(`System:AppData/${packageName}/${path}`);
    return JSON.parse(jsonTxt);
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
    return Promise.all([
      this.loadFileJson(DataLoaderFiles._files.example).then((value) => {
        return value as JsonType;
      }),
      this.loadFileJson(DataLoaderFiles._files.examplePDB).then((value) => {
        return value['pdb'];
      }),
      this.loadFileJson(DataLoaderFiles._files.realNums).then((value) => {
        return value as NumsType;
      }),
      this.loadFileJson(DataLoaderFiles._files.exampleOptm).then((value) => {
        return value['ptm_observed'] as ObsPtmType;
      }),
    ]);
  }
}
