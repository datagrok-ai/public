import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {encode as encodeToBase64, decode as decodeFromBase64} from 'uint8-to-base64';

import {
  catchToLog,
  CdrMapType,
  DataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType,
  ObsPtmType,
  PtmMapType,
} from './data-loader';
import {IAntigen, ICdr, IScheme, IVid, IVidObsPtm, MlbDatabase} from './mlb-database';

export class DataLoaderFiles extends DataLoader {
  private _startInit: number;
  private _serverListVersionDf: DG.DataFrame;

  private _files: { [name: string]: string } = {
    filterProps: 'properties.json',
    mutcodes: 'mutcodes.json',
    predictedPtmMap: 'ptm_map.json',
    predictedCdrMap: 'cdr_map.json',
    observedPtmMap: 'obs_ptm_map.json',
    observedCdrMap: 'obs_cdr_map.json',
    predictedPtm: 'ptm_in_cdr.d42',
    mlb: 'mlb.csv',
    tree: 'tree.csv',
    example: 'example.json',
    examplePDB: 'examplePDB.json',
    exampleOptm: 'exampleOptm.json',
    realNums: 'exampleNums.json',
  };

  private cache!: MlbDatabase;

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

  async init(startInit: number, serverListVersionDf: DG.DataFrame): Promise<void> {
    this._startInit = startInit;
    this._serverListVersionDf = serverListVersionDf;

    console.debug(`MLB: DataLoaderFiles.init(), ${this.fromStartInit()} s`);
    // Check files disabled while too much time consuming
    // await this.check_files(this._files);
    console.debug('MLB: DataLoaderFiles.init() check_files, ' +
      `${this.fromStartInit()} s`);

    await this.init2();

    console.debug('MLB: DataLoaderFiles.init() preload_data, ' +
      `${this.fromStartInit()} s`);
  }

  async init2(): Promise<void> {
    console.debug('MLB: DataLoaderFiles.init2() start, ' + `${this.fromStartInit()} s`);

    this.cache = new MlbDatabase();
    this.cache.init(this._serverListVersionDf);

    await Promise.all([
      //load numbering schemes
      this.cache.getData<IScheme, string[]>('scheme',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'listSchemes\': ',
          () => grok.functions.call(`${this._pName}:listSchemes`)),
        (dbRow: DG.Row) => ({scheme_id: dbRow.get('scheme_id'), scheme: dbRow.get('scheme')}),
        (objList: IScheme[]) => objList.map((obj: IScheme) => obj.scheme)
      ).then((value: string[]) => {
        console.debug('MLB: DataLoaderFiles.init2() set schemes, ' + `${this.fromStartInit()} s`);
        this._schemes = value;
      }),
      //load cdr definition list
      this.cache.getData<ICdr, string[]>('cdr',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'listCdrs\': ',
          () => grok.functions.call(`${this._pName}:listCdrs`)),
        (dbRow: DG.Row) => ({cdr_id: dbRow.get('cdr_id'), cdr: dbRow.get('cdr')}),
        (objList: ICdr[]) => objList.map((obj: ICdr) => obj.cdr)
      ).then((value: string[]) => {
        console.debug('MLB: DataLoaderFiles.init2() set cdrs, ' + `${this.fromStartInit()} s`);
        this._cdrs = value;
      }),
      //load antigen list
      this.cache.getData<IAntigen, DG.DataFrame>('antigen',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'listAntigens\': ',
          () => grok.functions.call(`${this._pName}:listAntigens`)),
        (dbRow: DG.Row) => Object.assign({},
          ...(['id', 'antigen', 'antigen_ncbi_id', 'antigen_gene_symbol']
            .map((fn: string) => ({[fn]: dbRow.get(fn)})))),
        (objList: IAntigen[]) => DG.DataFrame.fromObjects(objList)
      ).then((value: DG.DataFrame) => {
        console.debug('MLB: DataLoaderFiles.init2() set antigens, ' + `${this.fromStartInit()} s`);
        this._antigens = value;
      }),
      //load available Vids
      this.cache.getData<IVid, string[]>('vid',
        () => new Promise((resolve, reject) => {
          const df: DG.DataFrame = DG.DataFrame.fromObjects(
            [{v_id: 'VR000000008'}, {v_id: 'VR000000043'}, {v_id: 'VR000000044'}]);
          resolve(df);
        }),
        (dbRow: DG.Row) => Object.assign({},
          ...(['v_id'].map((fn: string) => ({[fn]: dbRow.get(fn)})))),
        (objList: IVid[]) => objList.map((obj) => obj.v_id)
      ).then((value: string[]) => {
        console.debug(`MLB: DataLoaderFiles.init() set vids, ${this.fromStartInit()} s`);
        this._vids = value;
      }),
      //load observed PTM data
      this.cache.getData<IVidObsPtm, string[]>('vidObsPtm',
        () => new Promise((resolve, reject) => {
          const df: DG.DataFrame = DG.DataFrame.fromObjects([{v_id: 'VR000000044'}]);
          resolve(df);
        }),
        (dbRow: DG.Row) => Object.assign({},
          ...(['v_id'].map((fn: string) => ({[fn]: dbRow.get(fn)})))),
        (objList: IVidObsPtm[]) => objList.map((obj: IVidObsPtm) => obj.v_id)
      ).then((value: string[]) => {
        console.debug(`MLB: DataLoaderFiles.init2() set obsPtmVids, ${this.fromStartInit()} s`);
        this._vidsObsPtm = value;
      }),
      //load properties data
      this.cache.getObject<FilterPropertiesType>(this._files.filterProps,
        async () => {
          const txt: string = await grok.dapi.files
            .readAsText(`System:AppData/${this._pName}/${this._files.filterProps}`);
          return JSON.parse(txt);
        })
        .then((value: FilterPropertiesType) => {
          console.debug(`MLB: DataLoaderFiles.init2() set filterProperties, ${this.fromStartInit()} s`);
          this._filterProperties = value;
        }),
      //load mutcodes data
      this.cache.getObject<MutcodesDataType>(this._files.mutcodes,
        async () => {
          const txt: string = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.mutcodes}`);
          return JSON.parse(txt);
        })
        .then((value: MutcodesDataType) => {
          console.debug(`MLB: DataLoaderFiles.init2() set mutcodes, ${this.fromStartInit()} s`);
          this._mutcodes = value;
        }),
      //load predicted ptm map data
      this.cache.getObject<PtmMapType>(this._files.predictedPtmMap,
        async () => {
          const txt: string = await grok.dapi.files
            .readAsText(`System:AppData/${this._pName}/${this._files.predictedPtmMap}`);
          return JSON.parse(txt);
        })
        .then((value: PtmMapType) => {
          console.debug(`MLB: DataLoaderFiles.init2() set predictedPtmMap, ${this.fromStartInit()} s`);
          this._predictedPtmMap = value;
        }),
      //load predicted cdr map data
      this.cache.getObject<CdrMapType>(this._files.predictedCdrMap,
        async () => {
          const txt: string = await grok.dapi.files
            .readAsText(`System:AppData/${this._pName}/${this._files.predictedCdrMap}`);
          return JSON.parse(txt);
        })
        .then((value: CdrMapType) => {
          console.debug(`MLB: DataLoaderFiles.init2() set predictedCdrMap, ${this.fromStartInit()} s`);
          this._predictedCdrMap = value;
        }),
      //load observed ptm map data
      this.cache.getObject<PtmMapType>(this._files.observedPtmMap,
        async () => {
          const jsonTxt: string = await grok.dapi.files
            .readAsText(`System:AppData/${this._pName}/${this._files.observedPtmMap}`);
          return JSON.parse(jsonTxt);
        })
        .then((value: PtmMapType) => {
          console.debug(`MLB: DataLoaderDb.init2() set observedPtmMap, ${this.fromStartInit()} s`);
          this._observedPtmMap = value;
        }),
      // load observed cdr map data
      this.cache.getObject<CdrMapType>(this._files.observedCdrMap,
        async () => {
          const jsonTxt: string = await grok.dapi.files
            .readAsText(`System:AppData/${this._pName}/${this._files.observedCdrMap}`);
          return JSON.parse(jsonTxt);
        })
        .then((value: CdrMapType) => {
          console.debug(`MLB: DataLoaderDb.init2() set observedCdrMap, ${this.fromStartInit()} s`);
          this._observedCdrMap = value;
        })
    ]);

    console.debug('MLB: DataLoaderFiles.init2() end, ' + `${this.fromStartInit()} s`);
  }

  async getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getMlbByAntigen`, {antigen: antigen});
    return df;
  }

  async getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getTreeByAntigen`, {antigen: antigen});
    return df;
  }

  async loadMlbDf(): Promise<DG.DataFrame> {
    const dfTxt: string = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.mlb}`);
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

  async loadTreeDf(): Promise<DG.DataFrame> {
    const dfTxt: string = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.tree}`);
    return DG.DataFrame.fromCsv(dfTxt);
  }

  private async loadFileJson(path: string): Promise<Object> {
    const jsonTxt = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${path}`);
    return JSON.parse(jsonTxt);
  }

  // -- PTM --

  async getPredictedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df = await grok.functions.call(`${this._pName}:getPredictedPtmByAntigen`, {antigen: antigen});
    return df;
  }

  async getObservedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df = await grok.functions.call(`${this._pName}:getObservedPtmByAntigen`, {antigen: antigen});
    return df;
  }

  // -- 3D --

  async load3D(vid: string): Promise<[JsonType, string, NumsType, ObsPtmType]> {
    return Promise.all([
      this.loadFileJson(this._files.example).then((value) => {
        return value as JsonType;
      }),
      this.loadFileJson(this._files.examplePDB).then((value) => {
        return value['pdb'];
      }),
      this.loadFileJson(this._files.realNums).then((value) => {
        return value as NumsType;
      }),
      this.loadFileJson(this._files.exampleOptm).then((value) => {
        return value['ptm_observed'] as ObsPtmType;
      }),
    ]);
  }
}
