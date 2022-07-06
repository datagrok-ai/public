import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {encode as encode_to_base64, decode as decode_from_base64} from 'uint8-to-base64';

import {
  catchToLog,
  CdrMapType,
  DataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType,
  ObsPtmType,
  PtmMapType
} from './data-loader';
import {IAntigen, ICdr, IVidObsPtm, IScheme, IVid, MlbDatabase} from './mlb-database';

export class DataLoaderDb extends DataLoader {
  private _startInit: number;
  private _serverListVersionDf: DG.DataFrame;

  private _files: { [name: string]: string } = {
    filterProps: 'properties.json',
    mutcodes: 'mutcodes.json',
    ptm_map: 'ptm_map.json',
    cdr_map: 'cdr_map.json',
    ptm_in_cdr: 'ptm_in_cdr.d42',
    h_out: 'h_out.csv',
    l_out: 'l_out.csv',
    tree: 'tree.csv',
  };

  private cache!: MlbDatabase;

  private _schemes: string[];
  private _cdrs: string[];
  private _antigens: DG.DataFrame;
  private _vids: string[];
  private _vidsObsPtm: string[];
  private _filterProperties: FilterPropertiesType;
  private _mutcodes: MutcodesDataType;
  private _ptmMap: PtmMapType;
  private _cdrMap: CdrMapType;
  private _refDf: DG.DataFrame;

  get schemes(): string[] { return this._schemes; }

  get cdrs(): string[] { return this._cdrs; }

  get antigens(): DG.DataFrame { return this._antigens; }

  get vids(): string[] { return this._vids; }

  get vidsObsPtm(): string[] { return this._vidsObsPtm; }

  get filterProperties(): FilterPropertiesType { return this._filterProperties; }

  get mutcodes(): MutcodesDataType { return this._mutcodes; }

  get ptmMap(): PtmMapType { return this._ptmMap; }

  get cdrMap(): CdrMapType { return this._cdrMap; }

  get refDf(): DG.DataFrame { return this._refDf; }

  private fromStartInit(): string {
    return ((Date.now() - this._startInit) / 1000).toString();
  }

  async init(startInit: number, serverListVersionDf: DG.DataFrame) {
    this._startInit = startInit;
    this._serverListVersionDf = serverListVersionDf;
    // Here we should load files from src/externalData
    // But if we will use require(), commit will fail
    console.debug(`MLB: DataLoaderDb.init(), ${this.fromStartInit()} s`);

    // Checking files is disabled because it takes too long
    // await this.check_files(this._files);
    console.debug(`MLB: DataLoaderDb.init() check_files, ${this.fromStartInit()} s`);

    await this.init2();

    await Promise.all([
      // catchToLog<Promise<DG.DataFrame>>('MLB database error \'listSchemes\': ',
      //   () => grok.functions.call(`${this._pName}:listSchemes`))
      //   .then((df: DG.DataFrame) => {
      //     console.debug(`MLB: DataLoaderDb.init() set schemes, ${this.fromStartInit()} s`);
      //     this._schemes = df.columns.byName('scheme').toList();
      //   }),
      // catchToLog<Promise<DG.DataFrame>>('MLB database error \'listCdrs\': ',
      //   () => grok.functions.call(`${this._pName}:listCdrs`))
      //   .then((df: DG.DataFrame) => {
      //     console.debug(`MLB: DataLoaderDb.init() set cdrs, ${this.fromStartInit()} s`);
      //     this._cdrs = df.columns.byName('cdr').toList();
      //   }),
      // catchToLog<Promise<DG.DataFrame>>('MLB database error \'listAntigens\': ',
      //   () => grok.functions.call(`${this._pName}:listAntigens`))
      //   .then((df: DG.DataFrame) => {
      //     console.debug(`MLB: DataLoaderDb.init() set antigens, ${this.fromStartInit()} s`);
      //     this._antigens = df;
      //   }),
      // grok.functions.call(`${this._pName}:getVids`).then(
      //   (df: DG.DataFrame) => {
      //     console.debug(`MLB: DataLoaderDb.init() set vids, ${this.fromStartInit()} s`);
      //     this._vids = df.columns.byIndex(0).toList();
      //   }),
      // grok.functions.call(`${this._pName}:getObservedPtmVids`).then(
      //   (df: DG.DataFrame) => {
      //     console.debug(`MLB: DataLoaderDb.init() set vidsObsPtm, ${this.fromStartInit()} s`);
      //     this._vidsObsPtm = df.columns.byIndex(0).toList();
      //   }),
      // grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.filterProps}`)
      //   .then((v) => {
      //     console.debug(`MLB: DataLoaderFiles.init() set filterProperties, ` +
      //       `${this.fromStartInit()} s`);
      //     this._filterProperties = JSON.parse(v);
      //   }),
      // // _package.files.readAsText(this._files.mutcodes)
      // grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.mutcodes}`)
      //   .then((v) => {
      //     console.debug('MLB: DataLoaderFiles.init() set mutcodes, ' +
      //       `${this.fromStartInit()} s`);
      //     this._mutcodes = JSON.parse(v);
      //   }),
      // //_package.files.readAsText(this._files.ptmMap)
      // grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.ptmMap}`)
      //   .then((v) => {
      //       console.debug('MLB: DataLoaderFiles.init() set ptmMap, ' +
      //         `${this.fromStartInit()} s`);
      //       this._ptmMap = JSON.parse(v);
      //     }
      //   ),
      // //_package.files.readAsText(this._files.cdrMap)
      // grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.cdrMap}`)
      //   .then((v) => {
      //     console.debug('MLB: DataLoaderFiles.init() set cdrMap, ' +
      //       `${this.fromStartInit()} s`);
      //     this._cdrMap = JSON.parse(v);
      //   }),
      // // _package.files.readBinaryDataFrames(this._files.ptmInCdr)
      // grok.dapi.files.readAsBytes(`System:AppData/${this._pName}/${this._files.ptmInCdr}`)
      //   .then((data) => {
      //     const dfList = DG.DataFrame.fromByteArray(data);
      //     console.debug('MLB: DataLoaderFiles.init() set refDf, ' +
      //       `${this.fromStartInit()} s`);
      //     this._refDf = dfList[0];
      //   }),
    ]);

    console.debug(`MLB: DataLoaderDb.init() preload_data, ${this.fromStartInit()} s`);
  }

  async init2(): Promise<void> {
    console.debug('MLB: DataLoaderDb.init2() start, ' + `${this.fromStartInit()} s`);

    const serverListVersionDf: DG.DataFrame = await grok.functions.call(`${this._pName}:getListVersion`);
    console.debug('MLB: DataLoaderDb.init2() serverListVersionDf, ' + `${this.fromStartInit()} s`);

    this.cache = new MlbDatabase();
    this.cache.init(this._serverListVersionDf);

    await Promise.all([
      this.cache.getData<IScheme, string[]>('scheme',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'listSchemes\': ',
          () => grok.functions.call(`${this._pName}:listSchemes`)),
        (dbRow: DG.Row) => ({scheme_id: dbRow.get('scheme_id'), scheme: dbRow.get('scheme')}),
        (objList: IScheme[]) => objList.map((obj: IScheme) => obj.scheme)
      ).then((value: string[]) => {
        console.debug('MLB: DataLoaderDb.init2() set schemes, ' + `${this.fromStartInit()} s`);
        this._schemes = value;
      }),
      this.cache.getData<ICdr, string[]>('cdr',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'listCdrs\': ',
          () => grok.functions.call(`${this._pName}:listCdrs`)),
        (dbRow: DG.Row) => ({cdr_id: dbRow.get('cdr_id'), cdr: dbRow.get('cdr')}),
        (objList: ICdr[]) => objList.map((obj: ICdr) => obj.cdr)
      ).then((value: string[]) => {
        console.debug('MLB: DataLoaderDb.init2() set cdrs, ' + `${this.fromStartInit()} s`);
        this._cdrs = value;
      }),
      this.cache.getData<IAntigen, DG.DataFrame>('antigen',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'listAntigens\': ',
          () => grok.functions.call(`${this._pName}:listAntigens`)),
        (dbRow: DG.Row) => Object.assign({},
          ...(['id', 'antigen', 'antigen_ncbi_id', 'antigen_gene_symbol'].map((fn: string) => ({[fn]: dbRow.get(fn)})))),
        (objList: IAntigen[]) => DG.DataFrame.fromObjects(objList)
      ).then((value: DG.DataFrame) => {
        console.debug('MLB: DataLoaderDb.init2() set antigens, ' + `${this.fromStartInit()} s`);
        this._antigens = value;
      }),
      this.cache.getData<IVid, string[]>('vid',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'getVids\': ',
          () => grok.functions.call(`${this._pName}:getVids`)),
        (dbRow: DG.Row) => Object.assign({},
          ...(['v_id',].map((fn: string) => ({[fn]: dbRow.get(fn)})))),
        (objList: IVid[]) => objList.map((obj: IVid) => obj.v_id)
      ).then((value: string[]) => {
        console.debug(`MLB: DataLoaderDb.init2() set vids, ${this.fromStartInit()} s`);
        this._vids = value;
      }),
      this.cache.getData<IVidObsPtm, string[]>('vid',
        () => catchToLog<Promise<DG.DataFrame>>(
          'MLB database error \'getObservedPtmVids\': ',
          () => grok.functions.call(`${this._pName}:getObservedPtmVids`)),
        (dbRow: DG.Row) => Object.assign({},
          ...(['v_id',].map((fn: string) => ({[fn]: dbRow.get(fn)})))),
        (objList: IVidObsPtm[]) => objList.map((obj: IVidObsPtm) => obj.v_id)
      ).then((value: string[]) => {
        console.debug(`MLB: DataLoaderDb.init2() set obsPtmVids, ${this.fromStartInit()} s`);
        this._vidsObsPtm = value;
      }),
      this.cache.getObject<FilterPropertiesType>(this._files.filterProps,
        async () => {
          const txt: string = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.filterProps}`);
          return JSON.parse(txt);
        })
        .then((value: FilterPropertiesType) => {
          console.debug(`MLB: DataLoaderDb.init2() set filterProperties, ${this.fromStartInit()} s`);
          this._filterProperties = value;
        }),
      this.cache.getObject<MutcodesDataType>(this._files.mutcodes,
        async () => {
          const txt: string = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.mutcodes}`);
          return JSON.parse(txt);
        })
        .then((value: MutcodesDataType) => {
          console.debug(`MLB: DataLoaderDb.init2() set mutcodes, ${this.fromStartInit()} s`);
          this._mutcodes = value;
        }),
      this.cache.getObject<PtmMapType>(this._files.ptmMap,
        async () => {
          const txt: string = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.ptmMap}`);
          return JSON.parse(txt);
        })
        .then((value: PtmMapType) => {
          console.debug(`MLB: DataLoaderDb.init2() set ptmMap, ${this.fromStartInit()} s`);
          this._ptmMap = value;
        }),
      this.cache.getObject<CdrMapType>(this._files.cdrMap,
        async () => {
          const txt: string = await grok.dapi.files.readAsText(`System:AppData/${this._pName}/${this._files.cdrMap}`);
          return JSON.parse(txt);
        })
        .then((value: CdrMapType) => {
          console.debug(`MLB: DataLoaderDb.init2() set cdrMap, ${this.fromStartInit()} s`);
          this._cdrMap = value;
        }),
      this.cache.getObject<string>(this._files.ptmInCdr,
        async () => {
          const t1: number = Date.now();
          const data: Uint8Array = await grok.dapi.files.readAsBytes(`System:AppData/${this._pName}/${this._files.ptmInCdr}`);
          const t2: number = Date.now();
          const txt_base64: string = encode_to_base64(data);
          const t3: number = Date.now();
          console.debug('MLB: DataLoaderDb.init2() refDf ' +
            `loading bytes ET ${((t2 - t1) / 1000).toString()} s, ` +
            `encode base64 ET ${((t3 - t2) / 1000).toString()} s`);
          return txt_base64;
        })
        .then((txt_base64: string) => {
          const t1: number = Date.now();
          const data: Uint8Array = decode_from_base64(txt_base64);
          const t2: number = Date.now();
          const df: DG.DataFrame = DG.DataFrame.fromByteArray(data);
          const t3: number = Date.now();
          console.debug('MLB: DataLoaderDb.init2() refDf ' +
            `decode base64 ET ${((t2 - t1) / 1000).toString()} s, ` +
            `fromByteArray ET ${((t3 - t2) / 1000).toString()} s`);
          console.debug(`MLB: DataLoaderDb.init2() set refDf, ${this.fromStartInit()} s`);
          this._refDf = df;
        }),
    ]);

    console.debug('MLB: DataLoaderDb.init2() end, ' + `${this.fromStartInit()} s`);
  }

  async getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getMlbByAntigen`, {antigen: antigen});
    return df;
  }

  async getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getTreeByAntigen`, {antigen: antigen});
    return df;
  }

  async loadJson(vid: string): Promise<JsonType> {
    if (!this.vids.includes(vid))
      return null;

    return JSON.parse((await grok.functions.call(`${this._pName}:getJsonByVid`, {vid: vid}))
      .columns.byIndex(0).get(0));
  }

  async loadPdb(vid: string): Promise<string> {
    if (!this.vids.includes(vid))
      return null;

    return (await grok.functions.call(`${this._pName}:getPdbByVid`, {vid: vid}))
      .columns.byIndex(0).get(0);
  }

  async loadRealNums(vid: string): Promise<NumsType> {
    if (!this.vids.includes(vid))
      return null;

    return JSON.parse((await grok.functions.call(`${this._pName}:getJsonComplementByVid`, {vid: vid}))
      .columns.byIndex(0).get(0));
  }

  async loadMlbDf(): Promise<DG.DataFrame> {
    const df = await grok.functions.call(`${this._pName}:GetMolecularLiabilityBrowser`);
    // 'ngl' column have been removed from query 2022-04
    df.columns.remove('ngl');
    return df;
  }

  async loadObsPtm(vid: string): Promise<ObsPtmType> {
    if (!this.vidsObsPtm.includes(vid))
      return null;

    return JSON.parse((await grok.functions.call(`${this._pName}:getJsonObsByVid`, {vid: vid}))
      .columns.byIndex(0).get(0));
  }

  async loadTreeDf(): Promise<DG.DataFrame> {
    const df_txt: string = await grok.dapi.files.readAsText(`System:Data/${this._pName}/${this._files.tree}`);
    return DG.DataFrame.fromCsv(df_txt);
  }
}
