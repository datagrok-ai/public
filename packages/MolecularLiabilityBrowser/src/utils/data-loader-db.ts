import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package, _startInit} from '../package';

import {
  catchToLog,
  CdrMapType,
  DataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType, NumsType,
  ObsPtmType,
  PtmMapType
} from './data-loader';

export class DataLoaderDb extends DataLoader {
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

  get antigens(): DG.DataFrame { return this.antigens; }

  get vids(): string[] { return this._vids; }

  get vidsObsPtm(): string[] { return this._vidsObsPtm; }

  get filterProperties(): FilterPropertiesType { return this._filterProperties; }

  get mutcodes(): MutcodesDataType { return this._mutcodes; }

  get ptmMap(): PtmMapType { return this._ptmMap; }

  get cdrMap(): CdrMapType { return this._cdrMap; }

  get refDf(): DG.DataFrame { return this._refDf; }

  async init() {
    // Here we should load files from src/externalData
    // But if we will use require(), commit will fail
    console.debug(`MLB: DataLoaderDb.init(), ${((Date.now() - _startInit) / 1000).toString()} s`);

    // Checking files is disabled because it takes too long
    // await this.check_files(this._files);
    console.debug(`MLB: DataLoaderDb.init() check_files, ${((Date.now() - _startInit) / 1000).toString()} s`);

    await Promise.all([
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listSchemes\': ',
        () => grok.functions.call(`${this._pName}:listSchemes`))
        .then((df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderDb.init() set schemes, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._schemes = df.columns.byName('scheme').toList();
        }),
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listCdrs\': ',
        () => grok.functions.call(`${this._pName}:listCdrs`))
        .then((df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderDb.init() set cdrs, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._cdrs = df.columns.byName('cdr').toList();
        }),
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listAntigens\': ',
        () => grok.functions.call(`${this._pName}:listAntigens`))
        .then((df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderDb.init() set antigens, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._antigens = df;
        }),
      grok.functions.call(`${this._pName}:getVids`).then(
        (df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderDb.init() set vids, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._vids = df.columns.byIndex(0).toList();
        }),
      grok.functions.call(`${this._pName}:getObservedPtmVids`).then(
        (df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderDb.init() set vidsObsPtm, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._vidsObsPtm = df.columns.byIndex(0).toList();
        }),
      _package.files.readAsText(this._files.filterProps).then(
        (v) => {
          console.debug(`MLB: DataLoaderDb.init() set filterProperties, ` +
            `${((Date.now() - _startInit) / 1000).toString()} s`);
          this._filterProperties = JSON.parse(v);
        }),
      _package.files.readAsText(this._files.mutcodes).then(
        (v) => {
          console.debug(`MLB: DataLoaderDb.init() set mutCodes, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._mutcodes = JSON.parse(v);
        }),
      _package.files.readAsText(this._files.ptm_map).then(
        (v) => {
          console.debug(`MLB: DataLoaderDb.init() set ptmMap, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._ptmMap = JSON.parse(v);
        }),
      _package.files.readAsText(this._files.cdr_map).then(
        (v) => {
          console.debug(`MLB: DataLoaderDb.init() set cdrMap, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._cdrMap = JSON.parse(v);
        }),
      _package.files.readBinaryDataFrames(this._files.ptm_in_cdr).then(
        (dfList) => {
          console.debug(`MLB: DataLoaderDb.init() set refDf, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._refDf = dfList[0];
        }),
    ]);

    console.debug(`MLB: DataLoaderDb.init() preload_data, ${((Date.now() - _startInit) / 1000).toString()} s`);
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

  async loadHChainDf(): Promise<DG.DataFrame> {
    throw new Error('Obsolete');
    // Could not find chains in old MLB
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.h_out));
  }

  async loadLChainDf(): Promise<DG.DataFrame> {
    throw new Error('Obsolete');
    // Could not find chains in old MLB
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.l_out));
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
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.tree));
  }
}
