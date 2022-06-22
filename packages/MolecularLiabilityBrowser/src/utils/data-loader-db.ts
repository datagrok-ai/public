import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';

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
  private _vids: string[];
  private _vidsObsPtm: string[];
  private _filterProperties: FilterPropertiesType;
  private _mutcodes: MutcodesDataType;
  private _ptmMap: PtmMapType;
  private _cdrMap: CdrMapType;
  private _refDf: DG.DataFrame;

  get schemes(): string[] { return this._schemes; }

  get cdrs(): string[] { return this._cdrs; }

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

    const t1 = Date.now();
    // Checking files is disabled because it takes too long
    // await this.check_files(this._files);
    const t2 = Date.now();

    //this._filterProperties = JSON.parse(await _package.files.readAsText(this._files.filterProps));
    //this._mutcodes = JSON.parse(await _package.files.readAsText(this._files.mutcodes));
    //this._ptmMap = JSON.parse(await _package.files.readAsText(this._files.ptm_map));
    //this._cdrMap = JSON.parse(await _package.files.readAsText(this._files.cdr_map));
    //this._refDf = (await _package.files.readBinaryDataFrames(this._files.ptm_in_cdr))[0];
    // this._realNums = JSON.parse(await _package.files.readAsText(this._files.realNums));

    await Promise.all([
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listSchemes\': ',
        () => grok.functions.call(`${this._pName}:listSchemes`))
        .then((df: DG.DataFrame) => { this._schemes = df.columns.byName('scheme').toList(); }),
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listCdrs\': ',
        () => grok.functions.call(`${this._pName}:listCdrs`))
        .then((df: DG.DataFrame) => { this._cdrs = df.columns.byName('cdr').toList(); }),
      grok.functions.call(`${this._pName}:getVids`).then(
        (df: DG.DataFrame) => {
          this._vids = df.columns.byIndex(0).toList();
        }),
      grok.functions.call(`${this._pName}:getObservedPtmVids`).then(
        (df: DG.DataFrame) => {
          this._vidsObsPtm = df.columns.byIndex(0).toList();
        }),
      _package.files.readAsText(this._files.filterProps).then(
        (v) => {
          this._filterProperties = JSON.parse(v);
        }),
      _package.files.readAsText(this._files.mutcodes).then(
        (v) => {
          this._mutcodes = JSON.parse(v);
        }),
      _package.files.readAsText(this._files.ptm_map).then(
        (v) => { this._ptmMap = JSON.parse(v); }),
      _package.files.readAsText(this._files.cdr_map).then(
        (v) => {
          this._cdrMap = JSON.parse(v);
        }),
      _package.files.readBinaryDataFrames(this._files.ptm_in_cdr).then(
        (dfList) => {
          this._refDf = dfList[0];
        }),
    ]);
    const t3 = Date.now();
    console.debug(`DataLoaderDb check_files ${((t2 - t1) / 1000).toString()} s`);
    console.debug(`DataLoaderDb preload_data ${((t3 - t2) / 1000).toString()} s`);
  }

  async listAntigens(): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>('MLB database access error: ', async () => {
      const df: DG.DataFrame = await grok.functions.call(`${this._pName}:listAntigens`);
      return df;
    });
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
    // Could not find chains in old MLB
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.h_out));
  }

  async loadLChainDf(): Promise<DG.DataFrame> {
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
