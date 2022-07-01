import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package, _startInit} from '../package';

import {
  CdrMapType,
  DataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType,
  ObsPtmType,
  PtmMapType,
  catchToLog,
} from './data-loader';

export class DataLoaderFiles extends DataLoader {
  private _files: { [name: string]: string } = {
    filterProps: 'properties.json',
    mutcodes: 'mutcodes.json',
    ptmMap: 'ptm_map.json',
    cdrMap: 'cdr_map.json',
    ptmInCdr: 'ptm_in_cdr.d42',
    hOut: 'h_out.csv',
    lOut: 'l_out.csv',
    mlb: 'mlb.csv',
    tree: 'tree.csv',
    example: 'example.json',
    examplePDB: 'examplePDB.json',
    exampleOptm: 'exampleOptm.json',
    realNums: 'exampleNums.json',
  };

  // private _localStorageKeys: { [name: string]: string } = Object.assign({}, DataLoaderFiles._files.);

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
  private _realNums: any;

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

  async init(): Promise<void> {
    Object.entries(this._files).map(([key, value]: [string, string]) => ({[key]: `MLB:{$value}`}));

    console.debug(`MLB: DataLoaderFiles.init(), ${((Date.now() - _startInit) / 1000).toString()} s`);
    // Check files disabled while too much time consuming
    // await this.check_files(this._files);
    console.debug(`MLB: DataLoaderFiles.init() check_files, ${((Date.now() - _startInit) / 1000).toString()} s`);

    await Promise.all([
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listSchemes\': ',
        () => grok.functions.call(`${this._pName}:listSchemes`))
        .then((df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderFiles.init() set schemes, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._schemes = df.columns.byName('scheme').toList();
        }),
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listCdrs\': ',
        () => grok.functions.call(`${this._pName}:listCdrs`))
        .then((df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderFiles.init() set cdrs, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._cdrs = df.columns.byName('cdr').toList();
        }),
      catchToLog<Promise<DG.DataFrame>>('MLB database error \'listAntigens\': ',
        () => grok.functions.call(`${this._pName}:listAntigens`))
        .then((df: DG.DataFrame) => {
          console.debug(`MLB: DataLoaderFiles.init() set antigens, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._antigens = df;
        }),
      new Promise((resolve: (value: unknown) => void, reject: (reason?: unknown) => void) => {
        console.debug(`MLB: DataLoaderFiles.init() set vids, ${((Date.now() - _startInit) / 1000).toString()} s`);
        this._vids = ['VR000000008', 'VR000000043', 'VR000000044'];
        resolve(true);
      }),
      new Promise((resolve: (value: unknown) => void, reject: (reason?: unknown) => void) => {
        this._vidsObsPtm = ['VR000000044'];
        resolve(true);
      }),
      _package.files.readAsText(this._files.filterProps).then(
        (v) => {
          console.debug(`MLB: DataLoaderFiles.init() set filterProperties, ` +
            `${((Date.now() - _startInit) / 1000).toString()} s`);
          this._filterProperties = JSON.parse(v);
        }
      ),
      _package.files.readAsText(this._files.mutcodes).then(
        (v) => {
          console.debug(`MLB: DataLoaderFiles.init() set mutcodes, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._mutcodes = JSON.parse(v);
        }
      ),
      _package.files.readAsText(this._files.ptmMap).then(
        (v) => {
          console.debug(`MLB: DataLoaderFiles.init() set ptmMap, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._ptmMap = JSON.parse(v);
        }
      ),
      _package.files.readAsText(this._files.cdrMap).then(
        (v) => {
          console.debug(`MLB: DataLoaderFiles.init() set cdrMap, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._cdrMap = JSON.parse(v);
        }
      ),
      // new Promise((resolve: (value: unknown) => void, reject: (reason?: unknown) => void) => {
      //   const dbr: IDBOpenDBRequest = indexedDB.open('MLB: ptm_in_cdr', 1);
      //
      //   dbr.
      // })
      _package.files.readBinaryDataFrames(this._files.ptmInCdr).then(
        (dfList) => {
          console.debug(`MLB: DataLoaderFiles.init() set refDf, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._refDf = dfList[0];
        }
      ),
      _package.files.readAsText(this._files.realNums).then(
        (v) => {
          console.debug(`MLB: DataLoaderFiles.init() set realNums, ${((Date.now() - _startInit) / 1000).toString()} s`);
          this._realNums = JSON.parse(v);
        }
      )]);

    console.debug(`MLB: DataLoaderFiles.init() preload_data, ${((Date.now() - _startInit) / 1000).toString()} s`);
  }

  async getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getMlbByAntigen`, {antigen: antigen});
    return df;
  }

  async getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getTreeByAntigen`, {antigen: antigen});
    return df;
  }

  async loadHChainDf(): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.hOut));
  }

  async loadLChainDf(): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.lOut));
  }

  async loadMlbDf(): Promise<DG.DataFrame> {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.mlb));

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
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.tree));
  }

  private async loadFileJson(path: string): Promise<Object> {
    return _package.files.readAsText(path)
      .then((data: string) => JSON.parse(data));
  }

  async loadJson(vid: string): Promise<JsonType> {
    // Always return example data due TwinPViewer inability to work with null data
    // if (!this.vids.includes(vid))
    //   return null;

    return (await this.loadFileJson(this._files.example)) as JsonType;
  }

  /** Load PDB structure data
   * @param {string} vid Molecule id
   */
  async loadPdb(vid: string): Promise<string> {
    // Always return example data due TwinPViewer inability to work with null data
    // if (!this.vids.includes(vid))
    //   return null;

    // TODO: Check for only allowed vid of example
    return (await this.loadFileJson(this._files.examplePDB))['pdb'];
  }

  async loadRealNums(vid: string): Promise<NumsType> {
    return (await this.loadFileJson(this._files.realNums)) as NumsType;
  }

  async loadObsPtm(vid: string): Promise<ObsPtmType> {
    if (!this.vidsObsPtm.includes(vid))
      return null;

    return (await this.loadFileJson(this._files.exampleOptm))['ptm_observed'] as ObsPtmType;
  }
}
