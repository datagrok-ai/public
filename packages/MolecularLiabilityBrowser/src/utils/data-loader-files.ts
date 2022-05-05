import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {
  CdrMapType,
  DataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType,
  ObsPtmType,
  PtmMapType
} from './data-loader';

export class DataLoaderFiles extends DataLoader {
  private _files = {
    filterProps: 'properties.json',
    mutcodes: 'mutcodes.json',
    ptm_map: 'ptm_map.json',
    cdr_map: 'cdr_map.json',
    ptm_in_cdr: 'ptm_in_cdr.d42',
    h_out: 'h_out.csv',
    l_out: 'l_out.csv',
    mlb: 'mlb.csv',
    tree: 'tree.csv',
    example: 'example.json',
    examplePDB: 'examplePDB.json',
    exampleOptm: 'exampleOptm.json',
    realNums: 'exampleNums.json',
  };
  private _filterProperties: FilterPropertiesType;
  private _mutcodes: MutcodesDataType;
  private _ptmMap: PtmMapType;
  private _cdrMap: CdrMapType;
  private _refDf: DG.DataFrame;
  private _realNums: any;

  get filterProperties(): FilterPropertiesType { return this._filterProperties; }

  get mutcodes(): MutcodesDataType { return this._mutcodes; }

  get ptmMap(): PtmMapType { return this._ptmMap; }

  get cdrMap(): CdrMapType { return this._cdrMap; }

  get refDf(): DG.DataFrame { return this._refDf; }

  get realNums(): NumsType { return this._realNums; }

  async init() {
    await this.check_files(this._files);

    this._filterProperties = JSON.parse(await _package.files.readAsText(this._files.filterProps));
    this._mutcodes = JSON.parse(await _package.files.readAsText(this._files.mutcodes));
    this._ptmMap = JSON.parse(await _package.files.readAsText(this._files.ptm_map));
    this._cdrMap = JSON.parse(await _package.files.readAsText(this._files.cdr_map));
    this._refDf = (await _package.files.readBinaryDataFrames(this._files.ptm_in_cdr))[0];
    this._realNums = JSON.parse(await _package.files.readAsText(this._files.realNums));
  }

  async getVids(): Promise<string[]> {
    const res: string[] = ['VR000000008', 'VR000000043', 'VR000000044'];
    return Promise.resolve(res);
  }

  async getObservedPtmVids(): Promise<string[]> {
    const res: string[] = ['VR000000044'];
    return Promise.resolve(res);
  }

  async load_hChainDf(): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.h_out));
  }

  async load_lChainDf(): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.l_out));
  }

  async load_mlbDf(): Promise<DG.DataFrame> {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.mlb));
    for (const column of df.columns)
      column.name = column.name.replaceAll('_', ' ');

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
    df.changeColumnType('antigen ncbi id', DG.COLUMN_TYPE.STRING);

    return df;
  }

  async load_treeDf(): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.tree));
  }

  private async load_file_json(path: string): Promise<Object> {
    return _package.files.readAsText(path)
      .then((data: string) => JSON.parse(data));
  }

  async load_example(vid: string): Promise<JsonType> {
    return this.load_file_json(this._files.example)
      .then((o) => <JsonType>o);
  }

  /** Load PDB structure data
   * @param {string} vid Molecule id
   */
  async load_pdb(vid: string): Promise<string> {
    // TODO: Check for only allowed vid of example
    return this.load_file_json(this._files.examplePDB)
      .then((o) => (o)['pdb']);
  }

  async load_obsPtm(vid: string): Promise<ObsPtmType> {
    return this.load_file_json(this._files.exampleOptm)
      .then((o) => <ObsPtmType>o['ptm_observed']);
  }
}
