import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {_package} from '../package';
import {
  CdrMapType,
  DataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType, NumsType,
  ObsPtmType,
  PtmMapType
} from './data-loader';

export class DataLoaderDb extends DataLoader {
  private _pName: string = 'MolecularLiabilityBrowser';

  private _files: { [name: string]: string } = {
    filterProps: 'properties.json',
    mutcodes: 'mutcodes.json',
    ptm_map: 'ptm_map.json',
    cdr_map: 'cdr_map.json',
    ptm_in_cdr: 'ptm_in_cdr.d42',
    h_out: 'h_out.csv',
    l_out: 'l_out.csv',
    tree: 'tree.csv',
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
    // Here we should load files from src/externalData
    // But if we will use require(), commit will fail

    // Checking files is disabled because it takes too long
    // await this.check_files(this._files);

    this._filterProperties = JSON.parse(await _package.files.readAsText(this._files.filterProps));
    this._mutcodes = JSON.parse(await _package.files.readAsText(this._files.mutcodes));
    this._ptmMap = JSON.parse(await _package.files.readAsText(this._files.ptm_map));
    this._cdrMap = JSON.parse(await _package.files.readAsText(this._files.cdr_map));
    this._refDf = await _package.files.readBinaryDataFrames(this._files.ptm_in_cdr)
      .then((dfList: DG.DataFrame[]) => dfList[0]);
    this._realNums = JSON.parse(await _package.files.readAsText(this._files.realNums));
  }

  async getVids(): Promise<string[]> {
    return await grok.functions.call(`${this._pName}:getVids`)
      .then((df: DG.DataFrame) => df.columns[0].toList());
  }

  async listAntigens(): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:listAntigens`);
    return df;
  }

  async getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getMlbByAntigen`, {antigen: antigen});
    return df;
  }

  async getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getTreeByAntigen`, {antigen: antigen});
    return df;
  }

  async getAnarci(scheme: string, chain: string, antigen: string): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await grok.functions.call(
      `${this._pName}:getAnarci_${scheme}_${chain}`, {antigen: antigen});
    return df;
  }

  async getObservedPtmVids(): Promise<string[]> {
    return await grok.functions.call(`${this._pName}:getObservedPtmVids`)
      .then((df: DG.DataFrame) => df.columns[0].toList());
  }

  async load_example(vid: string): Promise<JsonType> {
    return JSON.parse((await grok.functions.call(`${this._pName}:getJsonByVid`, {vid: vid})).columns[0].get(0));
  }

  async load_pdb(vid: string): Promise<string> {
    return (await grok.functions.call(`${this._pName}:getPdbByVid`, {vid: vid})).columns[0].get(0);
  }

  async load_hChainDf(): Promise<DG.DataFrame> {
    // Could not find chains in old MLB
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.h_out));
  }

  async load_lChainDf(): Promise<DG.DataFrame> {
    // Could not find chains in old MLB
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.l_out));
  }

  async load_mlbDf(): Promise<DG.DataFrame> {
    const df = await grok.functions.call(`${this._pName}:GetMolecularLiabilityBrowser`);
    // 'ngl' column have been removed from query 2022-04
    df.columns.remove('ngl');
    return df;
  }

  async load_obsPtm(vid: string): Promise<ObsPtmType> {
    return JSON.parse((await grok.functions.call(`${this._pName}:getJsonObsByVid`, {vid: vid})).columns[0].get(0));
  }

  async load_treeDf(): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await _package.files.readAsText(this._files.tree));
  }
}
