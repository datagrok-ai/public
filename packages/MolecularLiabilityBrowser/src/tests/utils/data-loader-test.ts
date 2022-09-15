import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {
  catchToLog,
  CdrMapType,
  DataLoader, FilesForDataLoader,
  FilterPropertiesType,
  JsonType,
  MutcodesDataType,
  NumsType, ObsPtmType,
  PtmMapType
} from '../../utils/data-loader';
import {ArgumentOutOfRangeError} from 'rxjs';
import {VdRegion} from '@datagrok-libraries/bio/src/vd-regions';
import {MlbDatabase} from '../../utils/mlb-database';
import {packageName} from '../../package';
import {DataLoaderDb} from '../../utils/data-loader-db';

export class DataLoaderTest extends DataLoader {
  static Files = class {
    static layoutSchemeCdr = 'layoutSchemeCdr';
    static mlbByAntigen = 'mlbByAntigen';
    static treeByAntigen: 'treeByAntigen';
    static predictedPtmByAntigen = 'predictedPtmByAntigen';
    static observedPtmByAntigen = 'observedPtmByAntigen';

    static load3dExample = 'load3dExample';
    static load3dExamplePdb = 'load3dExamplePdb';
    static load3dExampleOptm = 'load3dExampleOptm';
    static load3dExampleRealNums = 'load3dExampleRealNums';
  };
  static _files = {
    [FilesForDataLoader.Files.filterProps]: 'tests/properties.json',
    [FilesForDataLoader.Files.mutcodes]: 'tests/mutcodes.json',
    [FilesForDataLoader.Files.predictedPtmMap]: 'tests/ptm_map.json',
    [FilesForDataLoader.Files.predictedCdrMap]: 'tests/cdr_map.json',
    [FilesForDataLoader.Files.observedPtmMap]: 'tests/obs_ptm_map.json',
    [FilesForDataLoader.Files.observedCdrMap]: 'tests/obs_cdr_map.json',

    [DataLoaderTest.Files.layoutSchemeCdr]: 'tests/layoutSchemeCdr.csv',
    [DataLoaderTest.Files.mlbByAntigen]: 'tests/mlbByAntigen.csv',
    [DataLoaderTest.Files.treeByAntigen]: 'tests/treeByAntigen.csv',
    [DataLoaderTest.Files.predictedPtmByAntigen]: 'tests/predictedPtmByAntigen.csv',
    [DataLoaderTest.Files.observedPtmByAntigen]: 'tests/observedPtmByAntigen.csv',

    [DataLoaderTest.Files.load3dExample]: 'tests/load3D/example.json',
    [DataLoaderTest.Files.load3dExamplePdb]: 'tests/load3D/examplePDB.json',
    [DataLoaderTest.Files.load3dExampleOptm]: 'tests/load3D/exampleOptm.json',
    [DataLoaderTest.Files.load3dExampleRealNums]: 'tests/load3D/exampleNums.json',

    anarciAhoHeavyByAntigen: 'tests/anarciAhoHeavyByAntigen.csv',
    anarciAhoLightByAntigen: 'tests/anarciAhoLightByAntigen.csv',
    anarciImgtHeavyByAntigen: 'tests/anarciImgtHeavyByAntigen.csv',
    anarciImgtLightByAntigen: 'tests/anarciImgtLightByAntigen.csv',
  };

  private readonly _cache: MlbDatabase;
  private readonly _dlFiles: FilesForDataLoader;

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

  static fnFunc(file: string): string {
    return DataLoaderTest._files[file];
  }

  constructor() {
    super();

    this._cache = new MlbDatabase(null);
    this._dlFiles = new FilesForDataLoader(this._cache, DataLoaderTest.fnFunc);
  }

  async init(startInit: number): Promise<void> {
    //load numbering schemes
    this._schemes = ['imgt', 'aho'];
    //load cdr definition list
    this._cdrs = ['chothia', 'aroop'];

    //load antigen list
    this._antigens = DG.DataFrame.fromCsv(
      `id,antigen,antigen_ncbi_id,antigen_gene_symbol
1,B17W38,113223191,LOC113223191
252,IAPW8,3556,IL1RAP
`);

    // load available vids
    this._vids = ['VR000000008', 'VR000000043', 'VR000000044'];

    // load observed PTM data
    this._vidsObsPtm = ['VR000000044'];

    await Promise.all([
      this._dlFiles.getFilterProperties().then((value) => { this._filterProperties = value; }),
      this._dlFiles.getMutcodes().then((value) => { this._mutcodes = value; }),
      this._dlFiles.getPredictedPtmMap().then((value) => { this._predictedPtmMap = value; }),
      this._dlFiles.getPredictedCdrMap().then((value) => { this._predictedCdrMap = value; }),
      this._dlFiles.getObservedPtmMap().then((value) => { this._observedPtmMap = value; }),
      this._dlFiles.getObservedCdrMap().then((value) => { this._observedCdrMap = value; }),
    ]);
    let k = 11;
  }

  async getLayoutBySchemeCdr(scheme: string, cdr: string): Promise<VdRegion[]> {
    return catchToLog<Promise<VdRegion[]>>(
      'MLB: DataLoaderTest.getLayoutBySchemeCdr()',
      async () => {
        const csv: string = await this._dlFiles.readAsText(DataLoaderTest.Files.layoutSchemeCdr);
        const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
        const schemeCol: DG.Column<string> = df.getCol('scheme');
        const cdrCol: DG.Column<string> = df.getCol('cdr');
        const dfRes: DG.DataFrame = df.clone(DG.BitSet.create(df.rowCount, (rowI: number) => {
          return schemeCol.get(rowI) == scheme && cdrCol.get(rowI) == cdr;
        }));
        // dfRes.changeColumnType('position_start_name', 'string');
        // dfRes.changeColumnType('position_end_name', 'string');
        return DataLoader.DataFrameToVdRegionList(dfRes);
      });
  }

  async getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: DataLoaderTest.getMlbByAntigen()',
      async () => {
        const csv: string = await this._dlFiles.readAsText(DataLoaderTest.Files.mlbByAntigen);
        const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
        const antigenCol: DG.Column<string> = df.getCol('antigen');
        const dfRes: DG.DataFrame = df.clone(DG.BitSet.create(df.rowCount, (rowI: number) => {
          return antigenCol.get(rowI) == antigen;
        }));
        return dfRes;
      });
  }

  getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: DataLoaderTest.getTreeByAntigen()',
      async () => {
        const csv: string = await this._dlFiles.readAsText(DataLoaderTest.Files.treeByAntigen);
        const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
        const antigenCol: DG.Column<string> = df.getCol('antigen');
        const dfRes: DG.DataFrame = df.clone(DG.BitSet.create(df.rowCount, (rowI) => {
          return antigenCol.get(rowI) == antigen;
        }));
        return dfRes;
      });
  }

  getAnarci(scheme: string, chain: string, antigen: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: DataLoaderTest.getAnarci()',
      async () => {
        const scheme2: string = scheme.charAt(0).toUpperCase() + scheme.slice(1);
        const chain2: string = chain.charAt(0).toUpperCase() + chain.slice(1);
        const file = `anarci${scheme2}${chain2}ByAntigen`;
        const csv: string = await this._dlFiles.readAsText(file);
        const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
        const antigenCol: DG.Column<string> = df.getCol('antigen');
        const dfRes: DG.DataFrame = df.clone(DG.BitSet.create(df.rowCount, (rowI) => {
          return antigenCol.get(rowI) == antigen;
        }));
        return dfRes;
      });
  }

  /** deprecated */
  loadMlbDf(): Promise<DG.DataFrame> {
    throw new Error('Not implemented');
  }

  /** deprecated */
  loadTreeDf(): Promise<DG.DataFrame> {
    throw new Error('Not implemented');
  }

  // -- PTM --

  getObservedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: DataLoaderTest.getObservedPtmByAntigen()',
      async () => {
        const csv: string = await this._dlFiles.readAsText(DataLoaderTest.Files.observedPtmByAntigen);
        const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
        const antigenCol: DG.Column<string> = df.getCol('antigen');
        const dfRes: DG.DataFrame = df.clone(DG.BitSet.create(df.rowCount, (rowI) => {
          return antigenCol.get(rowI) == antigen;
        }));
        return dfRes;
      });
  }

  getPredictedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: DataLoaderTest.getPredictedPtmByAntigen()',
      async () => {
        const csv: string = await this._dlFiles.readAsText(DataLoaderTest.Files.predictedPtmByAntigen);
        const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
        const antigenCol: DG.Column<string> = df.getCol('antigen');
        const dfRes: DG.DataFrame = df.clone(DG.BitSet.create(df.rowCount, (rowI) => {
          return antigenCol.get(rowI) == antigen;
        }));
        return dfRes;
      });
  }

  // -- 3D --

  private async loadFileJson(path: string): Promise<Object> {
    const jsonTxt = await grok.dapi.files.readAsText(`System:AppData/${packageName}/${path}`);
    return JSON.parse(jsonTxt);
  }

  load3D(vid: string): Promise<[JsonType, string, NumsType, ObsPtmType]> {
    return catchToLog<Promise<[JsonType, string, NumsType, ObsPtmType]>>(
      'MLB: DataLoaderTest.load3D()',
      async () => {
        return Promise.all([
          this._dlFiles.readAsText(DataLoaderTest.Files.load3dExample).then((value) => {
            return JSON.parse(value) as JsonType;
          }),
          this._dlFiles.readAsText(DataLoaderTest.Files.load3dExamplePdb).then((value) => {
            return JSON.parse(value)['pdb'];
          }),
          this._dlFiles.readAsText(DataLoaderTest.Files.load3dExampleRealNums).then((value) => {
            return JSON.parse(value) as NumsType;
          }),
          this._dlFiles.readAsText(DataLoaderTest.Files.load3dExampleOptm).then((value) => {
            return JSON.parse(value)['ptm_observed'] as ObsPtmType;
          }),
        ]);
      });
  }
}