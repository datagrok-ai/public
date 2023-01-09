import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {VdRegion} from '@datagrok-libraries/bio/src/vd-regions';
import {_package, packageName} from '../package';
import {MlbDatabase} from './mlb-database';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

export type DataQueryDict = { [name: string]: DG.DataQuery };

export type PointType = [/** position */ number, /** probability */ number];

export interface FilterPropertiesType {
  source: string;
  names: string[];
  yellowLeft: number[];
  yellowRight: number[];
  redLeft: number[];
  redRight: number[];
  plotsX: number[][];
  plotsY: number[][];
  tooltips: string[];
}

export type MutcodesDataType = { [name: string]: string };

export type PtmMapType = { [key: string]: string };
export type CdrMapType = { [key: string]: string };

/**
 * TODO: descriptive name
 */
export interface JsonType {
  heavy_seq: string;
  light_seq: string;

  /** Heavy chain name/key of data section in 'ptm_predictions', 'parapred_predictions' */
  heavy_chain: string;
  /** Light chain name/key of data section in 'ptm_predictions', 'parapred_predictions' */
  light_chain: string;

  ptm_predictions: { [chain: string]: { [ptm: string]: PointType[] } };
  parapred_predictions: { [chain: string]: { [pos: string]: number } };

  cdr_ranges: { [cdr: string]: PointType[] };

  map_H: number[],
  map_L: number[],
}

export interface PdbType {
}

/**
 *
 */
export interface NumsType {
  heavy_numbering: string;
  light_numbering: string;
}

export interface ObsPtmType {
  H: any;
  L: any;
}

export interface IFeature {
  category: string;
  type: string;
  start: number;
  end: number;
  text: string;
}

export interface ColorSchemeType {
  colBackground: string;
  colHeavyChain: string;
  colLightChain: string;
  colCdr: string;
  colPara: string;
  colHighlight: string;
  colHighlightCdr: string;
  colParatopesLow: string;
  colParatopesHigh: string;
}

export interface PvizType {
  SeqEntry: { new(args: { sequence: string }): SeqEntry };
  DASReader: Function;
  FastaReader: Function;
  FeatureManager: {};
  IconFactory: {};
  SeqEntryAnnotInteractiveView: any;
  SeqEntryFastaView: {};
  FeatureDisplayer: FeatureDisplayer;
  OneLiner: Function;
}

export interface PvizParamType {
  seq: { [chain: string]: string };
  cdrMap: { [chain: string]: PvizCdrType };
  parMap: { [chain: string]: PvizParType };
  denMap: { [chain: string]: PvizDenType };
  ptmMap: { [chain: string]: PvizPtmType };
  motMap: { [chain: string]: PvizMotType };
  obsMap: { [chain: string]: PvizObsType };
}

export interface PvizCdrType {
  cdrFeatureMap: CdrFeatureType[];
}

export interface CdrFeatureType extends IFeature {
  improbable: true;
}

export interface PvizParType {
  parFeatureMap: ParFeatureType[];
  parColorObj: { P: string[] };
  parElObj: string[];
  parProbObj: number[];
}

export interface ParFeatureType extends IFeature {
  improbable: boolean;
}

export interface PvizDenType {
  denFeatureMap: DenFeatureType [];
  denColorObj: { D: string[] };
  denElObj: number[];
  denProbObj: number[];
  denPtmArr: [string, number][][];
}

export interface DenFeatureType extends IFeature {
  groupSet: string;
  improbable: boolean;
}

export class PvizPtmType {
  ptmFeatureMap: PtmFeatureType [];
  ptmColorObj: { [ptm: string]: string[] };
  ptmElObj: { [ptm: string]: number[] };
  ptmProbObj: Object;
}

export interface PtmFeatureType extends IFeature {
  groupSet: string;
  improbable: boolean;
}

export interface PvizMotType {
  motFeatureMap: MotFeatureType[];
  motColorObj: { [code: string]: string[] };
  motElObj: { [code: string]: number[] };
  motProbObj: Object;
}

export interface MotFeatureType extends IFeature {

}

export interface PvizObsType {
  obsFeatureMap: ObsFeatureType[];
  obsColorObj: { [ptm: string]: string[] };
  obsElObj: { [ptm: string]: number[] };
  obsProbObj: [string[], number[]][];
}

export interface ObsFeatureType extends IFeature {
  groupSet: string;
  improbable: boolean;
}

export interface SeqEntry {
  addFeatures: Function;
}

export interface FeatureDisplayer {
  addMouseoverCallback: Function;
  addClickCallback: Function;
}

export interface MotMapType {
  motFeatureMap: string[];
  motColorObj: {};
  motElObj: Object;
  motProbObj: Object;
}

/** Handles error and console.debug ET (elapsed time) */
export function catchToLog<T>(prefix: string, func: () => T): T {
  const t1: number = Date.now();
  try {
    const res: T = func();

    if (res instanceof Promise) {
      return res.catch((ex) => {
        console.error(prefix + ', ' + ex.toString());
        throw (ex);
      }).then((obj) => {
        const t2: number = Date.now();
        console.debug(prefix + `, ET: ${((t2 - t1) / 1000)} s`);

        return obj;
      }) as unknown as T;
    } else {
      const t2: number = Date.now();
      console.debug(prefix + `, ET: ${((t2 - t1) / 1000)} s`);

      return res;
    }
  } catch (err: any) {
    const errStr = errorToConsole(err);
    console.error(prefix + ', ' + errStr);
    throw (err);
  }
}

export async function initPackageQueries(packageName: string): Promise<DataQueryDict> {
  const queryList = await catchToLog<Promise<DG.DataQuery[]>>(
    `initPackageQueries( "${packageName}" )`,
    () => {
      const res = grok.dapi.queries
        .include('params,connection')
        .filter(`package.name = "${packageName}"`).list();
      return res;
    });

  const queries: DataQueryDict = Object.assign({}, ...queryList.map((q: DG.DataQuery) => {
    const qName = q.name.charAt(0).toLowerCase() + q.name.slice(1);
    return {[qName]: q};
  }));

  return queries;
}

export enum DataLoaderType {
  Unknown = 'unknown',
  Files = 'Files',
  Database = 'Database',
  Test = 'Test',
}

type FileNameFunc = (name: string) => string;

export class FilesForDataLoader {
  static Files = class {
    static filterProps = 'filterProps';
    static mutcodes = 'mutcodes';
    static predictedPtmMap = 'predictedPtmMap';
    static predictedCdrMap = 'predictedCdrMap';
    static observedPtmMap = 'observedPtmMap';
    static observedCdrMap = 'observedCdrMap';
  };

  static _files = {
    [FilesForDataLoader.Files.filterProps]: 'properties.json',
    [FilesForDataLoader.Files.mutcodes]: 'mutcodes.json',
    [FilesForDataLoader.Files.predictedPtmMap]: 'ptm_map.json',
    [FilesForDataLoader.Files.predictedCdrMap]: 'cdr_map.json',
    [FilesForDataLoader.Files.observedPtmMap]: 'obs_ptm_map.json',
    [FilesForDataLoader.Files.observedCdrMap]: 'obs_cdr_map.json',
    // [FilesForDataLoader.Files.predictedPtm] : 'ptm_in_cdr.d42',
  };
  private readonly _cache: MlbDatabase;
  private readonly _fnFunc: FileNameFunc;

  /** Default function to get file name of the file */
  static fnFuncDefault(file: string): string {
    return FilesForDataLoader._files[file];
  }

  constructor(cache: MlbDatabase, fnFunc: FileNameFunc = FilesForDataLoader.fnFuncDefault) {
    this._cache = cache;
    this._fnFunc = fnFunc;
  }

  async readAsText(file: string) {
    const txt: string = await grok.dapi.files
      .readAsText(`System:AppData/${packageName}/${this._fnFunc(file)}`);
    return txt;
  }

  async getFilterProperties(): Promise<FilterPropertiesType> {
    return this._cache.getObject<FilterPropertiesType>(
      FilesForDataLoader.Files.filterProps,
      () => catchToLog<Promise<FilterPropertiesType>>(
        'MLB: FilesForDataLoader.getFilterProperties()',
        async () => {
          const txt: string = await this.readAsText(FilesForDataLoader.Files.filterProps);
          return JSON.parse(txt);
        }));
    // .then((value: FilterPropertiesType) => {
    //   console.debug(`MLB: DataLoaderDb.init2() set filterProperties, ${this.fromStartInit()} s`);
    //   this._filterProperties = value;
    // })
  }

  async getMutcodes(): Promise<MutcodesDataType> {
    return this._cache.getObject<MutcodesDataType>(
      FilesForDataLoader.Files.mutcodes,
      () => catchToLog<Promise<MutcodesDataType>>(
        'MLB: FilesForDataLoader.getMutcodes()',
        async () => {
          const txt: string = await this.readAsText(FilesForDataLoader.Files.mutcodes);
          return JSON.parse(txt);
        }));
    // .then((value: MutcodesDataType) => {
    //   console.debug(`MLB: DataLoaderDb.init2() set mutcodes, ${this.fromStartInit()} s`);
    //   this._mutcodes = value;
    // })
  }

  async getPredictedPtmMap(): Promise<PtmMapType> {
    return this._cache.getObject<PtmMapType>(
      FilesForDataLoader.Files.predictedPtmMap,
      () => catchToLog<Promise<PtmMapType>>(
        'MLB: FilesForDataLoader.getPredictedPtmMap()',
        async () => {
          const jsonTxt: string = await this.readAsText(FilesForDataLoader.Files.predictedPtmMap);
          return JSON.parse(jsonTxt);
        }));
    // .then((value: PtmMapType) => {
    //   console.debug(`MLB: DataLoaderDb.init2() set predictedPtmMap, ${this.fromStartInit()} s`);
    //   this._predictedPtmMap = value;
    // })
  }

  async getPredictedCdrMap(): Promise<CdrMapType> {
    return this._cache.getObject<CdrMapType>(
      FilesForDataLoader.Files.predictedCdrMap,
      () => catchToLog<Promise<CdrMapType>>(
        'MLB: FilesForDataLoader.getPredictedCdrMap()',
        async () => {
          const jsonTxt: string = await this.readAsText(FilesForDataLoader.Files.predictedCdrMap);
          return JSON.parse(jsonTxt);
        }));
    // .then((value: CdrMapType) => {
    //   console.debug(`MLB: DataLoaderDb.init2() set predictedCdrMap, ${this.fromStartInit()} s`);
    //   this._predictedCdrMap = value;
    // }),
  }

  async getObservedPtmMap(): Promise<PtmMapType> {
    return this._cache.getObject<PtmMapType>(
      FilesForDataLoader.Files.observedPtmMap,
      () => catchToLog<Promise<PtmMapType>>(
        'MLB: FilesForDataLoader.getObservedPtmMap()',
        async () => {
          const jsonTxt: string = await this.readAsText(FilesForDataLoader.Files.observedPtmMap);
          return JSON.parse(jsonTxt);
        }));
    // .then((value: PtmMapType) => {
    //   console.debug(`MLB: DataLoaderDb.init2() set observedPtmMap, ${this.fromStartInit()} s`);
    //   this._observedPtmMap = value;
    // })
  }

  async getObservedCdrMap(): Promise<CdrMapType> {
    return this._cache.getObject<CdrMapType>(
      FilesForDataLoader.Files.observedCdrMap,
      () => catchToLog<Promise<CdrMapType>>(
        'MLB: FilesForDataLoader.getObservedCdrMap()',
        async () => {
          const jsonTxt: string = await this.readAsText(FilesForDataLoader.Files.observedCdrMap);
          return JSON.parse(jsonTxt);
        }));
    // .then((value: CdrMapType) => {
    //   console.debug(`MLB: DataLoaderDb.init2() set observedCdrMap, ${this.fromStartInit()} s`);
    //   this._observedCdrMap = value;
    // })
  }
}

export class QueriesForDataLoader {
  private _mlbQueries: DataQueryDict;
  private _cache: MlbDatabase;

  constructor(cache: MlbDatabase, mlbQueries: DataQueryDict) {
    this._mlbQueries = mlbQueries;

    this._cache = cache;
  }

  // -- Lists --

  async listSchemes(): Promise<string[]> {
    // load numbering schemes
    const df: DG.DataFrame = await this._cache.getDataFrame('scheme',
      () => catchToLog<Promise<DG.DataFrame>>(
        'MLB: QueriesForDataLoader.listSchemes()',
        async () => {
          const funcCall: DG.FuncCall = await this._mlbQueries['listSchemes'].prepare().call();
          const df: DG.DataFrame = funcCall.getOutputParamValue();
          return df;
        }));
    return df.getCol('scheme').toList();
  }

  async listCdrs(): Promise<string[]> {
    //load cdr definition list
    const df: DG.DataFrame = await this._cache.getDataFrame('cdr',
      () => catchToLog<Promise<DG.DataFrame>>(
        'MLB: QueriesForDataLoader.listCdrs()',
        async () => {
          const funcCall: DG.FuncCall = await this._mlbQueries['listCdrs'].prepare().call();
          const df: DG.DataFrame = funcCall.getOutputParamValue();
          return df;
        }));
    return df.getCol('cdr').toList();
  };

  async listAntigens(): Promise<DG.DataFrame> {
    const df: DG.DataFrame = await this._cache.getDataFrame('antigen',
      () => catchToLog<Promise<DG.DataFrame>>(
        'MLB: QueriesForDataLoader.listAntigens()',
        async () => {
          const funcCall: DG.FuncCall = await this._mlbQueries['listAntigens'].prepare().call();
          const df: DG.DataFrame = funcCall.getOutputParamValue();
          return df;
        }));
    return df;
  }

  async getVids(): Promise<string[]> {
    const df: DG.DataFrame = await this._cache.getDataFrame('vid',
      () => catchToLog<Promise<DG.DataFrame>>(
        'MLB: QueriesForDataLoader.getVids()',
        async () => {
          const funcCall: DG.FuncCall = await this._mlbQueries['getVids'].prepare().call();
          const df: DG.DataFrame = funcCall.getOutputParamValue();
          return df;
        }));
    return df.getCol('v_id').toList();
  }

  async getObservedPtmVids(): Promise<string[]> {
    const df: DG.DataFrame = await this._cache.getDataFrame('vidObsPtm',
      () => catchToLog<Promise<DG.DataFrame>>(
        'MLB: QueriesForDataLoader.getObservedPtmVids()',
        async () => {
          const funcCall: DG.FuncCall = await this._mlbQueries['getObservedPtmVids'].prepare().call();
          const df: DG.DataFrame = funcCall.getOutputParamValue();
          return df;
        }));
    return df.getCol('v_id').toList();
  }

  // -- --

  async getLayoutBySchemeCdr(scheme: string, cdr: string): Promise<VdRegion[]> {
    return catchToLog<Promise<VdRegion[]>>(
      'MLB: QueriesForDataLoader.getLayoutBySchemeCdr()',
      async () => {
        // const df: DG.DataFrame = await grok.functions
        //   .call(`${this._pName}:getLayoutBySchemeCdr`, {scheme: scheme, cdr: cdr}) as DG.DataFrame;
        const funcCall: DG.FuncCall = await this._mlbQueries['getLayoutBySchemeCdr']
          .prepare({scheme: scheme, cdr: cdr}).call();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        return DataLoader.dfToVdRegionList(df);
      });
  }

  async getMlbByAntigen(antigen: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: QueriesForDataLoader.getMlbByAntigen()',
      async () => {
        // const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getMlbByAntigen`, {antigen: antigen});
        const funcCall: DG.FuncCall = await this._mlbQueries['getMlbByAntigen']
          .prepare({antigen: antigen}).call();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        return df;
      });
  }

  async getTreeByAntigen(antigen: string): Promise<DG.DataFrame> {
    // const df: DG.DataFrame = await grok.functions.call(`${this._pName}:getTreeByAntigen`, {antigen: antigen});
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: QueriesForDataLoader.getTreeByAntigen()',
      async () => {
        const funcCall: DG.FuncCall = await this._mlbQueries['getTreeByAntigen']
          .prepare({antigen: antigen}).call();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        return df;
      });
  }

  async getAnarci(scheme: string, chain: string, antigen: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: QueriesForDataLoader.getAnarci()',
      async () => {
        // There is a problem with using underscore symbols in query names.
        const scheme2: string = scheme.charAt(0).toUpperCase() + scheme.slice(1);
        const chain2: string = chain.charAt(0).toUpperCase() + chain.slice(1);
        // const df: DG.DataFrame = await grok.functions
        //   .call(`${packageName}:getAnarci${scheme2}${chain2}`, {antigen: antigen}) as DG.DataFrame;
        const funcCall: DG.FuncCall = await this._mlbQueries[`getAnarci${scheme2}${chain2}`]
          .prepare({antigen: antigen}).call();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        return df;
      });
  }

  async loadMlbDf(): Promise<DG.DataFrame> {
    // const df = await grok.functions.call(`${this._pName}:GetMolecularLiabilityBrowser`);
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: QueriesForDataLoader.loadMlbDf() getMolecularLiabilityBrowser',
      async () => {
        const funcCall: DG.FuncCall = await this._mlbQueries['getMolecularLiabilityBrowser'].prepare();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        // 'ngl' column have been removed from query 2022-04
        df.columns.remove('ngl');
        return df;
      });
  }

  async getPredictedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    // const df = await grok.functions.call(`${this._pName}:getPredictedPtmByAntigen`, {antigen: antigen});
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: QueriesForDataLoader.getPredictedPtmByAntigen()',
      async () => {
        const funcCall: DG.FuncCall = await this._mlbQueries['getPredictedPtmByAntigen']
          .prepare({antigen: antigen}).call();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        return df;
      });
  }

  async getObservedPtmByAntigen(antigen: string): Promise<DG.DataFrame> {
    // const df = await grok.functions.call(`${this._pName}:getObservedPtmByAntigen`, {antigen: antigen});
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: QueriesForDataLoader.getObservedPtmByAntigen()',
      async () => {
        const funcCall: DG.FuncCall = await this._mlbQueries['getObservedPtmByAntigen']
          .prepare({antigen: antigen}).call();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        return df;
      });
  }

  async load3D(vid: string): Promise<DG.DataFrame> {
    return catchToLog<Promise<DG.DataFrame>>(
      'MLB: QueriesForDataLoader.get3D()', async () => {
        const funcCall: DG.FuncCall = await this._mlbQueries['get3D'].prepare({vid: vid}).call();
        const df: DG.DataFrame = funcCall.getOutputParamValue();
        return df;
      });
  }
}

export abstract class DataLoader {
  abstract get schemes(): string[];

  abstract get cdrs(): string[];

  abstract get antigens(): DG.DataFrame;

  abstract get vids(): string[];

  abstract get vidsObsPtm(): string[];

  /** Properties for filters
   */
  abstract get filterProperties(): FilterPropertiesType;

  abstract get mutcodes(): MutcodesDataType;

  abstract get predictedPtmMap(): PtmMapType;

  abstract get predictedCdrMap(): CdrMapType;

  abstract get observedPtmMap(): PtmMapType;

  abstract get observedCdrMap(): CdrMapType;

  abstract init(startInit: number): Promise<void>;

  protected async check_files(files: { [name: string]: string }): Promise<void> {
    // for (const filePath of Object.values(files)) {
    //   if (!(await _package.files.exists(filePath)))
    //     fileErrors.push(filePath);
    // }
    const fileNames = Object.values(files);
    const fileErrors: string[] = (await Promise.all(
      fileNames.map((filePath) => { return grok.dapi.files.exists(`System:AppData/${packageName}/${filePath}`); })
    ).then((exists) => fileNames.filter((v, i) => !exists[i])));

    if (fileErrors.length > 0)
      throw new Error(`Files errors:\n ${fileErrors.join('\n')}`);
  }

  abstract getLayoutBySchemeCdr(scheme: string, cdr: string): Promise<VdRegion[]>;

  abstract getMlbByAntigen(antigen: string): Promise<DG.DataFrame>;

  abstract getTreeByAntigen(antigen: string): Promise<DG.DataFrame>;

  abstract getAnarci(scheme: string, chain: string, antigen: string): Promise<DG.DataFrame>;

  // deprecated to load all data at once
  // /**
  //  * Heavy chain calculated data
  //  */
  // abstract loadHChainDf(): Promise<DG.DataFrame>;
  //
  // /**
  //  * Light chain calculated data
  //  */
  // abstract loadLChainDf(): Promise<DG.DataFrame>;

  /** deprecated */
  abstract loadMlbDf(): Promise<DG.DataFrame>;

  /** deprecated */
  abstract loadTreeDf(): Promise<DG.DataFrame>;

  // -- PTM --
  abstract getPredictedPtmByAntigen(antigen: string): Promise<DG.DataFrame>;

  abstract getObservedPtmByAntigen(antigen: string): Promise<DG.DataFrame>;

  // -- 3D --

  // abstract loadJson(vid: string): Promise<JsonType>;
  //
  // abstract loadPdb(vid: string): Promise<string>;
  //
  // abstract loadRealNums(vid: string): Promise<NumsType>;
  //
  // /** Get post observable translational modifications data for 'v_id' */
  // abstract loadObsPtm(vid: string): Promise<ObsPtmType>;

  abstract load3D(vid: string): Promise<[JsonType, string, NumsType, ObsPtmType]>;

  // -- Routines --

  static dfToVdRegionList(df: DG.DataFrame): VdRegion[] {
    const regionListRes: VdRegion[] = [];
    for (const dfRow of df.rows) {
      const region = new VdRegion(
        dfRow.get('type'), dfRow.get('name'), dfRow.get('chain'), dfRow.get('order'),
        dfRow.get('position_start_name').toString(), dfRow.get('position_end_name').toString());
      regionListRes.push(region);
    }
    return regionListRes;
  }
}
