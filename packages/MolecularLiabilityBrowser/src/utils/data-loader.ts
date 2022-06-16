import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

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

export interface MutcodesDataType {
}

export type PtmMapType = { [key: string]: string };
export type CdrMapType = { [key: string]: string };

/**
 * TODO: descriptive name
 */
export interface JsonType {
  heavy_seq: string;
  light_seq: string;
  /** Heavy chain name/key of data section in 'ptm_predictions', 'parapred_predictions'
   */
  heavy_chain: string;
  /** Light chain name/key of data section in 'ptm_predictions', 'parapred_predictions'
   */
  light_chain: string;
  ptm_predictions: any;
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

export function catchToLog<T>(prefix: string, func: () => T): T {
  try {
    const res: T = func();
    if (res instanceof Promise) {
      return res.catch((ex) => {
        console.error(prefix + ex.toString());
        throw (ex);
      }) as unknown as T;
    } else {
      return res;
    }
  } catch (ex: any) {
    console.error(prefix + ex.toString());
    throw (ex);
  }
}

export enum DataLoaderType {
  Unknown = 'unknown',
  Files = 'Files',
  Database = 'Database',
}

export abstract class DataLoader {
  /** Properties for filters
   */
  abstract get filterProperties(): FilterPropertiesType;

  abstract get mutcodes(): MutcodesDataType;

  abstract get ptmMap(): PtmMapType;

  abstract get cdrMap(): CdrMapType;

  abstract get refDf(): DG.DataFrame;

  abstract get realNums(): NumsType;

  abstract init();

  protected async check_files(files: { [name: string]: string }): Promise<void> {
    // for (const filePath of Object.values(files)) {
    //   if (!(await _package.files.exists(filePath)))
    //     fileErrors.push(filePath);
    // }
    const fileNames = Object.values(files);
    const fileErrors: string[] = (await Promise.all(
      fileNames.map((filePath) => { return _package.files.exists(filePath); })
    ).then((exists) => fileNames.filter((v, i) => !exists[i])));

    if (fileErrors.length > 0)
      throw new Error(`Files errors:\n ${fileErrors.join('\n')}`);
  }

  abstract getVids(): Promise<string[]>;

  abstract listAntigens(): Promise<DG.DataFrame>;

  abstract getMlbByAntigen(antigen: string): Promise<DG.DataFrame>;

  abstract getTreeByAntigen(antigen: string): Promise<DG.DataFrame>;

  abstract getAnarci(scheme: string, chain: string, antigen: string): Promise<DG.DataFrame>;

  abstract getObservedPtmVids(): Promise<string[]>;

  /**
   * Heavy chain calculated data
   */
  abstract loadHChainDf(): Promise<DG.DataFrame>;

  /**
   * Light chain calculated data
   */
  abstract loadLChainDf(): Promise<DG.DataFrame>;

  /**
   * TODO: Some description of data structure purpose
   */
  abstract loadMlbDf(): Promise<DG.DataFrame>;

  abstract loadTreeDf(): Promise<DG.DataFrame>;

  abstract loadExample(vid: string): Promise<JsonType>;

  /**
   */
  abstract loadPdb(vid: string): Promise<string>;

  /**
   * Get post observable translational modifications data for 'v_id'
   */
  abstract loadObsPtm(vid: string): Promise<ObsPtmType>;
}
