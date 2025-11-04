/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {HelmAtom, HelmType, IMonomerColors, IWebEditorMonomer,
  MonomerSetType, MonomerType, PolymerType} from '../helm/types';
import {Observable} from 'rxjs';
import {
  HELM_REQUIRED_FIELD as REQ,
  HELM_RGROUP_FIELDS as RGP, HELM_OPTIONAL_FIELDS as OPT,
} from '../utils/const';
import {errInfo} from '../utils/err-info';

export type RGroup = {
  [RGP.CAP_GROUP_SMILES]: string,
  [RGP.ALTERNATE_ID]: string,
  [RGP.CAP_GROUP_NAME]: string,
  [RGP.LABEL]: string,
}

/** https://github.com/PistoiaHELM/HELMMonomerSets/blob/master/HELMmonomerSchema.json */
export type Monomer = {
  [REQ.SYMBOL]: string,
  [REQ.NAME]: string,
  [REQ.MOLFILE]: string,
  [REQ.AUTHOR]: string,
  [REQ.ID]: number,
  [REQ.RGROUPS]: RGroup[],
  [REQ.SMILES]: string,
  [REQ.POLYMER_TYPE]: PolymerType,
  [REQ.MONOMER_TYPE]: MonomerType,
  [REQ.CREATE_DATE]: string | null,
  [OPT.NATURAL_ANALOG]?: string,
  [OPT.META]?: { [property: string]: any },

  lib?: IMonomerLibBase,
  wem?: IWebEditorMonomer,
};

export interface IMonomerLibBase {
  get source(): string;
  get isEmpty(): boolean;
  get onChanged(): Observable<any>;

  getMonomerSymbolsByType(polymerType: PolymerType): string[];

  /* Adds placeholder empty missing monomer to lib */
  addMissingMonomer(polymerType: PolymerType, monomerSymbol: string): Monomer;

  getMonomer(polymerType: PolymerType | null, monomerSymbol: string): Monomer | null;

  /** HELMWebEditor expects null for HelmTypes.LINKER and R-Group count != 2 */
  getWebEditorMonomer(a: HelmAtom | HelmType, symbol?: string): IWebEditorMonomer | null;

  /** Get R groups from smiles and returns as object. CCC[OH:1] would return {[R1]: 'OH'} */
  getRS(smiles: string): { [r: string]: string };

  getTooltip(biotype: HelmType, monomerSymbol: string): HTMLElement;

  getMonomerColors(biotype: HelmType, symbol: string): IMonomerColors;

  getMonomerTextColor(biotype: HelmType, symbol: string): string;
}


export interface IMonomerSetPlaceholder {
  get symbol(): string;
  get polymerType(): PolymerType;
  get monomerType(): MonomerType;
  get monomerLinks(): IMonomerLinkData[];

  get monomers(): Monomer[];
}

export interface IMonomerSet {
  get source(): string | undefined;
  get error(): string | undefined;

  get description(): string;
  get placeholders(): IMonomerSetPlaceholder[];
}

export type MonomerLibSummaryType = { [polymerType: string]: number };

export type MonomerLibData = { [polymerType: string]: { [symbol: string]: Monomer } };

export interface IMonomerLib extends IMonomerLibBase {
  get error(): string | undefined;
  get duplicateMonomers(): { [polymerType: string]: { [monomerSymbol: string]: Monomer[] } };

  getMonomerMolsByPolymerType(polymerType: PolymerType): { [monomerSymbol: string]: string } | null;
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: PolymerType, element?: string): string[];
  getPolymerTypes(): PolymerType[];
  toJSON(): Monomer[];
  toString(): string;

  /** Summary string with lib monomer count by type
   * @deprecated Keep for backward compatibility */
  getSummary(): string;

  /** Summary with lib monomer count by type */
  getSummaryObj(): MonomerLibSummaryType;

  /** Gets dataframe with columns 'polymerType', 'count'. */
  getSummaryDf(): DG.DataFrame;

  // For monomer palettes
  getMonomerSet(biotype: HelmType): MonomerSetType | null;

  override(overrideData: MonomerLibData, source: string): IMonomerLibBase;
}

export interface IMonomerLinkData {
  source: string;
  symbol: string;
  notes: string;
}

export interface IMonomerLink {
  get lib(): IMonomerLib;
  get symbol(): string;
}

export const DEFAULT_FILES_LIB_PROVIDER_NAME = 'Files';
export interface IMonomerLibProvider {
    name: string;
    listLibraries(): Promise<string[]>;
    refreshLists(): Promise<void>;
    get loadPromise(): Promise<void>;
    loadLibraries(names: string[]): Promise<MonomerLibData[]>;
    deleteLibrary(name: string): Promise<void>;
    addOrUpdateLibraryString(name: string, contentString: string): Promise<void>;
    addOrUpdateLibrary(libraryName: string, monomers: Monomer[]): Promise<void>;
    updateOrAddMonomersInLibrary(libraryName: string, monomers: Monomer[]): Promise<void>;
    getSingleLibrary(name: string): Promise<MonomerLibData | null>;
    getLibraryAsString(libName: string): Promise<string>;
    listSets(): Promise<string[]>;
    loadSets(names: string[]): Promise<IMonomerSet[]>;
    deleteSet(name: string): Promise<void>;
    addOrUpdateSetString(name: string, contentString: string): Promise<void>;
    addOrUpdateSet(setName: string, monomerSet: IMonomerSet): Promise<void>;
    deleteMonomersFromLibrary(libraryName: string, monomers: ({polymerType: PolymerType, symbol: string}[])): Promise<void>;
    get onChanged(): Observable<void>;
}

export async function findProviderWithLibraryName(providers: IMonomerLibProvider[], libraryName: string) {
  for (const provider of providers) {
    const libNames = await provider.listLibraries();
    if (libNames.includes(libraryName) || libNames.includes(`${libraryName}.json`))
      return provider;
  }
  return null;
}

export interface IMonomerLibHelper {

  get providersDataChanged(): Observable<void>;

  /** Returns list of monomer library providers. can be files, DBs, etc */
  getProviders(): Promise<IMonomerLibProvider[]>;

  /** Ensures files are loaded and validated, throws error after timeout */
  awaitLoaded(timeout?: number): Promise<void>;

  /** Singleton monomer library collected from various sources */
  getMonomerLib(): IMonomerLib;

  getAvaliableLibraryNames(): Promise<string[]>;

  refreshValidLibraryLists(): Promise<void>;

  getAvailableLibrariesPerProvider(): Promise<Map<string, string[]>>;

  /**  @deprecated Use {@link getMonomerLib} */
  getBioLib(): IMonomerLib;

  /** Singleton monomer set collected from various sources */
  getMonomerSets(): IMonomerSet;

  /** (Re)Loads libraries based on settings in user storage {@link LIB_STORAGE_NAME} to singleton.
   * @param {boolean} reload Clean {@link monomerLib} before load libraries [false]
   */
  loadMonomerLib(reload?: boolean): Promise<void>;

  /** @deprecated Use {@link loadMonomerLib} */
  loadLibraries(reload?: boolean): Promise<void>;

  /** (Re)loads monomer sets based on settings in user storage {@link SETS_STORAGE_NAME} to singleton.
   * @param {boolean} reload Clean {@link monomerSets} before load sets [false]
   */
  loadMonomerSets(reload?: boolean): Promise<void>;

  /** Read single library from given provider */
  readSingleLibrary(providerName: string, libraryName: string): Promise<IMonomerLib | null>;

  /** Reads single library by finding suitable provider itself. Condition must be met that there are no overlapping library names */
  readSingleLibraryByName(libraryName: string): Promise<IMonomerLib | null>;

  readSingleSet(providerName: string, setName: string): Promise<IMonomerSet | null>;
  // -- Settings --
  refreshValidLibraryLists(): Promise<void>;
  /** Changes user lib settings. */
  loadMonomerLibForTests(): Promise<void>;

  get fileUploadRequested(): Observable<void>;
  requestFileUpload(): void;
  get librarySelectionChanged(): Observable<void>;
  notifyLibrarySelectionChanged(): void;

  /** Use default file provider to read monomer library. just a shorthand */
  readLibraryFromFilePath(filePath: string): Promise<IMonomerLib>

}


export async function getMonomerLibHelper(): Promise<IMonomerLibHelper> {
  const funcList = DG.Func.find({package: 'Bio', name: 'getMonomerLibHelper'});
  if (funcList.length === 0)
    throw new Error('Package "Bio" must be installed for MonomerLibHelper.');

  const res: IMonomerLibHelper = (await funcList[0].prepare().call()).getOutputParamValue() as IMonomerLibHelper;
  return res;
}


export class MonomerSetPlaceholder implements IMonomerSetPlaceholder {
  public readonly monomers: Monomer[];

  public readonly error: string | null = null;

  constructor(
    private readonly monomerLib: IMonomerLib,
    public symbol: string,
    public readonly polymerType: PolymerType,
    public readonly monomerType: MonomerType,
    public readonly monomerLinks: IMonomerLinkData[],
  ) {
    try {
      this.monomers = this.monomerLinks.map((mLink) => {
        const resM = this.monomerLib.getMonomer(this.polymerType, mLink.symbol);
        if (!resM)
          throw new Error('Monomer not found: ');
        if (resM.lib?.source != mLink.source)
          throw new Error(`Monomer '${symbol}' found in different library.`);
        return resM;
      });
    } catch (err: any) {
      const [errMsg, _errStack] = errInfo(err);
      this.error = errMsg;
      this.monomers = [];
    }
  }
}

export class MonomerSet implements IMonomerSet {
  public constructor(
    public readonly description: string,
    public placeholders: IMonomerSetPlaceholder[],
    public readonly source: string | undefined = undefined,
    public readonly error: string | undefined = undefined,
  ) {}

  updateSets(setList: IMonomerSet[], reload: boolean = false): void {
    if (reload)
      this.placeholders = [];
    for (const _set of setList)
      if (!_set.error) this._updateSetInt(_set);

    // TODO: File onChanged
  }

  private _updateSetInt(set: IMonomerSet): void {
    for (const setPh of set.placeholders)
      this.placeholders.push(setPh);
  }
}
