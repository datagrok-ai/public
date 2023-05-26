// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {Observable, Subject} from 'rxjs';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {
  createJsonMonomerLibFromSdf,
  IMonomerLibHelper
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {HELM_REQUIRED_FIELDS as REQ, HELM_OPTIONAL_FIELDS as OPT} from '@datagrok-libraries/bio/src/utils/const';

const HELM_REQUIRED_FIELDS_ARRAY = [
  REQ.SYMBOL, REQ.NAME, REQ.MOLFILE, REQ.AUTHOR, REQ.ID,
  REQ.RGROUPS, REQ.SMILES, REQ.POLYMER_TYPE, REQ.MONOMER_TYPE, REQ.CREATE_DATE
] as const;

const HELM_OPTIONAL_FIELDS_ARRAY = [OPT.NATURAL_ANALOG, OPT.META] as const;
// -- Monomer libraries --
export const LIB_STORAGE_NAME = 'Libraries';
export const LIB_PATH = 'System:AppData/Bio/libraries/';
export const LIB_DEFAULT: { [fileName: string]: string } = {'HELMCoreLibrary.json': 'HELMCoreLibrary.json'};

/** Type for user settings of monomer library set to use. */
export type LibSettings = {
  exclude: string[],
}

export async function getLibFileNameList(): Promise<string[]> {
  const res: string[] = (await grok.dapi.files.list(`${LIB_PATH}`, false, ''))
    .map((it) => it.fileName);
  return res;
}

export async function getUserLibSettings(): Promise<LibSettings> {
  const resStr: string = await grok.dapi.userDataStorage.getValue(LIB_STORAGE_NAME, 'Settings', true);
  const res: LibSettings = resStr ? JSON.parse(resStr) : {exclude: []};

  // Fix empty object returned in case there is no settings stored for user
  res.exclude = res.exclude instanceof Array ? res.exclude : [];

  return res;
}

export async function setUserLibSetting(value: LibSettings): Promise<void> {
  await grok.dapi.userDataStorage.postValue(LIB_STORAGE_NAME, 'Settings', JSON.stringify(value), true);
}

export class MonomerLib implements IMonomerLib {
  private _monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
  private _onChanged = new Subject<any>();

  constructor(monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } }) {
    this._monomers = monomers;
  }

  getMonomer(polymerType: string, monomerSymbol: string): Monomer | null {
    if (polymerType in this._monomers! && monomerSymbol in this._monomers![polymerType])
      return this._monomers![polymerType][monomerSymbol];
    else
      return null;
  }

  getPolymerTypes(): string[] {
    return Object.keys(this._monomers);
  }

  getMonomerMolsByPolymerType(polymerType: string): { [monomerSymbol: string]: string } {
    const res: { [monomerSymbol: string]: string } = {};

    Object.keys(this._monomers[polymerType]).forEach((monomerSymbol) => {
      res[monomerSymbol] = this._monomers[polymerType][monomerSymbol].molfile;
    });

    return res;
  }

  getMonomerSymbolsByType(polymerType: string): string[] {
    return Object.keys(this._monomers[polymerType]);
  }

  get onChanged(): Observable<any> {
    return this._onChanged;
  }

  private _updateInt(lib: IMonomerLib): void {
    const typesNew = lib.getPolymerTypes();
    const types = this.getPolymerTypes();

    typesNew.forEach((type) => {
      //could possibly rewrite -> TODO: check duplicated monomer symbol

      if (!types.includes(type))
        this._monomers![type] = {};

      const monomers = lib.getMonomerSymbolsByType(type);
      monomers.forEach((monomerSymbol) => {
        this._monomers[type][monomerSymbol] = lib.getMonomer(type, monomerSymbol)!;
      });
    });
  }

  public update(lib: IMonomerLib): void {
    this._updateInt(lib);
    this._onChanged.next();
  }

  public updateLibs(libList: IMonomerLib[], reload: boolean = false): void {
    if (reload) this._monomers = {};
    for (const lib of libList) this._updateInt(lib);
    this._onChanged.next();
  }

  public clear(): void {
    this._monomers = {};
    this._onChanged.next();
  }
}

export class MonomerLibHelper implements IMonomerLibHelper {
  private readonly _monomerLib: MonomerLib = new MonomerLib({});

  /** Protect constructor to prevent multiple instantiation. */
  protected constructor() {}

  /** Singleton monomer library */
  getBioLib(): IMonomerLib {
    return this._monomerLib;
  }

  private loadLibrariesPromise: Promise<void> = Promise.resolve();

  /** Loads libraries based on settings in user storage {@link LIB_STORAGE_NAME}
   * @param {boolean} reload Clean {@link monomerLib} before load libraries [false]
   */
  async loadLibraries(reload: boolean = false): Promise<void> {
    return this.loadLibrariesPromise = this.loadLibrariesPromise.then(async () => {
      const [libFileNameList, settings]: [string[], LibSettings] = await Promise.all([
        getLibFileNameList(),
        getUserLibSettings()
      ]);
      const libs: IMonomerLib[] = await Promise.all(libFileNameList
        .filter((libFileName) => !settings.exclude.includes(libFileName))
        .map((libFileName) => {
          //TODO handle whether files are in place
          return this.readLibrary(LIB_PATH, libFileName);
        }));
      this._monomerLib.updateLibs(libs, reload);
    });
  }

  /** Reads library from file shares, handles .json and .sdf */
  async readLibrary(path: string, fileName: string): Promise<IMonomerLib> {
    let rawLibData: any[] = [];
    let file;
    let dfSdf;
    const fileSource = new DG.FileSource(path);
    if (fileName.endsWith('.sdf')) {
      const funcList: DG.Func[] = DG.Func.find({package: 'Chem', name: 'importSdf'});
      if (funcList.length === 1) {
        file = await fileSource.readAsBytes(fileName);
        dfSdf = await grok.functions.call('Chem:importSdf', {bytes: file});
        rawLibData = createJsonMonomerLibFromSdf(dfSdf[0]);
      } else {
        grok.shell.warning('Chem package is not installed');
      }
    } else {
      const file = await fileSource.readAsText(fileName);
      rawLibData = JSON.parse(file);
    }

    const monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
    const polymerTypes: string[] = [];
    rawLibData.forEach((monomer) => {
      if (!polymerTypes.includes(monomer[REQ.POLYMER_TYPE])) {
        monomers[monomer[REQ.POLYMER_TYPE]] = {};
        polymerTypes.push(monomer[REQ.POLYMER_TYPE]);
      }
      monomers[monomer[REQ.POLYMER_TYPE]][monomer[REQ.SYMBOL]] = monomer as Monomer;
    });

    return new MonomerLib(monomers);
  }

  // -- Instance singleton --
  private static _instance: MonomerLibHelper | null = null;

  public static get instance(): MonomerLibHelper {
    if (!MonomerLibHelper._instance) MonomerLibHelper._instance = new MonomerLibHelper();
    return MonomerLibHelper._instance;
  }
}
