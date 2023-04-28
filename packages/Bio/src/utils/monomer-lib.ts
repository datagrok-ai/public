// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {Observable, Subject} from 'rxjs';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {
  createJsonMonomerLibFromSdf,
  expectedMonomerData,
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

export class MonomerLib implements IMonomerLib {
  private _monomers: { [type: string]: { [monomerName: string]: Monomer } } = {};
  private _onChanged = new Subject<any>();

  constructor(monomers: { [type: string]: { [monomerName: string]: Monomer } }) {
    this._monomers = monomers;
  }

  getMonomer(monomerType: string, monomerName: string): Monomer | null {
    if (monomerType in this._monomers! && monomerName in this._monomers![monomerType])
      return this._monomers![monomerType][monomerName];
    else
      return null;
  }

  getTypes(): string[] {
    return Object.keys(this._monomers);
  }

  getMonomerMolsByType(type: string): { [monomerName: string]: string } {
    const res: { [monomerName: string]: string } = {};

    Object.keys(this._monomers[type]).forEach((monomerName) => {
      res[monomerName] = this._monomers[type][monomerName].molfile;
    });

    return res;
  }

  getMonomerNamesByType(type: string): string[] {
    return Object.keys(this._monomers[type]);
  }

  get onChanged(): Observable<any> {
    return this._onChanged;
  }

  private _updateInt(lib: IMonomerLib): void {
    const typesNew = lib.getTypes();
    const types = this.getTypes();

    typesNew.forEach((type) => {
      //could possibly rewrite -> TODO: check duplicated monomer symbol

      if (!types.includes(type))
        this._monomers![type] = {};

      const monomers = lib.getMonomerNamesByType(type);
      monomers.forEach((monomerName) => {
        this._monomers[type][monomerName] = lib.getMonomer(type, monomerName)!;
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
      const userLibrariesSettings: string[] = Object.keys(await grok.dapi.userDataStorage.get(LIB_STORAGE_NAME, true));
      const libs: IMonomerLib[] = await Promise.all(userLibrariesSettings.map((libFileName) => {
        //TODO handle whether files are in place
        return this.readLibrary(LIB_PATH, libFileName);
      }));
      this._monomerLib.updateLibs(libs, reload);
    });
  }

  /** Reads library from file shares, handles .json and .sdf */
  async readLibrary(path: string, fileName: string): Promise<IMonomerLib> {
    let data: any[] = [];
    let file;
    let dfSdf;
    const fileSource = new DG.FileSource(path);
    if (fileName.endsWith('.sdf')) {
      const funcList: DG.Func[] = DG.Func.find({package: 'Chem', name: 'importSdf'});
      if (funcList.length === 1) {
        file = await fileSource.readAsBytes(fileName);
        dfSdf = await grok.functions.call('Chem:importSdf', {bytes: file});
        data = createJsonMonomerLibFromSdf(dfSdf[0]);
      } else {
        grok.shell.warning('Chem package is not installed');
      }
    } else {
      const file = await fileSource.readAsText(fileName);
      data = JSON.parse(file);
    }

    const monomers: { [type: string]: { [name: string]: Monomer } } = {};
    const types: string[] = [];
    data.forEach((monomer) => {
      const monomerAdd: {[key: string]: any} = {};
      for (const field of HELM_REQUIRED_FIELDS_ARRAY)
        monomerAdd[field] = monomer[field];
      for (const optField of HELM_OPTIONAL_FIELDS_ARRAY) {
        const value = monomer[optField];
        if (value !== undefined)
          monomerAdd[optField] = value;
      }

      Object.keys(monomer).forEach((prop) => {
        if (!expectedMonomerData.includes(prop))
          monomerAdd[prop] = monomer[prop];
      });

      if (!types.includes(monomer[REQ.POLYMER_TYPE])) {
        monomers[monomer[REQ.POLYMER_TYPE]] = {};
        types.push(monomer[REQ.POLYMER_TYPE]);
      }

      monomers[monomer[REQ.POLYMER_TYPE]][monomer[REQ.NAME]] = monomerAdd as Monomer;
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
