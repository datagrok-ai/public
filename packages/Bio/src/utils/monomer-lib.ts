import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable, Subject} from 'rxjs';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {
  LibSettings, getUserLibSettings, setUserLibSettings, LIB_PATH
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {PolyToolMonomerLibHandler} from '@datagrok-libraries/bio/src/utils/poly-tool/monomer-lib-handler';
import {
  createJsonMonomerLibFromSdf,
  IMonomerLibHelper,
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  HELM_REQUIRED_FIELDS as REQ, HELM_OPTIONAL_FIELDS as OPT, HELM_POLYMER_TYPE
} from '@datagrok-libraries/bio/src/utils/const';

import {_package} from '../package';


export async function getLibFileNameList(): Promise<string[]> {
  // list files recursively because permissions are available for folders only
  const res: string[] = await Promise.all((await grok.dapi.files.list(LIB_PATH, true, ''))
    .map(async (it) => {
      // Get relative path (to LIB_PATH)
      return it.fullPath.substring(LIB_PATH.length);
    }));
  return res;
}

export async function manageFiles() {
  const a = ui.dialog({title: 'Manage files'})
    //@ts-ignore
    .add(ui.fileBrowser({path: 'System:AppData/Bio/libraries'}).root)
    .addButton('OK', () => a.close())
    .show();
}

export async function getLibraryPanelUI(): Promise<DG.Widget> {
  //@ts-ignore
  const filesButton: HTMLButtonElement = ui.button('Manage', manageFiles);
  const inputsForm: HTMLDivElement = ui.inputs([]);
  const libFileNameList: string[] = await getLibFileNameList();

  const settings = await getUserLibSettings();

  for (const libFileName of libFileNameList) {
    const libInput: DG.InputBase<boolean | null> = ui.boolInput(libFileName, !settings.exclude.includes(libFileName),
      () => {
        if (libInput.value == true) {
          // Checked library remove from excluded list
          settings.exclude = settings.exclude.filter((l) => l != libFileName);
        } else {
          // Unchecked library add to excluded list
          if (!settings.exclude.includes(libFileName)) settings.exclude.push(libFileName);
        }
        setUserLibSettings(settings).then(async () => {
          await MonomerLibHelper.instance.loadLibraries(true); // from libraryPanel()
          grok.shell.info('Monomer library user settings saved.');
        });
      });
    inputsForm.append(libInput.root);
  }
  return new DG.Widget(ui.divV([inputsForm, ui.div(filesButton)]));
}

export class MonomerLib implements IMonomerLib {
  public readonly error: string | undefined;

  private _monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
  private _onChanged = new Subject<any>();

  constructor(monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } }, error?: string) {
    this._monomers = monomers;
    this.error = error;
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

    Object.keys(this._monomers[polymerType] ?? {}).forEach((monomerSymbol) => {
      res[monomerSymbol] = this._monomers[polymerType][monomerSymbol].molfile;
    });

    return res;
  }

  getMonomerSymbolsByType(polymerType: string): string[] {
    return Object.keys(this._monomers[polymerType]);
  }

  /** Get a list of monomers with specified element attached to specified
   * R-group
   * WARNING: RGroup numbering starts from 1, not 0*/
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: string, element?: string): string[] {
    const monomerSymbols = this.getMonomerSymbolsByType(polymerType);
    let monomers = monomerSymbols.map((sym) => this.getMonomer(polymerType, sym));
    monomers = monomers.filter((el) => el !== null);
    if (monomers.length === 0)
      return [];

    function findAllIndices<T>(arr: T[], element: T): number[] {
      return arr.map((value, index) => (value === element ? index : -1))
        .filter((index) => index !== -1);
    }

    monomers = monomers.filter((monomer) => {
      if (!monomer?.rgroups)
        return false;
      let criterion = monomer?.rgroups.length >= rGroupNumber;
      const molfileHandler = MolfileHandler.getInstance(monomer.molfile);
      const rGroupIndices = findAllIndices(molfileHandler.atomTypes, 'R#');
      criterion &&= true;
      return criterion;
    });
    return monomers.map((monomer) => monomer?.symbol!);
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
    for (const lib of libList) if (!lib.error) this._updateInt(lib);
    this._onChanged.next();
  }

  public clear(): void {
    this._monomers = {};
    this._onChanged.next();
  }
}

type MonomerLibWindowType = Window & { $monomerLibHelper: MonomerLibHelper };
declare const window: MonomerLibWindowType;

export class MonomerLibHelper implements IMonomerLibHelper {
  private readonly _monomerLib: MonomerLib = new MonomerLib({});

  /** Protect constructor to prevent multiple instantiation. */
  protected constructor() {}

  /** Singleton monomer library
   * @return {MonomerLibHelper} MonomerLibHelper instance
   */
  getBioLib(): IMonomerLib {
    return this._monomerLib;
  }

  /** Allows syncing with managing settings/loading libraries */
  public loadLibrariesPromise: Promise<void> = Promise.resolve();

  /** Loads libraries based on settings in user storage {@link LIB_STORAGE_NAME}
   * @param {boolean} reload Clean {@link monomerLib} before load libraries [false]
   */
  async loadLibraries(reload: boolean = false): Promise<void> {
    return this.loadLibrariesPromise = this.loadLibrariesPromise.then(async () => {
      // This function is not allowed to throw any exception,
      // because it will prevent further handling monomer library settings
      // through blocking this.loadLibrariesPromise
      try {
        const [libFileNameList, settings]: [string[], LibSettings] = await Promise.all([
          getLibFileNameList(),
          getUserLibSettings(),
        ]);
        const filteredLibFnList = libFileNameList
          .filter((libFileName) => !settings.exclude.includes(libFileName))
          .filter((libFileName) => settings.explicit.length > 0 ? settings.explicit.includes(libFileName) : true);
        const libs: IMonomerLib[] = await Promise.all(filteredLibFnList
          .map((libFileName) => {
            //TODO handle whether files are in place
            return this.readLibrary(LIB_PATH, libFileName).catch((err: any) => {
              const errMsg: string = `Loading monomers from '${libFileName}' error: ` +
                `${err instanceof Error ? err.message : err.toString()}`;
              return new MonomerLib({}, errMsg);
            });
          }));
        this._monomerLib.updateLibs(libs, reload);
      } catch (err: any) {
        const errMsg: string = 'Loading monomer libraries error: ' +
          `${err instanceof Error ? err.message : err.toString()}`;
        grok.shell.warning(errMsg);

        const errStack = err instanceof Error ? err.stack : undefined;
        _package.logger.error(errMsg, undefined, errStack);
      }
    });
  }

  /** Reads library from file shares, handles .json and .sdf
   * @param {string} path Path to library file
   * @param {string} fileName Name of library file
   * @return {Promise<IMonomerLib>} Promise of IMonomerLib
   */
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
    } else if (fileName.endsWith('.json')) {
      const file = await fileSource.readAsText(fileName);
      rawLibData = JSON.parse(file);
    } else if (fileName.endsWith('.csv')) {
      // todo: replace by DataFrame's method after update of js-api
      function toJson(df: DG.DataFrame): any[] {
        return Array.from({length: df.rowCount}, (_, idx) =>
          df.columns.names().reduce((entry: { [key: string]: any }, colName) => {
            entry[colName] = df.get(colName, idx);
            return entry;
          }, {})
        );
      }

      const df = await fileSource.readCsv(fileName);
      const json = toJson(df);
      const polyToolMonomerLib = new PolyToolMonomerLibHandler(json);
      if (polyToolMonomerLib.isValid())
        rawLibData = polyToolMonomerLib.getJsonMonomerLib();
      else
        throw new Error('Invalid format of CSV monomer lib');
    } else {
      throw new Error('Monomer library of unknown file format, supported formats: SDF, JSON, CSV');
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
  public static get instance(): MonomerLibHelper {
    let res = window.$monomerLibHelper;
    if (!res) res = window.$monomerLibHelper = new MonomerLibHelper();
    return res;
  }
}
