/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/utils/src/test';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {IMonomerLib, IMonomerSet} from '@datagrok-libraries/bio/src/types';
import {
  getUserLibSettings, setUserLibSettings,
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  IMonomerLibFileEventManager, IMonomerLibHelper,
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {MonomerLib} from './monomer-lib';
import {MonomerSet} from './monomer-set';
import {MonomerLibFileManager} from './library-file-manager/file-manager';
import {MonomerLibFileEventManager} from './library-file-manager/event-manager';
import {LIB_PATH, LIB_SETTINGS_FOR_TESTS, SETS_PATH} from './consts';

import {_package} from '../../package';

type MonomerLibWindowType = Window & { $monomerLibHelperPromise?: Promise<MonomerLibManager> };
declare const window: MonomerLibWindowType;

/** Singleton wrapper for MonomerLib, provides API for managing libraries on
 * the platform  */
export class MonomerLibManager implements IMonomerLibHelper {
  private readonly _monomerLib = new MonomerLib({}, 'MAIN');
  private readonly _monomerSets = new MonomerSet('MAIN', []);
  private _initialLoadCompleted: boolean = false;
  public get initialLoadCompleted(): boolean { return this._initialLoadCompleted; }
  private _eventManager: MonomerLibFileEventManager;

  public get eventManager(): IMonomerLibFileEventManager { return this._eventManager; }

  public async awaitLoaded(timeout: number = Infinity): Promise<void> {
    timeout = timeout === Infinity ? 1 * 60000 : timeout; // Limit the max lib loaded timeout
    return await Promise.race([
      (async () => {
        const fileManager = await this.getFileManager();
        await fileManager.filesPromise;
        return true;
      })(),
      (async () => {
        await delay(timeout);
        return false;
      })(),
    ]).then((res) => {
      if (!res)
        throw new Error(`Loading monomer libraries timeout ${timeout} ms.`);
    });
  }

  /** Protect constructor to prevent multiple instantiation. */
  protected constructor(
    private readonly logger: ILogger,
  ) {}

  private static objCounter: number = -1;
  private readonly objId: number = (() => {
    if (++MonomerLibManager.objCounter > 0)
      throw new Error('MonomerLibManager MUST be a singleton.');
    return MonomerLibManager.objCounter;
  })();

  protected toLog(): string {
    return `MonomerLibManager<${this.objId}>`;
  }

  /** Singleton monomer library
   * @return {MonomerLibManager} MonomerLibHelper instance
   */
  getMonomerLib(): IMonomerLib {
    return this._monomerLib;
  }

  /**  @deprecated Use {@link v} */
  getBioLib(): IMonomerLib {
    return this.getMonomerLib();
  }

  getMonomerSets(): IMonomerSet {
    return this._monomerSets;
  }

  /** Object containing symbols for each type of polymer where duplicate monomers
   * are found in different libs (based on symbol as key) */
  get duplicateMonomers() {
    return this._monomerLib.duplicateMonomers;
  }

  /** Returns true if all duplicate monomers are assigned preferences (which library they are coming from) */
  get duplicatesHandled() {
    return this._monomerLib.duplicatesHandled;
  }

  assignDuplicatePreferances(settings: UserLibSettings) {
    this._monomerLib.assignDuplicatePreferences(settings);
  }

  /** Instance promise of {@link getFileManager} */
  private _fileManagerPromise?: Promise<MonomerLibFileManager>;

  async getFileManager(): Promise<MonomerLibFileManager> {
    if (this._fileManagerPromise === undefined) {
      this._fileManagerPromise = (async () => {
        const fileManager: MonomerLibFileManager =
          await MonomerLibFileManager.create(this, this._eventManager, this.logger);
        await fileManager.initializedPromise;
        return fileManager;
      })();
    }
    return this._fileManagerPromise;
  }

  /** Allows syncing with managing settings/loading libraries */
  public loadLibrariesPromise: Promise<void> = Promise.resolve();

  /** Loads libraries based on settings in user storage {@link LIB_STORAGE_NAME}
   * @param {boolean} reload Clean {@link monomerLib} before load libraries [false]
   */
  async loadMonomerLib(reload: boolean = false): Promise<void> {
    const logPrefix = `${this.toLog()}.loadMonomerLib()`;
    this.logger.debug(`${logPrefix}, start`);
    this.loadLibrariesPromise = this.loadLibrariesPromise.then(async () => {
      this.logger.debug(`${logPrefix}, IN`);

      const pi = DG.TaskBarProgressIndicator.create('Loading monomers ...');
      try {
        const [[libFileNameList,], settings]: [[string[],], UserLibSettings] =
          await Promise.all([
            await this.getFileManager().then((fileManager) => {
              return [fileManager.getValidLibraryPaths(),];
            }),
            getUserLibSettings(),
          ]);

        const filteredLibFnList = libFileNameList.filter((libFileName) => {
          const isFileIncluded = !settings.exclude.includes(libFileName);
          const isExplicit = settings.explicit.length === 0 || settings.explicit.includes(libFileName);

          return isFileIncluded && isExplicit;
        });

        let completedFileCount: number = 0;
        const allCount: number = filteredLibFnList.length;
        const [libs,]: [IMonomerLib[],] = await Promise.all([
          Promise.all(filteredLibFnList.map((libFileName) => {
            return this.readLibrary(LIB_PATH, libFileName)
              .catch((err: any) => {
                const errMsg: string = `Loading monomers from '${libFileName}' error: ` +
                  `${err instanceof Error ? err.message : err.toString()}`;
                return new MonomerLib({}, libFileName, errMsg);
              })
              .finally(() => {
                pi.update(Math.round(100 * (++completedFileCount) / allCount),
                  `Loading monomers ${completedFileCount}/${allCount}`);
              });
          })),]);
        this._monomerLib.updateLibs(libs, reload);
        this._initialLoadCompleted = true;
      } catch (err: any) {
        // WARNING: This function is not allowed to throw any exception,
        // because it will prevent further handling monomer library settings
        // through blocking this.loadLibrariesPromise

        const errMsg: string = 'Loading monomer libraries error: ' +
          `${err instanceof Error ? err.message : err.toString()}`;
        grok.shell.warning(errMsg);

        const errStack = err instanceof Error ? err.stack : undefined;
        this.logger.error(errMsg, undefined, errStack);
      } finally {
        pi.close();
        this.logger.debug(`${logPrefix}, OUT`);
      }
    });
    this.logger.debug(`${logPrefix}, end`);
    return this.loadLibrariesPromise;
  }

  /** @deprecated Use {@link loadMonomerLib} */
  async loadLibraries(reload?: boolean): Promise<void> {
    return this.loadMonomerLib(reload);
  }

  private loadSetsPromise: Promise<void> = Promise.resolve();

  async loadMonomerSets(reload: boolean = false): Promise<void> {
    const logPrefix = `${this.toLog()}.loadMonomerSets()`;

    this.logger.debug(`${logPrefix}, start`);
    this.loadSetsPromise = this.loadSetsPromise.then(async () => {
      this.logger.debug(`${logPrefix}, IN`);
      const pi = DG.TaskBarProgressIndicator.create(`Loading monomer sets ...`);
      try {
        const [[setFileNameList,]]: [[string[],],] = await Promise.all([
          await this.getFileManager().then((fileManager) => {
            return [fileManager.getValidSetPaths(),];
          })
        ]);

        // TODO: Filter for settings
        const filteredSetFnList = setFileNameList.filter((setFileName) => true);

        let completedFileCount: number = 0;
        const allCount: number = filteredSetFnList.length;
        const [sets,]: [IMonomerSet[],] = await Promise.all([
          Promise.all(filteredSetFnList.map((setFileName) => {
            // TODO: handle whether files are in place
            return this.readSet(SETS_PATH, setFileName)
              .catch((err: any) => {
                const errMsg: string = `Loading monomer sets from '${setFileName}' error: ` +
                  `${err instanceof Error ? err.message : err.toString()}`;
                return new MonomerSet('Broken monomer set', [], setFileName, errMsg);
              })
              .finally(() => {
                pi.update(Math.round(100 * (++completedFileCount) / allCount),
                  `Loading monomers ${completedFileCount}/${allCount}`);
              });
          })),]);
        this._monomerSets.updateSets(sets);
      } catch (err: any) {
        const errMsg: string = 'Loading monomer sets error: ' +
          `${err instanceof Error ? err.message : err.toString()}`;
        grok.shell.warning(errMsg);

        const errStack = err instanceof Error ? err.stack : undefined;
        this.logger.error(errMsg, undefined, errStack);
      } finally {
        pi.close();
        this.logger.debug(`${logPrefix}, OUT`);
      }
    });
    this.logger.debug(`${logPrefix}, end`);
    return this.loadSetsPromise;
  }

  /** Reads library from file shares, handles .json and .sdf
   * @param {string} path Path to library file
   * @param {string} fileName Name of library file
   * @return {Promise<IMonomerLib>} Promise of IMonomerLib
   */
  async readLibrary(path: string, fileName: string): Promise<IMonomerLib> {
    const fileManager = await this.getFileManager();
    const lib: IMonomerLib = await fileManager.loadLibraryFromFile(path, fileName);
    return lib;
  }

  async readSet(path: string, fileName: string): Promise<IMonomerSet> {
    const fileManager = await this.getFileManager();
    const set: IMonomerSet = await fileManager.loadSetFromFile(this._monomerLib, path, fileName);
    return set;
  }

  /** Reset user settings to the specified library. WARNING: clears user * settings */
  public async selectSpecifiedLibraries(libFileNameList: string[]): Promise<void> {
    const invalidNames = await this.getInvalidFileNames(libFileNameList);
    if (invalidNames.length > 0)
      throw new Error(`Cannot select libraries ${invalidNames}: no such library in the list`);
    const settings = await getUserLibSettings();
    settings.exclude = ((await this.getFileManager()).getValidLibraryPaths())
      .filter((fileName) => !libFileNameList.includes(fileName));
    await setUserLibSettings(settings);
  }

  private async getInvalidFileNames(libFileNameList: string[]): Promise<string[]> {
    const availableFileNames = (await this.getFileManager()).getValidLibraryPaths();
    const invalidNames = libFileNameList.filter((fileName) => !availableFileNames.includes(fileName));
    return invalidNames;
  }

  // -- Settings --

  /** Changes userLibSettings set only HELMCoreLibrary.json, polytool-lib.json */
  async loadMonomerLibForTests(): Promise<void> {
    await setUserLibSettings(LIB_SETTINGS_FOR_TESTS);
    await this.awaitLoaded(10000);
    await this.loadMonomerLib(true); // load default libraries
  }

  // -- Instance singleton --
  public static async getInstance(): Promise<MonomerLibManager> {
    let res = window.$monomerLibHelperPromise;
    if (res == undefined) {
      res = window.$monomerLibHelperPromise = (async () => {
        const instance = new MonomerLibManager(_package.logger);
        instance._eventManager = MonomerLibFileEventManager.getInstance();
        return instance;
      })();
    }
    return res;
  }
}
