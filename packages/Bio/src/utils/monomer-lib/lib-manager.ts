/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/utils/src/test';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {
  getUserLibSettings, setUserLibSettings, LIB_PATH
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  IMonomerLibFileEventManager, IMonomerLibHelper,
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {MonomerLib} from './monomer-lib';
import {MonomerLibFileManager} from './library-file-manager/file-manager';
import {MonomerLibFileEventManager} from './library-file-manager/event-manager';

import {_package} from '../../package';

type MonomerLibWindowType = Window & { $monomerLibHelperPromise?: Promise<MonomerLibManager> };
declare const window: MonomerLibWindowType;

/** Singleton wrapper for MonomerLib, provides API for managing libraries on
 * the platform  */
export class MonomerLibManager implements IMonomerLibHelper {
  private readonly _monomerLib = new MonomerLib({}, 'MAIN');

  private _eventManager: MonomerLibFileEventManager;

  public get eventManager(): IMonomerLibFileEventManager { return this._eventManager; }

  public async awaitLoaded(timeout: number = 3000): Promise<void> {
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
        throw new Error(`Loading monomer libraries is timeout ${timeout} ms.`);
    });
  }

  /** Protect constructor to prevent multiple instantiation. */
  protected constructor(
    private readonly logger: ILogger,
  ) {}

  /** Singleton monomer library
   * @return {MonomerLibManager} MonomerLibHelper instance
   */
  getBioLib(): IMonomerLib {
    return this._monomerLib;
  }

  /** Instance promise of {@link getFileManager} */
  private _fileManagerPromise?: Promise<MonomerLibFileManager>;

  async getFileManager(): Promise<MonomerLibFileManager> {
    if (this._fileManagerPromise === undefined) {
      this._fileManagerPromise = (async () => {
        const fileManager: MonomerLibFileManager =
          await MonomerLibFileManager.create(this, this._eventManager, this.logger);
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
  async loadLibraries(reload: boolean = false): Promise<void> {
    return this.loadLibrariesPromise = this.loadLibrariesPromise.then(async () => {
      // WARNING: This function is not allowed to throw any exception,
      // because it will prevent further handling monomer library settings
      // through blocking this.loadLibrariesPromise
      const pi = DG.TaskBarProgressIndicator.create('Loading monomers ...');
      try {
        const [libFileNameList, settings]: [string[], UserLibSettings] = await Promise.all([
          (await this.getFileManager()).getValidLibraryPaths(),
          getUserLibSettings(),
        ]);

        const filteredLibFnList = libFileNameList.filter((libFileName) => {
          const isFileIncluded = !settings.exclude.includes(libFileName);
          const isExplicit = settings.explicit.length === 0 || settings.explicit.includes(libFileName);

          return isFileIncluded && isExplicit;
        });

        let completedLibCount: number = 0;
        const libs: IMonomerLib[] = await Promise.all(filteredLibFnList
          .map((libFileName) => {
            //TODO handle whether files are in place
            return this.readLibrary(LIB_PATH, libFileName)
              .catch((err: any) => {
                const errMsg: string = `Loading monomers from '${libFileName}' error: ` +
                  `${err instanceof Error ? err.message : err.toString()}`;
                return new MonomerLib({}, libFileName, errMsg);
              }).finally(() => {
                pi.update(Math.round(100 * (++completedLibCount) / filteredLibFnList.length),
                  `Loading monomer libs ${completedLibCount}/${filteredLibFnList.length}`);
              });
          }));
        this._monomerLib.updateLibs(libs, reload);
      } catch (err: any) {
        const errMsg: string = 'Loading monomer libraries error: ' +
          `${err instanceof Error ? err.message : err.toString()}`;
        grok.shell.warning(errMsg);

        const errStack = err instanceof Error ? err.stack : undefined;
        _package.logger.error(errMsg, undefined, errStack);
      } finally {
        pi.close();
      }
    });
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

  // -- Instance singleton --
  public static async getInstance(): Promise<MonomerLibManager> {
    let res = window.$monomerLibHelperPromise;
    if (res === undefined) {
      res = window.$monomerLibHelperPromise = (async () => {
        const instance = new MonomerLibManager(_package.logger);
        instance._eventManager = MonomerLibFileEventManager.getInstance();
        return instance;
      })();
    }
    return res;
  }
}
