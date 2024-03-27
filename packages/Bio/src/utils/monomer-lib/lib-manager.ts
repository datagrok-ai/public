/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib} from '@datagrok-libraries/bio/src/types/index';
import {
  getUserLibSettings, setUserLibSettings, LIB_PATH
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  IMonomerLibHelper,
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {MonomerLib} from './monomer-lib';
import {MonomerLibFileManager} from './library-file-manager/file-manager';
import {MonomerLibFileEventManager} from './library-file-manager/event-manager';
import {_package} from '../../package';

type MonomerLibWindowType = Window & { $monomerLibHelper: MonomerLibManager };
declare const window: MonomerLibWindowType;

export async function getLibFileNameList(): Promise<string[]> {
  const fileEventManager = MonomerLibFileEventManager.getInstance();
  const fileManager = await MonomerLibFileManager.getInstance(fileEventManager);
  return fileManager.getValidLibraryPaths();
}

/** Singleton wrapper for MonomerLib, provides API for managing libraries on
 * the platform  */
export class MonomerLibManager implements IMonomerLibHelper {
  private readonly _monomerLib = new MonomerLib({});

  /** Protect constructor to prevent multiple instantiation. */
  protected constructor() {}

  /** Singleton monomer library
   * @return {MonomerLibManager} MonomerLibHelper instance
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
      // WARNING: This function is not allowed to throw any exception,
      // because it will prevent further handling monomer library settings
      // through blocking this.loadLibrariesPromise
      try {
        const [libFileNameList, settings]: [string[], UserLibSettings] = await Promise.all([
          getLibFileNameList(),
          getUserLibSettings(),
        ]);

        const filteredLibFnList = libFileNameList.filter((libFileName) => {
          const isFileIncluded = !settings.exclude.includes(libFileName);
          const isExplicit = settings.explicit.length === 0 || settings.explicit.includes(libFileName);

          return isFileIncluded && isExplicit;
        });

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
    const eventManager = MonomerLibFileEventManager.getInstance();
    const libFileManager = await MonomerLibFileManager.getInstance(eventManager);
    const lib: IMonomerLib = await libFileManager.loadLibraryFromFile(path, fileName);
    return lib;
  }

  /** Reset user settings to the specified library. WARNING: clears user * settings */
  public async selectSpecifiedLibraries(libFileNameList: string[]): Promise<void> {
    const invalidNames = await this.getInvalidFileNames(libFileNameList);
    if (invalidNames.length > 0)
      throw new Error(`Cannot select libraries ${invalidNames}: no such library in the list`);
    const settings = await getUserLibSettings();
    settings.exclude = (await getLibFileNameList()).filter((fileName) => !libFileNameList.includes(fileName));
    await setUserLibSettings(settings);
  }

  private async getInvalidFileNames(libFileNameList: string[]): Promise<string[]> {
    const availableFileNames = await getLibFileNameList();
    const invalidNames = libFileNameList.filter((fileName) => !availableFileNames.includes(fileName));
    return invalidNames;
  }

  // -- Instance singleton --
  public static get instance(): MonomerLibManager {
    if (!window.$monomerLibHelper) window.$monomerLibHelper = new MonomerLibManager();
    return window.$monomerLibHelper;
  }
}
