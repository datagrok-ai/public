/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/test/src/test';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {DEFAULT_FILES_LIB_PROVIDER_NAME, findProviderWithLibraryName, IMonomerLib, IMonomerSet} from '@datagrok-libraries/bio/src/types/monomer-library';
import {
  getUserLibSettings, setUserLibSettings,
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {MonomerLib} from './monomer-lib';
import {MonomerSet} from '@datagrok-libraries/bio/src/types/monomer-library';
import {LIB_SETTINGS_FOR_TESTS} from './consts';

import {_package} from '../../package';
import {IMonomerLibHelper, IMonomerLibProvider} from '@datagrok-libraries/bio/src/types/monomer-library';
import {merge, Observable, Subject} from 'rxjs';
import {MonomerLibFromFilesProvider} from './library-file-manager/monomers-lib-provider';

type MonomerLibWindowType = Window & { $monomerLibHelperPromise?: Promise<MonomerLibManager> };
declare const window: MonomerLibWindowType;

/** Singleton wrapper for MonomerLib, provides API for managing libraries on
 * the platform  */
export class MonomerLibManager implements IMonomerLibHelper {
  private readonly _monomerLib = new MonomerLib({}, 'MAIN');
  private readonly _monomerSets = new MonomerSet('MAIN', []);
  private _initialLoadCompleted: boolean = false;
  public get initialLoadCompleted(): boolean { return this._initialLoadCompleted; }

  private _providersDataChanged = new Subject<void>();
  public get providersDataChanged(): Observable<void> {
    return this._providersDataChanged;
  }

  private _fileUploadRequested = new Subject<void>();
  public get fileUploadRequested(): Observable<void> {
    return this._fileUploadRequested;
  }
  public requestFileUpload(): void {
    this._fileUploadRequested.next();
  }

  private _selectionChanged = new Subject<void>();
  public get librarySelectionChanged(): Observable<void> {
    return this._selectionChanged;
  }
  public notifyLibrarySelectionChanged(): void {
    this._selectionChanged.next();
  }

  async getAvaliableLibraryNames(refresh: boolean = false): Promise<string[]> {
    const providers = await this.getProviders();
    const libNames: string[] = [];
    for (const provider of providers) {
      if (refresh)
        await provider.refreshLists();
      const names = await provider.listLibraries();
      libNames.push(...names);
    }
    return libNames;
  }

  async refreshValidLibraryLists(): Promise<void> {
    const providers = await this.getProviders();
    await Promise.all(providers.map(async (provider) => provider.refreshLists()));
  }

  async getAvailableLibrariesPerProvider(): Promise<Map<string, string[]>> {
    const providers = await this.getProviders();
    const libNamesMap: Map<string, string[]> = new Map();
    for (const provider of providers) {
      const names = await provider.listLibraries();
      libNamesMap.set(provider.name, names);
    }
    return libNamesMap;
  }


  public async awaitLoaded(timeout: number = Infinity): Promise<void> {
    timeout = timeout === Infinity ? 1 * 60000 : timeout; // Limit the max lib loaded timeout
    return await Promise.race([
      (async () => {
        const providers = await this.getProviders();
        await Promise.all(providers.map(async (provider) => provider.loadPromise));
        await this.loadLibrariesPromise;
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

  private _monomerLibProviders: IMonomerLibProvider[] | null = null;
  public async getProviders(): Promise<IMonomerLibProvider[]> {
    if (this._monomerLibProviders == null) {
      const providerFuncs = DG.Func.find({meta: {role: 'monomer-lib-provider'}});
      this._monomerLibProviders = await Promise.all(providerFuncs.map(async (func) => {
        return (await func.apply({})) as IMonomerLibProvider;
      }));
      // hard coded files provider from bio package
      this._monomerLibProviders.push(await MonomerLibFromFilesProvider.getInstance(this, _package.logger));

      // eslint-disable-next-line rxjs/no-ignored-subscription --- needed for the lifetime of app
      DG.debounce(merge(...this._monomerLibProviders.map((p) => p.onChanged)), 200).subscribe(() => {
        this._providersDataChanged.next();
      });
    }
    return this._monomerLibProviders;
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


  /** Allows syncing with managing settings/loading libraries */
  private loadLibrariesPromise: Promise<void> = Promise.resolve();

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
        // const [[libFileNameList,], settings]: [[string[],], UserLibSettings] =
        //   await Promise.all([
        //     await this.getFileManager().then((fileManager) => {
        //       return [fileManager.getValidLibraryPaths(),];
        //     }),
        //     getUserLibSettings(),
        //   ]);
        const settings: UserLibSettings = await getUserLibSettings();
        const providers = await this.getProviders();
        const libNames = await Promise.all(providers.map(async (provider) => provider.listLibraries()));
        libNames.forEach((names, i) => {
          libNames[i] = names.filter((name) => {
            const isFileIncluded = !settings.exclude.includes(name);
            const isExplicit = (settings.explicit?.length ?? 0) === 0 || settings.explicit.includes(name);
            return isFileIncluded && isExplicit;
          });
        });

        let completedFileCount: number = 0;
        const libs: IMonomerLib[] = [];
        for (let i = 0; i < providers.length; i++) {
          const provider = providers[i];
          const names = libNames[i];
          try {
            const libDatas = await provider.loadLibraries(names);
            for (let j = 0; j < libDatas.length; j++)
              libs.push(new MonomerLib(libDatas[j], names[j]));
          } catch (err: any) {
            grok.shell.warning(`Loading monomer libraries from provider '${provider.name}' failed: ` +
            ``);
            this.logger.error(`Loading monomer libraries from provider '${provider.name}' failed: ` +
            `${err instanceof Error ? err.message : err.toString()}`);
          }
          pi.update(Math.round(100 * (++completedFileCount) / providers.length),
            `Loading monomers ${completedFileCount}/${providers.length}`);
        }


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

  async loadMonomerSets(_reload: boolean = false): Promise<void> {
    const logPrefix = `${this.toLog()}.loadMonomerSets()`;

    this.logger.debug(`${logPrefix}, start`);
    this.loadSetsPromise = this.loadSetsPromise.then(async () => {
      this.logger.debug(`${logPrefix}, IN`);
      const pi = DG.TaskBarProgressIndicator.create(`Loading monomer sets ...`);
      try {
        // const [[setFileNameList,]]: [[string[],],] = await Promise.all([
        //   await this.getFileManager().then((fileManager) => {
        //     return [fileManager.getValidSetPaths(),];
        //   })
        // ]);

        // // TODO: Filter for settings
        // const filteredSetFnList = setFileNameList.filter((setFileName) => true);

        // let completedFileCount: number = 0;
        // const allCount: number = filteredSetFnList.length;
        // const [sets,]: [IMonomerSet[],] = await Promise.all([
        //   Promise.all(filteredSetFnList.map((setFileName) => {
        //     // TODO: handle whether files are in place
        //     return this.readSet(SETS_PATH, setFileName)
        //       .catch((err: any) => {
        //         const errMsg: string = `Loading monomer sets from '${setFileName}' error: ` +
        //           `${err instanceof Error ? err.message : err.toString()}`;
        //         return new MonomerSet('Broken monomer set', [], setFileName, errMsg);
        //       })
        //       .finally(() => {
        //         pi.update(Math.round(100 * (++completedFileCount) / allCount),
        //           `Loading monomers ${completedFileCount}/${allCount}`);
        //       });
        //   })),]);
        const providers = await this.getProviders();
        const setNames = await Promise.all(providers.map(async (provider) => provider.listSets()));
        const sets: IMonomerSet[] = [];
        let completedFileCount: number = 0;
        for (let i = 0; i < providers.length; i++) {
          const provider = providers[i];
          const names = setNames[i];
          try {
            const setDatas = await provider.loadSets(names);
            for (let j = 0; j < setDatas.length; j++)
              sets.push(setDatas[j]);
          } catch (err: any) {
            grok.shell.warning(`Loading monomer sets from provider '${provider.name}' failed: ` +
            ``);
            this.logger.error(`Loading monomer sets from provider '${provider.name}' failed: ` +
            `${err instanceof Error ? err.message : err.toString()}`);
          }
          pi.update(Math.round(100 * (++completedFileCount) / providers.length),
            `Loading monomer sets ${completedFileCount}/${providers.length}`);
        }


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

  // /** Reads library from file shares, handles .json and .sdf
  //  * @param {string} path Path to library file
  //  * @param {string} fileName Name of library file
  //  * @return {Promise<IMonomerLib>} Promise of IMonomerLib
  //  */
  // async readLibrary(path: string, fileName: string): Promise<IMonomerLib> {
  //   const fileManager = await this.getFileManager();
  //   const lib: IMonomerLib = await fileManager.loadLibraryFromFile(path, fileName);
  //   return lib;
  // }

  // async readSet(path: string, fileName: string): Promise<IMonomerSet> {
  //   const fileManager = await this.getFileManager();
  //   const set: IMonomerSet = await fileManager.loadSetFromFile(this._monomerLib, path, fileName);
  //   return set;
  // }
  // -- Settings --

  /** Changes userLibSettings set only HELMCoreLibrary.json, polytool-lib.json */
  async loadMonomerLibForTests(): Promise<void> {
    await setUserLibSettings(LIB_SETTINGS_FOR_TESTS);
    await this.awaitLoaded(10000);
    await this.loadMonomerLib(true); // load default libraries
  }

  async readSingleLibrary(providerName: string, libraryName: string): Promise<IMonomerLib | null> {
    const providers = await this.getProviders();
    const provider = providers.find((p) => p.name === providerName);
    if (!provider) {
      this.logger.error(`Provider '${providerName}' not found.`);
      return null;
    }
    try {
      const libDatas = await provider.loadLibraries([libraryName]);
      if (libDatas.length === 0) {
        this.logger.error(`Library '${libraryName}' not found in provider '${providerName}'.`);
        return null;
      }
      return new MonomerLib(libDatas[0], libraryName);
    } catch (err: any) {
      this.logger.error(`Loading monomer library '${libraryName}' from provider '${provider.name}' failed: ` +
      `${err instanceof Error ? err.message : err.toString()}`);
      return null;
    }
  }


  async readLibraryFromFilePath(filePath: string): Promise<IMonomerLib> {
    const provider: MonomerLibFromFilesProvider = (await this.getProviders()).find((p) => p.name === DEFAULT_FILES_LIB_PROVIDER_NAME) as MonomerLibFromFilesProvider;
    if (!provider) {
      // Just in case, but this should never happen
      throw new Error(`Provider '${DEFAULT_FILES_LIB_PROVIDER_NAME}' not found.`);
    }

    const res = await provider.getSingleLibraryWithFullPath(filePath);
    if (!res)
      throw new Error(`Library at path '${filePath}' not found.`);
    return new MonomerLib(res, filePath);
  }

  async readSingleLibraryByName(libraryName: string): Promise<IMonomerLib | null> {
    const providers = await this.getProviders();
    const provider = await findProviderWithLibraryName(providers, libraryName);
    if (!provider) {
      this.logger.error(`No provider found for library '${libraryName}'.`);
      return null;
    }
    return this.readSingleLibrary(provider.name, libraryName);
  }

  async readSingleSet(providerName: string, setName: string): Promise<IMonomerSet | null> {
    const providers = await this.getProviders();
    const provider = providers.find((p) => p.name === providerName);
    if (!provider) {
      this.logger.error(`Provider '${providerName}' not found.`);
      return null;
    }
    try {
      const setDatas = await provider.loadSets([setName]);
      if (setDatas.length === 0) {
        this.logger.error(`Set '${setName}' not found in provider '${providerName}'.`);
        return null;
      }
      return setDatas[0];
    } catch (err: any) {
      this.logger.error(`Loading monomer set '${setName}' from provider '${provider.name}' failed: ` +
      `${err instanceof Error ? err.message : err.toString()}`);
      return null;
    }
  }

  // -- Instance singleton --
  public static async getInstance(): Promise<MonomerLibManager> {
    let res = window.$monomerLibHelperPromise;
    if (res == undefined) {
      res = window.$monomerLibHelperPromise = (async () => {
        const instance = new MonomerLibManager(_package.logger);
        return instance;
      })();
    }
    return res;
  }
}
