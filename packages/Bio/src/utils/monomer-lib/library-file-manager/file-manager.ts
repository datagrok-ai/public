/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {JSONSchemaType} from 'ajv';

import {MonomerType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLib, IMonomerLinkData, IMonomerSet, Monomer} from '@datagrok-libraries/bio/src/types';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {
  HELM_REQUIRED_FIELD as REQ,
} from '@datagrok-libraries/bio/src/utils/const';
import {
  IMonomerLibHelper, IMonomerLibFileManager
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {MonomerLib} from '../monomer-lib';
import {MonomerLibFileEventManager} from './event-manager';
import {MonomerLibFileValidator} from './file-validator';
import {MonomerSet, MonomerSetPlaceholder} from '../monomer-set';
import {HELM_JSON_SCHEMA_PATH, LIB_PATH, SETS_PATH} from '../consts';

import {_package} from '../../../package';


/** Singleton for adding, validation and reading of monomer library files.
 * All files **must** be aligned to the HELM standard before adding. */
export class MonomerLibFileManager implements IMonomerLibFileManager {
  public filesPromise: Promise<void> = Promise.resolve();
  public initializedPromise: Promise<void> = Promise.resolve();

  private constructor(
    private readonly fileValidator: MonomerLibFileValidator,
    private readonly libHelper: IMonomerLibHelper,
    public readonly eventManager: MonomerLibFileEventManager,
    private readonly logger: ILogger,
  ) {
    // these both are behavioral subjects, i.e. they emit their value when subscribed.
    // initial creation/request from bio package on awaiting files promise makes no sense,
    // until the subscription fires first time
    let resolveFilesPromise: () => void;
    let initialized = false;
    this.initializedPromise =
      Promise.race([DG.delay(1000), new Promise<void>((resolve) => resolveFilesPromise = resolve)]);
    const _libSub = this.eventManager.updateValidLibraryFileListRequested$.subscribe(() => {
      this.updateValidLibList().then(() => {});
      if (!initialized) {
        initialized = true;
        resolveFilesPromise();
      }
    });
    const _setSub = this.eventManager.updateValidSetFileListRequested$.subscribe(() => {
      this.updateValidSetList().then(() => {});
    });
  }

  private static objCounter: number = -1;
  private readonly objId: number = ++MonomerLibFileManager.objCounter;

  protected toLog(): string {
    return `MonomerLibFileManager<${this.objId}>`;
  }

  /** For internal use only, get {@link IMonomerLibHelper.getFileManager} */
  public static async create(
    libHelper: IMonomerLibHelper, eventManager: MonomerLibFileEventManager, logger: ILogger,
  ): Promise<MonomerLibFileManager> {
    const helmSchemaRaw = await grok.dapi.files.readAsText(HELM_JSON_SCHEMA_PATH);
    const helmSchema = JSON.parse(helmSchemaRaw) as JSONSchemaType<any>;

    const fileValidator = new MonomerLibFileValidator(helmSchema);
    return new MonomerLibFileManager(fileValidator, libHelper, eventManager, logger);
  }

  /** Add standard .json monomer library  */
  async addLibraryFile(fileContent: string, fileName: string, reload = true): Promise<void> {
    try {
      const alreadyFileExists = await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
      if (alreadyFileExists) {
        grok.shell.error(`File ${fileName} already exists`);
        return;
      }

      await this.validateAgainstHELM(fileContent, fileName);
      await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
      await this.updateValidLibList();
      const fileExists = await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
      if (!fileExists)
        grok.shell.error(`Failed to add ${fileName} library`);
      else
        grok.shell.info(`Added ${fileName} HELM library`);
    } catch (e) {
      console.error(e);
      grok.shell.error(`Failed to add ${fileName} library`);
    }
  }

  async deleteLibraryFile(fileName: string): Promise<void> {
    try {
      await grok.dapi.files.delete(LIB_PATH + `${fileName}`);
      await this.updateValidLibList();
      grok.shell.info(`Deleted ${fileName} library`);
    } catch (e) {
      console.error(e);
      const fileExists = await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
      if (fileExists)
        grok.shell.error(`Failed to delete ${fileName} library`);
      else
        grok.shell.warning(`File ${fileName} already deleted, refresh the list`);
    }
  }

  async loadLibraryFromFile(path: string, fileName: string): Promise<IMonomerLib> {
    let rawLibData: any[] = [];
    const fileSource = new DG.FileSource(path);
    const file = await fileSource.readAsText(fileName);
    rawLibData = JSON.parse(file);
    const monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
    const polymerTypes: string[] = [];
    rawLibData.forEach((monomer) => {
      const polymerType = monomer[REQ.POLYMER_TYPE];
      const monomerSymbol = monomer[REQ.SYMBOL];
      if (!polymerTypes.includes(polymerType)) {
        monomers[polymerType] = {};
        polymerTypes.push(polymerType);
      }
      monomers[polymerType][monomerSymbol] = monomer as Monomer;
    });

    return new MonomerLib(monomers, fileName);
  }

  async loadSetFromFile(monomerLib: IMonomerLib, path: string, fileName: string): Promise<IMonomerSet> {
    let raw: any = {};
    const fileSource = new DG.FileSource(path);
    const content = await fileSource.readAsText(fileName);
    raw = JSON.parse(content);

    const description = raw['description'];
    const placeholders = Object.entries(raw['placeholders']).map(([k, v]: [string, any]) => {
      const placeholderSymbol = k;
      const polymerType = v['polymerType'] as PolymerType;
      const monomerType = v['monomerType'] as MonomerType;
      const monomerLinks = v['set'] as IMonomerLinkData[];

      return new MonomerSetPlaceholder(monomerLib, placeholderSymbol, polymerType, monomerType, monomerLinks);
    });

    return new MonomerSet(description, placeholders);
  }

  getValidLibraryPaths(): string[] {
    return this.eventManager.getValidLibPathList();
  }

  getValidSetPaths(): string[] {
    return this.eventManager.getValidSetPathList();
  }

  // TODO: remove after adding init from user data storage
  // WARNING: a temporary solution
  async getValidLibraryPathsAsynchronously(): Promise<string[]> {
    return await this.eventManager.getValidLibraryPathsAsynchronously();
  }

  private async updateValidLibList(): Promise<void> {
    const logPrefix: string = `${this.toLog()}.updateValidLibList()`;
    this.logger.debug(`${logPrefix}, start`);
    this.filesPromise = this.filesPromise.then(async () => {
      this.logger.debug(`${logPrefix}, IN`);
      const invalidFiles = [] as string[];
      // console.log(`files before validation:`, this.libraryEventManager.getValidFilesPathList());
      const libPathList = await this.getLibFileListAtLocation();

      if (!this.libListHasChanged(libPathList)) {
        this.logger.debug(`${logPrefix}, end, not changed`);
        return;
      }

      for (const path of libPathList) {
        if (!path.endsWith('.json')) {
          invalidFiles.push(path);
          continue;
        }

        const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${path}`);
        if (!this.isValidHELMLibrary(fileContent, path))
          invalidFiles.push(path);
      }

      const validLibPathList = libPathList.filter((path) => !invalidFiles.includes(path));

      if (this.libListHasChanged(validLibPathList)) {
        this.eventManager.changeValidLibPathList(validLibPathList);
        await this.libHelper.loadMonomerLib(true);
      }
      // console.log(`files after validation:`, this.libraryEventManager.getValidFilesPathList());

      if (validLibPathList.some((el) => !el.endsWith('.json')))
        this.logger.warning(`Wrong validation: ${validLibPathList}`);

      if (invalidFiles.length > 0) {
        const message = `Invalid monomer library files in ${LIB_PATH}` +
          `, consider fixing or removing them: ${invalidFiles.join(', ')}`;

        this.logger.warning(message);
        // grok.shell.warning(message);
      }
      this.logger.debug(`${logPrefix}, OUT`);
    });
    this.logger.debug(`${logPrefix}, end`);
    return this.filesPromise;
  }

  private async updateValidSetList(): Promise<void> {
    const logPrefix: string = `${this.toLog()}.updateValidSetList()`;
    _package.logger.debug(`${logPrefix}, start`);
    this.filesPromise = this.filesPromise.then(async () => {
      _package.logger.debug(`${logPrefix}, IN`);
      const invalidFiles = [] as string[];
      const setPathList: string[] = await this.getSetFileListAtLocation();

      if (!this.setListHasChanged(setPathList)) {
        _package.logger.debug(`${logPrefix}, end, not changed`);
        return;
      }

      for (const path of setPathList) {
        if (!path.endsWith('.json')) {
          invalidFiles.push(path);
          continue;
        }

        const fileContent = await grok.dapi.files.readAsText(SETS_PATH + `${path}`);
        // TODO: Validate monomer set
        // if (!this.isValidMonomerSet(fileContent, path))
        //   invalidFiles.push(path);
      }

      const validSetPathList = setPathList.filter((path) => !invalidFiles.includes(path));
      if (this.setListHasChanged(validSetPathList)) {
        this.eventManager.changeValidSetPathList(validSetPathList);
        this.libHelper.loadMonomerSets(true);
      }

      _package.logger.debug(`${logPrefix}, OUT`);
    });
    _package.logger.debug(`${logPrefix}, end`);
    return this.filesPromise;
  }

  private libListHasChanged(newList: string[]): boolean {
    const currentList = this.eventManager.getValidLibPathList();
    return newList.length !== currentList.length || newList.some((el, i) => el !== currentList[i]);
  }

  private setListHasChanged(newList: string[]): boolean {
    const currentList = this.eventManager.getValidSetPathList();
    return newList.length !== currentList.length || newList.some((el, i) => el !== currentList[i]);
  }

  private async validateAgainstHELM(fileContent: string, fileName: string): Promise<void> {
    const isValid = this.isValidHELMLibrary(fileContent, fileName);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private isValidHELMLibrary(fileContent: string, fileName: string): boolean {
    return this.fileValidator.validateFile(fileContent, fileName);
  }

  /** Get relative paths for files in LIB_PATH, SET_PATH  */
  private async getLibFileListAtLocation(): Promise<string[]> {
    const logPrefix = `${this.toLog()}.getLibFileListAtLocation()`;
    this.logger.debug(`${logPrefix}, start`);

    const libPaths = await grok.dapi.files.list(LIB_PATH)
      .then((l) => l.filter((f) => f.isFile).map((fi) => fi.fullPath));

    const checkForUi = false;
    const existingLibPaths: string[] = [];
    if (checkForUi) {
      // WARNING: an extra sanity check,
      // caused by unexpected behavior of grok.dapi.files.list() when it returns non-existent paths
      for (const path of libPaths) {
        const exists = await grok.dapi.files.exists(path);
        if (exists)
          existingLibPaths.push(path);
      }
    } else
      existingLibPaths.push(...libPaths);

    return existingLibPaths.map((p) => /* relative to LIB_PATH */ p.substring(LIB_PATH.length));
  }

  private async getSetFileListAtLocation(): Promise<string[]> {
    const logPrefix = `${this.toLog()}.getSetFileListAtLocation()`;
    this.logger.debug(`${logPrefix}, start`);

    const setPaths = await grok.dapi.files.list(SETS_PATH)
      .then((l) => l.map((fi) => fi.fullPath));

    const checkForUi = false;
    const existingSetPaths: string[] = [];
    if (checkForUi) {
      for (const path of setPaths) {
        const exists = await grok.dapi.files.exists(path);
        if (exists)
          existingSetPaths.push(path);
      }
    } else
      existingSetPaths.push(...setPaths);

    return existingSetPaths.map((p) => /* relative to SET_PATH */ p.substring(SETS_PATH.length));
  }
}
