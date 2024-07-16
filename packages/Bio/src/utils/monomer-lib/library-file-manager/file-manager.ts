/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {JSONSchemaType} from 'ajv';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {
  HELM_REQUIRED_FIELD as REQ,
} from '@datagrok-libraries/bio/src/utils/const';
import {
  IMonomerLibHelper, IMonomerLibFileManager
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {MonomerLib} from '../monomer-lib';
import {HELM_JSON_SCHEMA_PATH} from './consts';
import {MonomerLibFileEventManager} from './event-manager';
import {MonomerLibFileValidator} from './file-validator';

import {_package} from '../../../package';


/** Singleton for adding, validation and reading of monomer library files.
 * All files **must** be aligned to the HELM standard before adding. */
export class MonomerLibFileManager implements IMonomerLibFileManager {
  public filesPromise: Promise<void> = Promise.resolve();

  private constructor(
    private readonly fileValidator: MonomerLibFileValidator,
    private readonly libHelper: IMonomerLibHelper,
    public readonly eventManager: MonomerLibFileEventManager,
    private readonly logger: ILogger,
  ) {
    this.eventManager.updateValidLibraryFileListRequested$.subscribe(async () => {
      await this.updateValidLibraryList();
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
  async addLibraryFile(fileContent: string, fileName: string): Promise<void> {
    try {
      const alreadyFileExists = await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
      if (alreadyFileExists) {
        grok.shell.error(`File ${fileName} already exists`);
        return;
      }

      await this.validateAgainstHELM(fileContent, fileName);
      await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
      await this.updateValidLibraryList();
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
      await this.updateValidLibraryList();
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

  getValidLibraryPaths(): string[] {
    return this.eventManager.getValidFilesPathList();
  }

  // TODO: remove after adding init from user data storage
  // WARNING: a temporary solution
  async getValidLibraryPathsAsynchronously(): Promise<string[]> {
    return await this.eventManager.getValidLibraryPathsAsynchronously();
  }

  private async updateValidLibraryList(): Promise<void> {
    const logPrefix: string = `${this.toLog()}.updateValidLibraryList()`;
    this.logger.debug(`${logPrefix}, start`);
    this.filesPromise = this.filesPromise.then(async () => {
      this.logger.debug(`${logPrefix}, IN`);
      const invalidFiles = [] as string[];
      // console.log(`files before validation:`, this.libraryEventManager.getValidFilesPathList());
      const filePaths = await this.getFilePathsAtDefaultLocation();

      if (!this.fileListHasChanged(filePaths)) {
        this.logger.debug(`${logPrefix}, end, not changed`);
        return;
      }

      for (const path of filePaths) {
        if (!path.endsWith('.json')) {
          invalidFiles.push(path);
          continue;
        }

        const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${path}`);
        if (!this.isValidHELMLibrary(fileContent, path))
          invalidFiles.push(path);
      }

      const validLibraryPaths = filePaths.filter((path) => !invalidFiles.includes(path));

      if (this.fileListHasChanged(validLibraryPaths)) {
        this.eventManager.changeValidFilesPathList(validLibraryPaths);
        this.libHelper.loadLibraries(true);
      }
      // console.log(`files after validation:`, this.libraryEventManager.getValidFilesPathList());

      if (validLibraryPaths.some((el) => !el.endsWith('.json')))
        this.logger.warning(`Wrong validation: ${validLibraryPaths}`);

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

  private fileListHasChanged(newList: string[]): boolean {
    const currentList = this.eventManager.getValidFilesPathList();
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

  /** Get relative paths for files in LIB_PATH  */
  private async getFilePathsAtDefaultLocation(): Promise<string[]> {
    const logPrefix = `${this.toLog()}.getFilePathsAtDefaultLocation()`;
    this.logger.debug(`${logPrefix}, start`);
    const list = await grok.dapi.files.list(LIB_PATH);
    const paths = list.map((fileInfo) => {
      return fileInfo.fullPath;
    });

    const checkForUi = false;
    const existingPaths = [] as string[];
    if (checkForUi) {
      // WARNING: an extra sanity check,
      // caused by unexpected behavior of grok.dapi.files.list() when it returns non-existent paths
      for (const path of paths) {
        const exists = await grok.dapi.files.exists(path);
        if (exists)
          existingPaths.push(path);
      }
    } else
      existingPaths.push(...paths);

    return existingPaths.map((path) => {
      // Get relative path (to LIB_PATH)
      return path.substring(LIB_PATH.length);
    });
  }
}
