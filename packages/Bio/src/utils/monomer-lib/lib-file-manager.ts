/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLib} from './monomer-lib';
import {
  HELM_REQUIRED_FIELD as REQ, HELM_REQUIRED_FIELDS,
  HELM_VALUE_TYPE,
  HELM_FIELD_TYPE
} from '@datagrok-libraries/bio/src/utils/const';
import {JSONSchemaType} from 'ajv';
import {HELM_JSON_SCHEMA_PATH} from './consts';

import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';
import {MonomerLibFileValidator} from './monomer-lib-file-validator';

/** Singleton for adding, validation and reading of monomer library files.
 * All files **must** be aligned to the HELM standard before adding. */
export class MonomerLibFileManager {
  private constructor(
    private fileValidator: MonomerLibFileValidator
  ) {
    this._libraryFilesUpdateSubject$.pipe(
      debounceTime(3000)
    ).subscribe(async () => {
      await this.updateValidLibraryFileList();
    });
  }

  private getValidFilesList(): string[] {
    return this._libraryFilesUpdateSubject$.getValue();
  }

  private _libraryFilesUpdateSubject$ = new rxjs.BehaviorSubject<string[]>([]);

  private static instance: MonomerLibFileManager | undefined;

  static async getInstance(): Promise<MonomerLibFileManager> {
    if (MonomerLibFileManager.instance === undefined) {
      const helmSchemaRaw = await grok.dapi.files.readAsText(HELM_JSON_SCHEMA_PATH);
      const helmSchema = JSON.parse(helmSchemaRaw) as JSONSchemaType<any>;

      const fileValidator = new MonomerLibFileValidator(helmSchema);
      MonomerLibFileManager.instance = new MonomerLibFileManager(fileValidator);

      await MonomerLibFileManager.instance.initialize();
    }

    return MonomerLibFileManager.instance;
  }

  get debouncedMonomerLibFileListChange$(): rxjs.Observable<string[]> {
    return this._libraryFilesUpdateSubject$.pipe(
      debounceTime(1000)
    );
  }

  private async initialize(): Promise<void> {
    await this.updateValidLibraryFileList();
  }

  /** Add standard .json monomer library  */
  async addLibraryFile(fileContent: string, fileName: string): Promise<void> {
    if (await this.fileExists(fileName)) {
      grok.shell.error(`File ${fileName} already exists`);
      return;
    }

    await this.validateLibraryFileAgainstHELM(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    await this.updateValidLibraryFileList();
    const fileExists = await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
    if (!fileExists)
      grok.shell.error(`Failed to add ${fileName} library`);
    else
      grok.shell.info(`Added ${fileName} HELM library`);
  }

  private async fileExists(fileName: string): Promise<boolean> {
    return await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
  }

  async deleteLibraryFile(fileName: string): Promise<void> {
    try {
      await grok.dapi.files.delete(LIB_PATH + `${fileName}`);
      await this.updateValidLibraryFileList();
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
      if (!polymerTypes.includes(monomer[REQ.POLYMER_TYPE])) {
        monomers[monomer[REQ.POLYMER_TYPE]] = {};
        polymerTypes.push(monomer[REQ.POLYMER_TYPE]);
      }
      monomers[monomer[REQ.POLYMER_TYPE]][monomer[REQ.SYMBOL]] = monomer as Monomer;
    });

    return new MonomerLib(monomers);
  }

  getRelativePathsOfValidLibraryFiles(): string[] {
    return this.getValidFilesList();
  }

  /** Necessary to prevent sync errors  */
  async refreshLibraryFilePaths(): Promise<void> {
    await this.updateValidLibraryFileList();
  }

  private async updateValidLibraryFileList(): Promise<void> {
    const invalidFiles = [] as string[];
    // todo: remove after debugging
    console.log(`files before validation:`, this.getValidFilesList());
    const filePaths = await this.getFilePathsAtDefaultLocation();

    if (!this.fileListHasChanged(filePaths))
      return;

    for (const path of filePaths) {
      if (!path.endsWith('.json')) {
        invalidFiles.push(path);
        continue;
      }

      const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${path}`);
      if (!this.isValidHELMFormatLibrary(fileContent, path))
        invalidFiles.push(path);
    }

    const validLibraryPaths = filePaths.filter((path) => !invalidFiles.includes(path));

    if (this.fileListHasChanged(validLibraryPaths))
      this._libraryFilesUpdateSubject$.next(validLibraryPaths);
    console.log(`files after validation:`, this.getValidFilesList());

    // todo: remove after debugging
    if (validLibraryPaths.some((el) => !el.endsWith('.json')))
      console.warn(`Wrong validation: ${validLibraryPaths}`);

    if (invalidFiles.length > 0) {
      const message = `Invalid monomer library files in ${LIB_PATH}` +
      `, consider fixing or removing them: ${invalidFiles.join(', ')}`;

      console.warn(message);
      grok.shell.warning(message);
    }
  }

  private fileListHasChanged(newList: string[]): boolean {
    const currentList = this.getValidFilesList();
    return newList.length !== currentList.length || newList.some((el, i) => el !== currentList[i]);
  }

  private async validateLibraryFileAgainstHELM(fileContent: string, fileName: string): Promise<void> {
    const isValid = this.isValidHELMFormatLibrary(fileContent, fileName);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private isValidHELMFormatLibrary(fileContent: string, fileName: string): boolean {
    return this.fileValidator.validateFile(fileContent, fileName);
  }

  /** Get relative paths for files in LIB_PATH  */
  private async getFilePathsAtDefaultLocation(): Promise<string[]> {
    const list = await grok.dapi.files.list(LIB_PATH);
    const paths = list.map((fileInfo) => {
      return fileInfo.fullPath;
    });

    // WARNING: an extra sanity check,
    // caused by unexpected behavior of grok.dapi.files.list() when it returns non-existent paths
    const existingPaths = [] as string[];
    for (const path of paths) {
      const exists = await grok.dapi.files.exists(path);
      if (exists)
        existingPaths.push(path);
    }

    return existingPaths.map((path) => {
      // Get relative path (to LIB_PATH)
      return path.substring(LIB_PATH.length);
    });
  }
}
