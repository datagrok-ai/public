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

import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';

const STRING_TYPE = 'string';

/** Singleton for adding, validation and reading of monomer library files.
 * All files **must** be aligned to the HELM standard before adding. */
export class MonomerLibFileManager {
  private constructor() {
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
      MonomerLibFileManager.instance = new MonomerLibFileManager();
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
      if (!this.isValidHELMFormatLibrary(fileContent))
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
    const isValid = this.isValidHELMFormatLibrary(fileContent);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  /** The file **must** strictly satisfy HELM standard */
  private isValidHELMFormatLibrary(fileContent: string): boolean {
    return new MonomerLibFileValidator().validate(fileContent);
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

class MonomerLibFileValidator {
  validate(fileContent: string): boolean {
    return this.validateHELMFormat(fileContent);
  }

  private validateHELMFormat(fileContent: string): boolean {
    let jsonContent: any[];
    try {
      jsonContent = JSON.parse(fileContent);
    } catch (e) {
      console.error('Bio: Monomer Library File Manager: Invalid JSON format:', e);
      return false;
    }

    if (!Array.isArray(jsonContent))
      return false;

    return jsonContent.every((monomer) => this.isValidHELMMonomerObject(monomer));
  }

  private isValidHELMMonomerObject(monomer: any): boolean {
    for (const field of HELM_REQUIRED_FIELDS) {
      const fieldType = HELM_FIELD_TYPE[field];

      if (!monomer.hasOwnProperty(field))
        return false;

      if (field.toLowerCase() === REQ.RGROUPS.toLowerCase() && !this.isValidRGroupsField(monomer[field]))
        return false;

      if (typeof fieldType === STRING_TYPE && !this.matchesValueType(monomer[field], fieldType as string))
        return false;
    }
    return true;
  }

  private isValidRGroupsField(rgroups: any[]): boolean {
    if (!Array.isArray(rgroups)) return false;

    return rgroups.every((rgroup) => {
      const fieldType = HELM_FIELD_TYPE[REQ.RGROUPS] as any;
      const itemsType = fieldType.itemsType as Record<string, string>;
      return Object.entries(itemsType).every(([field, type]) => {
        const hasField = rgroup.hasOwnProperty(field);
        const matchesType = this.matchesValueType(rgroup[field], type);

        return hasField && matchesType;
      });
    });
  }

  private matchesValueType(value: any, typeInfo: string): boolean {
    switch (typeInfo) {
      case HELM_VALUE_TYPE.STRING:
        return typeof value === STRING_TYPE;
      case HELM_VALUE_TYPE.STRING_OR_NULL:
        return typeof value === STRING_TYPE || value === null;
      case HELM_VALUE_TYPE.INTEGER:
        return Number.isInteger(value);
      case HELM_VALUE_TYPE.ARRAY:
        return Array.isArray(value);
      default:
        return false;
    }
  }
}
