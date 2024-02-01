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

/** Singleton for adding, validation and reading of monomer library files.
 * All files **must** be aligned to the HELM standard before adding. */
export class MonomerLibFileManager {
  private constructor() {
    DG.debounce<void>(this._monomerLibFileListChange$, 3000).subscribe(async () => {
      await this.updateFilePaths();
    });
  }

  private validFiles: string[] = [];

  private _monomerLibFileListChange$ = new rxjs.Subject<void>();

  private static instance: MonomerLibFileManager | undefined;

  static async getInstance(): Promise<MonomerLibFileManager> {
    if (MonomerLibFileManager.instance === undefined) {
      MonomerLibFileManager.instance = new MonomerLibFileManager();
      await MonomerLibFileManager.instance.init();
    }
    return MonomerLibFileManager.instance;
  }

  get monomerLibFileListChange$(): rxjs.Observable<void> {
    return this._monomerLibFileListChange$.asObservable();
  }

  private async init(): Promise<void> {
    await this.validateAllFiles();
  }

  /** Add standard .json monomer library  */
  async addLibFile(fileContent: string, fileName: string): Promise<void> {
    if (await this.fileExists(fileName)) {
      grok.shell.error(`File ${fileName} already exists`);
      return;
    }

    await this.validateFile(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    await this.validateAllFiles();
    this._monomerLibFileListChange$.next();
    const fileExists = await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
    if (!fileExists)
      grok.shell.error(`Failed to add ${fileName} library`);
    else
      grok.shell.info(`Added ${fileName} HELM library`);
  }

  private async fileExists(fileName: string): Promise<boolean> {
    return await grok.dapi.files.exists(LIB_PATH + `${fileName}`);
  }

  async deleteLibFile(fileName: string): Promise<void> {
    try {
      await grok.dapi.files.delete(LIB_PATH + `${fileName}`);
      await this.validateAllFiles();
      this._monomerLibFileListChange$.next();
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

  async readLibraryFile(path: string, fileName: string): Promise<IMonomerLib> {
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

  getRelativePathsOfValidFiles(): string[] {
    return this.validFiles;
  }

  /** Necessary to prevent sync errors  */
  async refreshValidFilePaths(): Promise<void> {
    await this.updateFilePaths();
  }

  private async updateFilePaths(): Promise<void> {
    const invalidFiles = [] as string[];
    // todo: remove after debugging
    console.log(`files before validation:`, this.validFiles);
    const filePaths = await this.getFilePaths();
    for (const path of filePaths) {
      if (!path.endsWith('.json')) {
        invalidFiles.push(path);
        continue;
      }
      const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${path}`);
      if (!this.isValid(fileContent))
        invalidFiles.push(path);
    }
    this.validFiles = filePaths.filter((path) => !invalidFiles.includes(path));

    // todo: remove after debugging
    if (this.validFiles.some((el) => !el.endsWith('.json')))
      console.warn(`Wrong validation: ${this.validFiles}`);

    if (invalidFiles.length > 0) {
      const message = `Invalid monomer library files in ${LIB_PATH}` +
      `, consider fixing or removing them: ${invalidFiles.join(', ')}`;
      console.warn(message);
      grok.shell.warning(message);
    }
  }

  private async validateFile(fileContent: string, fileName: string): Promise<void> {
    const isValid = this.isValid(fileContent);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private async validateAllFiles(): Promise<void> {
    await this.updateFilePaths();
  }

  /** The file **must** strictly satisfy HELM standard */
  private isValid(fileContent: string): boolean {
    return new MonomerLibFileValidator().validate(fileContent);
  }

  /** Get relative paths for files in LIB_PATH  */
  private async getFilePaths(): Promise<string[]> {
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
    let jsonContent: any[];
    try {
      jsonContent = JSON.parse(fileContent);
    } catch (e) {
      console.error('Bio: Monomer Library File Manager: Invalid JSON format:', e);
      return false;
    }

    if (!Array.isArray(jsonContent))
      return false;

    return jsonContent.every((monomer) => this.validateMonomer(monomer));
  }

  private validateMonomer(monomer: any): boolean {
    for (const field of HELM_REQUIRED_FIELDS) {
      const fieldType = HELM_FIELD_TYPE[field];

      if (!monomer.hasOwnProperty(field))
      {
        console.log('1', monomer);
        console.log(field);
        return false;
      }

      if (field.toLowerCase() === REQ.RGROUPS.toLowerCase() && !this.validateRGroups(monomer[field]))
      {
        console.log('2', monomer);
        console.log(field);
        return false;
      }

      if (typeof fieldType === 'string' && !this.matchesType(monomer[field], fieldType as string))
      {
        console.log('3', monomer);
        console.log(field);
        return false;
      }
    }
    return true;
  }

  private validateRGroups(rgroups: any[]): boolean {
    if (!Array.isArray(rgroups)) return false;

    return rgroups.every((rgroup) => {
      const fieldType = HELM_FIELD_TYPE[REQ.RGROUPS] as any;
      const itemsType = fieldType.itemsType as Record<string, string>;
      return Object.entries(itemsType).every(([field, type]) => {
        // WARNING: toLowerCase is necessary because HELMCoreLibrary has "capGroupSMILES" and "capGroupSmiles"
        const hasField = Object.keys(rgroup).map((key) => key.toLowerCase()).some((key) => key === field.toLowerCase());

        const result = hasField && this.matchesType(rgroup[field], type);
        console.log(rgroup, field, type, result);
        return result;
      });
    });
  }

  private matchesType(value: any, typeInfo: string): boolean {
    switch (typeInfo) {
      case HELM_VALUE_TYPE.STRING:
        return typeof value === 'string';
      case HELM_VALUE_TYPE.STRING_OR_NULL:
        return typeof value === 'string' || value === null;
      case HELM_VALUE_TYPE.INTEGER:
        return Number.isInteger(value);
      case HELM_VALUE_TYPE.ARRAY:
        return Array.isArray(value);
      default:
        return false;
    }
  }
}
