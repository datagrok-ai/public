/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLib} from './monomer-lib';
import {HELM_REQUIRED_FIELDS as REQ} from '@datagrok-libraries/bio/src/utils/const';

import * as rxjs from 'rxjs';

/** Singleton for adding, validation and reading of monomer library files.
 * All files **must** be aligned to HELM standard before adding. */
export class MonomerLibFileManager {
  private constructor() {
    DG.debounce<void>(this._onFileListChange, 3000).subscribe(async () => {
      await this.updateFilePaths();
    });
  }

  private validFiles: string[] = [];

  /** Tracks invalid files accidentally added to the library folder manually */
  private invalidFiles: string[] = [];

  private _onFileListChange = new rxjs.Subject<void>();

  private static instance: MonomerLibFileManager | undefined;

  static async getInstance(): Promise<MonomerLibFileManager> {
    if (MonomerLibFileManager.instance === undefined) {
      MonomerLibFileManager.instance = new MonomerLibFileManager();
      MonomerLibFileManager.instance.init();
    }
    return MonomerLibFileManager.instance;
  }

  get onFileListChange(): rxjs.Observable<void> {
    return this._onFileListChange.asObservable();
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
    this._onFileListChange.next();
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
      this._onFileListChange.next();
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
    this.invalidFiles = [];
    const filePaths = await this.getFilePaths();
    for (const path of filePaths) {
      if (!path.endsWith('.json')) {
        this.invalidFiles.push(path);
        continue;
      }
      const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${path}`);
      if (!this.isValid(fileContent))
        this.invalidFiles.push(path);
    }
    this.validFiles = filePaths.filter((path) => !this.invalidFiles.includes(path));
  }

  private async validateFile(fileContent: string, fileName: string): Promise<void> {
    const isValid = this.isValid(fileContent);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private async validateAllFiles(): Promise<void> {
    await this.updateFilePaths();

    if (this.invalidFiles.length > 0) {
      console.warn(
        `Invalid monomer library files in ${LIB_PATH}, consider fixing or removing them: ${this.invalidFiles}`);
      grok.shell.warning(
        `The following files in ${LIB_PATH} do not satisfy` +
        ` HELM standard and are not displayed in the list of available libraries:` +
        ` ${this.invalidFiles.join(', ')}` +
        `. Fix or delete them.`
      );
    }
  }

  /** The file **must** strictly satisfy HELM standard */
  private isValid(fileContent: string): boolean {
    const jsonContent = JSON.parse(fileContent);
    if (!Array.isArray(jsonContent))
      return false;
    const requiredFields = [
      REQ.AUTHOR, REQ.CREATE_DATE, REQ.ID, REQ.MOLFILE,
      REQ.MONOMER_TYPE, REQ.NAME, REQ.POLYMER_TYPE, REQ.RGROUPS, REQ.SMILES, REQ.SYMBOL
    ];
    const result = jsonContent.every((monomer: any) => {
      return requiredFields.every((field) => {
        return Object.keys(monomer).includes(field);
      });
    });
    return result;
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
