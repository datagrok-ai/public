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
    this._onFileListChange.subscribe(async () => {
      await this.validateAllFiles();
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

  get onFileListChange(): rxjs.Subject<void> {
    return this.onFileListChange;
  }

  private init(): void {
    this._onFileListChange.next();
  }

  /** Add standard .json monomer library  */
  async addLibFile(fileContent: string, fileName: string): Promise<void> {
    this._onFileListChange.next();
    await this.validateFile(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    this._onFileListChange.next();
    grok.shell.info(`Added ${fileName} HELM library`);
  }

  async deleteLibFile(fileName: string): Promise<void> {
    grok.dapi.files.delete(LIB_PATH + `${fileName}`);
    grok.shell.warning(`Deleted ${fileName} library`);
    this._onFileListChange.next();
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

  private async validateFile(fileContent: string, fileName: string): Promise<void> {
    const isValid = this.isValid(fileContent);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private async validateAllFiles(): Promise<void> {
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
    // valid files are difference of filenames in the list and the invalid
    this.validFiles = filePaths.filter((path) => !this.invalidFiles.includes(path));

    if (this.invalidFiles.length > 0) {
      console.log(`invalid files: ${this.invalidFiles}`);
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
    // check if jsonContent is an array
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
    return list.map((fileInfo) => {
      // Get relative path (to LIB_PATH)
      return fileInfo.fullPath.substring(LIB_PATH.length);
    });
  }
}
