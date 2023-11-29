/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLib} from './monomer-lib';
import {HELM_REQUIRED_FIELDS as REQ} from '@datagrok-libraries/bio/src/utils/const';


/** Singleton for adding, validation and reading of monomer library files.
 * All files **must** be aligned to HELM standard before adding. */
export class MonomerLibFileManager {
  private constructor() { }

  private static instance: MonomerLibFileManager | undefined;

  static async getInstance(): Promise<MonomerLibFileManager> {
    if (MonomerLibFileManager.instance === undefined) {
      MonomerLibFileManager.instance = new MonomerLibFileManager();
      await MonomerLibFileManager.instance.init();
    }
    return MonomerLibFileManager.instance;
  }

  private async init(): Promise<void> {
    this.validateAllFiles();
  }

  /** Add standard .json monomer library  */
  async addStandardLibFile(fileContent: string, fileName: string): Promise<void> {
    await this.validate(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    grok.shell.info(`Added ${fileName} HELM library`);
  }

  async deleteLibFile(fileName: string): Promise<void> {
    grok.dapi.files.delete(LIB_PATH + `${fileName}`);
    grok.shell.warning(`Deleted ${fileName} library`);
    await this.validateAllFiles();
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

  private async validate(fileContent: string, fileName: string): Promise<void> {
    await this.validateAllFiles();
    const isValid = this.isValid(fileContent);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private async validateAllFiles(): Promise<void> {
    const list = await grok.dapi.files.list(LIB_PATH);
    const invalidFiles: string[] = [];
    for (const file of list) {
      if (!file.name.endsWith('.json')) {
        invalidFiles.push(file.name);
        continue;
      }
      const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${file.name}`);
      if (!this.isValid(fileContent))
        invalidFiles.push(file.name);
    }

    if (invalidFiles.length > 0)
      grok.shell.warning(`The following files in ${LIB_PATH} do not satisfy HELM standard}, delete or fix them: ${invalidFiles.join(', ')}`);
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
}
