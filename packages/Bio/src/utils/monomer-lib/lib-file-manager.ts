/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';


/** Singleton for adding, validation and reading of monomer library files.
 * All files stored under LIB_PATH directory must be in .json format and satisfy HELM standard.
 * All libraries  that are not in .json format or do not satisfy HELM standard are considered as custom libraries.
 * They must be brought to the standard format before adding to LIB_PATH. */
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

  /** Transform non-standad monomer librarieies to standard format */
  async addCustomLibFile(fileContent: string, fileName: string): Promise<void> {
    await this.validate(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    grok.shell.info(`Added ${fileName} HELM library`);
  }

  async deleteLibFile(fileName: string): Promise<void> {
    grok.dapi.files.delete(LIB_PATH + `${fileName}`);
    grok.shell.warning(`Deleted ${fileName} library`);
    await this.validateAllFiles();
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
      const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${file.name}`);
      if (!this.isValid(fileContent))
        invalidFiles.push(file.name);
    }

    if (invalidFiles.length > 0)
      grok.shell.warning(`${invalidFiles.join(', ')} do not satisfy HELM standard}, consider deleting/fixing them`);
  }

  private isValid(fileContent: string): boolean {
    if (fileContent.length > 0)
      return true;
    return false;
  }
}
