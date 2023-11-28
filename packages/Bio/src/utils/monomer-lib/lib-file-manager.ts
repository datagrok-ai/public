/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';


/** Class for adding, validation and (todo) reading of monomer library files  */
export class MonomerLibFileManager {
  constructor() { }

  /** Add standard .json monomer library  */
  async addStandardLibFile(fileContent: string, fileName: string): Promise<void> {
    this.validate(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    grok.shell.info(`Library ${fileName} added`);
  }

  /** Transform non-standad monomer librarieies to standard format */
  async addCustomLibFile(fileContent: string, fileName: string): Promise<void> {
    this.validate(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    grok.shell.info(`Library ${fileName} added`);
  }

  async deleteLibFile(fileName: string): Promise<void> {
    grok.dapi.files.delete(LIB_PATH + `${fileName}`);
    grok.shell.warning(`Library ${fileName} deleted`);
  }

  private validate(fileContent: string, fileName: string): void {
    const isValid = this.isValid(fileContent);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private isValid(fileContent: string): boolean {
    if (fileContent.length > 0)
      return true;
    return false;
  }
}
