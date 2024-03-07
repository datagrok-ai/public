/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {createJsonMonomerLibFromSdf} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

interface CustomMonomerLibHandler {
  getHelmLibContent(): Promise<string>;
}

/** Casts specific SDF files as HELM-standard monomer libraries */
export class SdfMonomerLibHandler implements CustomMonomerLibHandler {
  constructor(private fileName: string, private fileContent: string) {
    this.validateChemInstallation();
    this.validate();
  }

  async getHelmLibContent(): Promise<string> {
    const bytes = this.transformToBytes();
    const dfSdf = await grok.functions.call('Chem:importSdf', {bytes});
    const jsonContent = createJsonMonomerLibFromSdf(dfSdf[0]);
    return JSON.stringify(jsonContent);
  }

  private transformToBytes(): Uint8Array {
    return new TextEncoder().encode(this.fileContent);
  }

  private validateChemInstallation(): void {
    const funcList: DG.Func[] = DG.Func.find({package: 'Chem', name: 'importSdf'});
    if (funcList.length === 1)
      return;
    throw new Error('MonomerLibFileManager: Chem package is not installed, cannot convert SDF to JSON');
  }

  private validate(): void {
    if (!this.fileName.endsWith('.sdf'))
      throw new Error(`File ${this.fileName} is not an SDF file`);
  }
}
