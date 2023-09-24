/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class HelmToMolfileConverter {
  constructor(private helmColumn: DG.Column<string> ) {
    this.helmColumn = helmColumn;
  }

  async convert(): Promise<DG.Column<string>> {
    const polymerMolfileCol: DG.Column<string> = await this.getPolymerMolfiles();
    return polymerMolfileCol;
  }

  async getPolymerMolfiles(): Promise<DG.Column<string>> {
    const polymerMolfileCol: DG.Column<string> = await grok.functions.call('HELM:getMolfiles', {col: this.helmColumn});
    return polymerMolfileCol;
  }
}
