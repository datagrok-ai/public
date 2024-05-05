import * as DG from 'datagrok-api/dg';

import {PolyToolMonomerLibHandler} from './monomer-lib-handler';

/** Handler of custom monomer libs for PolyTool  */
export class PolyToolCsvLibHandler {
  constructor(private fileName: string, private fileContent: string) {
    this.validateFileType();
    const df = DG.DataFrame.fromCsv(this.fileContent);
    const json = this.toJson(df);
    this.polyToolMonomerLib = new PolyToolMonomerLibHandler(json);
    this.validateContent();
  }

  private polyToolMonomerLib: PolyToolMonomerLibHandler;

  async getJson(): Promise<any> {
    const rawLibData = this.polyToolMonomerLib.getJsonMonomerLib();
    return rawLibData;
  }

  private toJson(df: DG.DataFrame): any[] {
    return Array.from({length: df.rowCount}, (_, idx) =>
      df.columns.names().reduce((entry: { [key: string]: any }, colName) => {
        entry[colName] = df.get(colName, idx);
        return entry;
      }, {})
    );
  }

  private validateFileType(): void {
    if (!this.fileName.endsWith('.csv'))
      throw new Error(`File ${this.fileName} is not an CSV file`);
  }

  private validateContent(): void {
    if (!this.polyToolMonomerLib.isValid())
      throw new Error('Invalid format of CSV monomer lib');
  }
}
