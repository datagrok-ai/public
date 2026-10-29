import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {readJmp} from './formats/jmp/jmp-parser';
import {JmpFormatError, JmpTable} from './formats/jmp/jmp-types';
import {jmpTableToDataFrame} from './formats/jmp/jmp-to-dataframe';
import {createJmpErrorView, createJmpTableView, warnSkipped} from './formats/jmp/jmp-viewer';

export const _package = new DG.Package();

export * from './package.g';

export class JmpPackageFunctions {
  @grok.decorators.fileHandler({ext: 'jmp', description: 'Opens JMP data table'})
  static async importJmp(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array,
  ): Promise<DG.DataFrame[]> {
    const table = await readJmp(new Uint8Array(bytes));
    warnSkipped(table);
    return [jmpTableToDataFrame(table)];
  }

  /** Browse previews a file with this on a click and keeps the view on a double-click; a file from disk
   * (dropped, or picked in File > Open) comes here too, since the shell prefers a file viewer to the
   * handler. Either way the file opens as its table, with its saved analyses; one that isn't a readable
   * JMP table shows why. */
  @grok.decorators.fileViewer({fileViewer: 'jmp'})
  static async previewJmp(file: DG.FileInfo): Promise<DG.ViewBase> {
    let table: JmpTable;
    try {
      table = await readJmp(await file.readAsBytes());
    } catch (e) {
      if (e instanceof JmpFormatError)
        return createJmpErrorView(file.name, e.message);
      throw e;
    }
    return createJmpTableView(table, file.name);
  }
}
