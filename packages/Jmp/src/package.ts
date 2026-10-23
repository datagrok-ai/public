import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {readJmp} from './formats/jmp/jmp-parser';
import {JmpFormatError, JmpTable} from './formats/jmp/jmp-types';
import {addSavedAnalyses, jmpTableToDataFrame, savedAnalyses, skippedColumns} from './formats/jmp/jmp-to-dataframe';
import {createJmpErrorView, createJmpPreview} from './formats/jmp/jmp-viewer';

export const _package = new DG.Package();

export * from './package.g';

const list = (names: string[]): string => `${names.slice(0, 5).join(', ')}${names.length > 5 ? ', …' : ''}`;

/** Says which columns and saved analyses didn't come across, so nothing disappears silently. */
function warnSkipped(table: JmpTable): void {
  const names = skippedColumns(table).map((s) => s.name);
  if (names.length > 0)
    grok.shell.warning(`JMP: ${names.length} column${names.length === 1 ? '' : 's'} not imported: ${list(names)}`);
  const analyses = savedAnalyses(table).filter((a) => a.viewers.length === 0).map((a) => a.name);
  if (analyses.length > 0) {
    grok.shell.warning(`JMP: ${analyses.length} saved ${analyses.length === 1 ? 'analysis' : 'analyses'} ` +
      `not converted: ${list(analyses)}`);
  }
}

export class JmpPackageFunctions {
  @grok.decorators.fileHandler({ext: 'jmp', description: 'Opens JMP data table'})
  static async importJmp(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array,
  ): Promise<DG.DataFrame[]> {
    const table = await readJmp(new Uint8Array(bytes));
    warnSkipped(table);
    return [jmpTableToDataFrame(table)];
  }

  /** A file in the Browse tree gets the preview. A file from disk (dropped, or picked in File > Open)
   * comes here too, since the shell prefers a file viewer to the handler; it has no connection, and
   * opens straight as a table, as every other table format does, with its saved analyses. A file that
   * isn't a readable JMP table shows why, in either case. */
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
    if (file.connection == null) {
      warnSkipped(table);
      const tv = DG.TableView.create(jmpTableToDataFrame(table, file.name.replace(/\.jmp$/i, '')), false);
      addSavedAnalyses(tv, savedAnalyses(table));
      return tv;
    }
    return createJmpPreview(table, file.name);
  }
}
