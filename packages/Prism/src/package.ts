import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parsePrismFile} from './prism/prism-parser';
import {prismSheetToDataFrame, prismAnalysisToDataFrame} from './prism/prism-to-dataframe';
import {hasXYData, prismToFitDataFrame} from './prism/prism-curves';
import {buildPrismView} from './prism/prism-viewer';


export const _package = new DG.Package();

export * from './package.g';


export class PrismPackageFunctions {
  @grok.decorators.fileHandler({ext: 'prism', description: 'Opens GraphPad Prism file'})
  static async importPrism(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array
  ): Promise<DG.DataFrame[]> {
    const prismFile = await parsePrismFile(new Uint8Array(bytes));
    const results: DG.DataFrame[] = [];

    for (const sheet of prismFile.sheets)
      results.push(prismSheetToDataFrame(sheet));

    // For XY sheets, produce a curves DataFrame with fit column
    const xySheets = prismFile.sheets.filter((s) => hasXYData(s));
    if (xySheets.length > 0)
      results.push(prismToFitDataFrame(xySheets));

    for (const analysis of prismFile.analyses)
      results.push(prismAnalysisToDataFrame(analysis));

    return results;
  }

  @grok.decorators.fileViewer({fileViewer: 'prism'})
  static async previewPrism(file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.name;

    const bytes = await file.readAsBytes();
    const prismFile = await parsePrismFile(bytes);

    if (prismFile.sheets.length === 0 && prismFile.analyses.length === 0) {
      view.append(ui.divText('No data found in the .prism file.'));
      return view;
    }

    const content = buildPrismView(prismFile);
    content.style.width = '100%';
    content.style.height = '100%';
    view.append(content);

    return view;
  }
}
