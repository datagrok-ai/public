import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parseTwbFile} from './tableau/tableau-parser';
import {twbDatasourceToDataFrame} from './tableau/tableau-to-dataframe';
import {buildTableauView, openAllDatasources} from './tableau/tableau-viewer';


export const _package = new DG.Package();

export * from './package.g';


export class TableauPackageFunctions {
  @grok.decorators.fileHandler({ext: 'twb', description: 'Opens Tableau workbook file'})
  static importTwb(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array
  ): DG.DataFrame[] | void {
    try {
      const text = new TextDecoder('utf-8').decode(new Uint8Array(bytes));
      const twbFile = parseTwbFile(text);
      const results: DG.DataFrame[] = [];
      for (const ds of twbFile.datasources)
        results.push(twbDatasourceToDataFrame(ds));
      return results;
    }
    catch (e: any) {
      grok.shell.warning('Tableau file is not supported or malformed');
      grok.shell.error(e);
    }
  }

  @grok.decorators.fileViewer({fileViewer: 'twb'})
  static async previewTwb(file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.name;

    const bytes = await file.readAsBytes();
    const text = new TextDecoder('utf-8').decode(bytes);
    const twbFile = parseTwbFile(text);

    if (twbFile.datasources.length === 0) {
      view.append(ui.divText('No data found in the .twb file.'));
      return view;
    }

    view.setRibbonPanels([[ui.bigButton('Open All', () => openAllDatasources(twbFile))]]);

    const content = buildTableauView(twbFile);
    content.style.width = '100%';
    content.style.height = '100%';
    view.append(content);

    return view;
  }
}
