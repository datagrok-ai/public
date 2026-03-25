import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parseMwx} from './formats/mwx/mwx-parser';
import {mwxWorksheetToDataFrame} from './formats/mwx/mwx-to-dataframe';
import {buildMwxView} from './formats/mwx/mwx-viewer';
import {parseMpx} from './formats/mpx/mpx-parser';
import {mpxProjectToDataFrames} from './formats/mpx/mpx-to-dataframe';
import {buildMpxView, openMpxProject} from './formats/mpx/mpx-viewer';

export const _package = new DG.Package();

export * from './package.g';


export class MinitabPackageFunctions {
  @grok.decorators.fileHandler({ext: 'mwx', description: 'Opens Minitab Worksheet file'})
  static async importMwx(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array
  ): Promise<DG.DataFrame[]> {
    const ws = await parseMwx(new Uint8Array(bytes));
    return [mwxWorksheetToDataFrame(ws)];
  }

  @grok.decorators.fileViewer({fileViewer: 'mwx'})
  static async previewMwx(file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.name;

    const bytes = await file.readAsBytes();
    const ws = await parseMwx(bytes);

    if (ws.columns.length === 0) {
      view.append(ui.divText('No data found in the .mwx file.'));
      return view;
    }

    const content = buildMwxView(ws);
    content.style.width = '100%';
    content.style.height = '100%';
    view.append(content);

    return view;
  }

  @grok.decorators.fileHandler({ext: 'mpx', description: 'Opens Minitab Project file'})
  static async importMpx(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array
  ): Promise<DG.DataFrame[]> {
    const project = await parseMpx(new Uint8Array(bytes));
    return mpxProjectToDataFrames(project);
  }

  @grok.decorators.fileViewer({fileViewer: 'mpx'})
  static async previewMpx(file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.name;

    const bytes = await file.readAsBytes();
    const project = await parseMpx(bytes);

    if (project.worksheets.length === 0) {
      view.append(ui.divText('No worksheets found in the .mpx file.'));
      return view;
    }

    const content = buildMpxView(project, () => openMpxProject(project));
    content.style.width = '100%';
    content.style.height = '100%';
    view.append(content);

    return view;
  }
}
