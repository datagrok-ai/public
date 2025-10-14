/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { renderExcalidraw } from './excal-component';
import './excal.css';

export * from './package.g';
export const _package = new DG.Package();
export class PackageFunctions {
  @grok.decorators.func({
    'name': 'info'
  })
  static info() {
    grok.shell.info(_package.webRoot);
  }


  @grok.decorators.app({
    'browsePath': 'Misc',
    'icon': 'images/excal.png',
    'outputs': [
      {
        'name': 'v',
        'type': 'view'
      }
    ],
    'name': 'Excalidraw'
  })
  static excalidrawApp() {
    const div = renderExcalidraw();
    const v= DG.View.fromRoot(div.actualRoot);
    v.name = 'Excalidraw';
    return v;
  }

  @grok.decorators.fileViewer({
    'fileViewer': 'excalidraw',
    'meta': {
      'ext': 'excalidraw'
    },
    'name': 'excalfileViewer',
    'description': 'Excalidraw viewer'
  })
  static async excalfileViewer(fileContent: DG.FileInfo) : Promise<DG.ViewBase | void> {
  
    if (!fileContent) {
      grok.shell.info('No file content provided');
      return;
    }
    const fileC = await fileContent.readAsString();
    if (!fileC) {
      grok.shell.info('No file content provided');
      return;
    }
    const json = JSON.parse(fileC);
    if (!json) {
      grok.shell.info('No JSON content provided');
      return;
    }
    const div = renderExcalidraw();
    const v = DG.View.fromRoot(div.actualRoot);
    v.name = fileContent.name ?? 'Excalidraw';
    setTimeout(() => {
      const api = div.getApi();
      api.updateScene(json);

    }, 500);
    return v;
  }
}