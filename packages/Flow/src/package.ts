/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

import {FuncFlowView} from './funcflow-view';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

export class PackageFunctions {
  @grok.decorators.app({
    name: 'Flow',
    description: 'Interactive function chain designer',
    tags: ['app'],
  })
  static funcflowApp(@grok.decorators.param({options: {metaUrl: true, optional: true}}) path?: string): DG.ViewBase {
    const url = new URL(window.location.href);
    const params = url.searchParams;
    console.log(params);
    setTimeout(() => {
      grok.shell.windows.showToolbox = false;
      grok.shell.windows.showContextPanel = false;
      grok.shell.windows.showHelp = false;
    }, 200);
    return new FuncFlowView();
  }

 @grok.decorators.fileViewer({fileViewer: 'ffjson'})
  static viewFuncFlow(file: DG.FileInfo): DG.ViewBase {
    const view = new FuncFlowView();
    file.readAsString().then((json) => view.loadFromJson(json));
    setTimeout(() => {
      grok.shell.windows.showToolbox = false;
      grok.shell.windows.showContextPanel = false;
      grok.shell.windows.showHelp = false;
    }, 200);
    return view;
  }
}
