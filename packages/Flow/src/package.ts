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
      grok.shell.windows.showBrowse = false;
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      grok.shell.windows.showBrowse = true;
      grok.shell.windows.showToolbox = true;
    }, 200);
    return new FuncFlowView();
  }

  @grok.decorators.fileViewer({fileViewer: 'ffjson'})
  static viewFuncFlow(file: DG.FileInfo): DG.ViewBase {
    const view = new FuncFlowView();
    file.readAsString().then((json) => view.loadFromJson(json));
    setTimeout(() => {
      grok.shell.windows.showToolbox = false;
      grok.shell.windows.showBrowse = false;
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      grok.shell.windows.showBrowse = true;
      grok.shell.windows.showToolbox = true;
    }, 200);
    return view;
  }

  /** Builds a flow from a table-creation script (the function-call cascade
   *  Datagrok records for reproducibly-created tables, used by data sync)
   *  and opens it in the Flow editor. */
  @grok.decorators.func({
    name: 'flowFromCreationScript',
    description: 'Builds a flow diagram from a table creation script and opens it in the Flow editor',
  })
  static async flowFromCreationScript(script: string): Promise<DG.ViewBase> {
    const view = new FuncFlowView();
    await view.loadFromCreationScript(script);
    return view;
  }

  @grok.decorators.func({
    name: 'openCreationScriptFlowDialog',
    meta: {role: 'creationScriptEditor'},
  })
  static async openCreationScriptFlowDialog(script: string, show: boolean = true): Promise<DG.Dialog> {
    const view = new FuncFlowView();
    view.name = `Creation Script`;
    try {
      await view.loadFromCreationScript(script);
    } catch (e) {
      grok.shell.error(`Failed to load flow from creation script`);
      console.error(e);
    }
    const d = ui.dialog({title: 'Creation Script Flow'})
      .add(view.root)
      .addButton('Open In Editor', () => {
        grok.shell.addView(view);
        d.close();
      });
    if (show)
      d.show({resizable: true, width: 800, height: 600});
    return d;
  }
}
