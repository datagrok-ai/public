/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions{
  @grok.decorators.app({
    icon: 'images/signals-icon.png',
    browsePath: 'Chem',
    name: 'Signals',
  })
  static async signalsApp(): Promise<DG.ViewBase> {
    const appHeader = u2.appHeader({
      iconPath: _package.webRoot + '/images/signals-icon.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/Signals/README.md',
      description: '- Connect to the Revvity Signals ELN.\n' +
        '- Analyze assay data.\n' +
        '- Find contextual information on molecules.\n'
    });

    const view = DG.View.fromRoot(appHeader);
    view.name = 'Signals';
    return view;
  }

  @grok.decorators.appTreeBrowser({app: 'Signals'})
  static async signalsAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    treeNode.group('Projects');
    treeNode.group('Experiments');
    treeNode.group('Equipment');
    treeNode.group('Requests');
    treeNode.group('Favorites');
  }

  @grok.decorators.func({name: 'Signals'})
  static async saveToSignalsEln(
    data: string) {
  
    grok.shell.info('Successfully saved to Signals ELN.');
  }
}