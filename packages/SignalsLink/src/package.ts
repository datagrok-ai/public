/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AppHeader} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions{
  @grok.decorators.app({
    icon: 'images/signals-icon.png',
    browsePath: 'Chem',
    name: 'Signals',
  })
  static async signalsApp(): Promise<DG.ViewBase> {
    const view = DG.View.create();
    const appHeader = ViewHost.of(view).build(() => new AppHeader({
      icon: _package.webRoot + '/images/signals-icon.png',
      learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/Signals/README.md',
      items: [
        'Connect to the Revvity Signals ELN.',
        'Analyze assay data.',
        'Find contextual information on molecules.',
      ],
    }));

    view.root.append(appHeader.root);
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