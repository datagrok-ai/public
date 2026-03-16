/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {getKnimeClient} from './knime-client-factory';
import {getOrRegisterFunc} from './function-registry';
import '../css/knime-link.css';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.app({
    browsePath: 'Compute',
    name: 'KNIME',
  })
  static async knimeLinkApp(): Promise<DG.ViewBase> {
    const appHeader = u2.appHeader({
      iconPath: _package.webRoot + '/images/knime.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/KnimeLink/README.md',
      description: '- Connect to your KNIME Business Hub.\n' +
        '- Browse available deployments.\n' +
        '- Pass tables, parameters, and files as input.\n' +
        '- Execute workflows and view results as DataFrames.\n',
    });

    const view = DG.View.fromRoot(appHeader);
    view.name = 'KNIME';
    return view;
  }

  @grok.decorators.appTreeBrowser({app: 'KNIME'})
  static async knimeLinkAppTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
    let client;
    try {
      client = getKnimeClient();
    }
    catch (e: any) {
      grok.shell.warning('Set KNIME base URL in package settings (Manage > Plugins > KnimeLink)');
      return;
    }

    const setBreadcrumbs = (path: string[], tree: DG.TreeViewGroup) => {
      const breadcrumbs = ui.breadcrumbs(path);
      breadcrumbs.onPathClick.subscribe(async (value: string[]) => {
        const actualItem = value[value.length - 1];
        if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
          return;
        if (actualItem === 'KNIME')
          tree.currentItem = tree;
      });
      if (grok.shell.v) {
        if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') {
          const homeIcon = ui.iconFA('home', () => {
            grok.shell.v.close();
            grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME);
          });
          breadcrumbs.root.firstElementChild!.replaceWith(homeIcon);
        }
        const viewNameRoot = grok.shell.v.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
        if (viewNameRoot) {
          viewNameRoot.textContent = '';
          viewNameRoot.appendChild(breadcrumbs.root);
        }
      }
    };

    treeNode.items.length = 0;
    try {
      const deployments = await client.listDeployments('rest');
      if (deployments.length === 0) {
        treeNode.item('No deployments found');
        return;
      }
      for (const dep of deployments) {
        const node = treeNode.item(dep.name);
        node.onSelected.subscribe(async () => {
          if (!node.value) {
            setBreadcrumbs(['Home', 'KNIME', dep.name], treeNode);
            const func = await getOrRegisterFunc(dep, client);
            node.value = func;
          }
          const objHandler = DG.ObjectHandler.forEntity(node.value);
          if (objHandler)
            grok.shell.addPreview(await (objHandler.renderPreview(node.value)))
        });
      }
    }
    catch (e: any) {
      treeNode.item(`Error loading deployments: ${e?.message ?? e}`);
    }
  }
}