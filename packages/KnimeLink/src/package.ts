/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {getKnimeClient} from './knime-client-factory';
import {getOrRegisterFunc} from './function-registry';
import {loadCachedEntries, registerFromCache, refreshAndUpdateCache} from './function-cache';
import {IKnimeClient} from './knime-client';
import type {KnimeDeployment} from './types';
import '../css/knime-link.css';

export * from './package.g';
export const _package = new DG.Package();

const TEST_PREFIX = 'datagrok_test_';

export class PackageFunctions {
  @grok.decorators.autostart({description: 'KnimeLink function registration'})
  static async knimeLinkAutostart(): Promise<void> {
    let client;
    try {
      client = getKnimeClient();
    }
    catch {
      return;
    }

    const cached = loadCachedEntries();
    if (cached.length > 0)
      registerFromCache(cached, client);

    try {
      await refreshAndUpdateCache(client);
    }
    catch (e: any) {
      grok.shell.error(`KnimeLink: background cache refresh failed: ${e}`);
    }
  }

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
    const errors: string[] = [];
    try {
      const deployments = await client.listDeployments('rest');
      const testDeployments: typeof deployments = [];
      const addDeploymentNode = (dep: typeof deployments[0], parent: DG.TreeViewGroup) => {
        const node = parent.item(dep.name);
        node.onSelected.subscribe(async () => {
          if (!node.value) {
            setBreadcrumbs(['Home', 'KNIME', dep.name], treeNode);
            const func = await getOrRegisterFunc(dep, client);
            node.value = func;
          }
          const objHandler = DG.ObjectHandler.forEntity(node.value);
          if (objHandler)
            grok.shell.preview = await objHandler.renderPreview(node.value);
          showWorkflowImage(dep, client);
        });
      };
      for (const dep of deployments) {
        try {
          if (dep.name.toLowerCase().startsWith(TEST_PREFIX))
            testDeployments.push(dep);
          else
            addDeploymentNode(dep, treeNode);
        }
        catch (e: any) {
          errors.push(`${dep.name}: ${e?.message ?? e}`);
        }
      }
      if (testDeployments.length > 0) {
        const testGroup = treeNode.group('Test workflows', undefined, false);
        for (const dep of testDeployments) {
          try { addDeploymentNode(dep, testGroup); }
          catch (e: any) { errors.push(`${dep.name}: ${e?.message ?? e}`); }
        }
      }
      if (treeNode.items.length === 0 && errors.length === 0)
        treeNode.item('No deployments found');
    }
    catch (e: any) {
      errors.push(e?.message ?? e);
    }
    if (errors.length > 0)
      grok.shell.warning(`KNIME: Failed to load some workflows:\n${errors.join('\n')}`);
  }
}

async function showWorkflowImage(dep: KnimeDeployment, client: IKnimeClient): Promise<void> {
  if (!dep.workflowId)
    return;
  const imageUrl = await client.getWorkflowImageUrl(dep.workflowId);
  if (!imageUrl)
    return;
  const header = ui.h2(`${dep.name} workflow preview`);
  header.className = 'knime-preview-header';
  const imgContainer = ui.div([], 'knime-workflow-svg-container');
  const container = ui.div([header, imgContainer], 'knime-workflow-preview');
  grok.shell.o = container;

  ui.setUpdateIndicator(imgContainer, true, 'Loading workflow image...');
  const img = document.createElement('img');
  img.className = 'knime-workflow-svg';
  img.alt = `Workflow: ${dep.name}`;
  img.onload = () => {
    ui.setUpdateIndicator(imgContainer, false);
    imgContainer.appendChild(img);
  };
  img.onerror = () => {
    ui.setUpdateIndicator(imgContainer, false);
    imgContainer.appendChild(ui.divText('Failed to load workflow image'));
  };
  img.src = imageUrl;
}