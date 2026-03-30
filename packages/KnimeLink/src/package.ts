/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {getKnimeClient} from './knime-client-factory';
import {getOrRegisterFunc, sanitizeFuncName} from './function-registry';
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

    const addNode = (dep: KnimeDeployment, parent: DG.TreeViewGroup, func?: DG.Func): DG.TreeViewNode => {
      const node = parent.item(dep.name);
      if (func)
        node.value = func;
      node.onSelected.subscribe(async () => {
        setBreadcrumbs(['Home', 'KNIME', dep.name], treeNode);
        if (!node.value) {
          const f = await getOrRegisterFunc(dep, client!);
          node.value = f;
        }
        const objHandler = DG.ObjectHandler.forEntity(node.value);
        if (objHandler)
          grok.shell.preview = await objHandler.renderPreview(node.value);
        showWorkflowImage(dep, client!);
      });
      return node;
    };

    const removeStaleItems = (cached: Map<string, DG.TreeViewNode>, liveNames: Set<string>) => {
      for (const [name, node] of cached) {
        if (!liveNames.has(name))
          node.remove();
      }
    };

    const updateOrAddNodes = async (
      deployments: KnimeDeployment[], parent: DG.TreeViewGroup,
      cached: Map<string, DG.TreeViewNode>,
    ) => {
      for (const dep of deployments) {
        try {
          const existingNode = cached.get(dep.name);
          if (existingNode)
            existingNode.value = await getOrRegisterFunc(dep, client!);
          else
            addNode(dep, parent);
        }
        catch (e: any) {
          errors.push(`${dep.name}: ${e?.message ?? e}`);
        }
      }
    };

    // Phase 1: Create tree items from cached functions
    const cached = loadCachedEntries();
    const cachedItemsByName = new Map<string, DG.TreeViewNode>();
    const cachedTestItemsByName = new Map<string, DG.TreeViewNode>();
    let testGroup: DG.TreeViewGroup | null = null;

    const cachedTest: typeof cached = [];
    for (const entry of cached) {
      if (entry.deployment.name.toLowerCase().startsWith(TEST_PREFIX)) {
        cachedTest.push(entry);
        continue;
      }
      const funcName = sanitizeFuncName(entry.deployment.name, entry.deployment.id);
      const existing = DG.Func.find({package: 'KnimeLink', name: funcName});
      cachedItemsByName.set(entry.deployment.name,
        addNode(entry.deployment, treeNode, existing.length > 0 ? existing[0] : undefined));
    }
    if (cachedTest.length > 0) {
      testGroup = treeNode.group('Test workflows', undefined, false);
      for (const entry of cachedTest) {
        const funcName = sanitizeFuncName(entry.deployment.name, entry.deployment.id);
        const existing = DG.Func.find({package: 'KnimeLink', name: funcName});
        cachedTestItemsByName.set(entry.deployment.name,
          addNode(entry.deployment, testGroup, existing.length > 0 ? existing[0] : undefined));
      }
    }

    // Phase 2: Refresh from live API with progress
    const pi = DG.TaskBarProgressIndicator.create('Updating KNIME deployments...');
    try {
      const deployments = await client.listDeployments('rest');
      const liveNames = new Set(deployments.map((d) => d.name));

      removeStaleItems(cachedItemsByName, liveNames);
      removeStaleItems(cachedTestItemsByName, liveNames);

      const regularDeps: KnimeDeployment[] = [];
      const testDeps: KnimeDeployment[] = [];
      for (const dep of deployments) {
        if (dep.name.toLowerCase().startsWith(TEST_PREFIX))
          testDeps.push(dep);
        else
          regularDeps.push(dep);
      }

      await updateOrAddNodes(regularDeps, treeNode, cachedItemsByName);

      if (testDeps.length > 0) {
        if (!testGroup)
          testGroup = treeNode.group('Test workflows', undefined, false);
        await updateOrAddNodes(testDeps, testGroup, cachedTestItemsByName);
      }

      if (treeNode.items.length === 0 && errors.length === 0)
        treeNode.item('No deployments found');
    }
    catch (e: any) {
      errors.push(e?.message ?? e);
    }
    finally {
      pi.close();
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