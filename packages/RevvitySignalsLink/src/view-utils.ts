import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { awaitCheck } from '@datagrok-libraries/utils/src/test';
import { getRevvityLibraries } from './libraries';
import { initializeFilters } from './search-utils';
import { getViewNameByCompoundType } from './utils';
import { retrieveQueriesMap } from './compounds';


const REVVITY_SIGNALS_APP_PATH: string = 'apps/Revvitysignalslink';
export let initSearchQuery = undefined;

export function resetInitSearchQuery() {
  initSearchQuery = undefined;
} 

export function createPath(viewName?: string[]) {
  let path = `${REVVITY_SIGNALS_APP_PATH}`;
  if (viewName)
    path += `/${viewName.map((it) => encodeURIComponent(it)).join('/')}`;
  return path;
}

export function setBreadcrumbsInViewName(viewPath: string[], tree: DG.TreeViewGroup, view?: DG.ViewBase): void {
  const usedView = view ?? grok.shell.v;
  const path = ['Home', 'Revvity Signals', ...viewPath.filter((v) => v !== 'Home' && v !== 'Demo')];
  const breadcrumbs = ui.breadcrumbs(path);

  breadcrumbs.onPathClick.subscribe(async (value) => {
    const actualItem = value[value.length - 1];
    if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
      return;
    tree.currentItem = actualItem === 'Revvity Signals' ? tree : tree.items.find((item) => item.text === actualItem)!;
  });

  if (usedView) {
    if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') { // integrate it to the actual breadcrumbs element
      const homeIcon = ui.iconFA('home', () => {
        grok.shell.v.close();
        grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME);
      });
      breadcrumbs.root.firstElementChild!.replaceWith(homeIcon);
    }
    const viewNameRoot = usedView.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
    if (viewNameRoot) {
      viewNameRoot.textContent = '';
      viewNameRoot.appendChild(breadcrumbs.root);
    }
  }
}

export async function openRevvityNode(treeNode: DG.TreeViewGroup, libToOpen: string, typeToOpen?: string) {
  //need to wait for tree to become available
  await awaitCheck(() => treeNode.items.find((node) => node.text.toLowerCase() === libToOpen.toLowerCase()) !== undefined,
    `Libraries haven't been loaded in 10 seconds`, 10000);
  const libNode = treeNode.items.find((node) => node.text.toLowerCase() === libToOpen.toLowerCase()) as DG.TreeViewGroup;
  libNode.expanded = true;
  if (typeToOpen) {
    const compoundType = getViewNameByCompoundType(typeToOpen);
    //need to wait for types to become available after lib node was expanded
    await awaitCheck(() => libNode.items.find((node) => node.text.toLowerCase() === compoundType.toLowerCase()) !== undefined,
      `Library types haven't been loaded in 10 seconds`, 10000);
    treeNode.items.find((node) => node.text.toLowerCase() === compoundType.toLowerCase())?.root.click();
  }
}

export async function handleInitialURL(treeNode: DG.TreeViewGroup, url: string) {

  if (url.startsWith('/'))
    url = url.slice(1);
  const componentsArr = url.split('/');
  if (componentsArr.length) {
    const libs = await getRevvityLibraries();
    //check for library
    const libIdx = libs.findIndex((it) => it.name.toLocaleLowerCase() === componentsArr[0].toLowerCase())
    if (libIdx === -1) {
      grok.shell.error(`Library ${componentsArr[0]} doesn't exist`);
      return;
    }
    //in case path contains only library
    if (componentsArr.length < 2) {
      openRevvityNode(treeNode, componentsArr[0]);
      return;
    }
    //check for type
    const typeIdx = libs[libIdx].types.findIndex((it) => it.name.toLowerCase() === componentsArr[1].toLowerCase());
    if (typeIdx === -1) {
      grok.shell.error(`Type ${componentsArr[1]} doesn't exist`);
      openRevvityNode(treeNode, componentsArr[0]);
      return;
    }

    //create initial serach condition from url
    if (componentsArr.length === 4 && componentsArr[2] === 'search') {
        const condition = JSON.parse(componentsArr[3]);
        initSearchQuery = condition;
    }

    openRevvityNode(treeNode, componentsArr[0], componentsArr[1]);
  }
}


export async function createViewFromPreDefinedQuery(treeNode: DG.TreeViewGroup, query: string, name: string, libName: string, compoundType: string) {
  //!!!!!!!!!!TODO: Change to .then()
  const df = DG.DataFrame.create();
  const tv = grok.shell.addTablePreview(df);
  tv.name = name;
  tv.path = createPath([libName, compoundType]);
  if (!initSearchQuery) {
    ui.setUpdateIndicator(tv.root, true, `Loading ${name}...`);
    grok.functions.call('RevvitySignalsLink:searchEntitiesWithStructures', {
      query: query,
      params: '{}'
    }).then((res: DG.DataFrame) => {
      ui.setUpdateIndicator(tv.root, false);
      tv.dataFrame = res;
      setBreadcrumbsInViewName([libName, compoundType], treeNode);
      const filtersDiv = ui.div([]);
      initializeFilters(tv, filtersDiv, libName, compoundType);
    })
      .catch((e: any) => {
        grok.shell.error(e?.message ?? e);
        ui.setUpdateIndicator(tv.root, false);
      });
  }
}