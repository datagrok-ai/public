import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { awaitCheck } from '@datagrok-libraries/utils/src/test';
import { getRevvityLibraries } from './libraries';
import { initializeFilters, SAVED_SEARCH_STORAGE } from './search-utils';
import { getCompoundTypeByViewName, getViewNameByCompoundType } from './utils';
import { retrieveQueriesMap } from './compounds';
import { ComplexCondition } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { RevvityUser } from './revvity-api';
import { USER_FIELDS } from './constants';



const REVVITY_SIGNALS_APP_PATH: string = 'apps/Revvitysignalslink';

let openedView: DG.TableView | DG.ViewBase | null = null;

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

export async function openRevvityNode(treeNode: DG.TreeViewGroup, nodesToExpand: string[], nodeToSelectName: string,
  libToOpen?: string, typeToOpen?: string, initialSearchQuery?: ComplexCondition, isSavedSearch?: boolean) {
  let lastExpandedNode = null;
  for (let nodeName of nodesToExpand) {
    try {
      await awaitCheck(() => treeNode.items.find((node) => node.text.toLowerCase() === nodeName.toLowerCase()) !== undefined,
        `${nodeName} haven't been loaded in 10 seconds`, 10000);
    } catch(e: any) {
      grok.shell.error(e?.message ?? e);
      return;
    }
    lastExpandedNode = treeNode.items.find((node) => node.text.toLowerCase() === nodeName.toLowerCase()) as DG.TreeViewGroup;
    lastExpandedNode.expanded = true;
  }
  const path = nodesToExpand.concat(nodeToSelectName)
    .concat(initialSearchQuery ? ['search', JSON.stringify(initialSearchQuery)] : []);
  if (libToOpen && typeToOpen) {
    createViewFromPreDefinedQuery(treeNode, path, libToOpen, typeToOpen, initialSearchQuery, isSavedSearch);

    // //selecting node manually
    // if (lastExpandedNode) {
    //   await awaitCheck(() => lastExpandedNode.items.find((node) => node.text.toLowerCase() === nodeToSelectName.toLowerCase()) !== undefined,
    //     `${nodeToSelectName} types haven't been loaded in 10 seconds`, 10000);
    //   const nodeToSelect = treeNode.items.find((node) => node.text.toLowerCase() === nodeToSelectName.toLowerCase());
    //   nodeToSelect?.root.classList.add('d4-tree-view-node-selected');
    // }
  }
}

export async function handleInitialURL(treeNode: DG.TreeViewGroup, url: string) {

  const tv = grok.shell.addTablePreview(DG.DataFrame.create());
  ui.setUpdateIndicator(tv.root, true, 'Loading Revvity Signals...');
  if (url.startsWith('/'))
    url = url.slice(1);
  const componentsArr = url.split('/');
  if (componentsArr.length) {
    let libs;
    try {
      libs = await getRevvityLibraries();
    } catch(e: any) {
      ui.setUpdateIndicator(tv.root, false);
      grok.shell.error(`Revvity libraries haven't been loaded: ${e?.message ?? e}`);
      return;
    }

    let idx = 0;
    const nodesToExpand = [];
    if (componentsArr[idx] === 'saved searches') {
      nodesToExpand.push(componentsArr[idx]);
      idx++;
    }
    if (componentsArr.length > idx) {
      //collect library name
      const libName = componentsArr[idx];
      //check if library exists
      const libIdx = libs.findIndex((it) => it.name.toLocaleLowerCase() === componentsArr[idx].toLowerCase())
      if (libIdx === -1) {
        grok.shell.error(`Library ${componentsArr[idx]} doesn't exist`);
        openRevvityNode(treeNode, nodesToExpand, nodesToExpand[nodesToExpand.length - 1]);
        return;
      }
      nodesToExpand.push(componentsArr[idx]);
      idx++;

      if (componentsArr.length > idx) {
        //collect type name
        const entityType = getCompoundTypeByViewName(componentsArr[idx]);
        //check if entity type exists
        const typeIdx = libs[libIdx].types.findIndex((it) => it.name.toLowerCase() === entityType.toLowerCase());
        if (typeIdx === -1) {
          grok.shell.error(`Type ${entityType} doesn't exist`);
          openRevvityNode(treeNode, nodesToExpand, componentsArr[idx], libName);
          return;
        }

        //in case of saved searches
        if (nodesToExpand[0] === 'saved searches' && componentsArr.length > idx + 1) {
          nodesToExpand.push(componentsArr[idx]);
          openRevvityNode(treeNode, nodesToExpand, componentsArr[idx + 1], libName,
            entityType, undefined, true);
          return;
        }

        //in case of search query, 'search' should be with index 2, components array should contain 4 items 
        if (componentsArr[idx + 1] === 'search' && componentsArr.length === idx + 3) {
          const condition = JSON.parse(componentsArr[idx + 2]);
          openRevvityNode(treeNode, nodesToExpand, componentsArr[idx], libName,
            entityType, condition, false);
            return;
        }
        openRevvityNode(treeNode, nodesToExpand, componentsArr[idx], libName, entityType);
      }
      openRevvityNode(treeNode, nodesToExpand, componentsArr[idx], libName);
    }
  }
}

export function createViewForExpandabelNode(viewName: string,
  getElement:(root: HTMLElement, libName?: string, typeName?: string) => Promise<void>, libName?: string, typeName?: string) {
  openedView?.close();
  const div = ui.div();
  openedView = grok.shell.addPreview(DG.View.fromRoot(div));
  const name = typeName ? getViewNameByCompoundType(typeName) : libName ?? viewName;
  openedView.name = name.charAt(0).toUpperCase() + name.slice(1);
  getElement(div, libName, typeName);
}

export async function createViewFromPreDefinedQuery(treeNode: DG.TreeViewGroup, path: string[], libName: string,
  compoundType: string, initialSearchQuery?: ComplexCondition, isSavedSearch?: boolean): Promise<void> {
  let name = path[path.length - 1];
  if (!isSavedSearch)
    name = getViewNameByCompoundType(compoundType);
  if (openedView && openedView.name.toLowerCase() === name.toLowerCase() && ("dockNode" in openedView && openedView.dockNode.parent)) {
    grok.shell.v = openedView;
    return;
  }
  openedView?.close();

  const df = DG.DataFrame.create();
  const tv = grok.shell.addTablePreview(df);
  openedView = tv;
  openedView.name = name.charAt(0).toUpperCase() + name.slice(1);
  openedView.path = createPath(path);

  if (isSavedSearch && !initialSearchQuery) {
    const savedSearchesStr = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, `${libName}|${compoundType}`) || '{}';
    const savedSearches: { [key: string]: string } = JSON.parse(savedSearchesStr);
    initialSearchQuery = JSON.parse(savedSearches[name]);
  }

  ui.setUpdateIndicator(openedView.root, true, `Loading ${name}...`);


  const initFilters = () => {
    ui.setUpdateIndicator(openedView!.root, false);
    setBreadcrumbsInViewName([libName, compoundType], treeNode);
    const filtersDiv = ui.divV([], 'revvity-signals-filter-panel');
    initializeFilters(treeNode, tv, filtersDiv, libName, compoundType, initialSearchQuery, !isSavedSearch);
  }

  if (!initialSearchQuery) {
    grok.functions.call('RevvitySignalsLink:searchEntitiesWithStructures', {
      query: JSON.stringify(retrieveQueriesMap[compoundType]),
      params: '{}'
    }).then((res: DG.DataFrame) => {
      (openedView! as DG.TableView).dataFrame = res;
      setUserColumnsStyle(openedView!  as DG.TableView);
      initFilters();
    })
      .catch((e: any) => {
        grok.shell.error(e?.message ?? e);
        ui.setUpdateIndicator(openedView!.root, false);
      });
  } else
    initFilters();
}

export function setUserColumnsStyle(tv: DG.TableView) {
  USER_FIELDS.forEach((it) => {
    const col = tv.grid.col(it);
    if (col) {
      col.cellType = 'html';
      col.width = 100;
    }
  });
}