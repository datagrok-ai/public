import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { awaitCheck, delay } from '@datagrok-libraries/test/src/test';
import { getRevvityLibraries } from './libraries';
import { filterProperties, getPropertiesForLibAndEntityType, initializeFilters, runSearchQuery, SAVED_SEARCH_STORAGE } from './search-utils';
import { getCompoundTypeByViewName, getViewNameByCompoundType } from './utils';
import { getConditionForLibAndType, retrieveQueriesMap } from './compounds';
import { ComplexCondition, Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { RevvityUser } from './revvity-api';
import { LAYOUT_STORAGE, USER_FIELDS } from './constants';
import { funcs } from './package-api';
import { _package } from './package';
import { applyRevvityLayout, saveRevvityLayout } from './layout';


const REVVITY_SIGNALS_APP_PATH: string = 'apps/Revvitysignalslink';

export let openedView: DG.TableView | DG.ViewBase | null = null;

export async function updateView(tv: DG.TableView, df: DG.DataFrame, compoundType: string, libName: string, libId?: string, filtersDiv?: HTMLDivElement) {
  tv.dataFrame = df;

  setColumnsFormat(tv);

  //look for friendly names and set
  if (!libId) {
    const libs = await getRevvityLibraries();
    const lib = libs.filter((it) => it.name.toLowerCase() === libName.toLowerCase());
    if (lib.length)
      libId = lib[0].id;
  }
  if (libId)
    await setColumnsFriendlyNames(libId, compoundType, tv.dataFrame);

  applyRevvityLayout(tv.grid, `${libName}|${compoundType}`);
  const columns = tv.dataFrame.columns.names();
  if (columns.length)
    tv.dataFrame.currentCell = tv.dataFrame.cell(-1, tv.dataFrame.col(columns[0])!.name);

  DG.debounce(tv.grid.onPropertyValueChanged, 5000).subscribe(() => {
    saveRevvityLayout(tv.grid, libName, compoundType)
  });
  DG.debounce(tv.dataFrame.onMetadataChanged, 5000).subscribe(() => {
    saveRevvityLayout(tv.grid, libName, compoundType)
  });
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

export async function openRevvityNode(treeNode: DG.TreeViewGroup, nodesToExpand: string[], nodeToSelectName: string,
  libToOpen?: string, typeToOpen?: string, initialSearchQuery?: ComplexCondition, isSavedSearch?: boolean) {
  let lastExpandedNode = null;
  for (let nodeName of nodesToExpand) {
    try {
      await awaitCheck(() => treeNode.items.find((node) => node.text.toLowerCase() === nodeName.toLowerCase()) !== undefined,
        `${nodeName} haven't been loaded in 10 seconds`, 10000);
    } catch (e: any) {
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

  if (url.startsWith('/'))
    url = url.slice(1);
  const componentsArr = url.split('/');
  if (componentsArr.length) {
    let libs;
    try {
      libs = await getRevvityLibraries();
    } catch (e: any) {
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
  getElement: (root: HTMLElement, libName?: string, typeName?: string) => Promise<HTMLDivElement>, libName?: string, typeName?: string) {
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
  openedView?.close();

  const df = DG.DataFrame.create();
  const tv = grok.shell.addTablePreview(df);
  const icon = ui.iconImage('', `${_package.webRoot}img/revvity.png`);
  //tv.setIcon(icon); TODO!!! Uncomment when setIcon method is available in js-api

  openedView = tv;
  openedView.name = name.charAt(0).toUpperCase() + name.slice(1);
  openedView.path = createPath(path);
  const filtersDiv = ui.divV([], 'revvity-signals-filter-panel');

  if (isSavedSearch && !initialSearchQuery) {
    const savedSearchesStr = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, `${libName}|${compoundType}`) || '{}';
    const savedSearches: { [key: string]: string } = JSON.parse(savedSearchesStr);
    initialSearchQuery = JSON.parse(savedSearches[name]);
  }

  ui.setUpdateIndicator(openedView.root, true, `Loading ${name}...`);

  const initFilters = () => {
    ui.setUpdateIndicator(openedView!.root, false);
    setBreadcrumbsInViewName([libName, compoundType], treeNode);
    initializeFilters(treeNode, tv, filtersDiv, libName, compoundType, initialSearchQuery, !isSavedSearch);
  }

  if (!initialSearchQuery) {
    getRevvityLibraries().then((libs) => {
      if (!libs.length) {
        grok.shell.warning(`No materials libraries found`);
        return;
      }
      const libId = libs.filter((it) => it.name.toLowerCase() === libName.toLowerCase());
      if (!libId.length) {
        grok.shell.warning(`Library ${libName} not found`);
        return;
      }
      initialSearchQuery = {
        logicalOperator: Operators.Logical.and,
        conditions: []
      }
      initFilters();
    })
  } else
    initFilters();
}

export function setColumnsFormat(tv: DG.TableView) {
  const colNames = tv.dataFrame.columns.names();
  for (const colName of colNames) {
    const dfColType = tv.dataFrame.col(colName)!.type;
    const format = dfColType === DG.TYPE.DATE_TIME ? 'M/d/yyyy' : dfColType === DG.TYPE.FLOAT ? 'two digits after comma' : '';
    if (format) {
      const gridCol = tv.grid.col(colName)!
      gridCol.format = format;
    }
  }
}

export async function setColumnsFriendlyNames(libId: string, entityType: string, dataFrame: DG.DataFrame) {
  const colNames = dataFrame.columns.names();
  const props = await getPropertiesForLibAndEntityType(libId, entityType);
  for (const colName of colNames) {
    const prop = props.filter((it) => it.name === colName);
    if (prop.length && prop[0].friendlyName)
      dataFrame.col(colName)!.setTag('friendlyName', prop[0].friendlyName);
  }
}