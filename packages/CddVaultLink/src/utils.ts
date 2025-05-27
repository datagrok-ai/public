
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {queryExportStatus, queryExportResult, ExportStatus, ApiResponse, Batch, Project, Vault, Molecule} from "./cdd-vault-api";
import { ALL_TABS, EXPANDABLE_TABS } from './constants';
import { awaitCheck } from '@datagrok-libraries/utils/src/test';

export const CDD_HOST = 'https://app.collaborativedrug.com/';
export const PREVIEW_ROW_NUM = 100;
const CDD_VAULT_APP_PATH: string = 'apps/Cddvaultlink';

export type CDDVaultStats = {
  name: string,
  projects?: string,
  molecules?: number,
  protocols?: number,
  batches?: number,
  collections?: number,
}

export async function getAsyncResultsAsDf(vaultId: number, exportId: number, timeoutMinutes: number, sdf: boolean): Promise<DG.DataFrame> {
  try {
    const resultResponse = await getAsyncResults(vaultId, exportId, timeoutMinutes, sdf);
    let df = DG.DataFrame.create();
    if (sdf) {
      const dfs = await grok.functions.call('Chem:importSdf', { bytes: resultResponse.data });
      if (dfs.length)
        df = dfs[0];
    } else {
      if (resultResponse.data?.objects) {
        prepareDataForDf(resultResponse.data.objects as any[]);
        df = DG.DataFrame.fromObjects(resultResponse.data.objects)!;
      }
    }
    return df;
  } catch (e: any) {
    throw e;
  }
}

const EXCLUDE_FIELDS = ['udfs', 'source_files']; 

export function prepareDataForDf(objects: any[]) {
  for (let i = 0; i < objects.length; i++) {
    EXCLUDE_FIELDS.forEach((key) => delete objects[i][key]);
    if (objects[i]['batches']) {
      objects[i]['molecule_batch_identifiers'] = (objects[i]['batches'] as Batch[]).map((it) => it.molecule_batch_identifier).filter((it) => it !== undefined);
      delete objects[i]['batches'];
    }
    if (objects[i]['projects'])
      objects[i]['projects'] = (objects[i]['projects'] as Project[]).map((it) => it.name);
    if (objects[i]['molecule_fields']) {
      Object.keys(objects[i]['molecule_fields']).forEach((key: string) => objects[i][key] = objects[i]['molecule_fields'][key]);
      delete objects[i]['molecule_fields'];
    }
  }
}

export async function getAsyncResults(vaultId: number, exportId: number, timeoutMinutes: number, text: boolean): Promise<ApiResponse<any>> {

  const timeoutMs = timeoutMinutes * 60 * 1000;
  const startTime = Date.now();

  while (Date.now() - startTime < timeoutMs) {
    const statusResponse = await queryExportStatus(vaultId, exportId);

    if (statusResponse.error)
      throw statusResponse.error;

    const status = statusResponse.data?.status;
    if (status === "finished") {
      const resultResponse = await queryExportResult(vaultId, exportId, text);
      if (resultResponse.error)
        throw resultResponse.error
      return resultResponse;
    }

    // Wait for 2 seconds before next check
    await new Promise(resolve => setTimeout(resolve, 2000));
  }

  throw `Export timed out after ${timeoutMinutes} minutes`;
}

export function getExportId(exportResponse: ApiResponse<ExportStatus>): number {
  if (exportResponse.error)
    throw exportResponse.error;

  const exportId = exportResponse.data?.id;
  if (!exportId)
    throw 'No export ID received';

  return exportId;
}

export async function createMoleculesDfFromObjects(vaultId: number, objects?: any[]) {
  if (!objects)
    return DG.DataFrame.create();
  prepareDataForDf(objects as any[]);
  const df = DG.DataFrame.fromObjects(objects)!;
  if (!df)
    return DG.DataFrame.create();
  createLinksFromIds(vaultId, df);
  reorderColummns(df);
  await grok.data.detectSemanticTypes(df);
  return df;
}

export async function createLinksFromIds(vaultId: number, df: DG.DataFrame) {
    const idCol = df.col('id');
    if (idCol) {
        const linkIdsCol = DG.Column.string('id', df.rowCount).init((i) => {
            const id = idCol.get(i);
            return `[${id}](${`${CDD_HOST}vaults/${vaultId}/molecules/${id}/`})`;
        });
        df.columns.replace(idCol, linkIdsCol);
    }
}

export async function reorderColummns(df: DG.DataFrame) {
  const colNames = df.columns.names();
  const firstColumns = ['id', 'name', 'smiles'];
  const newColOrder = [];
  for (const colName of firstColumns) {
    const index = colNames.indexOf(colName);
    if (index > -1) {
      colNames.splice(index, 1);
      newColOrder.push(colName);
    }
  }
  df.columns.setOrder(newColOrder.concat(colNames));
}

export function paramsStringFromObj(params: any): string {
    let paramsStr = '';
    const paramNames = Object.keys(params);
    for (let i = 0; i < paramNames.length; i++) {
        const paramVal = params[paramNames[i]];
        if (paramVal) {
            paramsStr += paramsStr === '' ? `?${paramNames[i]}=${paramVal}` : `&${paramNames[i]}=${paramVal}`;
        }
    }
    return paramsStr;
}

export function createObjectViewer(obj: any, title: string = 'Object Viewer', additionalHeaderEl?: HTMLElement): HTMLElement {
  // Helper function to determine if a value is a dictionary/object
  function isDictionary(value: any): boolean {
    return value !== null && 
           typeof value === 'object' && 
           !Array.isArray(value) && 
           !(value instanceof Date);
  }

  // Helper function to determine if a value is a simple array
  function isSimpleArray(value: any): boolean {
    return Array.isArray(value) && 
           value.every(item => 
             typeof item === 'string' || 
             typeof item === 'number' || 
             typeof item === 'boolean'
           );
  }

  // Helper function to get pane name for array items
  function getPaneName(item: any, index: number, parentName: string): string {
    if (item.name) return item.name;
    if (item.id) return item.id.toString();
    if (parentName === 'protocol_statistics' && item.readout_definition?.name)
      return item.readout_definition.name;
    return `Item ${index}`;
  }

  // Helper function to create a view for a specific value
  function createValueView(value: any, parentName: string): HTMLElement {
    if (value === null || value === undefined) {
      return ui.divText('null');
    }

    if (isSimpleArray(value)) {
      return ui.divText(value.join(', '));
    }

    if (Array.isArray(value)) {
      const accordion = ui.accordion();
      value.forEach((item, index) => {
        if (isDictionary(item) || Array.isArray(item)) {
          accordion.addPane(getPaneName(item, index, parentName), () => createValueView(item, parentName));
        } else {
          accordion.addPane(`Item ${index}`, () => ui.divText(String(item)));
        }
      });
      return accordion.root;
    }

    if (isDictionary(value)) {
      const simpleProperties: { [key: string]: any } = {};
      const complexProperties: { [key: string]: any } = {};

      // Separate simple and complex properties
      for (const [key, val] of Object.entries(value)) {
        if (isDictionary(val) || Array.isArray(val)) {
          complexProperties[key] = val;
        } else {
          simpleProperties[key] = val;
        }
      }

      // Create container for both simple properties table and complex properties accordion
      const container = ui.divV([]);

      // Add simple properties table if any exist
      if (Object.keys(simpleProperties).length > 0) {
        container.appendChild(ui.tableFromMap(simpleProperties));
      }

      // Add complex properties as nested accordions if any exist
      if (Object.keys(complexProperties).length > 0) {
        const accordion = ui.accordion();
        for (const [key, val] of Object.entries(complexProperties)) {
          accordion.addPane(key, () => createValueView(val, parentName));
        }
        container.appendChild(accordion.root);
      }

      return container;
    }

    if (value instanceof Date) {
      return ui.divText(value.toISOString());
    }

    return ui.divText(String(value));
  }

  // Create dictionaries for simple and complex properties
  const simpleProperties: { [key: string]: any } = {};
  const complexProperties: { [key: string]: any } = {};

  // Separate simple and complex properties
  for (const [key, value] of Object.entries(obj)) {
    if (isDictionary(value) || Array.isArray(value)) {
      complexProperties[key] = value;
    } else {
      simpleProperties[key] = value;
    }
  }

  const header = ui.divH([ui.h2(title)]);
  if (additionalHeaderEl)
    header.append(additionalHeaderEl);
  // Create the main container
  const container = ui.divV([header]);

  // Add simple properties table if any exist
  if (Object.keys(simpleProperties).length > 0) {
    container.appendChild(ui.tableFromMap(simpleProperties));
  }

  // Add complex properties as accordion if any exist
  if (Object.keys(complexProperties).length > 0) {
    const accordion = ui.accordion();
    for (const [key, value] of Object.entries(complexProperties)) {
      accordion.addPane(key, () => createValueView(value, key));
    }
    container.appendChild(accordion.root);
  }

  return container;
}

export async function handleInitialURL(treeNode: DG.TreeViewGroup, url: URL) {

  const currentTabs = url.pathname.includes(`${CDD_VAULT_APP_PATH}`) ?
    url.pathname.replace(`${CDD_VAULT_APP_PATH}`, ``).replace(/^\/+/, '').split('/').map((it) => decodeURIComponent(it)) : [];
  const currentVault = currentTabs.length > 0 ? currentTabs[0] : undefined;
  const currentView = currentTabs.length > 1 ? currentTabs[1] : undefined;
  const currentSubView = currentTabs.length > 2 ? currentTabs[2] : undefined;
  openCddNode(treeNode, currentVault, currentView, currentSubView);
}

export async function openCddNode(treeNode: DG.TreeViewGroup, currentVault?: string, currentView?: string, currentSubView?: string) {
  if (currentVault) {
    //need to wait for tree to become available
    await awaitCheck(() => treeNode.items.find((node) => node.text === currentVault) !== undefined, `CDD tabs haven't been loaded in 10 seconds`, 10000);
    const vault = treeNode.items.find((node) => node.text === currentVault) as DG.TreeViewGroup;
    //look for current vault from URL and open it
    if (vault) {
      currentView ? vault.expanded = true : vault.root.click();
      //look for current vault from URL and open it
      if (currentView) {
        if (!ALL_TABS.includes(currentView)) {
          vault.root.click();
          grok.shell.warning(`${currentView} section doesn't exist`);
          return;
        }
        try {
          await awaitCheck(() => treeNode.items.find((node) => node.text === currentView) !== undefined, `CDD tabs haven't been loaded in 5 seconds`, 5000);
          if (EXPANDABLE_TABS.includes(currentView)) {
            const tab = treeNode.items.find((node) => node.text === currentView) as DG.TreeViewGroup;
            currentSubView ? tab.expanded = true : tab.root.click();
            //look for inner category from URL and open it
            if (currentSubView) {
              await awaitCheck(() => tab.items.length > 0, `CDD tabs haven't been loaded in 10 seconds`, 10000);
              const innerView = treeNode.items.find((node) => node.text === currentSubView);
              if (!innerView) {
                tab.root.click();
                grok.shell.warning(`${currentSubView} id doesn't exist in ${currentView}`);
                return;
              }
              innerView!.root.click();
            }

          } else
            treeNode.items.find((node) => node.text === currentView)!.root.click();
        } catch(e) {
          console.log(e);
        }
      }
    }
    else
      grok.shell.warning(`Vault ${currentVault} doesn't exist`);
  }
}

export function createNestedCDDNode(openedViews: {[key: string]: DG.View}, items: any[] | null, nodeName: string, vaultNode: DG.TreeViewGroup,
  getItemsFunsName: string, getItemsFuncParams: any, treeNode: DG.TreeViewGroup, vault: Vault,
  onItemSelected: (item: any) => Promise<void>) {
  const nestedNode = vaultNode.group(nodeName, null, false);
  const fullNodeName = `${vault.id}|${nodeName}`;
  const loadData = async () => {
    if (!items) {
      try {
        const itemsStr = await grok.functions.call(getItemsFunsName, getItemsFuncParams);
        items = itemsStr !== '' ? JSON.parse(itemsStr) as any[] : [];
      } catch (e: any) {
        grok.shell.error(e?.message ?? e);
        ui.setUpdateIndicator(openedViews[fullNodeName].root, false);
        ui.empty(openedViews[fullNodeName].root);
        openedViews[fullNodeName].root.append(ui.divText(`Error`));
      }
    }
  }
  nestedNode.onSelected.subscribe(async () => {
    openedViews[fullNodeName]?.close();
    openedViews[fullNodeName] = DG.View.create();
    openedViews[fullNodeName].path = createPath(vault.name, [nodeName]);
    grok.shell.addPreview(openedViews[fullNodeName]);
    ui.setUpdateIndicator(openedViews[fullNodeName].root, true, `Loading ${nodeName}...`);
    openedViews[fullNodeName].name = nodeName;
    setBreadcrumbsInViewName([vault.name, nodeName], treeNode);
    await loadData();
    const links = items!.length ? createLinks(nodeName, items!.map((it) => it.name), treeNode, openedViews[fullNodeName]) :
      ui.divText(`No ${nodeName} found in the vault`);
    ui.empty(openedViews[fullNodeName].root);
    openedViews[fullNodeName].append(links);
    ui.setUpdateIndicator(openedViews[fullNodeName].root, false);
    if (!nestedNode.expanded)
      nestedNode.expanded = true;
  });

  nestedNode.onNodeExpanding.subscribe(async () => {
    await loadData();
    for (const item of items!) {
      const itemNode = nestedNode.item(item.name);
      itemNode.onSelected.subscribe(async () => {
        await onItemSelected(item);
        setBreadcrumbsInViewName([nodeName, item.name], treeNode);
      });
    }
  });
}

export function createPath(vaultName: string, viewName?: string[]) {
  let path = `${CDD_VAULT_APP_PATH}/`;
  path += encodeURIComponent(vaultName);
  if (viewName)
    path += `/${viewName.map((it) => encodeURIComponent(it)).join('/')}`;
  return path;
}

function updateView(views: {[key: string]: DG.ViewBase}, viewName: string[], viewFullName: string,
  vaultName: string, treeNode: DG.TreeViewGroup, progressMessage: string, res?: DG.DataFrame) {
  views[viewFullName]?.close();
  views[viewFullName] = res ? grok.shell.addTablePreview(res) : grok.shell.addPreview(DG.View.create());
  views[viewFullName].name = viewName[viewName.length - 1];
  views[viewFullName].path = createPath(vaultName, viewName);
  setBreadcrumbsInViewName([vaultName].concat(viewName), treeNode, views[viewFullName]);
  if (!res)
    ui.setUpdateIndicator(views[viewFullName].root, true, progressMessage);
}

export async function createCDDTableView(vaultId: number, views: {[key: string]: DG.View}, viewName: string[], progressMessage: string, funcName: string,
  funcParams: {[key: string]: any}, vaultName: string, treeNode: DG.TreeViewGroup) {
  const viewFullName = `${vaultId}|${viewName.join('|')}`;
  try {
    updateView(views, viewName, viewFullName, vaultName, treeNode, progressMessage);
    const df: DG.DataFrame = await grok.functions.call(funcName, funcParams);
    updateView(views, viewName, viewFullName, vaultName, treeNode, progressMessage, df);
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
    updateView(views, viewName, viewFullName, vaultName, treeNode, progressMessage, DG.DataFrame.create());
  }
}

export async function createCDDTableViewWithPreview(vaultId: number, views: {[key: string]: DG.ViewBase}, viewName: string[], progressMessage: string, syncfuncName: string,
  syncfuncParams: {[key: string]: any}, asyncfuncName: string, asyncfuncParams: {[key: string]: any},
  vaultName: string, treeNode: DG.TreeViewGroup) {

  const handleError = (e: any) => {
    errors++;
    grok.shell.error(e?.message ?? e);
    if (errors === 2) {
      progressBar.close();
      //if both sync and async requests are failed - showing empty dataframe
      updateView(views, viewName, viewFullName, vaultName, treeNode, progressMessage, DG.DataFrame.create());
    }
  }

  const viewFullName = `${vaultId}|${viewName.join('|')}`;
  let asyncRequestCompleted = false;
  let errors = 0;
  updateView(views, viewName, viewFullName, vaultName, treeNode, progressMessage);

  //run sync function with offset and create a preview
  grok.functions.call(syncfuncName, syncfuncParams).then(async (df: DG.DataFrame) => {
    //in case asyn request is completed before sync - returning from function
    if (asyncRequestCompleted)
      return;
    if (df.rowCount < PREVIEW_ROW_NUM) {
      asyncRequestCompleted = true; //we will not need async request results
      progressBar.close();
    }
    updateView(views, viewName, viewFullName, vaultName, treeNode, progressMessage, df);
  }).catch((e) => {
    handleError(e);
    if (!asyncRequestCompleted)
      grok.shell.warning(`Request for first ${PREVIEW_ROW_NUM} rows failed. Waiting for full results`);
  });

  //reset tableView with asynchronously received results
  grok.functions.call(asyncfuncName, asyncfuncParams).then(async (df: DG.DataFrame) => {
    if (asyncRequestCompleted)
      return;
    updateView(views, viewName, viewFullName, vaultName, treeNode, progressMessage, df);
    progressBar.close();
  }).catch((e) => {
    handleError(e);
    if (asyncRequestCompleted)
      grok.shell.warning(`Loaded results only for first ${PREVIEW_ROW_NUM} rows. LOading full results failed`);
  }).finally(() => {
    asyncRequestCompleted = true;
  });
  const progressBar = DG.TaskBarProgressIndicator.create(`Loading ${viewName[viewName.length - 1]}...`);
}

export function createLinks(header: string, nodeNames: string[], tree: DG.TreeViewGroup, view: DG.ViewBase): HTMLDivElement {
  const table = ui.table(nodeNames, (item) =>
  ([ui.link(item, () => {
    view.close();
    tree.currentItem = tree.items.find((it) => it.text === item)!
  }, 'Click to open')]), [header]);
  return table;
}

export function setBreadcrumbsInViewName(viewPath: string[], tree: DG.TreeViewGroup, view?: DG.ViewBase): void {
  const usedView = view ?? grok.shell.v;
  const path = ['Home', 'CDD Vault', ...viewPath.filter((v) => v !== 'Home' && v !== 'Demo')];
  const breadcrumbs = ui.breadcrumbs(path);

  breadcrumbs.onPathClick.subscribe(async (value) => {
    const actualItem = value[value.length - 1];
    if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
      return;
    tree.currentItem = actualItem === 'CDD Vault' ? tree : tree.items.find((item) => item.text === actualItem)!;
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

export function createCDDContextPanel(obj: Molecule | Batch, vaultId?: number): HTMLElement {
  const keys = Object.keys(obj);
  const resDiv = ui.divV([], 'cdd-context-panel');
  const accordions = ui.divV([]);
  const dictForTableView: {[key: string]: any} = {};
  for(const key of keys) {
    if (key === 'batches') {
      const batchesAcc = ui.accordion(key);
      const batches = (obj as Molecule)[key] as Batch[];
      batchesAcc.addPane(key, () => {
        const innerbatchesAcc = ui.accordion();
        for (const batch of batches) {
          innerbatchesAcc.addPane(batch.id.toString(), () => createCDDContextPanel(batch));
        }
        return innerbatchesAcc.root;
      });
      accordions.append(batchesAcc.root);
    } else if (key === 'collections' || key === 'projects' || key === 'source_files' || key === 'batch_fields') {
      const acc = ui.accordion(key);
      acc.addPane(key, () => {
        const div = ui.divV([]);
        ((obj as any)[key] as any[]).forEach((it) => div.append(ui.tableFromMap(it, true)));
        return div;
      });
      accordions.append(acc.root);
    } else if (key === 'molecule_fields' || key === 'udfs' || key === 'stoichiometry') {
      const acc = ui.accordion(key);
      acc.addPane(key, () => ui.tableFromMap((obj as any)[key] as any, true));
      accordions.append(acc.root);
    } else {
      const initialValue = (obj as any)[key];
      let value: any | null = null;
      if (key === 'id' && vaultId)
        value = ui.link(initialValue, () => window.open(`${CDD_HOST}vaults/${vaultId}/molecules/${initialValue}/`));
      else
        value = initialValue instanceof Array ? initialValue.join(', ') : initialValue;
      dictForTableView[key] = value;
    }
  }
  resDiv.append(ui.tableFromMap(dictForTableView, true));
  resDiv.append(accordions);

  return resDiv;
}


export function createInitialSatistics(statsDiv: HTMLDivElement) {
  
  grok.functions.call('CDDVaultLink:getVaults').then(async (res: string) => {
    const stats: CDDVaultStats[] = [];
    if (!res)
      return;
    const vaults = JSON.parse(res) as Vault[];
    for (const vault of vaults) {
      try {
        const resStr = await grok.functions.call('CDDVaultLink:getVaultStats', {vaultId: vault.id, vaultName: vault.name});
        stats.push(JSON.parse(resStr));

      } catch (e: any) {
        grok.shell.error(`Cannot get statistics for vault ${vault.name}: ${e?.message ?? e}`);
        continue;
      }
    }

    const table = ui.table(stats, (info) =>
      ([ui.link(info.name, () => {
        const cddNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('CDD Vault');
        cddNode.expanded = true;
        openCddNode(cddNode, info.name);
      }),
        info.projects ?? '',
        info.molecules ?? '',
        info.protocols ?? '',
        info.batches ?? '',
        info.collections ?? '',
      ]),
    ['Vault', 'Projects', 'Molecules', 'Protocols', 'Batches', 'Collections']);
    statsDiv.append(table);
    ui.setUpdateIndicator(statsDiv, false);

  }).catch((e: any) => {
    grok.shell.error(e?.message ?? e);
  })
}
