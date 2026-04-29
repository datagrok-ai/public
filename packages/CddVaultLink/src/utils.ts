/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {queryExportStatus, queryExportResult, ExportStatus, ApiResponse,
  Batch, Project, Vault, Molecule, CDD_HOST} from './cdd-vault-api';
import {ALL_TABS, BATCHES_TAB, CDDVaultSearchType, COLLECTIONS_TAB, EXPANDABLE_TABS,
  MOLECULES_TAB, PROTOCOLS_TAB, SAVED_SEARCHES_TAB, SEARCH_TAB} from './constants';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';
import {SearchEditor} from './search-function-editor';
import {funcs} from './package-api';

export const PREVIEW_ROW_NUM = 100;
const CDD_VAULT_APP_PATH: string = 'apps/Cddvaultlink';
const PROTOCOL_STATISTICS_KEY = 'protocol_statistics';

let openedView: DG.ViewBase | null = null;

export type CDDVaultStats = {
  name: string,
  projects?: string,
  molecules?: number,
  protocols?: number,
  batches?: number,
  collections?: number,
}

export async function getAsyncResultsAsDf(vaultId: number, exportId: number,
  timeoutMinutes: number, sdf: boolean): Promise<DG.DataFrame> {
  const resultResponse = await getAsyncResults(vaultId, exportId, timeoutMinutes, sdf);
  if (sdf) {
    const dfs = await grok.functions.call('Chem:importSdf', {bytes: resultResponse.data});
    return dfs.length ? dfs[0] : DG.DataFrame.create();
  }
  if (resultResponse.data?.objects) {
    prepareDataForDf(resultResponse.data.objects as any[]);
    return DG.DataFrame.fromObjects(resultResponse.data.objects) ?? DG.DataFrame.create();
  }
  return DG.DataFrame.create();
}

const EXCLUDE_FIELDS = ['udfs', 'source_files'];

export function prepareDataForDf(objects: any[]) {
  for (let i = 0; i < objects.length; i++) {
    EXCLUDE_FIELDS.forEach((key) => delete objects[i][key]);
    if (objects[i]['batches']) {
      objects[i]['molecule_batch_identifiers'] = (objects[i]['batches'] as Batch[])
        .map((it) => it.molecule_batch_identifier).filter((it) => it !== undefined);
      delete objects[i]['batches'];
    }
    if (objects[i]['molecule'] && typeof objects[i]['molecule'] === 'object') {
      const mol = objects[i]['molecule'] as Molecule;
      objects[i]['molecule_id'] = mol.id;
      objects[i]['molecule_smiles'] = mol.smiles;
      delete objects[i]['molecule'];
    }
    if (objects[i]['projects'])
      objects[i]['projects'] = (objects[i]['projects'] as Project[]).map((it) => it.name);
    if (objects[i]['molecule_fields']) {
      Object.keys(objects[i]['molecule_fields'])
        .forEach((key: string) => objects[i][key] = objects[i]['molecule_fields'][key]);
      delete objects[i]['molecule_fields'];
    }
  }
}

export async function getAsyncResults(vaultId: number, exportId: number,
  timeoutMinutes: number, binary: boolean): Promise<ApiResponse<any>> {
  const timeoutMs = timeoutMinutes * 60 * 1000;
  const startTime = Date.now();

  while (Date.now() - startTime < timeoutMs) {
    const statusResponse = await queryExportStatus(vaultId, exportId);

    if (statusResponse.error)
      throw statusResponse.error;

    const status = statusResponse.data?.status;
    if (status === 'finished') {
      const resultResponse = await queryExportResult(vaultId, exportId, binary);
      if (resultResponse.error)
        throw resultResponse.error;
      return resultResponse;
    }

    // Wait for 2 seconds before next check
    await new Promise((resolve) => setTimeout(resolve, 2000));
  }

  throw new Error(`Export timed out after ${timeoutMinutes} minutes`);
}

/** Run a CDD async-export query end-to-end: start it, poll until done, return the raw ApiResponse.
 *  Throws on any error from `startQuery` or during polling (same as calling the parts manually). */
export async function runAsyncExport(vaultId: number,
  startQuery: () => Promise<ApiResponse<ExportStatus>>,
  timeoutMinutes: number, binary: boolean = false): Promise<ApiResponse<any>> {
  const exportResponse = await startQuery();
  const exportId = getExportId(exportResponse);
  return await getAsyncResults(vaultId, exportId, timeoutMinutes, binary);
}

/** Same as runAsyncExport, but converts the result into a DataFrame (SDF or JSON objects). */
export async function runAsyncExportAsDf(vaultId: number,
  startQuery: () => Promise<ApiResponse<ExportStatus>>,
  timeoutMinutes: number, sdf: boolean = false): Promise<DG.DataFrame> {
  const exportResponse = await startQuery();
  const exportId = getExportId(exportResponse);
  return await getAsyncResultsAsDf(vaultId, exportId, timeoutMinutes, sdf);
}

export function getExportId(exportResponse: ApiResponse<ExportStatus>): number {
  if (exportResponse.error)
    throw exportResponse.error;

  const exportId = exportResponse.data?.id;
  if (!exportId)
    throw new Error('No export ID received');

  return exportId;
}

/** Shared df-from-CDD-objects pipeline: flatten via prepareDataForDf, build a DF,
 *  apply an optional per-endpoint transform, reorder canonical columns, detect sem types. */
export async function createCddDfFromObjects(objects: any[] | undefined,
  postProcess?: (df: DG.DataFrame) => void): Promise<DG.DataFrame> {
  if (!objects)
    return DG.DataFrame.create();
  prepareDataForDf(objects);
  const df = DG.DataFrame.fromObjects(objects) ?? DG.DataFrame.create();
  postProcess?.(df);
  reorderColumns(df);
  await grok.data.detectSemanticTypes(df);
  return df;
}

export async function createMoleculesDfFromObjects(vaultId: number, objects?: any[]): Promise<DG.DataFrame> {
  return createCddDfFromObjects(objects, (df) => createLinksFromIds(vaultId, df));
}

export function createMoleculeIdLinks(vaultId: number, df: DG.DataFrame): void {
  const idCol = df.col('molecule_id');
  if (idCol) {
    const linkCol = DG.Column.string('molecule_id', df.rowCount).init((i) => {
      const id = idCol.get(i);
      return id != null ? `[${id}](${`${CDD_HOST}/vaults/${vaultId}/molecules/${id}/`})` : '';
    });
    df.columns.replace(idCol, linkCol);
  }
}

export async function createBatchesDfFromObjects(vaultId: number, objects?: any[]): Promise<DG.DataFrame> {
  return createCddDfFromObjects(objects, (df) => createMoleculeIdLinks(vaultId, df));
}

export function createLinksFromIds(vaultId: number, df: DG.DataFrame): void {
  const idCol = df.col('id');
  if (idCol) {
    const linkIdsCol = DG.Column.string('id', df.rowCount).init((i) => {
      const id = idCol.get(i);
      return id != null ? `[${id}](${`${CDD_HOST}/vaults/${vaultId}/molecules/${id}/`})` : '';
    });
    df.columns.replace(idCol, linkIdsCol);
  }
}

export function reorderColumns(df: DG.DataFrame): void {
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

export function paramsStringFromObj(params: Record<string, any>): string {
  const usp = new URLSearchParams();
  for (const [key, val] of Object.entries(params)) {
    if (val !== undefined && val !== null)
      usp.append(key, String(val));
  }
  const str = usp.toString();
  return str ? `?${str}` : '';
}

export function createObjectViewer(obj: any, title: string = 'Object Viewer',
  additionalHeaderEl?: HTMLElement): HTMLElement {
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
           value.every((item) =>
             typeof item === 'string' ||
             typeof item === 'number' ||
             typeof item === 'boolean',
           );
  }

  // Helper function to get pane name for array items
  function getPaneName(item: any, index: number, parentName: string): string {
    if (item.name) return item.name;
    if (item.id) return item.id.toString();
    if (parentName === PROTOCOL_STATISTICS_KEY && item.readout_definition?.name)
      return item.readout_definition.name;
    return `Item ${index}`;
  }

  // Helper function to create a view for a specific value
  function createValueView(value: any, parentName: string): HTMLElement {
    if (value === null || value === undefined)
      return ui.divText('null');


    if (isSimpleArray(value))
      return ui.divText(value.join(', '));


    if (Array.isArray(value)) {
      const accordion = ui.accordion();
      value.forEach((item, index) => {
        if (isDictionary(item) || Array.isArray(item))
          accordion.addPane(getPaneName(item, index, parentName), () => createValueView(item, parentName));
        else
          accordion.addPane(`Item ${index}`, () => ui.divText(String(item)));
      });
      return accordion.root;
    }

    if (isDictionary(value)) {
      const simpleProperties: { [key: string]: any } = {};
      const complexProperties: { [key: string]: any } = {};

      // Separate simple and complex properties
      for (const [key, val] of Object.entries(value)) {
        if (isDictionary(val) || Array.isArray(val))
          complexProperties[key] = val;
        else
          simpleProperties[key] = val;
      }

      // Create container for both simple properties table and complex properties accordion
      const container = ui.divV([]);

      // Add simple properties table if any exist
      if (Object.keys(simpleProperties).length > 0)
        container.appendChild(ui.tableFromMap(simpleProperties));


      // Add complex properties as nested accordions if any exist
      if (Object.keys(complexProperties).length > 0) {
        const accordion = ui.accordion();
        for (const [key, val] of Object.entries(complexProperties))
          accordion.addPane(key, () => createValueView(val, parentName));

        container.appendChild(accordion.root);
      }

      return container;
    }

    if (value instanceof Date)
      return ui.divText(value.toISOString());


    return ui.divText(String(value));
  }

  // Create dictionaries for simple and complex properties
  const simpleProperties: { [key: string]: any } = {};
  const complexProperties: { [key: string]: any } = {};

  // Separate simple and complex properties
  for (const [key, value] of Object.entries(obj)) {
    if (isDictionary(value) || Array.isArray(value))
      complexProperties[key] = value;
    else
      simpleProperties[key] = value;
  }

  const header = ui.divH([ui.h2(title)]);
  if (additionalHeaderEl)
    header.append(additionalHeaderEl);
  // Create the main container
  const container = ui.divV([header]);

  // Add simple properties table if any exist
  if (Object.keys(simpleProperties).length > 0)
    container.appendChild(ui.tableFromMap(simpleProperties));


  // Add complex properties as accordion if any exist
  if (Object.keys(complexProperties).length > 0) {
    const accordion = ui.accordion();
    for (const [key, value] of Object.entries(complexProperties))
      accordion.addPane(key, () => createValueView(value, key));

    container.appendChild(accordion.root);
  }

  return container;
}

export async function handleInitialURL(treeNode: DG.TreeViewGroup, url: URL) {
  const currentTabs = url.pathname.includes(`${CDD_VAULT_APP_PATH}`) ?
    url.pathname.replace(`${CDD_VAULT_APP_PATH}`, ``).replace(/^\/+/, '')
      .split('/').map((it) => decodeURIComponent(it)) : [];
  const currentVault = currentTabs.length > 0 ? currentTabs[0] : undefined;
  const currentView = currentTabs.length > 1 ? currentTabs[1] : undefined;
  const currentSubView = currentTabs.length > 2 ? currentTabs[2] : undefined;
  openCddNode(treeNode, currentVault, currentView, currentSubView);
}

export async function openCddNode(treeNode: DG.TreeViewGroup, currentVault?: string,
  currentView?: string, currentSubView?: string) {
  if (currentVault) {
    //need to wait for tree to become available
    await awaitCheck(() => treeNode.items.find((node) => node.text === currentVault) !== undefined,
      `CDD tabs haven't been loaded in 10 seconds`, 10000);
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
          await awaitCheck(() => treeNode.items.find((node) => node.text === currentView) !== undefined,
            `CDD tabs haven't been loaded in 5 seconds`, 5000);
          if (EXPANDABLE_TABS.includes(currentView)) {
            const tab = treeNode.items.find((node) => node.text === currentView) as DG.TreeViewGroup;
            currentSubView ? tab.expanded = true : tab.root.click();
            //look for inner category from URL and open it
            if (currentSubView) {
              await awaitCheck(() => tab.items.length > 0, `CDD tabs haven't been loaded in 10 seconds`, 10000);
              const innerView = tab.items.find((node) => node.text === currentSubView);
              if (!innerView) {
                tab.root.click();
                grok.shell.warning(`${currentSubView} id doesn't exist in ${currentView}`);
                return;
              }
              innerView!.root.click();
            }
          } else
            treeNode.items.find((node) => node.text === currentView)!.root.click();
        } catch (e: any) {
          grok.shell.error(`Failed to open CDD view: ${e?.message ?? e}`);
        }
      }
    } else
      grok.shell.warning(`Vault ${currentVault} doesn't exist`);
  }
}

export function createNestedCDDNode(_initialItems: any[] | null, nodeName: string, vaultNode: DG.TreeViewGroup,
  getItemsFuncName: string, getItemsFuncParams: any, treeNode: DG.TreeViewGroup, vault: Vault,
  onItemSelected: (item: any) => Promise<void>) {
  const nestedNode = vaultNode.group(nodeName, null, false);
  let items: any[] | null = _initialItems;
  let childItemsBuilt = false;
  const loadData = async (): Promise<boolean> => {
    if (items) return true;
    try {
      const itemsStr = await grok.functions.call(getItemsFuncName, getItemsFuncParams);
      items = itemsStr !== '' ? JSON.parse(itemsStr) as any[] : [];
      return true;
    } catch (e: any) {
      grok.shell.error(e?.message ?? e);
      if (openedView) {
        ui.setUpdateIndicator(openedView.root, false);
        ui.empty(openedView.root);
        openedView.root.append(ui.divText(`Error`));
      }
      return false;
    }
  };
  const buildChildItems = () => {
    if (childItemsBuilt || !items) return;
    for (const item of items) {
      const itemNode = nestedNode.item(item.name);
      itemNode.onSelected.subscribe(async () => {
        await onItemSelected(item);
        setBreadcrumbsInViewName([nodeName, item.name], treeNode);
      });
    }
    childItemsBuilt = true;
  };
  nestedNode.onSelected.subscribe(async () => {
    openedView?.close();
    openedView = DG.View.create();
    openedView.path = createPath(vault.name, [nodeName]);
    grok.shell.addPreview(openedView);
    ui.setUpdateIndicator(openedView.root, true, `Loading ${nodeName}...`);
    openedView.name = nodeName;
    setBreadcrumbsInViewName([vault.name, nodeName], treeNode);
    const ok = await loadData();
    if (!ok) return;
    const links = items!.length ? createLinks(nodeName, items!.map((it) => it.name), treeNode, openedView) :
      ui.divText(`No ${nodeName} found in the vault`);
    ui.empty(openedView.root);
    openedView.append(links);
    ui.setUpdateIndicator(openedView.root, false);
    if (!nestedNode.expanded)
      nestedNode.expanded = true;
  });

  nestedNode.onNodeExpanding.subscribe(async () => {
    const ok = await loadData();
    if (!ok) return;
    buildChildItems();
  });
}

export function createPath(vaultName: string, viewName?: string[]) {
  let path = `${CDD_VAULT_APP_PATH}/`;
  path += encodeURIComponent(vaultName);
  if (viewName)
    path += `/${viewName.map((it) => encodeURIComponent(it)).join('/')}`;
  return path;
}

function updateView(viewName: string[], vaultName: string, treeNode: DG.TreeViewGroup,
  progressMessage: string, res?: DG.DataFrame) {
  openedView?.close();
  openedView = res ? grok.shell.addTablePreview(res) : grok.shell.addPreview(DG.View.create());
  if (res)
    adjustIdColumnWidth(openedView as DG.TableView);
  openedView.name = viewName[viewName.length - 1];
  openedView.path = createPath(vaultName, viewName);
  setBreadcrumbsInViewName([vaultName].concat(viewName), treeNode, openedView);
  if (!res)
    ui.setUpdateIndicator(openedView.root, true, progressMessage);
}

/**
 * Opens a tab and loads data using a preview-then-load-all pattern.
 *
 * Flow:
 *   1. Call `syncFuncName` with `page_size: PREVIEW_ROW_NUM` — show the preview df immediately.
 *   2. If the preview returned fewer than PREVIEW_ROW_NUM rows, the full dataset is already loaded — done.
 *   3. Otherwise attach a ribbon: "Showing first N rows" + [Load all] button.
 *   4. On [Load all] click: run `asyncFuncName` in background with a DG.TaskBarProgressIndicator.
 *      When it resolves, swap the view's DataFrame and update the ribbon to "Showing all N rows".
 *   5. If the user navigates to another tab before the async completes, the result is discarded
 *      (we guard with the module-level `openedView` — the swap only happens if the view is still live).
 *
 * Pass `asyncFuncName: null` for tabs without a sync/async pair (e.g. Saved Search results):
 * the function will just run syncFuncName as a single call with no ribbon.
 */
export async function createCDDTableView(viewName: string[], progressMessage: string, syncFuncName: string,
  syncFuncParams: {[key: string]: any}, asyncFuncName: string | null, asyncFuncParams: {[key: string]: any} | null,
  vault: Vault, treeNode: DG.TreeViewGroup, addFilters?: boolean) {
  try {
    updateView(viewName, vault.name, treeNode, progressMessage);
    const viewToken = openedView; // captured to detect tab switches

    const df: DG.DataFrame = await grok.functions.call(syncFuncName, syncFuncParams);

    // User navigated away while sync was loading — drop result silently.
    if (openedView !== viewToken)
      return;

    updateView(viewName, vault.name, treeNode, progressMessage, df);
    if (addFilters && openedView)
      initializeFilters(openedView as DG.TableView, vault);

    // No async pairing, or preview already contains everything — done.
    if (!asyncFuncName || df.rowCount < PREVIEW_ROW_NUM)
      return;

    attachLoadAllRibbon(openedView as DG.TableView, viewName, asyncFuncName, asyncFuncParams ?? {});
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
    updateView(viewName, vault.name, treeNode, progressMessage, DG.DataFrame.create());
  }
}

/** Adds a "Showing first N rows / Load all" ribbon row to a TableView. Triggers the async full-fetch on click. */
function attachLoadAllRibbon(tv: DG.TableView, viewName: string[], asyncFuncName: string,
  asyncFuncParams: {[key: string]: any}) {
  const ribbonViewToken = tv; // used to check the view is still the live one at resolve time

  const info = ui.divText(`Showing first ${PREVIEW_ROW_NUM} rows`,
    {style: {alignSelf: 'center', marginRight: '8px', pointerEvents: 'none', cursor: 'default'}});
  const loadAllButton = ui.button('Load all', async () => {
    loadAllButton.disabled = true;
    info.textContent = 'Loading all rows...';
    const progressBar = DG.TaskBarProgressIndicator.create(`Loading all ${viewName[viewName.length - 1]}...`);
    try {
      const fullDf: DG.DataFrame = await grok.functions.call(asyncFuncName, asyncFuncParams);

      // User navigated away — discard, matches the "cancel on tab switch" contract.
      if (openedView !== ribbonViewToken) {
        progressBar.close();
        return;
      }

      tv.dataFrame = fullDf;
      adjustIdColumnWidth(tv);
      info.textContent = `Showing all ${fullDf.rowCount} rows`;
      loadAllButton.style.display = 'none';
    } catch (e: any) {
      grok.shell.error(e?.message ?? e);
      info.textContent = `Showing first ${PREVIEW_ROW_NUM} rows`;
      loadAllButton.disabled = false;
    } finally {
      progressBar.close();
    }
  });

  // Compose with any ribbon already set (e.g. the Filters button). Append as an extra row.
  const existing = tv.getRibbonPanels();
  tv.setRibbonPanels([...existing, [info, loadAllButton]]);
}

export async function initializeFilters(tv: DG.TableView, vault: Vault) {
  const filtersDiv = ui.divV([]);
  //create filters icon
  const externalFilterIcon = document.createElement('div');
  externalFilterIcon.innerHTML = `
<svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
  <!-- Funnel Body -->
  <path d="M4 4H15L10 10V17L8 19V10L4 4Z" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round"/>
  
  <!-- The arrow group is rotated -45 degrees around point (16, 13) -->
  <g transform="rotate(-45, 16, 13)">
    <!-- Long Horizontal Arrow Line (starts at the new x=16, y=13) -->
    <path d="M16 13H24" stroke="currentColor" stroke-width="1" stroke-linecap="round"/>
    <!-- Arrowhead (adjusted for the new starting point) -->
    <path d="M21 10L24 13L21 16" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round"/>
  </g>
</svg>`;
  externalFilterIcon.className = 'cdd-filters-button-icon';
  externalFilterIcon.onclick = () => {
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.3);
    externalFilterIcon.classList.remove('cdd-filters-button-icon-show');
  };
  ui.tooltip.bind(externalFilterIcon, 'Add CDD filters');

  const filtersButton = ui.div(externalFilterIcon);

  tv.setRibbonPanels([[filtersButton]]);
  //create filters panel
  tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.3);
  tv.dockManager.onClosed.subscribe((_el: any) => {
    externalFilterIcon.classList.add('cdd-filters-button-icon-show');
  });

  const funcEditor = new SearchEditor(vault.id);
  const acc = funcEditor.getEditor();

  const runSearch = async () => {
    ui.setUpdateIndicator(tv.grid.root, true);
    funcEditor.saveLastSearch();
    const params = funcEditor.getParams();
    try {
      const df = await funcs.cDDVaultSearchAsync(vault.id, params.structure ?? null,
        params.structure_search_type ?? null, params.structure_similarity_threshold ?? null,
        params.protocol ?? null, params.run ?? null);
      if (df) {
        const protocol = params.protocol ? `, protocol: ${params.protocol}` : '';
        const run = params.run ? `, run: ${params.run}` : '';
        const search = params.structure ? `, ${params.structure_search_type}${params.structure_search_type === CDDVaultSearchType.SIMILARITY ?
          `:${params.structure_similarity_threshold}` : ''} search for ${params.structure}` : '';

        df!.name = `Vault: ${vault.id}${protocol}${run}${search}`;
        tv.dataFrame = df;
        adjustIdColumnWidth(tv);

        // Search returns the full result — drop the preview's "Load all" row and show the final row count.
        const info = ui.divText(`Showing ${df.rowCount} rows`,
          {style: {alignSelf: 'center', marginRight: '8px', pointerEvents: 'none', cursor: 'default'}});
        tv.setRibbonPanels([[filtersButton], [info]]);
      }
    } finally {
      ui.setUpdateIndicator(tv.grid.root, false);
    }
  };
  const runSearchButton = ui.button('Search', runSearch);

  const resetIcon = ui.iconFA('redo', async () => {
    await funcEditor.reset();
    ui.setUpdateIndicator(tv.grid.root, true);
    try {
      const df: DG.DataFrame = await funcs.getMolecules(vault.id, '');
      if (df) {
        df.name = `Vault: ${vault.id}`;
        tv.dataFrame = df;
        adjustIdColumnWidth(tv);

        if (df.rowCount < PREVIEW_ROW_NUM) {
          const info = ui.divText(`Showing all ${df.rowCount} rows`,
            {style: {alignSelf: 'center', marginRight: '8px', pointerEvents: 'none', cursor: 'default'}});
          tv.setRibbonPanels([[filtersButton], [info]]);
        } else {
          tv.setRibbonPanels([[filtersButton]]);
          attachLoadAllRibbon(tv, [MOLECULES_TAB], 'CDDVaultLink:getMoleculesAsync',
            {vaultId: vault.id, moleculesIds: '', timeoutMinutes: 5});
        }
      }
    } catch (e: any) {
      grok.shell.error(e?.message ?? e);
    } finally {
      ui.setUpdateIndicator(tv.grid.root, false);
    }
  }, 'Reset search');

  filtersDiv.append(acc);
  filtersDiv.append(ui.divH([runSearchButton, resetIcon],
    {style: {paddingLeft: '4px', alignItems: 'center', gap: '8px'}}));

  const viewToken = tv;
  funcEditor.initComplete.then(() => {
    // Guard against view being closed before init resolves.
    if (openedView !== viewToken) return;
    if (funcEditor.hasRestoredSearch)
      runSearch();
  });
}

export function createLinks(header: string, nodeNames: string[], tree: DG.TreeViewGroup, view: DG.ViewBase): HTMLDivElement {
  const table = ui.table(nodeNames, (item) =>
    ([ui.link(item, () => {
      view.close();
      const target = tree.items.find((it) => it.text === item);
      if (target) tree.currentItem = target;
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
    if (actualItem === 'CDD Vault')
      tree.currentItem = tree;
    else {
      const target = tree.items.find((item) => item.text === actualItem);
      if (target) tree.currentItem = target;
    }
  });

  if (usedView) {
    if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') { // integrate it to the actual breadcrumbs element
      const homeIcon = ui.iconFA('home', () => {
        grok.shell.v.close();
        grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME);
      }, 'Home');
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
  for (const key of keys) {
    if (key === 'batches') {
      const batchesAcc = ui.accordion(key);
      const batches = (obj as Molecule)[key] as Batch[];
      batchesAcc.addPane(key, () => {
        const innerBatchesAcc = ui.accordion();
        for (const batch of batches)
          innerBatchesAcc.addPane(batch.id.toString(), () => createCDDContextPanel(batch));

        return innerBatchesAcc.root;
      });
      accordions.append(batchesAcc.root);
    } else if (key === 'collections' || key === 'projects' || key === 'source_files') {
      const acc = ui.accordion(key);
      acc.addPane(key, () => {
        const div = ui.divV([]);
        ((obj as any)[key] as any[]).forEach((it) => div.append(ui.tableFromMap(it, true)));
        return div;
      });
      accordions.append(acc.root);
    } else if (key === 'molecule_fields' || key === 'udfs' || key === 'stoichiometry' || key === 'batch_fields') {
      const acc = ui.accordion(key);
      acc.addPane(key, () => ui.tableFromMap((obj as any)[key] as any, true));
      accordions.append(acc.root);
    } else {
      const initialValue = (obj as any)[key];
      let value: any | null = null;
      if (key === 'id' && vaultId)
        value = ui.link(initialValue, () => window.open(`${CDD_HOST}/vaults/${vaultId}/molecules/${initialValue}/`));
      else
        value = initialValue instanceof Array ? initialValue.join(', ') : initialValue;
      dictForTableView[key] = value;
    }
  }
  resDiv.append(ui.tableFromMap(dictForTableView, true));
  resDiv.append(accordions);

  return resDiv;
}


export function createInitialStatistics(statsDiv: HTMLDivElement) {
  funcs.getVaults().then(async (res: string) => {
    const stats: CDDVaultStats[] = [];
    if (!res)
      return;
    const vaults = JSON.parse(res) as Vault[];
    for (const vault of vaults) {
      try {
        const resStr = await funcs.getVaultStats(vault.id, vault.name);
        stats.push(JSON.parse(resStr));
      } catch (e: any) {
        grok.shell.error(`Cannot get statistics for vault ${vault.name}: ${e?.message ?? e}`);
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
  });
}

export function createVaultNode(vault: Vault, treeNode: DG.TreeViewGroup) {
  openedView?.close();
  openedView = DG.View.create();
  openedView.name = vault.name;
  const tabs = createLinks(vault.name, [PROTOCOLS_TAB, SAVED_SEARCHES_TAB, COLLECTIONS_TAB, MOLECULES_TAB, BATCHES_TAB, SEARCH_TAB], treeNode, openedView);
  openedView.append(tabs);
  grok.shell.addPreview(openedView);
  setBreadcrumbsInViewName([vault.name], treeNode, openedView);
  openedView.path = createPath(vault.name);
}

export function addNodeWithEmptyResults(name: string, warningMessage?: string) {
  openedView?.close();
  const df = DG.DataFrame.create();
  df.name = name;
  openedView = grok.shell.addTablePreview(df);
  if (warningMessage)
    grok.shell.warning(warningMessage);
}

async function adjustIdColumnWidth(tv: DG.TableView) {
  await DG.delay(100);
  const idCol = tv.grid.col('id');
  if (idCol)
    idCol.width = 100;
}

