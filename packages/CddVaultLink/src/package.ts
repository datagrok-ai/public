/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import {MoleculeFieldSearch, getVaults, MoleculeQueryParams, queryMolecules, queryReadoutRows, Molecule, Batch, querySavedSearches, SavedSearch, querySavedSearchById, queryExportStatus, queryExportResult,
  queryMoleculesAsync, queryReadoutRowsAsync, ApiResponse, MoleculesQueryResult, ProtocolQueryResult, Protocol, queryProtocolsAsync, Vault, Collection,
  queryCollectionsAsync} from "./cdd-vault-api";
import { CDDVaultSearchType } from './constants';
import '../css/cdd-vault.css';
import { SeachEditor } from './search-function-editor';
import { CDD_HOST, createLinksFromIds, createObjectViewer, getAsyncResults, getAsyncResultsAsDf, prepareDataForDf, reorderColummns } from './utils';

export const _package = new DG.Package();

export const PREVIEW_ROW_NUM = 100;

//tags: app
//name: CDD Vault
//meta.icon: images/cdd-icon-small.png
//output: view v
//meta.browsePath: Chem
export async function cddVaultApp(): Promise<DG.ViewBase> {

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/cdd-icon-big.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/CddVaultLink/README.md',
    description: '- Integrate with your CDD Vault.\n' +
      '- Analyze assay data.\n' +
      '- Find contextual information on molecules.\n' +
      '- Browse the vault content.\n'
  });
  const view = DG.View.fromRoot(appHeader);
  view.name = 'CDD Vault';
  return view;
}

//input: dynamic treeNode
//input: view browseView
export async function cddVaultAppTreeBrowser(treeNode: DG.TreeViewGroup) {
  const vaults = await getVaults();
  if (vaults.error) {
    grok.shell.error(vaults.error);
    return;
  }

  if(!vaults?.data?.length) {
    grok.shell.error(`No vaults found`);
    return;
  }

  for (const vault of vaults.data!) {
    //vault node
    const vaultNode = treeNode.group(vault.name);
    vaultNode.onSelected.subscribe(() => {
      const view = DG.View.create();
      view.name = vault.name;
      const tabs = createLinks(['Protocols', 'Saved searches', 'Collections', 'Molecules', 'Search'], treeNode, view);
      view.append(tabs);
      grok.shell.addPreview(view);
      setBreadcrumbsInViewName([vault.name], treeNode, view);
    });

    //protocols node

    let protocols: Protocol[] | null = null;
    createNestedCDDNode(protocols, 'Protocols', vaultNode, 'CDDVaultLink:getProtocolsAsync', { vaultId: vault.id, timeoutMinutes: 5 }, treeNode, vault,
      async (item: any) => {
        createCDDTableView(['Protocols', item.name], 'Waiting for molecules', 'CDDVaultLink:cDDVaultSearchAsync',
          {
            vaultId: vault.id, structure: '', structure_search_type: CDDVaultSearchType.SUBSTRUCTURE,
            structure_similarity_threshold: 0, protocol: item.id, run: undefined
          }, vault.name, treeNode);
        grok.shell.windows.context.visible = true;
        grok.shell.o = createObjectViewer(item, item.name);
      }
    );

    //saved searches node
    let savedSearches: SavedSearch[] | null = null; 
    createNestedCDDNode(savedSearches, 'Saved searches', vaultNode, 'CDDVaultLink:getSavedSearches', { vaultId: vault.id }, treeNode, vault,
      async (item: any) => {
        createCDDTableView(['Saved searches', item.name], `Waiting for ${item.name} results`, 'CDDVaultLink:getSavedSearchResults',
          { vaultId: vault.id, searchId: item.id, timeoutMinutes: 5}, vault.name, treeNode);
      }
    );

    //collections
    let collections: Collection[] | null = null;
    createNestedCDDNode(collections, 'Collections', vaultNode, 'CDDVaultLink:getCollectionsAsync', { vaultId: vault.id, timeoutMinutes: 5 }, treeNode, vault,
      async (item: any) => {
        //in case collection doesn't contain molecules - add empty tableView
        if (!item.molecules || !item.molecules.length) {
          const df = DG.DataFrame.create();
          df.name = item.name;
          grok.shell.addTablePreview(DG.DataFrame.create());
          grok.shell.warning(`No molecules found for ${item.name} collection`);
        }
        createCDDTableViewWithPreview(['Collections', item.name], `Waiting for ${item.name} results`,
          'CDDVaultLink:getMolecules',
          {
            vaultId: vault.id,
            moleculesIds: item.molecules.join(',')
          },
          'CDDVaultLink:getMoleculesAsync',
          { 
            vaultId: vault.id,
            moleculesIds: item.molecules.join(','),
            timeoutMinutes: 5
          }, vault.name, treeNode);
      }
    );

    //molecules node
    const moleculesNode = vaultNode.item('Molecules');
    moleculesNode.onSelected.subscribe(async (_) => {
      createCDDTableViewWithPreview(['Molecules'], 'Waiting for molecules', 'CDDVaultLink:getMolecules',
        {vaultId: vault.id, moleculesIds: ''}, 'CDDVaultLink:getMoleculesAsync', { vaultId: vault.id, moleculesIds: '', timeoutMinutes: 5}, vault.name, treeNode);
    });

    //search node
    const searchNode = vaultNode.item('Search');
    searchNode.onSelected.subscribe(() => {
      const view = DG.View.create();
      const funcEditor = new SeachEditor(vault.id);
      const acc = funcEditor.getEditor();
      let df: DG.DataFrame | null = null;
      const runButton = ui.bigButton('SEARCH', async () => {
        ui.setUpdateIndicator(gridDiv, true);
        const params = funcEditor.getParams();
        df = await grok.functions.call('CDDVaultLink:cDDVaultSearchAsync',
          {
            vaultId: vault.id, structure: params.structure, structure_search_type: params.structure_search_type,
            structure_similarity_threshold: params.structure_similarity_threshold, protocol: params.protocol, run: params.run 
          });
        ui.empty(gridDiv);
        if (df) {
          const protocol = params.protocol ? `, protocol: ${params.protocol}` : '';
          const run = params.run ? `, run: ${params.run}` : '';
          const search = params.structure ? `, ${params.structure_search_type}${params.structure_search_type === CDDVaultSearchType.SIMILARITY ?
            `:${params.structure_similarity_threshold}`: ''} search for ${params.structure}` : '';

          df!.name = `Vault: ${vault.id}${protocol}${run}${search}`;
          gridDiv.append(df.plot.grid().root);
        }
        ui.setUpdateIndicator(gridDiv, false);
      });
      const gridDiv = ui.div('', 'cdd-vault-search-res-div');
      runButton.classList.add('cdd-vault-run-search-button');

      const addToWorkspaceButton = ui.icons.add(() => {
        if (df)
          grok.shell.addTablePreview(df);
      }, 'Add results to workspace');
      view.setRibbonPanels([[addToWorkspaceButton]]);
      view.name = 'Search CDD Vault'
      view.root.append(ui.divV([
        acc,
        runButton,
        gridDiv
      ], {style: {height: '100%'}}));
      grok.shell.addPreview(view);
      setBreadcrumbsInViewName([vault.name, 'Search'], treeNode, view);
    });
   //TODO! unlock other tabs
   // vaultNode.group('Plates');
   // vaultNode.group('Assays');
  }
}

function createNestedCDDNode(items: any[] | null, nodeName: string, vaultNode: DG.TreeViewGroup,
  getItemsFunsName: string, getItemsFuncParams: any, treeNode: DG.TreeViewGroup, vault: Vault,
  onItemSelected: (item: any) => Promise<void>) {
  const nestedNode = vaultNode.group(nodeName, null, false);
  const loadData = async () => {
    if (!items) {
      const itemsStr = await grok.functions.call(getItemsFunsName, getItemsFuncParams);
      items = itemsStr !== '' ? JSON.parse(itemsStr) as any[] : [];
    }
  }
  nestedNode.onSelected.subscribe(async () => {
    nestedNode.expanded = true;
    await loadData();
    const view = DG.View.create();
    view.name = nodeName;
    const tabs = createLinks(items!.map((it) => it.name), treeNode, view);
    view.append(tabs);
    grok.shell.addPreview(view);
    setBreadcrumbsInViewName([vault.name, nodeName], treeNode);
  });

  nestedNode.onNodeExpanding.subscribe(async () => {
    await loadData();
    for (const item of items!) {
      const protocolItem = nestedNode.item(item.name);
      protocolItem.onSelected.subscribe(async () => {
        await onItemSelected(item);
        setBreadcrumbsInViewName([nodeName, item.name], treeNode);
      });
    }
  });
}

async function createCDDTableView(viewName: string[], progressMessage: string, funcName: string,
  funcParams: {[key: string]: any}, vaultName: string, treeNode: DG.TreeViewGroup) {
  const view = DG.View.create();
  view.name = viewName[viewName.length - 1];
  grok.shell.addPreview(view);
  ui.setUpdateIndicator(view.root, true, progressMessage);
  const df: DG.DataFrame = await grok.functions.call(funcName, funcParams);
  view.close();
  df.name = viewName[viewName.length - 1];
  const tv = grok.shell.addTablePreview(df);
  setBreadcrumbsInViewName([vaultName].concat(viewName), treeNode, tv);
}

async function createCDDTableViewWithPreview(viewName: string[], progressMessage: string, syncfuncName: string,
  syncfuncParams: {[key: string]: any}, asyncfuncName: string, asyncfuncParams: {[key: string]: any},
  vaultName: string, treeNode: DG.TreeViewGroup) {
  const view = DG.View.create();
  view.name = viewName[viewName.length - 1];
  grok.shell.addPreview(view);
  ui.setUpdateIndicator(view.root, true, progressMessage);
  let tv: DG.TableView | null = null;
  let updateDataWithAsyncResults = true;
  //run sync function with offset and create a preview
  grok.functions.call(syncfuncName, syncfuncParams).then((res: DG.DataFrame) => {
    if (!tv) {
      if (res.rowCount < PREVIEW_ROW_NUM) {
        updateDataWithAsyncResults = false;
        progressBar.close();
      }
      res.name = viewName[viewName.length - 1];
      tv = grok.shell.addTablePreview(res);
      view.close();
      setBreadcrumbsInViewName([vaultName].concat(viewName), treeNode, tv);
    }
  });
  //reset tableView with asynchronously received results
  grok.functions.call(asyncfuncName, asyncfuncParams).then((res: DG.DataFrame) => {
    if (!updateDataWithAsyncResults)
      return;
    if (tv)
      tv.close();
    else {
      view.close();
    }
    res.name = viewName[viewName.length - 1];
    tv = grok.shell.addTablePreview(res);
    setBreadcrumbsInViewName([vaultName].concat(viewName), treeNode, tv);
    progressBar.close();
  });
  const progressBar = DG.TaskBarProgressIndicator.create(`Loading ${viewName[viewName.length - 1]}...`);
}

function createLinks(nodeNames: string[], tree: DG.TreeViewGroup, view: DG.ViewBase): HTMLDivElement {
  const div = ui.divV([]);
  for (const name of nodeNames) {
    const button = ui.link(name, () => {
      view.close();
      tree.currentItem = tree.items.find((item) => item.text === name)!
    }, 'Click to open');
    button.classList.add('cdd-menu-point');
    div.append(button);
  }
  return div;
}

function setBreadcrumbsInViewName(viewPath: string[], tree: DG.TreeViewGroup, view?: DG.View): void {
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


//name: Databases | CDD Vault
//input: string mol {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molecule: string): DG.Widget {
  return DG.Widget.fromRoot(ui.wait(async () => {
    const vaults = await getVaults();
    if (vaults.error) {
      return ui.divText(vaults.error);
    }
    const vaultId = vaults.data![0].id;
    const cddMols = await queryMolecules(vaultId, { structure: molecule, structure_search_type: "exact"});

    if (!cddMols.data?.objects?.length)
      return ui.divText('Not found');

    return createCDDContextPanel(cddMols.data.objects[0], vaultId);
  }));
}

function createCDDContextPanel(obj: Molecule | Batch, vaultId?: number): HTMLElement {
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

//name: CDDVaultSearchEditor
//tags: editor
//input: funccall call
export async function CDDVaultSearchEditor(call: DG.FuncCall): Promise<void> {
  const vaults = await getVaults();
  if (vaults.error) {
    grok.shell.error(vaults.error);
    return;
  }
  const vaultId = vaults.data![0].id;
  const funcEditor = new SeachEditor(vaultId);
  const dialog = ui.dialog({title: 'CDD search'})
    .add(funcEditor.getEditor())
    .onOK(async () => {
      const params = funcEditor.getParams();
    });
  //dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
  dialog.show();
}

//name: CDD Vault search 2
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: string structure {category: Structure; nullable: true; semType: Molecule} [SMILES; cxsmiles or mol string]
//input: string structure_search_type {category: Structure; nullable: true; choices: ["exact", "similarity", "substructure"]} [SMILES, cxsmiles or mol string]
//input: double structure_similarity_threshold {category: Structure; nullable: true} [A number between 0 and 1]
//input: int protocol {category: Protocol; nullable: true} [Protocol id]
//input: int run {category: Protocol; nullable: true} [Specific run id]
//output: dataframe df
//editor: Cddvaultlink:CDDVaultSearchEditor
export async function cDDVaultSearch2(vaultId: number, structure?: string, structure_search_type?: CDDVaultSearchType,
  structure_similarity_threshold?: number, protocol?: number, run?: number): Promise<DG.DataFrame> {
  //collecting molecule ids according to protocol query params
  const molIds: number[] = [];
  if (protocol) {
    //TODO! Make async request and remove page size
    const readoutRowsRes = await queryReadoutRows(vaultId, {protocols: protocol.toString(), runs: run?.toString(), page_size: 1000});
    if (readoutRowsRes.error) {
      grok.shell.error(readoutRowsRes.error);
      return DG.DataFrame.create();
    }
    const readoutRows= readoutRowsRes.data?.objects;
    if (readoutRows) {
      for (const readoutRow of readoutRows)
        if (!molIds.includes(readoutRow.molecule))
          molIds.push(readoutRow.molecule)
    }      
  }
  const molQueryParams: MoleculeQueryParams = !structure ? {molecules: molIds.join(',')} :
    {structure: structure, structure_search_type: structure_search_type, structure_similarity_threshold: structure_similarity_threshold};

  //TODO! Make async request and remove page size
  molQueryParams.page_size = 1000;
  const cddMols = await queryMolecules(vaultId, molQueryParams);
  if (cddMols.error) {
    grok.shell.error(cddMols.error);
    return DG.DataFrame.create();
  }
  if (cddMols.data?.objects && cddMols.data?.objects.length) {
    //in case we had both protocol and structure conditions - combine results together
    const molsRes = protocol && structure ? cddMols.data!.objects!.filter((it) => molIds.includes(it.id)) : cddMols.data?.objects;
    prepareDataForDf(cddMols.data?.objects);
    const df = DG.DataFrame.fromObjects(molsRes)!;
    if (!df)
      return DG.DataFrame.create();
    createLinksFromIds(vaultId, df);
    reorderColummns(df);
    await grok.data.detectSemanticTypes(df);
    return df;
  }
  return DG.DataFrame.create();
}

//name: CDD Vault Search Async 
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: string structure {category: Structure; nullable: true; semType: Molecule} [SMILES; cxsmiles or mol string]
//input: string structure_search_type {category: Structure; nullable: true; choices: ["exact", "similarity", "substructure"]} [SMILES, cxsmiles or mol string]
//input: double structure_similarity_threshold {category: Structure; nullable: true} [A number between 0 and 1]
//input: int protocol {category: Protocol; nullable: true} [Protocol id]
//input: int run {category: Protocol; nullable: true} [Specific run id]
//output: dataframe df
//editor: Cddvaultlink:CDDVaultSearchEditor
export async function cDDVaultSearchAsync(vaultId: number, structure?: string, structure_search_type?: CDDVaultSearchType,
  structure_similarity_threshold?: number, protocol?: number, run?: number): Promise<DG.DataFrame> {
  //collecting molecule ids according to protocol query params
  let molIds: number[] = [];
  if (protocol) {
    const exportResponse = await queryReadoutRowsAsync(vaultId, {protocols: protocol.toString(), runs: run?.toString()});
    const readoutRowsRes = await getAsyncResults(vaultId, exportResponse, 5, false);
    if (!readoutRowsRes)
      return DG.DataFrame.create();
    const readoutRows = readoutRowsRes.data?.objects;
    if (readoutRows) {
      for (const readoutRow of readoutRows)
        if (!molIds.includes(readoutRow.molecule))
          molIds.push(readoutRow.molecule)
    }    
    
  }
  const molQueryParams: MoleculeQueryParams = !structure ? {molecules: molIds.join(',')} :
    {structure: structure, structure_search_type: structure_search_type, structure_similarity_threshold: structure_similarity_threshold};

  const moleculesExportResponse = await queryMoleculesAsync(vaultId, molQueryParams);
  const cddMols = await getAsyncResults(vaultId, moleculesExportResponse, 5, false) as ApiResponse<MoleculesQueryResult>;
  if (!cddMols)
    return DG.DataFrame.create();

  if (cddMols.data?.objects && cddMols.data?.objects.length) {
    //in case we had both protocol and structure conditions - combine results together
    const molsRes = protocol && structure ? cddMols.data!.objects!.filter((it) => molIds.includes(it.id)) : cddMols.data?.objects;
    prepareDataForDf(cddMols.data?.objects);
    const df = DG.DataFrame.fromObjects(molsRes)!;
    if (!df)
      return DG.DataFrame.create();
    createLinksFromIds(vaultId, df);
    reorderColummns(df);
    await grok.data.detectSemanticTypes(df);
    return df;
  }
  return DG.DataFrame.create();
}

//name: Get Molecules
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: string moleculesIds
//output: dataframe df
export async function getMolecules(vaultId: number, moleculesIds: string): Promise<DG.DataFrame> {
  const params: {[key: string]: any} = {page_size: PREVIEW_ROW_NUM};
  if (moleculesIds)
    params.molecules = moleculesIds;
  const molecules = await queryMolecules(vaultId, params); //TODO! Make async request and remove page size
  if (molecules.error) {
    grok.shell.error(molecules.error);
    return DG.DataFrame.create();
  }
  let df: DG.DataFrame | null = null;
  if (molecules.data?.objects) {
    prepareDataForDf(molecules.data?.objects as any[]);
    df = DG.DataFrame.fromObjects(molecules.data.objects)!;
  }
  if (!df)
    return DG.DataFrame.create();
  createLinksFromIds(vaultId, df);
  reorderColummns(df);
  await grok.data.detectSemanticTypes(df);
  return df;
}


//name: Get Molecules Async
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: string moleculesIds
//input: int timeoutMinutes
//output: dataframe df
export async function getMoleculesAsync(vaultId: number, moleculesIds: string, timeoutMinutes: number): Promise<DG.DataFrame> {
  const params: {[key: string]: any} = {};
  if (moleculesIds)
    params.molecules = moleculesIds;
  const exportResponse = await queryMoleculesAsync(vaultId, params);
  const df = await getAsyncResultsAsDf(vaultId, exportResponse, timeoutMinutes, false);
  createLinksFromIds(vaultId, df);
  reorderColummns(df);
  return df;
}

//name: Get Protocols Async
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: int timeoutMinutes
//output: string protocols
export async function getProtocolsAsync(vaultId: number, timeoutMinutes: number): Promise<string> {
  const exportResponse = await queryProtocolsAsync(vaultId);
  const protocols = await getAsyncResults(vaultId, exportResponse, timeoutMinutes, false);
  return protocols?.data?.objects ? JSON.stringify(protocols.data.objects) : '';
}

//name: Get Collections Async
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: int timeoutMinutes
//output: string protocols
export async function getCollectionsAsync(vaultId: number, timeoutMinutes: number): Promise<string> {
  const exportResponse = await queryCollectionsAsync(vaultId, {include_molecule_ids: true});
  const collections = await getAsyncResults(vaultId, exportResponse, timeoutMinutes, false);
  return collections?.data?.objects ? JSON.stringify(collections.data.objects) : '';
}

//name: Get Saved Searches
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//output: string result
export async function getSavedSearches(vaultId: number): Promise<string> {
  const savedSearches = await querySavedSearches(vaultId);
  if (savedSearches.error) {
    grok.shell.error(savedSearches.error);
    return '';
  }
  return savedSearches.data ? JSON.stringify(savedSearches.data) : '';
}

//name: Get Saved Search Results
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId
//input: int searchId
//input: int timeoutMinutes
//output: dataframe df
export async function getSavedSearchResults(vaultId: number, searchId: number, timeoutMinutes: number): Promise<DG.DataFrame> {
  const exportResponse = await querySavedSearchById(vaultId, searchId);
  const res = await getAsyncResultsAsDf(vaultId, exportResponse, timeoutMinutes, true);
  reorderColummns(res);
  return res;
}




//name: CDD Vault search
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: string molecules {category: General; nullable: true} [Comma separated list of ids] 
//input: string names {category: General; nullable: true} [Comma separated list of names/synonyms]
//input: bool include_original_structures {category: General; nullable: true} [If true, include the original user defined structure for each molecule]
//input: bool only_ids {category: General; nullable: true} [If true, only the Molecule IDs are returned, allowing for a smaller and faster response]
//input: bool only_batch_ids {category: General; nullable: true} [If true, the full Molecule details are still returned but the Batch-level information is left out of the JSON results. (Only the IDs of the Batches belonging to the Molecules are still included.)]
//input: string created_before {category: General; nullable: true} [ISO 8601 date] 
//input: string created_after {category: General; nullable: true} [ISO 8601 date] 
//input: string modified_before {category: General; nullable: true} [ISO 8601 date] 
//input: string modified_after {category: General; nullable: true} [ISO 8601 date]
//input: string batch_created_before {category: Batch fields; nullable: true} [ISO 8601 date. A molecule with any batch that has a creation date on or before the parameter will be included] 
//input: string batch_created_after {category: Batch fields; nullable: true} [ISO 8601 date. A molecule with any batch that has a creation date on or after the parameter will be included] 
//input: string batch_field_before_name {category: Batch fields; nullable: true} [Specifes a user-defined batch field for batch_field_before_date] 
//input: string batch_field_before_date {category: Batch fields; nullable: true} [ISO 8601 date. A molecule with any batch that has a batch_field_before_name value date on or before the parameter will be included] 
//input: string batch_field_after_name {category: Batch fields; nullable: true} [Specifes a user-defined batch field for batch_field_after_date] 
//input: string batch_field_after_date {category: Batch fields; nullable: true} [ISO 8601 date. A molecule with any batch that has a batch_field_after_name value date on or after the parameter will be included] 
//input: string projects {category: Projects; nullable: true} [Comma separated list of project ids]
//input: string data_sets {category: Datasets; nullable: true} [Comma separated list of dataset ids] 
//input: string structure {category: Structure; nullable: true; semType: Molecule} [SMILES; cxsmiles or mol string]
//input: string structure_search_type {category: Structure; nullable: true; choices: ["exact", "similarity", "substructure"]} [SMILES, cxsmiles or mol string]
//input: double structure_similarity_threshold {category: Structure; nullable: true} [A number between 0 and 1]
//input: string inchikey {category: Structure; nullable: true} [Use this parameter instead of the "structure" and "structure_search_type" parameters]
//input: list<string> molecule_fields {category: Filelds; nullable: true} [Use this parameter to limit the number of Molecule UDF Fields to return]
//input: list<string> batch_fields {category: Filelds; nullable: true} [Use this parameter to limit the number of Batch UDF Fields to return]
//input: list<string> fields_search {category: Molecules filelds search; nullable: true} [This parameter is used for searching across the custom user-defined Molecule fields created by your Vault Administrator]
//output: dataframe df
export async function cDDVaultSearch(vaultId: number, molecules: string, names: string, include_original_structures: boolean,
  only_ids: boolean, only_batch_ids: boolean, created_before: string, created_after: string, modified_before: string,
  modified_after: string, batch_created_before: string, batch_created_after: string, batch_field_before_name: string,
  batch_field_before_date: string, batch_field_after_name: string, batch_field_after_date: string, projects: string, data_sets: string,
  structure: string, structure_search_type: CDDVaultSearchType, structure_similarity_threshold: number, inchikey: string, 
  molecule_fields: string[], batch_fields: string[], fields_search: string[]
): Promise<DG.DataFrame> {
  const params: MoleculeQueryParams = {};
  if (molecules)
    params.molecules = molecules;
  if (names)
    params.names = names;
  params.include_original_structures = include_original_structures;
  params.only_ids = only_ids;
  params.only_batch_ids = only_batch_ids;
  if (created_before)
    params.created_before = created_before;
  if(created_after)
    params.created_after = created_after;
  if (modified_before)
    params.modified_before = modified_before;
  if(modified_after)
    params.modified_after = modified_after;
  if(batch_created_before)
    params.batch_created_before = batch_created_before;
  if (batch_created_after)
    params.batch_created_after = batch_created_after;
  if (batch_field_before_name)
    params.batch_field_before_name = batch_field_before_name;
  if (batch_field_before_date)
    params.batch_field_before_date = batch_field_before_date;
  if (batch_field_after_name)
    params.batch_field_after_name = batch_field_after_name;
  if (batch_field_after_date)
    params.batch_field_after_date = batch_field_after_date;
  if (projects)
    params.projects = projects;
  if (data_sets)
    params.data_sets = data_sets;
  if (structure)
    params.structure = structure;
  if (structure_search_type)
    params.structure_search_type = structure_search_type;
  if (structure_similarity_threshold)
    params.structure_similarity_threshold = structure_similarity_threshold;
  if (inchikey)
    params.inchikey = inchikey;
  if (molecule_fields)
    params.molecule_fields = molecule_fields;
  if (batch_fields)
    params.batch_fields = batch_fields;
  if (fields_search) {
    const fieldsValues: MoleculeFieldSearch[] = [];
    for (const val of fields_search) {
      const fieldParamJson = JSON.parse(val) as MoleculeFieldSearch;
      fieldsValues.push(fieldParamJson);
    }
    params.fields_search = fieldsValues;
  }


  const cddMols = await queryMolecules(vaultId, params);
  let df: DG.DataFrame | null = null;
  if (cddMols.data?.objects) {
    prepareDataForDf(cddMols.data?.objects as any[]);
    df = DG.DataFrame.fromObjects(cddMols.data.objects)!;
  }
  if (!df)
    return DG.DataFrame.create();
  if (params.fields_search) {
    for (const moleculeUDF of params.fields_search) {
      const colType = moleculeUDF.text_value ? DG.TYPE.STRING : moleculeUDF.float_value ? DG.TYPE.FLOAT : moleculeUDF.date_value ? DG.TYPE.DATE_TIME : null;
      if (colType) {
        const col = DG.Column.fromType(colType,  moleculeUDF.name, df.rowCount).init((i) => cddMols.data!.objects![i].molecule_fields![moleculeUDF.name]);
        df.columns.add(col);
      }
    }
  }
  await grok.data.detectSemanticTypes(df);
  return df;
}
