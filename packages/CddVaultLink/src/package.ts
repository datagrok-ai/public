/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import {MoleculeFieldSearch, queryVaults, MoleculeQueryParams, queryMolecules, queryReadoutRows, querySavedSearches, SavedSearch, querySavedSearchById,
  queryMoleculesAsync, queryReadoutRowsAsync, ApiResponse, MoleculesQueryResult, Protocol, queryProtocolsAsync, Collection,
  queryCollectionsAsync,
  queryBatchesAsync,
  queryProjects,
  Vault} from "./cdd-vault-api";
import { CDDVaultSearchType, COLLECTIONS_TAB, MOLECULES_TAB, PROTOCOLS_TAB, SAVED_SEARCHES_TAB, SEARCH_TAB } from './constants';
import '../css/cdd-vault.css';
import { SeachEditor } from './search-function-editor';
import { addNodeWithEmptyResults, CDDVaultStats, createCDDContextPanel, createCDDTableView, createCDDTableViewWithPreview, createInitialSatistics, createLinks,
  createLinksFromIds, createMoleculesDfFromObjects, createNestedCDDNode, createObjectViewer, createPath, createSearchNode, createVaultNode, getAsyncResults, getAsyncResultsAsDf,
  getExportId, handleInitialURL, prepareDataForDf, PREVIEW_ROW_NUM, reorderColummns, setBreadcrumbsInViewName } from './utils';

export const _package = new DG.Package();

//tags: app
//name: CDD Vault
//meta.icon: images/cdd-icon-small.png
//input: string path {meta.url: true; optional: true}
//input: string filter {optional: true}
//output: view v
//meta.browsePath: Chem
export async function cddVaultApp(path: string, filter: string): Promise<DG.ViewBase> {

  const initialUrl = new URL(window.location.href);

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/cdd-icon-big.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/CddVaultLink/README.md',
    description: '- Integrate with your CDD Vault.\n' +
      '- Analyze assay data.\n' +
      '- Find contextual information on molecules.\n' +
      '- Browse the vault content.\n'
  });

  const statsDiv = ui.div('', {style: {position: 'relative'}});
  ui.setUpdateIndicator(statsDiv, true, 'Loading vaults statistics...');
  const view = DG.View.fromRoot(ui.divV([appHeader, statsDiv]));
  view.name = 'CDD Vault';

  if (path) {
    const cddNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('CDD Vault');
    cddNode.expanded = true;
    handleInitialURL(cddNode, initialUrl);
  } else
    createInitialSatistics(statsDiv);


  return view;
}

//input: dynamic treeNode
//input: view browseView
export async function cddVaultAppTreeBrowser(treeNode: DG.TreeViewGroup) {
  try {
    const vaults = JSON.parse(await grok.functions.call('CDDVaultLink:getVaults')) as Vault[];

    for (const vault of vaults) {
      //vault node
      const vaultNode = treeNode.group(vault.name);
      vaultNode.onSelected.subscribe(() => {
        createVaultNode(vault, treeNode);
      });
  
      //protocols node
  
      let protocols: Protocol[] | null = null;
      createNestedCDDNode(protocols, PROTOCOLS_TAB, vaultNode, 'CDDVaultLink:getProtocolsAsync', { vaultId: vault.id, timeoutMinutes: 5 }, treeNode, vault,
        async (item: any) => {
          createCDDTableView([PROTOCOLS_TAB, item.name], 'Waiting for molecules', 'CDDVaultLink:cDDVaultSearchAsync',
            {
              vaultId: vault.id, structure: '', structure_search_type: CDDVaultSearchType.SUBSTRUCTURE,
              structure_similarity_threshold: 0, protocol: item.id, run: undefined
            }, vault.name, treeNode);
          grok.shell.windows.context.visible = true;
          grok.shell.o = createObjectViewer(item, item.name);
        }
      );
  
      //saved searches node - only asyn method is available (so createCDDTableViewWithPreview function is not applicable)
      let savedSearches: SavedSearch[] | null = null; 
      createNestedCDDNode(savedSearches, SAVED_SEARCHES_TAB, vaultNode, 'CDDVaultLink:getSavedSearches', { vaultId: vault.id }, treeNode, vault,
        async (item: any) => {
          createCDDTableView([SAVED_SEARCHES_TAB, item.name], `Waiting for ${item.name} results`, 'CDDVaultLink:getSavedSearchResults',
            { vaultId: vault.id, searchId: item.id, timeoutMinutes: 5}, vault.name, treeNode);
        }
      );
  
      //collections
      let collections: Collection[] | null = null;
      createNestedCDDNode(collections, COLLECTIONS_TAB, vaultNode, 'CDDVaultLink:getCollectionsAsync', { vaultId: vault.id, timeoutMinutes: 5 }, treeNode, vault,
        async (item: any) => {
          //in case collection doesn't contain molecules - add empty tableView
          if (!item.molecules || !item.molecules.length) {
            addNodeWithEmptyResults(item.name, `No molecules found for ${item.name} collection`);
          }
          createCDDTableViewWithPreview([COLLECTIONS_TAB, item.name], `Waiting for ${item.name} results`,
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
      const moleculesNode = vaultNode.item(MOLECULES_TAB);
      moleculesNode.onSelected.subscribe(async (_) => {
        createCDDTableViewWithPreview([MOLECULES_TAB], 'Waiting for molecules', 'CDDVaultLink:getMolecules',
          {vaultId: vault.id, moleculesIds: ''}, 'CDDVaultLink:getMoleculesAsync', { vaultId: vault.id, moleculesIds: '', timeoutMinutes: 5}, vault.name, treeNode);
      });
  
      //search node
      const searchNode = vaultNode.item(SEARCH_TAB);
      searchNode.onSelected.subscribe(() => {
        createSearchNode(vault, treeNode);
      });
     //TODO! unlock other tabs
     // vaultNode.group('Plates');
     // vaultNode.group('Assays');
    }
   // handleInitialURL(treeNode);
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
  }
}


//name: Databases | CDD Vault
//input: string mol {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molecule: string): DG.Widget {
  return DG.Widget.fromRoot(ui.wait(async () => {
    try {
      const vaults = JSON.parse(await grok.functions.call('CDDVaultLink:getVaults')) as Vault[];
      //looking for molecule in the first vault
      const vaultId = vaults[0].id;
      const cddMols = await queryMolecules(vaultId, { structure: molecule, structure_search_type: "exact"});
  
      if (!cddMols.data?.objects?.length)
        return ui.divText('Not found');
  
      return createCDDContextPanel(cddMols.data.objects[0], vaultId);
    } catch(e: any) {
      return ui.divText(e?.message ?? e);
    }
  }));
}

//name: CDDVaultSearchEditor
//tags: editor
//input: funccall call
export async function CDDVaultSearchEditor(call: DG.FuncCall): Promise<void> { //is not used at the moment
  try {
    const vaults = JSON.parse(await grok.functions.call('CDDVaultLink:getVaults')) as Vault[];
    const vaultId = vaults[0].id;
    const funcEditor = new SeachEditor(vaultId);
    const dialog = ui.dialog({title: 'CDD search'})
      .add(funcEditor.getEditor())
      .onOK(async () => {
        const params = funcEditor.getParams();
      });
    //dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
    dialog.show();
  } catch(e: any) {
    grok.shell.error(e?.message ?? e);
  }
}

//name: Get Vault Stats
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId
//input: string vaultName
//output: string stats
export async function getVaultStats(vaultId: number, vaultName: string): Promise<string> {
  const promises: Promise<any>[] = [];
  promises.push(queryProjects(vaultId));
  promises.push(getAsyncResults(vaultId, getExportId(await queryMoleculesAsync(vaultId, {only_ids: true})), 1, false));
  promises.push(getAsyncResults(vaultId, getExportId(await queryProtocolsAsync(vaultId, {only_ids: true})), 1, false));
  promises.push(getAsyncResults(vaultId, getExportId(await queryCollectionsAsync(vaultId, {only_ids: true})), 1, false));
  promises.push(getAsyncResults(vaultId, getExportId(await queryBatchesAsync(vaultId, {only_ids: true})), 1, false));

  const res = await Promise.all(promises);
  const stats: CDDVaultStats = {
    name: vaultName,
    projects: res[0].data?.length ? res[0].data?.map((it: DG.Project) => it.name).join(', ') : '',
    molecules: res[1].data?.count ?? undefined,
    protocols: res[2].data?.count ?? undefined,
    batches: res[3].data?.count ?? undefined,
    collections: res[4].data?.count ?? undefined,
  }
  return JSON.stringify(stats);
}

//name: Get Vaults
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//output: string vaults
export async function getVaults(): Promise<string> {
  const vaults = await queryVaults();
  if (vaults.error)
    throw vaults.error;
  if(!vaults?.data?.length)
    throw `No vaults found`;
  return JSON.stringify(vaults.data);
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
  const molecules = await queryMolecules(vaultId, params);
  if (molecules.error) {
    throw molecules.error;
  }
  const df = await createMoleculesDfFromObjects(vaultId, molecules.data?.objects as any[]);
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
  const exportId = getExportId(exportResponse);
  const df = await getAsyncResultsAsDf(vaultId, exportId, timeoutMinutes, false);
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
  const exportResponse = await queryProtocolsAsync(vaultId, {});
  const exportId = getExportId(exportResponse);
  const protocols = await getAsyncResults(vaultId, exportId, timeoutMinutes, false);
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
  const exportId = getExportId(exportResponse);
  const collections = await getAsyncResults(vaultId, exportId, timeoutMinutes, false);
  return collections?.data?.objects ? JSON.stringify(collections.data.objects) : '';
}

//name: Get Saved Searches
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//output: string result
export async function getSavedSearches(vaultId: number): Promise<string> {
  const savedSearches = await querySavedSearches(vaultId);
  if (savedSearches.error)
    throw savedSearches.error;
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
  const exportId = getExportId(exportResponse);
  const res = await getAsyncResultsAsDf(vaultId, exportId, timeoutMinutes, true);
  reorderColummns(res);
  return res;
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
    const exportId = getExportId(exportResponse);
    const readoutRowsRes = await getAsyncResults(vaultId, exportId, 5, false);
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
  const moleculesExportId = getExportId(moleculesExportResponse);
  const cddMols = await getAsyncResults(vaultId, moleculesExportId, 5, false) as ApiResponse<MoleculesQueryResult>;
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