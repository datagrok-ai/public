/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import {MoleculeFieldSearch, getVaults, MoleculeQueryParams, queryMolecules, queryReadoutRows, Molecule, Batch, querySavedSearches, SavedSearch, querySavedSearchById, queryExportStatus, queryExportResult, queryMoleculesAsync, queryReadoutRowsAsync, ApiResponse, MoleculesQueryResult} from "./cdd-vault-api";
import { CDDVaultSearchType } from './constants';
import '../css/cdd-vault.css';
import { SeachEditor } from './search-function-editor';
import { CDD_HOST, createLinksFromIds, getAsyncResults, getAsyncResultsAsDf } from './utils';

export const _package = new DG.Package();

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

  return DG.View.fromRoot(appHeader);
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
    const vaultNode = treeNode.group(vault.name);

    const moleculesNode = vaultNode.item('Molecules');
    moleculesNode.onSelected.subscribe(async (_) => {
      const view = DG.View.create({name: 'Molecules'});
      grok.shell.addPreview(view);
      ui.setUpdateIndicator(view.root, true, 'Waiting for molecules');
      const df = await grok.functions.call('CDDVaultLink:getMoleculesAsync', { vaultId: vault.id, timeoutMinutes: 5});
      view.close();
      df.name = 'Molecules';
      grok.shell.addTablePreview(df);   

    });

    const searchNode = vaultNode.item('Search');
    searchNode.onSelected.subscribe(() => {
      const view = DG.View.create();
      const funcEditor = new SeachEditor(vault.id);
      const acc = funcEditor.getEditor();
      const runButton = ui.bigButton('SEARCH', async () => {
        ui.setUpdateIndicator(gridDiv, true);
        const params = funcEditor.getParams();
        const df = await grok.functions.call('CDDVaultLink:cDDVaultSearchAsync',
          {
            vaultId: vault.id, structure: params.structure, structure_search_type: params.structure_search_type,
            structure_similarity_threshold: params.structure_similarity_threshold, protocol: params.protocol, run: params.run 
          });
        ui.empty(gridDiv);
        gridDiv.append(df.plot.grid().root);
        ui.setUpdateIndicator(gridDiv, false);
      });
      const gridDiv = ui.div('', 'cdd-vault-search-res-div');
      runButton.classList.add('cdd-vault-run-search-button');
      view.name = 'Search CDD Vault'
      view.root.append(ui.divV([
        acc,
        runButton,
        gridDiv
      ], {style: {height: '100%'}}));
      grok.shell.addPreview(view);
    });

    const savedSearchesNode = vaultNode.group('Saved searches', null, false);
    savedSearchesNode.onNodeExpanding.subscribe(async () => {
      const savedSearchesStr =  await grok.functions.call('CDDVaultLink:getSavedSearches', { vaultId: vault.id });
      if (savedSearchesStr !== '') {
        const savedSearchesObj = JSON.parse(savedSearchesStr) as SavedSearch[];
        for (const search of savedSearchesObj) {
          const searchItem = savedSearchesNode.item(search.name);
          searchItem.onSelected.subscribe(async () => {
            const view = DG.View.create();
            const gridDiv = ui.div('', 'cdd-vault-search-res-div');
            ui.setUpdateIndicator(gridDiv, true);
            view.append(gridDiv);
            grok.shell.addPreview(view);
            const df = await grok.functions.call('CDDVaultLink:getSavedSearchResults', { vaultId: vault.id, searchId: search.id, timeoutMinutes: 5});
            gridDiv.append(df.plot.grid().root);
            ui.setUpdateIndicator(gridDiv, false);
          });
        }
      }
    });
   //TODO! unlock other tabs
   // vaultNode.group('Protocols');
   // vaultNode.group('Plates');
   // vaultNode.group('Assays');
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
    const df = DG.DataFrame.fromObjects(molsRes)!;
    if (!df)
      return DG.DataFrame.create();
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
    const df = DG.DataFrame.fromObjects(molsRes)!;
    if (!df)
      return DG.DataFrame.create();
    await grok.data.detectSemanticTypes(df);
    return df;
  }
  return DG.DataFrame.create();
}

//name: Get Molecules
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//output: dataframe df
export async function getMolecules(vaultId: number): Promise<DG.DataFrame> {
  const molecules = await queryMolecules(vaultId, { page_size: 100 }); //TODO! Make async request and remove page size
  if (molecules.error) {
    grok.shell.error(molecules.error);
    return DG.DataFrame.create();
  }
  const df = DG.DataFrame.fromObjects(molecules.data!.objects!)!;
  if (!df)
    return DG.DataFrame.create();
  createLinksFromIds(vaultId, df);
  await grok.data.detectSemanticTypes(df);
  return df;
}


//name: Get Molecules Async
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: int vaultId {nullable: true}
//input: int timeoutMinutes
//output: dataframe df
export async function getMoleculesAsync(vaultId: number, timeoutMinutes: number): Promise<DG.DataFrame> {
  const exportResponse = await queryMoleculesAsync(vaultId, {});
  const df = await getAsyncResultsAsDf(vaultId, exportResponse, timeoutMinutes, false);
  createLinksFromIds(vaultId, df);
  return df;
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
  const df = DG.DataFrame.fromObjects(cddMols.data!.objects!)!;
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
