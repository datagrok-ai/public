/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from '@datagrok-libraries/utils/src/u2';
import { MolTrackDockerService } from './utils/moltrack-docker-service';
import { RegistrationView } from './utils/registration-tab';
import { createPath, registerAllData, registerAssayData, updateAllMolTrackSchemas } from './utils/utils';
import { EntityBaseView } from './utils/registration-entity-base';
import { MolTrackProp, Scope } from './utils/constants';
import { createSearchPanel } from './utils/search';
import { MolTrackEntityType, MolTrackProperty } from './utils/types';

export const _package = new DG.Package();

//tags: init
export async function init(): Promise<void> {
  try {
    await grok.functions.call('Chem:initChemAutostart');
  } catch (e) {}
  // This will be used for the updated docker setup later.
  const connection = await grok.dapi.connections.filter('name = "moltrack"').first();
  const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
  for (const query of queries)
    await (query.prepare()).call();
}

//name: initDB
export async function initDB() {
  await updateAllMolTrackSchemas();
  await registerAssayData();
  await registerAllData();
}

//tags: app
//name: MolTrack
//meta.icon: images/moltrack.png
//input: string path {meta.url: true; optional: true}
//output: view v
//meta.browsePath: Chem
export async function molTrackApp(path: string): Promise<DG.ViewBase> {
  const url = new URL(window.location.href);
  const corporateCompoundId = url.searchParams.get('corporate_compound_id');
  const corporateBatchId = url.searchParams.get('corporate_batch_id');

  if (corporateCompoundId)
    return await compoundView(corporateCompoundId);

  if (corporateBatchId)
    return await batchView(corporateBatchId);

  if (path && path.includes('Compound'))
    return initRegisterView('Compound');

  if (path && path.includes('Batch'))
    return initRegisterView('Batch');

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/moltrack.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/MolTrack/README.md',
    description: '- Chemical compound registration system\n' +
      '- Analyze assay data\n' +
      '- Find contextual information on molecules.\n',
  });

  return DG.View.fromRoot(ui.divV([
    appHeader,
    ui.wait(async () => (await grok.functions.call('MolTrack:getCompounds') as DG.DataFrame).plot.grid().root),
  ]));
}

async function buildRegistrationView({
  title,
  smiles,
  pathQueryParam,
  propsList,
  batchSection,
  path,
}: {
  title: string;
  smiles: string;
  pathQueryParam: string;
  propsList: any[];
  batchSection: boolean,
  path: string,
}) {
  const registrationView = new EntityBaseView(false);
  registrationView.initialSmiles = smiles;
  registrationView.singleRetrieved = true;
  registrationView.title = title;
  registrationView.isBatchSectionExpanded = batchSection;
  await registrationView.buildUIMethod();

  registrationView.formBackingObject = {};
  for (const prop of propsList) {
    const value = prop.value_string ?? prop.value_num ?? prop.value_datetime ?? prop.value_uuid;
    if (value !== null && value !== undefined)
      registrationView.formBackingObject[prop.name] = value;
  }

  for (const input of registrationView.inputs) {
    const propName = input.property.name;
    if (registrationView.formBackingObject.hasOwnProperty(propName))
      input.value = registrationView.formBackingObject[propName];
  }

  registrationView.show();
  registrationView.view.path = `${createPath(path)}?${pathQueryParam}`;

  return registrationView.view;
}

async function compoundView(corporateCompoundId: string) {
  const df = await retrieveEntity(Scope.COMPOUNDS);
  const idx = findEntityIndex(df!, 'properties', corporateCompoundId);
  const smiles = df!.get('canonical_smiles', idx);
  const propertiesJson = df!.get('properties', idx);
  const properties: any[] = JSON.parse(propertiesJson ?? '[]');

  return buildRegistrationView({
    title: `Compound: ${corporateCompoundId}`,
    smiles,
    pathQueryParam: `corporate_compound_id=${encodeURIComponent(corporateCompoundId)}`,
    propsList: properties,
    batchSection: false,
    path: 'Compound',
  });
}

async function batchView(corporateBatchId: string) {
  const df = await retrieveEntity(Scope.BATCHES);
  const idx = findEntityIndex(df!, 'properties', corporateBatchId);

  const compoundId = df!.get('compound_id', idx);
  const compoundInfo = await getCompoundById(compoundId);
  const smiles = compoundInfo?.get('canonical_smiles', 0);

  const batchProps: any[] = JSON.parse(df!.get('properties', idx) ?? '[]');
  const compoundProps: any[] = JSON.parse(compoundInfo?.get('properties', 0) ?? '[]');
  const combinedPropsMap = new Map<string, any>();
  for (const prop of [...compoundProps, ...batchProps])
    combinedPropsMap.set(prop.name, prop);

  return buildRegistrationView({
    title: `Batch: ${corporateBatchId}`,
    smiles,
    pathQueryParam: `corporate_batch_id=${encodeURIComponent(corporateBatchId)}`,
    propsList: Array.from(combinedPropsMap.values()),
    batchSection: true,
    path: 'Batch',
  });
}

function findEntityIndex(df: DG.DataFrame, colName: string, searchValue: string): number {
  for (let i = 0; i < df!.rowCount; i++) {
    if (df!.get(colName, i)?.includes(searchValue))
      return i;
  }
  throw new Error(`Entity with value '${searchValue}' not found in column '${colName}'`);
}

function initRegisterView(entity: 'Compound' | 'Batch') {
  const isBatch = entity === 'Batch';
  const view = new EntityBaseView(!isBatch);

  if (isBatch) {
    view.title = 'Register a new batch';
    view.isBatchSectionExpanded = true;
    view.path = 'Batch';
    view.buildUIMethod();
    view.view.name = 'Register a batch';
  } else {
    view.view.name = 'Register a compound';
    view.view.path = createPath('Compound');
  }

  view.show();
  return view.view;
}

//input: dynamic treeNode
//input: view browseView
export async function molTrackAppTreeBrowser(appNode: DG.TreeViewGroup, browseView: any) {
  function createRegisterNode(label: string, initView: () => void) {
    appNode.getOrCreateGroup('Register').item(label).onSelected.subscribe(initView);
  }

  function createRetrieveNode(scope: string) {
    const formattedScope = scope
      .toLowerCase()
      .replace(/_/g, ' ')
      .replace(/\b\w/g, c => c.toUpperCase());

    appNode.getOrCreateGroup('Search').item(formattedScope).onSelected.subscribe(async () => {
      const data = await grok.functions.call('MolTrack:retrieveEntity', { scope });
      const tv = grok.shell.addTablePreview(data);
      tv.name = formattedScope;
      createSearchPanel(tv, scope as Scope);
    });
  }

  createRegisterNode('Compound', () => initRegisterView('Compound'));
  createRegisterNode('Batch', () => initRegisterView('Batch'));
  createRegisterNode('Bulk...', () => new RegistrationView().show());

  const excludedScopes = [Scope.ASSAY_RUNS, Scope.ASSAY_RESULTS];
  
  Object.values(Scope)
    .filter(scope => !excludedScopes.includes(scope))
    .forEach(scope => createRetrieveNode(scope));
}

//name: checkMolTrackHealth
//description: Checks whether the MolTrack service is running and responsive
//output: string result
export async function checkMolTrackHealth(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.checkHealth();
}

//name: fetchCompoundProperties
//description: Retrieves all properties defined for the 'compound' scope
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
//output: string result
export async function fetchCompoundProperties(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.fetchCompoundProperties();
}

//name: fetchBatchProperties
//description: Retrieves all properties defined for the 'batch' scope
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
//output: string result
export async function fetchBatchProperties(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.fetchBatchProperties();
}

//name: fetchSchema
//description: Retrieves all dynamic fields
//output: string result
export async function fetchSchema(): Promise<string> {
  await MolTrackDockerService.init();
  const res = await MolTrackDockerService.fetchSchema();
  return JSON.stringify(res);
}

//name: fetchDirectSchema
//description: Retrieves all static fields
//output: string result
export async function fetchDirectSchema(): Promise<string> {
  await MolTrackDockerService.init();
  const res = await MolTrackDockerService.fetchDirectSchema();
  return JSON.stringify(res);
}

//name: updateMolTrackProperties
//input: string jsonPayload
//description: Registers compound properties in the MolTrack service based on the given JSON data
//output: string result
export async function updateMolTrackSchema(jsonPayload: string): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.updateSchema(jsonPayload);
}

//name: registerAssays
//input: string assayPayload
//output: string result
export async function registerAssays(assayPayload: string): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.registerAssay(assayPayload);
}

//name: registerBulk
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
//input: file csvFile
//input: string scope
//input: string mapping
//input: string errorHandling
//output: dataframe result
export async function registerBulk(
  csvFile: DG.FileInfo,
  scope: string,
  mapping: string,
  errorHandling: string,
): Promise<DG.DataFrame> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.registerBulk(csvFile, scope, mapping, errorHandling);
}

//name: searchTest
//input: string operator
export async function searchTest(operator: string) {
  await MolTrackDockerService.init();
  const query = {
    'level': 'compounds',
    'output': ['compounds.canonical_smiles'],
    'filter': {
      'field': 'compounds.structure',
      'operator': operator,
      'value': 'c1ccccc1',
    },
    'output_format': 'json',
  };
  return await MolTrackDockerService.search(query, 'compounds');
}

//name: retrieveEntity
//input: string scope
//output: dataframe result
export async function retrieveEntity(scope: string): Promise<DG.DataFrame | undefined> {
  await MolTrackDockerService.init();
  const resultJson = await MolTrackDockerService.retrieveEntity(scope);
  const flattenedRes = resultJson.map((item: any) => flattened(item));
  return DG.DataFrame.fromObjects(flattenedRes);
}

export async function getCompoundById(id: number): Promise<DG.DataFrame | undefined> {
  await MolTrackDockerService.init();
  const resultJson = await MolTrackDockerService.getCompoundById(id);
  const flattenedRes = [resultJson].map((item) => flattened(item));
  return DG.DataFrame.fromObjects(flattenedRes);
}

export async function checkCompoundExists(smiles: string): Promise<boolean> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.checkCompoundExists(smiles);
}

function flattened(item: any) {
  const row: any = {};
  for (const [key, value] of Object.entries(item)) {
    row[key] = (typeof value === 'object' && value !== null) ?
      JSON.stringify(value) :
      value;
  }
  return row;
}


//name: Databases | MolTrack
//input: string mol {semType: Molecule}
//tags: panel
//output: widget res
export async function getMoltrackPropPanel(mol: string): Promise<DG.Widget> {
  // TODO: Implement a server endpoint to fetch information by structure and update accordingly
  const compound = await getCompound(mol);
  if (compound.detail) {
    const registerButton = ui.bigButton('Register', async () => {
      const registrationView = new EntityBaseView(false);
      registrationView.initialSmiles = mol;
      await registrationView.buildUIMethod();
      registrationView.show();
    });
    registerButton.style.width = '80px';
    const div = ui.divV([ui.divText('No matches'), registerButton]);
    div.style.gap = '10px';
    return DG.Widget.fromRoot(div);
  }

  const acc = ui.accordion();

  const createTableFromObject = (obj: Record<string, unknown>, excludeKeys: string[] = []) => {
    const map: Record<string, unknown> = {};
    for (const [key, value] of Object.entries(obj)) {
      if (!excludeKeys.includes(key)) map[key] = value;
    }
    return ui.tableFromMap(map);
  };

  const extractPropValue = (prop: Record<string, any>) =>
    prop.value_uuid ?? prop.value_num ?? prop.value_string ?? prop.value_qualifier ?? prop.value_datetime ?? null;

  const createPropertiesTable = (properties: any[]) => {
    if (!properties || !properties.length) return ui.divText('No properties');
    const map: Record<string, unknown> = {};
    for (const prop of properties) map[prop.name] = extractPropValue(prop);
    return ui.tableFromMap(map);
  };

  const createNestedAccordion = (items: any[], idKey: string) => {
    if (!items || !items.length) return ui.divText('No items');
    const nested = ui.accordion();
    for (const item of items) {
      nested.addPane(item[idKey], () => createTableFromObject(item));
    }
    return nested.root;
  };

  const panes = [
    {
      title: 'Batches',
      content: () => createNestedAccordion(JSON.parse(compound.batches || '[]'), 'id'),
    },
    {
      title: 'Properties',
      content: () => createPropertiesTable(JSON.parse(compound.properties || '[]')),
    },
  ];

   // --- Create the corporate_compound_id link first ---
  const properties = JSON.parse(compound.properties || '[]');
  const corpProp = properties.find((p: any) => p.name === 'corporate_compound_id');
  const value = extractPropValue(corpProp);
  const corpLink = ui.link(value.toString() ?? '', async () => {
    await compoundView(value.toString());
  });

  const compoundContent = ui.div(createTableFromObject(compound, ['properties', 'batches']));

  for (const pane of panes) acc.addPane(pane.title, pane.content);

  return DG.Widget.fromRoot(ui.divV([corpLink, compoundContent, acc.root]));
}

export async function getCompound(mol: string) {
  await MolTrackDockerService.init();
  const query = {
    'level': 'compounds',
    'output': ['compounds.id'],
    'filter': {
      'field': 'compounds.canonical_smiles',
      'operator': 'CONTAINS',
      'value': mol,
      'threshold': 1,
    },
    'output_format': 'json',
  };
  
  // const searchResult = await MolTrackDockerService.search(query, 'compounds'); 
  // const id = searchResult.data[0]["compounds.id"];
  const compoundInfo = await getCompoundById(1);
  const compound = compoundInfo?.toJson()[0];
  return compound;
}