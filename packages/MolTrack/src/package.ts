/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from '@datagrok-libraries/utils/src/u2';
import { MolTrackDockerService } from './utils/moltrack-docker-service';
import { RegistrationView } from './utils/registration-tab';
import { createPath, registerAllData, registerAssayData, updateAllMolTrackSchemas } from './utils/utils';
import { RegistrationCompoundView } from './utils/registration-compound-tab';
import { RegistrationBatchView } from './utils/registration-batch-tab';
import { Scope } from './utils/constants';

export const _package = new DG.Package();

//tags: init
export function init(): void {
  setTimeout(async () => {
    // Run independent tasks in parallel
    await Promise.all([
      updateAllMolTrackSchemas(),
      registerAssayData(),
      registerAllData(),
    ]);
  }, 0);

  // This will be used for the updated docker setup later.
  // const connection = await grok.dapi.connections.filter('name = "moltrack"').first();
  // const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
  // for (const query of queries)
  //   await (query.prepare()).call();
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

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/cdd-icon-big.png',
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
}: {
  title: string;
  smiles: string;
  pathQueryParam: string;
  propsList: any[];
  batchSection: boolean
}) {
  const registrationView = new RegistrationCompoundView(false);
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
  registrationView.view.path = `${createPath('Compound')}?${pathQueryParam}`;

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
  });
}

function findEntityIndex(df: DG.DataFrame, colName: string, searchValue: string): number {
  for (let i = 0; i < df!.rowCount; i++) {
    if (df!.get(colName, i)?.includes(searchValue))
      return i;
  }
  throw new Error(`Entity with value '${searchValue}' not found in column '${colName}'`);
}

//input: dynamic treeNode
//input: view browseView
export async function molTrackAppTreeBrowser(appNode: DG.TreeViewGroup, browseView: any) {
  //search node
  const registerCompoundNode = appNode.getOrCreateGroup('Register').item('Compound');
  registerCompoundNode.onSelected.subscribe(() => {
    const registrationCompoundView = new RegistrationCompoundView();
    registrationCompoundView.show();
    registrationCompoundView.view.path = createPath('Compound');
  });
  const registerBatchNode = appNode.getOrCreateGroup('Register').item('Batch');
  registerBatchNode.onSelected.subscribe(() => {
    const registrationBatchView = new RegistrationBatchView();
    registrationBatchView.show();
  });
  const registerBulkNode = appNode.getOrCreateGroup('Register').item('Bulk...');
  registerBulkNode.onSelected.subscribe(() => {
    const registrationView = new RegistrationView();
    registrationView.show();
  });

  for (const scope of Object.values(Scope)) {
    const formattedScope = scope
      .toLowerCase()
      .replace(/_/g, ' ')
      .replace(/\b\w/g, (char) => char.toUpperCase());
    const retrieveNode = appNode.getOrCreateGroup('Retrieve').item(formattedScope);
    retrieveNode.onSelected.subscribe(async () => {
      const data = await grok.functions.call('MolTrack:retrieveEntity', { scope });
      grok.shell.addTablePreview(data);
    });
  }
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

//name: fetchProperties
//output: string result
export async function fetchProperties(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.fetchProperties();
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

function flattened(item: any) {
  const row: any = {};
  for (const [key, value] of Object.entries(item)) {
    row[key] = (typeof value === 'object' && value !== null) ?
      JSON.stringify(value) :
      value;
  }
  return row;
}
