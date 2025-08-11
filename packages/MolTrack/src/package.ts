/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from '@datagrok-libraries/utils/src/u2';
import { MolTrackDockerService } from './utils/moltrack-docker-service';
import { RegistrationView } from './utils/registration-tab';
import { registerAllData, registerAssayData, updateAllMolTrackSchemas } from './utils/utils';
import { RegistrationCompoundView } from './utils/registration-compound-tab';
import { RegistrationBatchView } from './utils/registration-batch-tab';
import { RegistrationAssayView } from './utils/registration-assay-tab';
import { Scope } from './utils/constants';
// import { RegistrationAssayRunView } from './utils/registration-assay-run-tab';
// import { RegistrationAssayResultView } from './utils/registration-assay-results-tab';

export const _package = new DG.Package();

//tags: init
export async function init(): Promise<void> {
  await updateAllMolTrackSchemas();
  await registerAssayData();
  await registerAllData();

  // This will be used for the updated docker setup later.
  // const connection = await grok.dapi.connections.filter('name = "moltrack"').first();
  // const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
  // for (const query of queries)
  //   await (query.prepare()).call();
}

//tags: app
//name: MolTrack
//meta.icon: images/moltrack.png
//output: view v
//meta.browsePath: Chem
export async function molTrackApp(): Promise<DG.ViewBase> {
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

//input: dynamic treeNode
//input: view browseView
export async function molTrackAppTreeBrowser(appNode: DG.TreeViewGroup, browseView: any) {
  //search node
  const registerCompoundNode = appNode.getOrCreateGroup('Register').item('Compound');
  registerCompoundNode.onSelected.subscribe(() => {
    const registrationCompoundView = new RegistrationCompoundView();
    registrationCompoundView.show();
  });
  const registerBatchNode = appNode.getOrCreateGroup('Register').item('Batch');
  registerBatchNode.onSelected.subscribe(() => {
    const registrationBatchView = new RegistrationBatchView();
    registrationBatchView.show();
  });
  const registerAssayNode = appNode.getOrCreateGroup('Register').item('Assay');
  registerAssayNode.onSelected.subscribe(() => {
    const registrationAssayView = new RegistrationAssayView();
    registrationAssayView.show();
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

  const flattened = resultJson.map((item: any) => {
    const row: any = {};
    for (const [key, value] of Object.entries(item)) {
      // If value is an object or array, stringify it
      if (typeof value === 'object' && value !== null) {
        row[key] = JSON.stringify(value);
      } else {
        row[key] = value;
      }
    }
    return row;
  });
  
  return DG.DataFrame.fromObjects(flattened);
}