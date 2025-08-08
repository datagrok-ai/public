/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from '@datagrok-libraries/utils/src/u2';
import { MolTrackDockerService } from './utils/moltrack-docker-service';
import { RegistrationView } from './utils/registration-tab';
import { RegistrationSingleView } from './utils/registration-single-tab'
import { RegistrationCompoundView } from './utils/registration-compound-tab';
import { RegistrationBatchView } from './utils/registration-batch-tab';
import { RegistrationAssayView } from './utils/registration-assay-tab';
import { RegistrationAssayRunView } from './utils/registration-assay-run-tab';
import { RegistrationAssayResultView } from './utils/registration-assay-results-tab';

export const _package = new DG.Package();

//tags: init
export async function init(): Promise<void> {
  const connection = await grok.dapi.connections.filter('name = "moltrack"').first();
  const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
  for (const query of queries)
    await (query.prepare()).call();
}

//tags: app
//name: MolTrack App
//meta.icon: images/cdd-icon-big.png
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
    ui.wait(async () => (await grok.functions.call('MolTrack:getCompounds') as DG.DataFrame).plot.grid().root)
  ]));
}

//input: dynamic treeNode
//input: view browseView
export async function molTrackAppTreeBrowser(appNode: DG.TreeViewGroup, browseView: any) {
  //search node
  const registerBulkNode = appNode.getOrCreateGroup("Register").item("Register bulk");
  registerBulkNode.onSelected.subscribe(() => {
    const registrationView = new RegistrationView();
    registrationView.show();
  });
  const registerCompoundNode = appNode.getOrCreateGroup("Register").item("Register compound");
  registerCompoundNode.onSelected.subscribe(() => {
    const registrationCompoundView = new RegistrationCompoundView();
    registrationCompoundView.show();
  });
  const registerBatchNode = appNode.getOrCreateGroup("Register").item("Register batch");
  registerBatchNode.onSelected.subscribe(() => {
    const registrationBatchView = new RegistrationBatchView();
    registrationBatchView.show();
  });
  const registerAssayNode = appNode.getOrCreateGroup("Register").item("Register assay");
  registerAssayNode.onSelected.subscribe(() => {
    const registrationAssayView = new RegistrationAssayView();
    registrationAssayView.show();
  });
  const registerAssayRunNode = appNode.getOrCreateGroup("Register").item("Register assay run");
  registerAssayRunNode.onSelected.subscribe(() => {
    const registrationAssayRunView = new RegistrationAssayRunView();
    registrationAssayRunView.show();
  });
    const registerAssayResultNode = appNode.getOrCreateGroup("Register").item("Register assay result");
  registerAssayResultNode.onSelected.subscribe(() => {
    const registrationAssayResultView = new RegistrationAssayResultView();
    registrationAssayResultView.show();
  });
}


//name: checkMolTrackHealth
//description: Checks whether the MolTrack service is running and responsive
//output: string result
export async function checkMolTrackHealth(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.checkHealth();
}

//name: fetchMolTrackProperties
//description: Retrieves all properties defined for the 'compound' scope
//output: string result
export async function fetchMolTrackProperties(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.fetchProperties();
}

//name: updateMolTrackProperties
//input: string jsonPayload
//description: Registers compound properties in the MolTrack service based on the given JSON data
//output: string result
export async function updateMolTrackSchema(jsonPayload: string): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.updateSchema(jsonPayload);
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
