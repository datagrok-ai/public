/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from "@datagrok-libraries/utils/src/u2";
import { createRegistrationNode, getMolTrackContainer } from './utils';
import { scopeToUrl } from './constants';
import '../css/moltrack.css';

export const _package = new DG.Package();

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
      '- Find contextual information on molecules.\n'
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
  const searchNode = appNode.item("Register");
  searchNode.onSelected.subscribe(() => {
    createRegistrationNode(appNode);
  });
  appNode.group('Protocols');
  appNode.group('Plates');
  appNode.group('Assays');
}


//name: checkMolTrackHealth
//description: Checks whether the MolTrack service is running and responsive
//output: string result
export async function checkMolTrackHealth(): Promise<string> {
  const container = await getMolTrackContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/v1/health', {
    method: 'GET',
  });
  return response.text();
}

//name: fetchMolTrackProperties
//description: Retrieves all properties defined for the 'compound' scope
//output: string result
export async function fetchMolTrackProperties(): Promise<string> {
  const container = await getMolTrackContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/v1/schema/compounds', {
    method: 'GET',
  });
  return response.text();
}

//name: updateMolTrackProperties
//input: string jsonPayload
//description: Registers compound properties in the MolTrack service based on the given JSON data
//output: string result
export async function updateMolTrackSchema(jsonPayload: string): Promise<string> {
  const container = await getMolTrackContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/v1/schema/', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: jsonPayload,
  });
  return response.text();
}



//name: registerBulk
//input: file csv_file
//input: string scope
//input: string mapping
//input: string errorHandling
//output: dataframe result
export async function registerBulk(csv_file: DG.FileInfo, scope: string, mapping: string, errorHandling: string): Promise<DG.DataFrame> {
  let resultJson ="";
  const formData = new FormData();
  const content = await csv_file.readAsBytes();
  const blob = new Blob([content], { type: 'text/csv' });
  const file = new File([blob], csv_file.fileName, { type: 'text/csv' });
  formData.append('csv_file', file, csv_file.name);
  formData.append('error_handling', errorHandling || 'reject_row');
  formData.append('mapping', mapping || '');
  formData.append('output_format', 'json');
  try {
    const container = await getMolTrackContainer();
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, scopeToUrl[scope], {
      method: 'POST',
      body: formData
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`MolTrack API error: ${response.status} ${errorText}`);
    }
    resultJson = await response.text();
  } catch (e) {
    grok.shell.error(String(e));
  }
  const json = JSON.stringify(JSON.parse(resultJson)['data']);
  const ret_value = DG.DataFrame.fromJson(json);
  return ret_value;
}