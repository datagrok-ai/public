/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from "@datagrok-libraries/utils/src/u2";
import { createRegistrationNode } from './utils';

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

async function getMoltrackContainer() {
  return await grok.dapi.docker.dockerContainers.filter('moltrack').first();
}

//name: checkMoltrackHealth
//description: Checks whether the Moltrack service is running and responsive
//output: string result
export async function checkMoltrackHealth(): Promise<string> {
  const container = await getMoltrackContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/v1/health', {
    method: 'GET',
  });
  return response.text();
}

//name: fetchMoltrackProperties
//description: Retrieves all properties defined for the 'compound' scope
//output: string result
export async function fetchMoltrackProperties(): Promise<string> {
  const container = await getMoltrackContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/v1/schema/compounds', {
    method: 'GET',
  });
  return response.text();
}

//name: updateMoltrackProperties
//input: string jsonPayload
//description: Registers compound properties in the Moltrack service based on the given JSON data
//output: string result
export async function updateMoltrackSchema(jsonPayload: string): Promise<string> {
  const container = await getMoltrackContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/v1/schema/', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: jsonPayload,
  });
  return response.text();
}
