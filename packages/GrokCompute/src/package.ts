/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {toDart} from "datagrok-api/dg";


export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}


let container: DG.DockerContainer;

async function getContainer() {
  if (!container)
    container = await grok.dapi.docker.dockerContainers.filter(_package.name).first();
  if (container.status !== 'started' && container.status !== 'checking')
    await grok.dapi.docker.dockerContainers.run(container.id, true);
  return container;
}

//name: getAllModelingEngines
//description: Gets all registered modeling engines with parameters
//output: map models
export async function getAllModelingEngines(): Promise<object> {
  const container = await getContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/modeling/engines');
  if (response.status !== 200)
    throw new Error(response.statusText);
  return toDart(await response.json());
}

//name: trainModel
//description: Train model
//input: string id
//input: string type
//input: string tableServerUrl
//input: string tableToken
//input: string predict
//input: map parameterValues
//output: blob result
export async function trainModel(id: string, type: string, tableServerUrl: string, tableToken: string,
                                 predict: string, parameterValues: {[_: string]: any}): Promise<Uint8Array> {
  const container = await getContainer();
  const uriParams = new URLSearchParams({
    'id': id,
    'type': type,
    'table_server_url': tableServerUrl,
    'table_token': tableToken,
    'predict': predict
  });
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id,
    '/modeling/train?' + uriParams, {method: 'POST', body: JSON.stringify(parameterValues),
      headers: {'Content-Type': 'application/json'}});
  if (response.status !== 200)
    throw new Error(response.statusText);
  return new Uint8Array(await response.arrayBuffer());
}

//name: applyModel
//description: Apply model
//input: string id
//input: string type
//input: blob modelBlob
//input: string tableServerUrl
//input: string tableToken
//output: list result
export async function applyModel(id: string, type: string, modelBlob: Uint8Array, tableServerUrl: string,
                                 tableToken: string): Promise<DG.Column[]> {
  const container = await getContainer();
  const uriParams = new URLSearchParams({
    'id': id,
    'type': type,
    'table_server_url': tableServerUrl,
    'table_token': tableToken
  });
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id,
    '/modeling/predict?' + uriParams, {method: 'POST', body: modelBlob,
      headers: {'Content-Type': 'application/octet-stream'}});
  if (response.status !== 200)
    throw new Error(response.statusText);
  const data = await response.json();
  return [DG.Column.fromStrings('outcome', Array.from(data['outcome'], (v: any, _) => v?.toString()))]
}
