import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {_package} from "../package";
import { gzip } from "pako";

let descriptorsCached: object[];

async function getContainer() {
  const container = await grok.dapi.docker.dockerContainers.filter('name = "chem-chem"').first();
  if (!container)
    throw new Error('Chem service is not available: "chem-chem" Docker container not found');
  if (container.status !== 'started' && container.status !== 'checking') {
    try {
      await grok.dapi.docker.dockerContainers.run(container.id, true);
    }
    catch (e) {
      _package.logger.warning(`Could not start the "chem-chem" container: ${e}`);
    }
  }
  return container;
}

export async function getDescriptorsTree(): Promise<object> {
  if (!descriptorsCached) {
    const container = await getContainer();
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/chem/descriptors/tree');
    if (response.status !== 200)
      throw await requestError(response, '/chem/descriptors/tree');
    descriptorsCached = await response.json();
  }
  return descriptorsCached;
}

export async function calculateDescriptors(molecules: DG.Column, descriptors: string[]): Promise<DG.Column[]> {
  const result: object = await gzipPostRequest({'molecules': molecules.toList(), 'descriptors': descriptors},
    '/chem/descriptors');
  const colArray = Array<DG.Column>(Object.entries(result).length);
  let counter = 0;
  for (const [key, value] of Object.entries(result)) {
    const column = DG.Column.fromList(value['type'], key, value['value']);
    colArray[counter] = column;
    counter++;
  }
  return colArray;
}

export async function smilesToCanonical(smiles: string[]): Promise<string[]> {
  return await gzipPostRequest(smiles, '/chem/molecules_to_canonical');
}

export async function smilesToInchi(smiles: string[]): Promise<string[]> {
  return await gzipPostRequest(smiles, '/chem/molecules_to_inchi');
}

export async function smilesToInchiKey(smiles: string[]): Promise<string[]> {
  return await gzipPostRequest(smiles, '/chem/molecules_to_inchi_key');
}

export async function inchiToSmiles(inchi: string[]): Promise<string[]> {
  return await gzipPostRequest(inchi, '/chem/inchi_to_smiles');
}

export async function inchiToInchiKey(inchi: string[]): Promise<string[]> {
  return await gzipPostRequest(inchi, '/chem/inchi_to_inchi_key');
}

async function gzipPostRequest(data: any, path: string): Promise<any> {
  const container = await getContainer();
  data = JSON.stringify(data);
  const len = data.length;
  // use approximate size of data to decide what compression level to use
  let body: Uint8Array = await gzip(data, {'level': len > 100000 && len < 10000000 ? 6 : len <= 100000 ? 1 : 9});
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, path,
    {method: 'POST', body: body, headers: {'Content-Encoding': 'gzip',  'Content-Type': 'application/json',
        'Content-Length': `${body.length}`}});
  if (response.status !== 200)
    throw await requestError(response, path);
  return await response.json();
}

async function requestError(response: Response, path: string): Promise<Error> {
  const body = await response.json().catch(() => null);
  const reason = body?.['datagrok-error'] ?? (response.statusText || '(no statusText)');
  return new Error(`chem-chem container request to ${path} failed: ${response.status} ${reason}`);
}
