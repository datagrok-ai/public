import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {_package} from "../package";
import { gzip } from "pako";

let chemContainer: DG.DockerContainer;

let descriptorsCached: object[];

async function getContainer() {
  if (!chemContainer)
    chemContainer = await grok.dapi.docker.dockerContainers.filter('name = "chem-chem"').first();
  if (chemContainer.status !== 'started' && chemContainer.status !== 'checking')
    await grok.dapi.docker.dockerContainers.run(chemContainer.id, true);
  return chemContainer;
}

export async function getDescriptorsTree(): Promise<object> {
  if (!descriptorsCached) {
    const container = await getContainer();
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/chem/descriptors/tree');
    if (response.status !== 200)
      throw new Error(response.statusText);
    descriptorsCached = await response.json();
  }
  return descriptorsCached;
}

export async function calculateDescriptors(molecules: DG.Column, descriptors: string[]): Promise<DG.Column[]> {
  const result: object = await gzipPostRequest({'molecules': molecules.toList(), 'descriptors': descriptors},
    '/chem/descriptors');
  const colArray = Array<DG.Column>(Object.entries(result).length);
  for (const [key, value] of Object.entries(result)) {
    const column = DG.Column.fromList(value['type'], key, value['value']);
    colArray.push(column);
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
    throw new Error(response.statusText);
  return await response.json();
}
