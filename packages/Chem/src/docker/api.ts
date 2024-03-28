import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {_package} from "../package";

let chemContainer: DG.DockerContainer;

let descriptorsCached: object[];

async function getContainer() {
  if (!chemContainer)
    chemContainer = await grok.dapi.docker.dockerContainers.filter(_package.name).first();
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

export async function calculateDescriptors(table: DG.DataFrame, molecules: DG.Column, descriptors: string[]): Promise<void> {
  const container = await getContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/chem/descriptors',
    {method: 'POST', body: JSON.stringify({'molecules': molecules.toList(), 'descriptors': descriptors})});
  if (response.status !== 200)
    throw new Error(response.statusText);
  const result: object = await response.json();
  for (const [key, value] of Object.entries(result)) {
    const column = DG.Column.fromList(value['type'], table.columns.getUnusedName(key), value['value']);
    table.columns.add(column);
  }
}

export async function smilesToCanonical(smiles: string[]): Promise<string[]> {
  return await convert(smiles, '/chem/molecules_to_canonical');
}

export async function smilesToInchi(smiles: string[]): Promise<string[]> {
  return await convert(smiles, '/chem/molecules_to_inchi');
}

export async function smilesToInchiKey(smiles: string[]): Promise<string[]> {
  return await convert(smiles, '/chem/molecules_to_inchi_key');
}

export async function inchiToSmiles(inchi: string[]): Promise<string[]> {
  return await convert(inchi, '/chem/inchi_to_smiles');
}

export async function inchiToInchiKey(inchi: string[]): Promise<string[]> {
  return await convert(inchi, '/chem/inchi_to_inchi_key');
}

async function convert(molecules: string[], path: string): Promise<string[]> {
  const container = await getContainer();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, path,
    {method: 'POST', body: JSON.stringify(molecules)});
  if (response.status !== 200)
    throw new Error(response.statusText);
  return await response.json();
}
