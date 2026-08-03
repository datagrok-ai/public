import * as grok from 'datagrok-api/grok';

export type pubChemSearchType = 'similarity' | 'substructure' | 'identity';
export type pubChemIdType = number | string;
export type anyObject = {[key: string]: any};
export type paramsType = {[key: string]: pubChemIdType};

export function urlParamsFromObject(obj: {[key: string]: pubChemIdType}): string {
  return new URLSearchParams(
    Object.entries(obj).map(([key, value]) => [key, String(value)])).toString();
}

export async function getSmiles(molString: string): Promise<string> {
  molString = await grok.functions.call('Chem:convertMolNotation',
    {molecule: molString, sourceNotation: 'unknown', targetNotation: 'smiles'});
  return molString;
}
