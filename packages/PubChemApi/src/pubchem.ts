import * as grok from 'datagrok-api/grok';

import {anyObject, getSmiles, paramsType, pubChemIdType, pubChemSearchType, urlParamsFromObject} from './utils';
import {delay} from '@datagrok-libraries/utils/src/test';
import {pubChemPug, pubChemRest} from './tests/const';

export async function similaritySearch(
  idType: string, id: pubChemIdType, params?: paramsType): Promise<anyObject[] | null> {
  params ??= {};
  const listId = await _asyncSearchId('similarity', idType, id, params);
  if (!listId)
    return null;

  const json: anyObject[] | anyObject | null = await _getListById(listId);
  return Array.isArray(json) ? json : null;
}

export async function identitySearch(
  idType: string, id: pubChemIdType, params?: anyObject): Promise<anyObject[] | null> {
  params ??= {};
  const listId = await _asyncSearchId('identity', idType, id, params);
  if (!listId)
    return null;

  const json: anyObject[] | anyObject | null = await _getListById(listId, [], {});
  return Array.isArray(json) ? json : null;
}

export async function substructureSearch(
  idType: string, id: pubChemIdType, params?: anyObject): Promise<anyObject[] | null> {
  params ??= {};
  const listId = await _asyncSearchId('substructure', idType, id, params);
  if (!listId)
    return null;

  const json: anyObject[] | anyObject | null = await _getListById(listId);
  return Array.isArray(json) ? json : null;
}

export async function smilesToPubChem(molString: string): Promise<number | null> {
  try {
    molString = await getSmiles(molString);
  } catch (e) {
    return null;
  }
  const s = await getBy('smiles', 'cids', molString);
  const cids = s.IdentifierList?.CID[0];
  return cids;
}

export async function getBy(
  idType: string, idTypeReturn: string, id: pubChemIdType, params?: anyObject): Promise<anyObject> {
  params ??= {};
  const url = `${pubChemPug}/compound/${idType}/${id}/${idTypeReturn}/JSON?${urlParamsFromObject(params)}`;
  const response = await grok.dapi.fetchProxy(url);
  const json = await response.json();
  return json;
}

export async function _getListById(
  listId: string, propertyList: string[] = ['CanonicalSMILES'], params?: paramsType): Promise<anyObject[] | null> {
  params ??= {};
  const properties = propertyList.length ? `/property/${propertyList.join(',')}` : '';
  const url =
    `${pubChemPug}/compound/listkey/${listId}${properties}/JSON?${urlParamsFromObject(params)}`;
  let json: anyObject;
  let maxRequests = 10;
  do {
    await delay(100);
    maxRequests--;
    const response = await grok.dapi.fetchProxy(url);
    json = await response.json();
  } while (json.hasOwnProperty('Waiting') && maxRequests > 0);

  return json.PropertyTable?.Properties ?? json.PC_Compounds ?? null;
}

export async function _asyncSearchId(
  searchType: pubChemSearchType, idType: string, id: pubChemIdType, params?: paramsType): Promise<string | null> {
  params ??= {};
  params.MaxRecords ??= 20;
  const url =
    `${pubChemPug}/compound/${searchType}/${idType}/${encodeURIComponent(id)}/JSON?${urlParamsFromObject(params)}`;
  const response = await grok.dapi.fetchProxy(url);
  const json: anyObject = await response.json();
  return json.Waiting?.ListKey ?? null;
}

export async function getCompoundInfo(id: pubChemIdType): Promise<anyObject> {
  const url = `${pubChemRest}/pug_view/data/compound/${id}/JSON`;
  const response = await grok.dapi.fetchProxy(url);
  return response.json();
}
