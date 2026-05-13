import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {anyObject, getSmiles, paramsType, pubChemIdType, pubChemSearchType, urlParamsFromObject} from './utils';
import {COLUMN_NAMES, pubChemPug, pubChemRest} from './constants';

// PubChem's current PUG REST spec only supports the async `fast*` search endpoints;
// the bare `similarity` / `identity` / `substructure` endpoints return 500 intermittently.
const SEARCH_ENDPOINTS: {[K in pubChemSearchType]: string} = {
  similarity: 'fastsimilarity_2d',
  identity: 'fastidentity',
  substructure: 'fastsubstructure',
};

// PC_Compound shape returned by fast* endpoints → flat widget shape ({CID, CanonicalSMILES})
// widget.ts reads CanonicalSMILES first, falls back to ConnectivitySMILES
function flattenPcCompound(compound: anyObject): anyObject {
  const cid = compound.id?.id?.cid;
  const props: anyObject[] = compound.props ?? [];
  const absolute = props.find((p) => p.urn?.label === 'SMILES' && p.urn?.name === 'Absolute');
  const connectivity = props.find((p) => p.urn?.label === 'SMILES' && p.urn?.name === 'Connectivity');
  return {
    [COLUMN_NAMES.CID]: cid,
    [COLUMN_NAMES.CANONICAL_SMILES]: absolute?.value?.sval ?? connectivity?.value?.sval,
    [COLUMN_NAMES.CONNECTIVITY_SMILES]: connectivity?.value?.sval,
  };
}

export async function search(
  searchType: pubChemSearchType, idType: string, id: pubChemIdType,
  params?: paramsType): Promise<anyObject[] | null> {
  params = {MaxRecords: 20, ...params};
  const endpoint = SEARCH_ENDPOINTS[searchType];
  const url =
    `${pubChemPug}/compound/${endpoint}/${idType}/${encodeURIComponent(id)}/JSON?${urlParamsFromObject(params)}`;
  const response = await grok.dapi.fetchProxy(url);
  const json: anyObject = await response.json();

  // fast* endpoints may respond as: Waiting.ListKey (async), {PC_Compounds: [...]}
  // (fastsimilarity_2d), or a bare PC_Compound[] array (fastidentity).
  let records: anyObject[] | null;
  if (json.Waiting?.ListKey) {
    // async path: identity wants full records (empty propertyList), others want CanonicalSMILES
    const propertyList = searchType === 'identity' ? [] : ['CanonicalSMILES'];
    const polled = await getListById(json.Waiting.ListKey, propertyList);
    records = Array.isArray(polled) ? polled : null;
  } else if (Array.isArray(json))
    records = json;
  else if (Array.isArray(json.PC_Compounds))
    records = json.PC_Compounds;
  else
    records = null;

  if (!records)
    return null;

  // identity widget reads PC_Compounds.props directly; similarity/substructure expect flat rows
  if (searchType === 'identity')
    return records;
  return records.map((r) => r.props ? flattenPcCompound(r) : r);
}

export async function similaritySearch(
  idType: string, id: pubChemIdType, params?: paramsType): Promise<anyObject[] | null> {
  return search('similarity', idType, id, params);
}

export async function identitySearch(
  idType: string, id: pubChemIdType, params?: paramsType): Promise<anyObject[] | null> {
  return search('identity', idType, id, params);
}

export async function substructureSearch(
  idType: string, id: pubChemIdType, params?: paramsType): Promise<anyObject[] | null> {
  return search('substructure', idType, id, params);
}

export async function smilesToPubChem(molString: string): Promise<number | null> {
  try {
    const smiles = await getSmiles(molString);
    const s = await getBy('smiles', 'cids', smiles);
    return s.IdentifierList?.CID?.[0] ?? null;
  } catch (e) {
    return null;
  }
}

export async function getBy(
  idType: string, idTypeReturn: string, id: pubChemIdType, params?: anyObject): Promise<anyObject> {
  params ??= {};
  const url =
    `${pubChemPug}/compound/${idType}/${encodeURIComponent(id)}/${idTypeReturn}/JSON?${urlParamsFromObject(params)}`;
  const response = await grok.dapi.fetchProxy(url);
  const json = await response.json();
  return json;
}

const LIST_POLL_MAX_REQUESTS = 30;
const LIST_POLL_DELAY_MS = 500;

async function getListById(
  listId: string, propertyList: string[] = ['CanonicalSMILES'], params?: paramsType): Promise<anyObject[] | null> {
  params ??= {};
  const properties = propertyList.length ? `/property/${propertyList.join(',')}` : '';
  const url =
    `${pubChemPug}/compound/listkey/${listId}${properties}/JSON?${urlParamsFromObject(params)}`;
  let json: anyObject;
  let maxRequests = LIST_POLL_MAX_REQUESTS;
  do {
    await DG.delay(LIST_POLL_DELAY_MS);
    maxRequests--;
    const response = await grok.dapi.fetchProxy(url);
    json = await response.json();
  } while ('Waiting' in json && maxRequests > 0);

  return json.PropertyTable?.Properties ?? json.PC_Compounds ?? null;
}

export async function getCompoundInfo(id: pubChemIdType): Promise<anyObject> {
  const url = `${pubChemRest}/pug_view/data/compound/${id}/JSON`;
  const response = await grok.dapi.fetchProxy(url);
  return await response.json();
}
