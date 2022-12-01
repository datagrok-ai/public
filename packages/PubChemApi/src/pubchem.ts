import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {anyObject, paramsType, pubChemIdType, pubChemSearchType, urlParamsFromObject} from './utils';
import {delay} from '@datagrok-libraries/utils/src/test';

const pubChemBaseURL = 'https://pubchem.ncbi.nlm.nih.gov';
const pubChemRest = `${pubChemBaseURL}/rest`;
const pubChemPug = `${pubChemRest}/pug`;

export function renderInfoValue(info: anyObject, refs: anyObject): HTMLElement {
  const infoValue: {[key: string]: any[]} = info['Value'];
  let infoValueList = Object.values(infoValue)[0];
  infoValueList = infoValueList.length ? Object.values(infoValueList[0]) : infoValueList;
  let text: string = infoValueList.length ? Object.values(infoValueList[0]).toString() : '';

  if (text && info['ValueUnit'])
    text += ` ${info['ValueUnit']}`;

  if (info['Table']) {
    const columnNames: string[] = info['Table']['ColumnName'];
    const rows: any[] = info['Table']['Row'] || [];

    return ui.table(rows, (row, _) => {
      row['Cell'].map((x: any) => renderInfoValue(x, refs)).toList();
    }, columnNames);
  }

  let result;
  if (infoValue['String'] && info['URL'])
    result = ui.div(text.replaceAll('<img src="/', `<img src="${pubChemBaseURL}`));
  else if (text && info['URL'])
    result = ui.link(text, info['URL']);
  else if (infoValue['StringList'])
    result = ui.list(infoValue['StringList']);
  else if (text)
    result = ui.divText(text);
  else if (info['ExternalDataMimeType'] === 'image/png' && info['ExternalDataURL'])
    result = ui.image(info['ExternalDataURL'], 200, 150);
  else
    result = null;


  if (result instanceof HTMLElement && info['ReferenceNumber']) {
    //@ts-ignore: api types are wrong
    ui.tooltip.bind(result, () => ui.tableFromMap(refs[info['ReferenceNumber']]));
  }

  return result ?? ui.divText('unknown');
}

export function renderSection(section: anyObject, refs: anyObject): HTMLDivElement {
  const content = ui.divV([]);

  const description = section['Description'];
  if (description)
    content.append(ui.div(description));


  const information: any[] = section['Information'];
  if (information) {
    const table = ui.table(information, (inf, _) => [inf['Description'] || inf['Name'], renderInfoValue(inf, refs)]);
    content.append(table);
  }

  const sections: anyObject[] = section['Section'];
  if (sections) {
    const acc = ui.accordion(`pubChem/${section['TOCHeading']}`);
    for (const section of sections) {
      const paneName = section['TOCHeading'];
      acc.addPane(paneName, () => renderSection(section, refs));
      if (section['Description']) {
        ui.tooltip.bind(ui.divText(paneName), () => ui.div(section['Description']));
      }
    }
    content.append(acc.root);
  }
  return content;
}

export async function init(id: pubChemIdType): Promise<HTMLElement> {
  if (id === '0' || id === null)
    return ui.div('Not found in PubChem.');

  const url = `${pubChemRest}/pug_view/data/compound/${id}/JSON`;
  const response = await grok.dapi.fetchProxy(url);
  const json = await response.json();
  const acc = ui.accordion('pubChem');
  acc.header = ui.label(`pubchem: ${id}`);

  const references: anyObject = {};
  const recordJson = json['Record'];
  for (const ref of recordJson['Reference'])
    references[ref['ReferenceNumber']] = ref;

  const sections = recordJson['Section'];
  for (const section of sections)
    acc.addPane(section['TOCHeading'], () => renderSection(section, references));

  return acc.root;
}


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

export async function smilesToPubChem(smiles: string) {
  const s = await getBy('smiles', 'cids', smiles);
  const cids = s['IdentifierList']['CID'][0];
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

  return json.PropertyTable?.Properties ?? json['PC_Compounds'] ?? null;
}

export async function _asyncSearchId(
  searchType: pubChemSearchType, idType: string, id: pubChemIdType, params?: paramsType): Promise<string | undefined> {
  params ??= {};
  params['MaxRecords'] ??= 20;
  const url =
    `${pubChemPug}/compound/${searchType}/${idType}/${encodeURIComponent(id)}/JSON?${urlParamsFromObject(params)}`;
  const response = await grok.dapi.fetchProxy(url);
  const json: anyObject = await response.json();
  return json.Waiting?.ListKey;
}
