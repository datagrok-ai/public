import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {anyObject, paramsType, pubChemIdType, pubChemSearchType, urlParamsFromObject} from './utils';

const pubChemBaseURL = 'https://pubchem.ncbi.nlm.nih.gov';
const pubChemRest = `${pubChemBaseURL}/rest`;
const pubChemPug = `${pubChemRest}/pug`;

export function renderInfoName(info: anyObject): string {
  return info['Description'] || info['Name'];
}

export function renderInfoValue(info: anyObject, refs: anyObject): HTMLElement {
  let text = info['StringValue'] ?? info['NumValue']?.toString() ?? info['DateValue']?.toString() ??
    info['BoolValue']?.toString() ?? null;

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
  if (info['StringValue'] && info['URL'])
    result = ui.div(text.replaceAll(`<img src="/', '<img src="${pubChemBaseURL}`));
  else if (text && info['URL'])
    result = ui.link(text, info['URL']);
  else if (info['StringValueList'])
    result = ui.list(info['StringValueList']);
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
    const table = ui.table(information, (info, _) => [renderInfoName(info), renderInfoValue(info, refs)]);
    content.append(table);
  }

  const sections: anyObject[] = section['Section'];
  if (sections) {
    const acc = ui.accordion(`pubChem/${section['TOCHeading']}`);
    for (const section of sections) {
      const paneName = section['TOCHeading'];
      const pane = acc.addPane(paneName, () => renderSection(section, refs));
      if (section['Description']) {
        //@ts-ignore: api types are wrong
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
  const response = await fetch(url);
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
  idType: string, id: pubChemIdType, params?: paramsType): Promise<DG.DataFrame | null> {
  params ??= {};
  const listId = await _asyncSearchId('similarity', idType, id, params);
  let json: object[];
  let maxRequests = 10;
  do {
    maxRequests--;
    json = await _getListById(listId);
  } while (typeof json === 'undefined' && maxRequests > 0);

  const df = DG.DataFrame.fromObjects(json);
  return df ?? null;
}

export async function identitySearch(
  idType: string, id: pubChemIdType, params?: anyObject): Promise<DG.DataFrame | null> {
  params ??= {};
  const listId = await _asyncSearchId('identity', idType, id, params);
  const json = await _getListById(listId, [], {});
  const df = DG.DataFrame.fromObjects(json);
  return df ?? null;
}

export async function substructureSearch(
  idType: string, id: pubChemIdType, params?: anyObject): Promise<DG.DataFrame | null> {
  params ??= {};
  const listId = await _asyncSearchId('substructure', idType, id, params);
  let json: object[];
  let maxRequests = 10;
  do {
    maxRequests--;
    json = await _getListById(listId);
  } while (typeof json === 'undefined' && maxRequests > 0);

  const df = DG.DataFrame.fromObjects(json);
  return df ?? null;
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
  const response = await fetch(url);
  const json = await response.json();
  return json;
}

export async function _getListById(
  listId: string, propertyList: string[] = ['CanonicalSMILES'], params?: paramsType): Promise<anyObject[]> {
  params ??= {};
  const properties = propertyList.length ? `/property/${propertyList.join(',')}` : '';
  const url =
    `${pubChemPug}/compound/listkey/${listId}${properties}/JSON?${urlParamsFromObject(params)}`;
  const response = await fetch(url);
  const json = await response.json();
  return json.PropertyTable?.Properties ?? json['PC_Compounds'] ?? await _getListById(listId, propertyList, params);
}

export async function _asyncSearchId(
  searchType: pubChemSearchType, idType: string, id: pubChemIdType, params?: paramsType): Promise<string> {
  params ??= {};
  params['MaxRecords'] ??= 20;
  const url =
    `${pubChemPug}/compound/${searchType}/${idType}/${encodeURIComponent(id)}/JSON?${urlParamsFromObject(params)}`;
  const response = await fetch(url);
  const json: anyObject = await response.json();
  return json['Waiting']['ListKey'];
}
