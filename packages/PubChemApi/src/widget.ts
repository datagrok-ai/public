import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {anyObject, getSmiles, pubChemIdType, pubChemSearchType} from './utils';
import {getCompoundInfo, identitySearch, similaritySearch, substructureSearch} from './pubchem';
import {COLUMN_NAMES} from './constants';

const WIDTH = 200;
const HEIGHT = 100;

const pubChemCompoundUrl = (cid: pubChemIdType): string =>
  `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`;

// Walk the pug_view Section tree by TOCHeading (e.g. ['Names and Identifiers', 'Molecular Formula']).
function sectionAt(record: anyObject, path: string[]): anyObject | undefined {
  let current: anyObject | undefined = record;
  for (const heading of path)
    current = current?.Section?.find((s: anyObject) => s.TOCHeading === heading);
  return current;
}

// Read the first stringifiable Value under a section path. Handles StringWithMarkup and Number[+Unit].
function readString(record: anyObject, path: string[]): string | null {
  const value = sectionAt(record, path)?.Information?.[0]?.Value;
  const str = value?.StringWithMarkup?.[0]?.String ?? value?.Number?.[0]?.toString();
  if (str == null)
    return null;
  return value.Unit ? `${str} ${value.Unit}` : str;
}

function extractHazardIcons(record: anyObject): HTMLElement[] {
  const markups = sectionAt(record, ['Primary Hazards'])
    ?.Information?.[0]?.Value?.StringWithMarkup?.[0]?.Markup ?? [];
  return markups
    .filter((m: anyObject) => m.Type === 'Icon' && m.URL)
    .map((m: anyObject) => {
      const img = ui.image(m.URL, 28, 28);
      if (m.Extra)
        ui.tooltip.bind(img, m.Extra);
      return img;
    });
}

// ---------- Top-level panel builder ----------

export async function buildInfoPanel(id: pubChemIdType | null | undefined): Promise<HTMLElement> {
  if (id == null || id === '0')
    return ui.div('Not found in PubChem');

  const json = await getCompoundInfo(id);
  const record = json?.Record;
  if (!record)
    return ui.div('Not found in PubChem');

  const cid = record.RecordNumber;

  const map: anyObject = {};
  if (record.RecordTitle)
    map['Name'] = record.RecordTitle;
  map['CID'] = String(cid);
  const addRow = (label: string, path: string[]) => {
    const val = readString(record, path);
    if (val)
      map[label] = val;
  };
  addRow('Formula', ['Names and Identifiers', 'Molecular Formula']);
  addRow('MW', ['Chemical and Physical Properties', 'Computed Properties', 'Molecular Weight']);
  addRow('CAS', ['Names and Identifiers', 'Other Identifiers', 'CAS']);
  addRow('IUPAC Name', ['Names and Identifiers', 'Computed Descriptors', 'IUPAC Name']);
  addRow('SMILES', ['Names and Identifiers', 'Computed Descriptors', 'SMILES']);
  addRow('InChIKey', ['Names and Identifiers', 'Computed Descriptors', 'InChIKey']);

  const root = ui.divV([]);
  const icons = extractHazardIcons(record);
  if (icons.length > 0) {
    const iconsRow = ui.divH(icons);
    iconsRow.style.marginBottom = '8px';
    root.append(iconsRow);
  }
  root.append(ui.tableFromMap(map));
  const footer = ui.link('View full record on PubChem ↗', pubChemCompoundUrl(cid));
  footer.style.display = 'inline-block';
  footer.style.marginTop = '8px';
  root.append(footer);
  return root;
}

export async function getSearchWidget(molString: string, searchType: pubChemSearchType): Promise<DG.Widget> {
  try {
    molString = await getSmiles(molString);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule string is malformed'));
  }
  const headerHost = ui.div();
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid chem-search-panel-wrapper');
  const widget = new DG.Widget(ui.divV([headerHost, compsHost]));
  const showNoMatches = () => {
    compsHost.firstChild?.remove();
    compsHost.appendChild(ui.divText('No matches'));
    return widget;
  };

  let moleculesJson: anyObject[] | null;
  switch (searchType) {
  case 'similarity':
    moleculesJson = await similaritySearch('smiles', molString);
    break;
  case 'substructure':
    moleculesJson = await substructureSearch('smiles', molString);
    break;
  case 'identity':
    moleculesJson = await identitySearch('smiles', molString);
    break;
  default:
    throw new Error(`PubChemSearch: Search type \`${searchType}\` not found`);
  }

  if (moleculesJson === null || moleculesJson.length === 0)
    return showNoMatches();

  if (searchType === 'identity') {
    const props: {value: anyObject, urn: anyObject}[] | undefined = moleculesJson[0]?.['props'];
    if (!props || props.length === 0)
      return showNoMatches();
    const result: anyObject = {};
    const bannedKeys = ['label', 'name', 'implementation', 'datatype'];

    for (const prop of props) {
      const urn = prop.urn;
      const label: string = urn.label;
      const name: string | undefined = urn.name;
      const value = Object.values(prop.value)[0];
      const boxedValue = ui.divText(value);

      for (const bannedKey of bannedKeys)
        delete urn[bannedKey];

      ui.tooltip.bind(boxedValue, () => ui.tableFromMap(urn));
      result[`${name ? name + ' ' : ''}${label}`] = boxedValue;
    }

    const resultMap = ui.tableFromMap(result);

    return new DG.Widget(resultMap);
  }
  const resultDf = DG.DataFrame.fromObjects(moleculesJson);
  if (!resultDf)
    throw new Error('Failed to build a dataframe from PubChem results');

  let similarStructures: DG.DataFrame | null = null;
  const moleculesCol = resultDf.col(COLUMN_NAMES.CANONICAL_SMILES) ?? resultDf.col(COLUMN_NAMES.CONNECTIVITY_SMILES);
  if (!moleculesCol)
    throw new Error(`Resulting table doesn't contain ${COLUMN_NAMES.CANONICAL_SMILES} or ${COLUMN_NAMES.CONNECTIVITY_SMILES} column`);
  let indexes = new Int32Array(0);
  let scoreCol: DG.Column<number> | null = null;
  let rowCount = resultDf.rowCount;
  if (searchType === 'similarity') {
    similarStructures = await grok.chem.findSimilar(moleculesCol, molString, {limit: 20, cutoff: 0.75});

    if (similarStructures === null || similarStructures.rowCount === 0)
      return showNoMatches();
    scoreCol = similarStructures.getCol(COLUMN_NAMES.SCORE);
    indexes = similarStructures.getCol(COLUMN_NAMES.INDEX).getRawData() as Int32Array;
    rowCount = similarStructures.rowCount;
  }

  const cidCol = resultDf.getCol(COLUMN_NAMES.CID);

  for (let idx = 0; idx < rowCount; idx++) {
    const piv = searchType === 'similarity' ? indexes[idx] : idx;
    const molHost = ui.divV([]);
    const res = grok.chem.drawMolecule(moleculesCol.get(piv), WIDTH, HEIGHT, true);
    molHost.append(res);
    if (searchType === 'similarity' && scoreCol !== null)
      molHost.append(ui.divText(`Score: ${scoreCol.get(idx)?.toFixed(2)}`));

    ui.tooltip.bind(molHost, () => ui.divText(`CID: ${cidCol.get(piv)}\nClick to open in PubChem`));
    molHost.addEventListener('click',
      () => window.open(pubChemCompoundUrl(cidCol.get(piv)), '_blank'));
    compsHost.appendChild(molHost);
  }

  headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
    resultDf.name = `PubChem ${searchType[0].toUpperCase() + searchType.slice(1)} Search`;
    grok.shell.addTableView(resultDf);
  }, 'Open compounds as table'));
  compsHost.style.overflowY = 'auto';
  if (compsHost.parentElement)
    compsHost.parentElement.style.width = 'auto';

  compsHost.firstChild?.remove();

  return widget;
}
