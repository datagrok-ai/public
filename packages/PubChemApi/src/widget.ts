import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {anyObject, pubChemSearchType} from './utils';
import {identitySearch, similaritySearch, substructureSearch} from './pubchem';

export async function getSearchWidget(molString: string, searchType: pubChemSearchType): Promise<DG.Widget> {
  const headerHost = ui.divH([]);
  const compsHost = ui.divH([ui.loader()]);
  const widget = new DG.Widget(compsHost);

  let json: anyObject[] | null;
  switch (searchType) {
  case 'similarity':
    json = await similaritySearch('smiles', molString);
    break;
  case 'substructure':
    json = await substructureSearch('smiles', molString);
    break;
  case 'identity':
    json = await identitySearch('smiles', molString);
    break;
  default:
    throw new Error(`DrugBankSearch: Search type \`${searchType}\` not found`);
  }

  if (json == null) {
    compsHost.firstChild?.remove();
    compsHost.appendChild(ui.divText('No matches'));
    return widget;
  }

  if (searchType == 'identity') {
    const props: {value: anyObject, urn: anyObject}[] = json[0]['props'];
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

  const table = DG.DataFrame.fromObjects(json);

  if (!table || table.filter.trueCount === 0) {
    compsHost.firstChild?.remove();
    compsHost.appendChild(ui.divText('No matches'));
    return widget;
  }

  const smilesCol = table.col('CanonicalSMILES');
  if (smilesCol !== null) {
    smilesCol.semType = 'Molecule';
    smilesCol.setTag('cell.renderer', 'Molecule');
    if (searchType === 'substructure') {
      smilesCol.temp['chem-scaffold-filter'] = molString;
      smilesCol.temp['chem-scaffold'] = molString;
    }
  }

  const grid = table.plot.grid();
  grid.columns.setOrder(['CanonicalSMILES', 'CID']);
  compsHost.appendChild(grid.root);
  headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
    table!.name = 'PubChem Similarity Search';
    grok.shell.addTableView(table!);
  }, 'Open compounds as table'));
  compsHost.style.overflowY = 'auto';
  grid.root.style.width = 'auto';
  compsHost.firstChild?.remove();

  if (compsHost.children.length === 0) {
    compsHost.appendChild(ui.divText('No matches'));
  }

  return widget;
}
