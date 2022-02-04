import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {pubChemSearchType} from './utils';
import {identitySearch, similaritySearch, substructureSearch} from './pubchem';

export async function pubChemSearchWidget(molString: string, searchType: pubChemSearchType): Promise<DG.Widget> {
  const headerHost = ui.divH([]);
  const compsHost = ui.divH([ui.loader()]);
  const widget = new DG.Widget(compsHost);

  let table: DG.DataFrame | null;
  switch (searchType) {
  case 'similarity':
    table = await similaritySearch('smiles', molString);
    break;
  case 'substructure':
    table = await substructureSearch('smiles', molString);
    break;
  case 'identity':
    table = await identitySearch('smiles', molString);
    break;
  default:
    throw new Error(`DrugBankSearch: Search type \`${searchType}\` not found`);
  }

  compsHost.firstChild?.remove();
  if (table === null || table.filter.trueCount === 0) {
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

  if (compsHost.children.length === 0)
    compsHost.appendChild(ui.divText('No matches'));

  return widget;
}
