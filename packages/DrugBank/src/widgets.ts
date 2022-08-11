import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {drugBankSimilaritySearch, drugBankSubstructureSearch} from './searches';

import * as OCL from 'openchemlib/full';
import {drugBankSearchTypes, getTooltip} from './utils';

export async function drugBankSearchWidget(
  molString: string, searchType: drugBankSearchTypes, dbdf: DG.DataFrame): Promise<DG.Widget> {
  const headerHost = ui.divH([]);
  const compsHost = ui.divH([]);
  const panel = ui.divV([compsHost]);

  let table: DG.DataFrame | null;
  switch (searchType) {
  case 'similarity':
    table = await drugBankSimilaritySearch(molString, 20, 0, dbdf);
    break;
  case 'substructure':
    table = await drugBankSubstructureSearch(molString, false, dbdf);
    break;
  default:
    throw new Error(`DrugBankSearch: No such search type ${searchType}`);
  }

  compsHost.firstChild?.remove();
  // compsHost.removeChild(compsHost.firstChild!);
  if (table === null || table.filter.trueCount === 0) {
    compsHost.appendChild(ui.divText('No matches'));
    return new DG.Widget(panel);
  }
  table.name = `DrugBank ${searchType === 'similarity' ? 'Similarity' : 'Substructure'} Search`;

  const bitsetIndexes = table.filter.getSelectedIndexes();
  const iterations = Math.min(bitsetIndexes.length, 20);
  const smilesCol: DG.Column<string> = table.getCol('SMILES');
  const idCol: DG.Column<string> = table.getCol('DRUGBANK_ID');
  const linkCol: DG.Column<string> = table.getCol('link');
  for (let n = 0; n < iterations; n++) {
    const piv = bitsetIndexes[n];
    const molHost = ui.canvas();
    const smiles = smilesCol.get(piv)!;
    const molecule = OCL.Molecule.fromSmiles(smiles);
    OCL.StructureView.drawMolecule(molHost, molecule, {'suppressChiralText': true});

    ui.tooltip.bind(molHost, () => getTooltip(idCol.get(piv)!));
    molHost.addEventListener('click', () => window.open(linkCol.get(piv)!, '_blank'));
    compsHost.appendChild(molHost);
  }
  headerHost.appendChild(
    ui.iconFA('arrow-square-down', () => grok.shell.addTableView(table!), 'Open compounds as table'));
  compsHost.style.overflowY = 'auto';

  return new DG.Widget(panel);
}

export function drugNameMoleculeWidget(id: string, dbdf: DG.DataFrame): string | null {
  const drugName = id.slice(3);
  for (var i = 0; i < dbdf.rowCount; i++) {
    if (dbdf.get('SYNONYMS', i).toLowerCase().includes(drugName))
      return dbdf.get('Smiles', i);
  }
  return null;
}
