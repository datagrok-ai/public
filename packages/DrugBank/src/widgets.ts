import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {drugBankSimilaritySearch, drugBankSubstructureSearch} from './searches';

import * as OCL from 'openchemlib/full';
import {drugBankSearchTypes} from './utils';

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

  const bitsetIndexes = table.filter.getSelectedIndexes();
  const iterations = Math.min(bitsetIndexes.length, 20);
  for (let n = 0; n < iterations; n++) {
    const piv = bitsetIndexes[n];
    const molHost = ui.canvas();
    const smiles = table.get('SMILES', piv);
    const molecule = OCL.Molecule.fromSmiles(smiles);
    OCL.StructureView.drawMolecule(molHost, molecule, {'suppressChiralText': true});

    //@ts-ignore: second argument expects string type, but works with callable too
    ui.tooltip.bind(molHost, () => getTooltip(table.get('DRUGBANK_ID', piv)));
    molHost.addEventListener('click', () => {
      window.open(table!.get('link', piv), '_blank');
    });
    compsHost.appendChild(molHost);
  }
  headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
    table!.name = 'DrugBank Similarity Search';
    grok.shell.addTableView(table!);
  }, 'Open compounds as table'));
  compsHost.style.overflowY = 'auto';

  return new DG.Widget(panel);
}
