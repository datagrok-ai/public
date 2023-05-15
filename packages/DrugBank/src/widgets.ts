import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as OCL from 'openchemlib/full';
import {findSimilar, searchSubstructure} from './searches';

const WIDTH = 200;
const HEIGHT = 100;

export async function searchWidget(molString: string, searchType: 'similarity' | 'substructure', dbdf: DG.DataFrame,
): Promise<DG.Widget> {
  const headerHost = ui.div();
  const compsHost = ui.divV([]);
  const panel = ui.divV([headerHost, compsHost]);

  let table: DG.DataFrame | null;
  try {
    switch (searchType) {
    case 'similarity':
      table = await findSimilar(molString, 20, 0, dbdf);
      break;
    case 'substructure':
      table = await searchSubstructure(molString, dbdf);
      break;
    default:
      throw new Error(`DrugBankSearch: No such search type ${searchType}`);
    }
  } catch (e) {
    return new DG.Widget(ui.divText('Error occurred during search. Molecule is possible malformed'));
  }

  compsHost.firstChild?.remove();
  if (table === null || table.filter.trueCount === 0) {
    compsHost.appendChild(ui.divText('No matches'));
    return new DG.Widget(panel);
  }
  table.name = `DrugBank ${searchType === 'similarity' ? 'Similarity' : 'Substructure'} Search`;

  const bitsetIndexes = table.filter.getSelectedIndexes();
  const iterations = Math.min(bitsetIndexes.length, 20);
  const moleculeCol: DG.Column<string> = table.getCol('molecule');
  const idCol: DG.Column<string> = table.getCol('DRUGBANK_ID');
  const nameCol: DG.Column<string> = table.getCol('COMMON_NAME');

  for (let n = 0; n < iterations; n++) {
    const piv = bitsetIndexes[n];
    const molfile = moleculeCol.get(piv)!;

    const molHost = ui.canvas(WIDTH, HEIGHT);
    molHost.classList.add('chem-canvas');
    const r = window.devicePixelRatio;
    molHost.width = WIDTH * r;
    molHost.height = HEIGHT * r;
    molHost.style.width = (WIDTH).toString() + 'px';
    molHost.style.height = (HEIGHT).toString() + 'px';
    const renderFunctions = DG.Func.find({meta: {chemRendererName: 'RDKit'}});
    if (renderFunctions.length > 0) {
    renderFunctions[0].apply().then((rendndererObj) => {
      rendndererObj.render(molHost.getContext('2d')!, 0, 0, WIDTH, HEIGHT,
        DG.GridCell.fromValue(molfile));
      });
    }
  
    ui.tooltip.bind(molHost, () => ui.divText(`Common name: ${nameCol.get(piv)!}\nClick to open in DrugBank Online`));
    molHost.addEventListener('click', () => window.open(`https://go.drugbank.com/drugs/${idCol.get(piv)}`, '_blank'));
    compsHost.appendChild(molHost);
  }
  headerHost.appendChild(
    ui.iconFA('arrow-square-down', () => grok.shell.addTableView(table!), 'Open compounds as table'));
  compsHost.style.overflowY = 'auto';

  return new DG.Widget(panel);
}

export function drugNameMoleculeConvert(id: string, dbdfRowCount: number, synonymsCol: DG.Column<string>,
  smilesCol: DG.Column<string>): string {
  const drugName = id.slice(3).toLowerCase();
  //TODO: consider using raw data instead of .get or column iterator (see ApiSamples)
  //TODO: benchmark it
  for (let i = 0; i < dbdfRowCount; i++) {
    const currentSynonym = synonymsCol.get(i)!.toLowerCase();
    //TODO: check why .includes & consider using hash-map
    if (currentSynonym.includes(drugName))
      return smilesCol.get(i)!;
  }
  return '';
}
