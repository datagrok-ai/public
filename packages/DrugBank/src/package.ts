/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as OCL from 'openchemlib/full';

export const _package = new DG.Package();

type searchResult = DG.DataFrame | null;

let dbdf: DG.DataFrame | null = null;

//tags: init
export async function initDrugBank() {
  dbdf = DG.DataFrame.fromCsv(await _package.files.readAsText('db.csv'));
}

export async function drugBankSubstructureSearch(mol: string, substructLibrary: boolean): Promise<searchResult> {
  const bitset = await grok.chem.substructureSearch(
    dbdf!.getCol('molecule'), mol, {'substructLibrary': substructLibrary});
  if (bitset === null) {
    return null;
  }

  dbdf!.filter.copyFrom(bitset);
  return dbdf!;
}

export async function dbSimilaritySearch(molecule: string, limit: number, cutoff: number): Promise<searchResult> {
  let searchdf = await grok.chem.findSimilar(dbdf!.getCol('molecule'), molecule, {'limit': limit, 'cutoff': cutoff});
  if (searchdf == null)
    return null;

  if (dbdf!.col('index') === null)
    (dbdf!.columns as DG.ColumnList).addNewInt('index').init((i) => i);
  return dbdf!.join(searchdf, ['index'], ['index'], ['molecule'], ['molecule'], DG.JOIN_TYPE.INNER, true);
}

function getTooltip(value: string) {
  let props = {
    'DRUGBANK_ID': value,
  };
  return ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]);
}

export async function dbSearchWidget(mol: string, searchType: 'similarity' | 'substructure'): Promise<DG.Widget> {
  const headerHost = ui.divH([]);
  const compsHost = ui.divH([]);
  const panel = ui.divV([compsHost]);

  let table: DG.DataFrame | null;
  switch (searchType) {
    case 'similarity':
      table = await dbSimilaritySearch(mol, 20, 0);
      break;
    case 'substructure':
      table = await drugBankSubstructureSearch(mol, false);
      break;
    default:
      throw `DrugBankSearch: No such search type ${searchType}`;
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
    const molfile = table.get('molecule', piv);
    const molecule = OCL.Molecule.fromMolfile(molfile);
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


//name: DrugBank Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function DrugBankSubstructureSearchPanel(mol: string) {
  return await dbSearchWidget(mol, 'substructure');
}

//name: DrugBank Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function DrugBankSimilaritySearchPanel(mol: string) {
  return await dbSearchWidget(mol, 'similarity');
}
