import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

const WIDTH = 200;
const HEIGHT = 100;

export async function chemblSubstructureSearch(mol: string): Promise<DG.DataFrame | null> {
  try {
    let df: DG.DataFrame = await grok.data.query(`${_package.name}:SubstructureSmile`, {'smile': mol});
    if (df.rowCount === 0)
      return null;
    df = df.clone(null, ['canonical_smiles', 'molecule_chembl_id']);

    return df;
  } catch (e: any) {
    console.error('In SubstructureSearch: ' + e.toString());
    throw e;
  }
}

export async function chemblSimilaritySearch(molecule: string): Promise<DG.DataFrame | null> {
  try {
    let df = await grok.data.query(`${_package.name}:SimilaritySmileScore`, {'smile': molecule, 'score': 40});
    if (df.rowCount === 0)
      return null;
    df = df.clone(null, ['canonical_smiles', 'molecule_chembl_id']);

    return df;
  } catch (e: any) {
    console.error('In SimilaritySearch: ' + e.toString());
    throw e;
  }
}

//name: ChEMBL Search Widget
//tags: widgets
//input: string mol {semType: Molecule}
//input: string searchType
//output: widget result
export function chemblSearchWidget(mol: string, substructure: boolean = false): DG.Widget {
  const headerHost = ui.divH([]);
  const compsHost = ui.divH([ui.loader(), headerHost]);
  const panel = ui.divV([compsHost]);
  const searchFunc = substructure ?
    async () => chemblSubstructureSearch(mol) :
    async () => chemblSimilaritySearch(mol);

  searchFunc().then((table: DG.DataFrame | null) => {
    compsHost.removeChild(compsHost.firstChild!);
    if (table === null) {
      compsHost.appendChild(ui.divText('No matches'));
      return;
    }

    const moleculeCol = table.getCol('canonical_smiles');
    const chemblIdCol = table.getCol('molecule_chembl_id');

    const molCount = Math.min(table.rowCount, 20);
    const r = window.devicePixelRatio;

    const renderFunctions = DG.Func.find({meta: {chemRendererName: 'RDKit'}});
    if (renderFunctions.length === 0)
      throw new Error('RDKit renderer is not available');

    for (let i = 0; i < molCount; i++) {
      const molHost = ui.canvas(WIDTH, HEIGHT);
      molHost.classList.add('chem-canvas');
      molHost.width = WIDTH * r;
      molHost.height = HEIGHT * r;
      molHost.style.width = (WIDTH).toString() + 'px';
      molHost.style.height = (HEIGHT).toString() + 'px';

      renderFunctions[0].apply().then((rendndererObj) => {
        rendndererObj.render(molHost.getContext('2d')!, 0, 0, WIDTH, HEIGHT, DG.GridCell.fromValue(moleculeCol.get(i)));
      });

      ui.tooltip.bind(molHost,
        () => ui.divText(`ChEMBL ID: ${chemblIdCol.get(i)}\nClick to open in ChEMBL Database`));
      molHost.addEventListener('click',
        () => window.open(`https://www.ebi.ac.uk/chembl/compound_report_card/${chemblIdCol.get(i)}`, '_blank'));
      compsHost.appendChild(molHost);
    }

    headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
      table.name = `"ChEMBL Similarity Search"`;
      grok.shell.addTableView(table);
    }, 'Open compounds as table'));
    compsHost.style.overflowY = 'auto';
  },
  )
    .catch((err: any) => {
      if (compsHost.children.length > 0)
        compsHost.removeChild(compsHost.firstChild!);

      const div = ui.divText('No matches');
      ui.tooltip.bind(div, `${err}`);
      compsHost.appendChild(div);
    });
  return new DG.Widget(panel);
}

//name Chembl Get by Id
//input string id
//output dataframe
export async function getById(id: string): Promise<DG.DataFrame | null> {
  if (!id.toLowerCase().startsWith('chembl'))
    id = `CHEMBL${id}`;

  try {
    return await grok.data.query(`${_package.name}:MoleculeJson`, {'molecule_chembl_id__exact': id});
  } catch (e: any) {
    console.error(e);
    return null;
  }
}

//name: Databases | ChEMBL | Substructure Search API
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSubstructureSearchPanel(mol: string): DG.Widget {
  return mol ? chemblSearchWidget(mol, true) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Databases | ChEMBL | Similarity Search API
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSimilaritySearchPanel(mol: string): DG.Widget {
  return mol ? chemblSearchWidget(mol) : new DG.Widget(ui.divText('SMILES is empty'));
}
