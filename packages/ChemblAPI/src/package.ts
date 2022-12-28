import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();


export async function chemblSubstructureSearch(mol: string): Promise<DG.DataFrame> {
  try {
    let df = await grok.data.query(`${_package.name}:SubstructureSmile`, {'smile': mol})
    df = df.clone(null, ['canonical_smiles', 'molecule_chembl_id']);
    if (df == null)
      return DG.DataFrame.create();

    return df;
  } catch (e: any) {
    console.error("In SubstructureSearch: " + e.toString());
    throw e;
  }
}

export async function chemblSimilaritySearch(molecule: string): Promise<DG.DataFrame> {
  try {
    let df = await grok.data.query(`${_package.name}:SimilaritySmileScore`, {'smile': molecule, 'score': 40})
    df = df.clone(null, ['canonical_smiles', 'molecule_chembl_id']);
    if (df == null)
      return DG.DataFrame.create();

    return df;
  } catch (e: any) {
    console.error("In SimilaritySearch: " + e.toString());
    throw e;
  }

}

//name: Chembl Search Widget
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

  searchFunc().then((t: any) => {
      compsHost.removeChild(compsHost.firstChild!);
      if (t == null) {

        compsHost.appendChild(ui.divText('No matches'));
        return;
      }
      t.col('canonical_smiles').semType = 'Molecule';
      t.col('canonical_smiles').setTag('cell.renderer', 'Molecule');


      const grid = t.plot.grid();
      const col = grid.columns.byName('molecule_chembl_id');
      col.cellType = 'html';
      grid.onCellPrepare(function (gc: DG.GridCell) {
        if (gc.isTableCell && gc.gridColumn.name === 'molecule_chembl_id') {
          const link = `https://www.ebi.ac.uk/chembl/compound_report_card/${gc.cell.value}/`;
          gc.style.element = ui.divV([
            ui.link(gc.cell.value, link, link)
          ]);
        }
      });
      compsHost.appendChild(grid.root)
      headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
        t.name = `"DrugBank Similarity Search"`;
        grok.shell.addTableView(t);
      }, 'Open compounds as table'));
      compsHost.style.overflowY = 'auto';
    }
  )
  .catch((err: any) => {
    if (compsHost.children.length > 0) {
      compsHost.removeChild(compsHost.firstChild!);
    }
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
  if (!id.toLowerCase().startsWith("chembl")) {
    id = "CHEMBL" + id
  }
  try {
    return await grok.data.query(`${_package.name}:MoleculeJson`, {'molecule_chembl_id__exact': id});
  } catch (e: any) {
    console.error(e);
    return null;
  }
}

//name: Chembl API | Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSubstructureSearchPanel(mol: string): DG.Widget {
  return mol ? chemblSearchWidget(mol, true) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Chembl API | Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSimilaritySearchPanel(mol: string): DG.Widget {
  return mol ? chemblSearchWidget(mol) : new DG.Widget(ui.divText('SMILES is empty'));
}
