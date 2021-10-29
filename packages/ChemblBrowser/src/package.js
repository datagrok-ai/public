/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

//name: ChemblSubstructureSearch
//input: string mol {semType: Molecule}
//output: dataframe dbResult
export async function
chemblSubstructureSearch(mol) {
    try {
        let df = await grok.data.query('Chemblapi:SubstructureSmile', {'smile': mol})
        df = df.clone(null, ['canonical_smiles', 'molecule_chembl_id']);
        if (df == null) {
            return DG.DataFrame.create();
        }
        return df;
    } catch (e) {
        console.error("In SubstructureSearch: " + e.toString());
        throw e;
    }
}

//name: ChemblSimilaritySearch
//input: string molecule {semType: Molecule}
//output: dataframe drbResult
export async function chemblSimilaritySearch(molecule) {
    try {
        let df = await grok.data.query('Chemblapi:SimilaritySmileScore', {'smile': molecule, 'score': 40})
        df = df.clone(null, ['canonical_smiles', 'molecule_chembl_id']);
        if (df == null) {
            return DG.DataFrame.create();
        }
        return df;

    } catch (e) {
        console.error(e);
        return null;
    }

}

//name: Chembl Search Widget
//tags: widgets
//input: string mol {semType: Molecule}
//input: string searchType
//output: widget result
export function ChemblSearchWidget(mol, searchType) {
    let headerHost = ui.divH([]);
    let compsHost = ui.divH([ui.loader(), headerHost]);
    let panel = ui.divV([compsHost]);
    let search = {
        'similarity': async () => chemblSimilaritySearch(mol),
        'substructure': async () => chemblSubstructureSearch(mol)
    }
    if (!searchType in search) {

        throw "DrugBankSearch: No such search type" + searchType;
    }

    search[searchType]().then(t => {
          compsHost.removeChild(compsHost.firstChild);
          if (t == null) {

              compsHost.appendChild(ui.divText('No matches'));
              return;
          }
          t.col('canonical_smiles').semType = 'Molecule';
          t.col('canonical_smiles').setTag('cell.renderer', 'Molecule');


          let grid = t.plot.grid();
          let col = grid.columns.byName('molecule_chembl_id');
          col.cellType = 'html';
          grid.onCellPrepare(function (gc) {
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
    .catch(err => {
        if (compsHost.children.length > 0) {
            compsHost.removeChild(compsHost.firstChild);
        }
        let div = ui.divText('No matches');
        ui.tooltip.bind(div, `${err}`);
        compsHost.appendChild(div);
    });
    return new DG.Widget(panel);
}

//name Chembl Get by Id
//input string id
//output dataframe
export async function getById(id) {
    if (!id.toLowerCase().startsWith("chembl")) {
        id = "CHEMBL" + id
    }
    try {
        return await grok.data.query(`${packageName}:MoleculeJson`, {'molecule_chembl_id__exact': id});
    } catch (e) {
        console.error(e);
        return null;
    }

}

//name: Chembl Substructure Search Widget
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function ChemblSubstructureSearchPanel(mol) {
    return ChemblSearchWidget(mol, 'substructure');
}

//name: Chembl Similarity Search Widget
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function ChemblSimilaritySearchPanel(mol) {
    return ChemblSearchWidget(mol, 'similarity');
}
