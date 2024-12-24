/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

const WIDTH = 200;
const HEIGHT = 100;

export function init() {
}


//name: Databases | SureChembl | Substructure Search
//tags: panel, widgets
//input: string molecule {semType: Molecule}
//output: widget result
//condition: true
export function SureChembl(molecule: string): DG.Widget {
  return molecule ? patentSearchBySubstructure(molecule) : new DG.Widget(ui.divText('SMILES is empty'));
}


type PatentInfo = {
  title: string, 
  id: string,
  language: string, 
  assign_applic: string,
  published: string
}


function patentSearchBySubstructure(molecule: string): DG.Widget {

  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid chem-search-panel-wrapper');
  compsHost.style.overflowY = 'scroll';
  const panel = ui.divV([compsHost]);

  SureChemblSubstructureSearch(molecule).then((table: DG.DataFrame | null) => {

    compsHost.removeChild(compsHost.firstChild!);
    if (table === null || table.rowCount === 0) {
      compsHost.appendChild(ui.divText('No matches'));
      return;
    }

    const numPatentsByMol: {[key: string]: PatentInfo[]} = {};

    for (let i =0; i < table.rowCount; i++) {
      const smiles: string = table.col('smiles')?.get(i);
      if (!numPatentsByMol[smiles]) {
        numPatentsByMol[smiles] = [];
      }
      const patentInfo = {
        title: table.col('title')?.get(i),
        id: table.col('doc_id')?.get(i),
        language: table.col('language')?.get(i),
        assign_applic: table.col('assign_applic')?.get(i),
        published: table.col('published')?.get(i),
      };
      numPatentsByMol[smiles].push(patentInfo);
    }

    Object.keys(numPatentsByMol).forEach((key: string) => {
      const molHost = ui.div();
      grok.functions.call('Chem:drawMolecule', {'molStr': key, 'w': WIDTH, 'h': HEIGHT, 'popupMenu': true})
        .then((res: HTMLElement) => {
          molHost.append(res);
        });
      const acc = ui.accordion();
      acc.root.style.paddingLeft = '25px';
      acc.addPane(`Patents found: ${numPatentsByMol[key].length}`, () => {
        const patentsPane = ui.divV([]);
        numPatentsByMol[key].forEach((patent) => {
          patentsPane.append(ui.tableFromMap(patent));
          patentsPane.append(ui.link('Go to patent page', () => window.open(`https://www.surechembl.org/patent/${patent.id}`, '_blank')))
        });
        return patentsPane;
      });
      compsHost.append(ui.divV([molHost, acc.root]));
    })

  }).catch((err: any) => {
    if (compsHost.children.length > 0)
      compsHost.removeChild(compsHost.firstChild!);

    const div = ui.divText('Error');
    grok.shell.error(err);
    ui.tooltip.bind(div, `${err}`);
    compsHost.appendChild(div);
  });
  return new DG.Widget(panel);
}


export async function SureChemblSubstructureSearch(molecule: string): Promise<DG.DataFrame | null> {
  try {
    const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
    const smarts = mol.get_smarts();
    mol?.delete();
    const df: DG.DataFrame | null =
      await grok.data.query(`${_package.name}:searchPatentBySubstructure`, {'pattern': smarts, 'maxMols': 10});
    return df;
  } catch (e: any) {
    console.error('In SubstructureSearch: ' + e.toString());
    throw e;
  }
}

export async function SureChemblSimilaritySearch(molecule: string): Promise<DG.DataFrame | null> {
  try {
    const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
    const smiles = mol.get_smiles();
    mol?.delete();
    const df: DG.DataFrame | null =
      await grok.data.query(`${_package.name}:patternSimilaritySearch`, {'pattern': smiles, 'maxRows': 100});
    return df;
  } catch (e: any) {
    console.error('In SimilaritySearch: ' + e.toString());
    throw e;
  }
}