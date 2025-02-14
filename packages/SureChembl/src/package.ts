/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/surechembl.css';

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

  const compsHost = ui.div([ui.loader()]);
  compsHost.style.overflowY = 'scroll';

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
      const patentId = table.col('doc_id')?.get(i);
      const patentInfo = {
        title: table.col('title')?.get(i),
        id: `[${patentId}](${`https://www.surechembl.org/patent/${patentId}`})`,
        language: table.col('language')?.get(i),
        assign_applic: table.col('assign_applic')?.get(i),
        published: table.col('published')?.get(i).toString(),
      };
      numPatentsByMol[smiles].push(patentInfo);
    }

    Object.keys(numPatentsByMol).forEach((key: string) => {
      const molHost = ui.div();
      grok.functions.call('Chem:drawMolecule', {'molStr': key, 'w': WIDTH, 'h': HEIGHT, 'popupMenu': false})
        .then((res: HTMLElement) => {
          res.style.float = 'left';
          molHost.append(res);
        });
      const acc = ui.accordion();
      acc.root.style.paddingLeft = '25px';
      const accPane = acc.addPane(`Patents found: ${numPatentsByMol[key].length}`, () => {
        const df = DG.DataFrame.fromObjects(numPatentsByMol[key])!;
        return df.plot.grid().root;
      });
      const addToWorkspaceButton = ui.icons.add(() => {
        const df = DG.DataFrame.fromObjects(numPatentsByMol[key])!;
        df.name = `Patents for ${key}`;
        grok.shell.addTableView(df);
      }, 'Add table to workspace');
      addToWorkspaceButton.classList.add('surechembl-add-patents-to-workspace-button');
      const accPaneHeader = accPane.root.getElementsByClassName('d4-accordion-pane-header');
      if (accPaneHeader.length)
        accPaneHeader[0].append(addToWorkspaceButton);
      compsHost.append(ui.divV([molHost, acc.root], {style: {paddingBottom: '15px'}}));
    });
  }).catch((err: any) => {
    if (compsHost.children.length > 0)
      compsHost.removeChild(compsHost.firstChild!);

    const div = ui.divText('Error');
    grok.shell.error(err);
    ui.tooltip.bind(div, `${err}`);
    compsHost.appendChild(div);
  });
  return new DG.Widget(compsHost);
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