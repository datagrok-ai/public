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
export function sureChemblSubstructureSearchWidget(molecule: string): DG.Widget {
  return molecule ? patentSearch(molecule, sureChemblSubstructureSearch) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Databases | SureChembl | Similarity Search
//tags: panel, widgets
//input: string molecule {semType: Molecule}
//output: widget result
//condition: true
export function sureChemblSimilaritySearchWidget(molecule: string): DG.Widget {
  return molecule ? patentSearch(molecule, sureChemblSimilaritySearch) : new DG.Widget(ui.divText('SMILES is empty'));
}


type PatentInfo = {
  title: string,
  id: string,
  language: string,
  assign_applic: string,
  published: string,
}


function patentSearch(molecule: string, searchFunction: (molecule: string) => Promise<DG.DataFrame | null>): DG.Widget {
  const compsHost = ui.div([ui.loader()]);
  compsHost.style.overflowY = 'scroll';

  searchFunction(molecule).then((table: DG.DataFrame | null) => {
    updateSearchPanel(table, compsHost);
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

function updateSearchPanel(table: DG.DataFrame | null, compsHost: HTMLDivElement) {
  compsHost.removeChild(compsHost.firstChild!);
  if (table === null || table.rowCount === 0) {
    compsHost.appendChild(ui.divText('No matches'));
    return;
  }

  const numPatentsByMol: {[key: string]: PatentInfo[]} = {};
  const similarities: {[key: string]: number} = {};
  const isSimilarity = table.col('similarity');

  for (let i =0; i < table.rowCount; i++) {
    const smiles: string = table.col('smiles')?.get(i);
    if (!numPatentsByMol[smiles])
      numPatentsByMol[smiles] = [];
    const patentId = table.col('doc_id')?.get(i);
    const patentInfo = {
      title: table.col('title')?.get(i),
      id: `[${patentId}](${`https://www.surechembl.org/patent/${patentId}`})`,
      language: table.col('language')?.get(i),
      assign_applic: table.col('assign_applic')?.get(i),
      published: table.col('published')?.get(i).toString(),
    };
    numPatentsByMol[smiles].push(patentInfo);
    if (isSimilarity) {
      if (!similarities[smiles])
        similarities[smiles] = table.col('similarity')?.get(i);
    }
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
    const molDiv = ui.divV([molHost], {style: {paddingBottom: '15px'}});
    if (isSimilarity)
      molDiv.prepend(ui.divText(similarities[key].toFixed(4).toString()));
    molDiv.append(acc.root);
    compsHost.append(molDiv);
  });
}


export async function sureChemblSubstructureSearch(molecule: string): Promise<DG.DataFrame | null> {
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

export async function sureChemblSimilaritySearch(molecule: string): Promise<DG.DataFrame | null> {
  try {
    const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
    const smiles = mol.get_smiles();
    mol?.delete();
    const df: DG.DataFrame | null =
      await grok.data.query(`${_package.name}:searchPatentBySimilarity`, {'pattern': smiles, 'maxMols': 10});
    return df;
  } catch (e: any) {
    console.error('In SimilaritySearch: ' + e.toString());
    throw e;
  }
}
