/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/surechembl.css';
import {downloadPatentDocuments} from './download-patents';
import {SearchType} from './constants';

export const _package = new DG.Package();

const WIDTH = 200;
const HEIGHT = 100;

const defaultMolLimit = 10;
const defaultSimilarityThreshold = 0.6;

export function init() {
}


//name: Databases | SureChEMBL | Substructure Search
//tags: panel, widgets
//input: string molecule {semType: Molecule}
//output: widget result
//condition: true
export function sureChemblSubstructureSearchWidget(molecule: string): DG.Widget {
  return molecule ? patentSearch(molecule, sureChemblSubstructureSearch) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Databases | SureChEMBL | Similarity Search
//tags: panel, widgets
//input: string molecule {semType: Molecule}
//output: widget result
//condition: true
export function sureChemblSimilaritySearchWidget(molecule: string): DG.Widget {
  return molecule ? patentSearch(molecule, sureChemblSimilaritySearch, true) : new DG.Widget(ui.divText('SMILES is empty'));
}


type PatentInfo = {
  title: string,
  id: string,
  language: string,
  assign_applic: string,
  published: string,
  smiles: string,
}


function patentSearch(molecule: string,
  searchFunction: (molecule: string, limit: number, threshold?: number) => Promise<DG.DataFrame | null>,
  isSimilarity?: boolean): DG.Widget {
  const widget = ui.divV([]);

  const compsHost = ui.div();
  compsHost.style.overflowY = 'scroll';

  const molLimit = ui.input.int('Molecules limit', {value: 10});
  DG.debounce(molLimit.onChanged, 1000).subscribe(() => {
    runSearch(molecule, searchFunction, compsHost, molLimit.value ?? defaultMolLimit,
      isSimilarity ? SearchType.similarity: SearchType.substructure, similarityThreshold?.value ?? undefined);
  });

  widget.append(molLimit.root);

  let similarityThreshold: DG.InputBase | null = null;
  if (isSimilarity) {
    similarityThreshold = ui.input.float('Similarity cutoff', {value: defaultSimilarityThreshold, min: 0, max: 1, step: 0.05, showSlider: true});
    DG.debounce(similarityThreshold.onChanged, 1000).subscribe(() => {
      runSearch(molecule, searchFunction, compsHost, molLimit.value ?? defaultMolLimit,
        isSimilarity ? SearchType.similarity: SearchType.substructure, similarityThreshold!.value ?? defaultSimilarityThreshold);
    });
    widget.append(similarityThreshold.root);
  }

  widget.append(compsHost);

  runSearch(molecule, searchFunction, compsHost, molLimit.value ?? defaultMolLimit, isSimilarity ? SearchType.similarity: SearchType.substructure,
    similarityThreshold && similarityThreshold.value ? similarityThreshold.value : defaultSimilarityThreshold );
  return new DG.Widget(widget);
}

//these two variables are required to terminate previous searches if new one has started while previous was running
let currentSimSearch: string = '';
let currentSubSearch: string = '';


function runSearch(molecule: string,
  searchFunction: (molecule: string, limit: number, threshold?: number) => Promise<DG.DataFrame | null>,
  compsHost: HTMLDivElement, limit: number, searchType: SearchType, threshold?: number) {
  const currentSearch = `${molecule}|${limit}|${threshold}`;
  if (searchType === SearchType.substructure) {
    if (currentSubSearch == currentSearch)
      return;
    else
      currentSubSearch = currentSearch;
  } else {
    if (currentSimSearch == currentSearch)
      return;
    else
      currentSimSearch = currentSearch;
  }

  ui.empty(compsHost);
  compsHost.append(ui.loader());

  searchFunction(molecule, limit, threshold)
    .then((table: DG.DataFrame | null) => {
      if ((searchType === SearchType.substructure && currentSearch == currentSubSearch) ||
        (searchType === SearchType.similarity && currentSearch == currentSimSearch))
        updateSearchPanel(table, compsHost);
    })
    .catch((err: any) => {
      if ((searchType === SearchType.substructure && currentSearch !== currentSubSearch) ||
        (searchType === SearchType.similarity && currentSearch !== currentSimSearch))
        return;

      if (compsHost.children.length > 0)
        compsHost.removeChild(compsHost.firstChild!);

      const div = ui.divText('Error');
      grok.shell.error(err);
      ui.tooltip.bind(div, `${err}`);
      compsHost.appendChild(div);
    });
}

function updateSearchPanel(table: DG.DataFrame | null, compsHost: HTMLDivElement) {
  compsHost.removeChild(compsHost.firstChild!);
  if (table === null || table.rowCount === 0) {
    compsHost.appendChild(ui.divText('No matches'));
    return;
  }

  const createAddToWorkspaceButton = (objects: any[], dfName: string, tooltip: string, className: string): HTMLElement => {
    const addToWorkspaceButton = ui.icons.add(async () => {
      const df = DG.DataFrame.fromObjects(objects)!;
      df.name = dfName;
      await grok.data.detectSemanticTypes(df);
      grok.shell.addTableView(df);
    }, tooltip);
    addToWorkspaceButton.classList.add(className);
    return addToWorkspaceButton;
  };

  const downloadPatents = (patentsIds: string[], className: string): HTMLElement => {
    const downloadPatentsButton = ui.iconFA('arrow-to-bottom', () => {
      grok.shell.info('Started downloading patents');
      downloadPatentDocuments(patentsIds).then((res: any) => {
        const infoDiv = ui.div([ui.divText('Patents download finished')]);
        const df = DG.DataFrame.fromObjects(res);
        if (df) {
          df.name = 'Patents download results';
          const openResultsButton = ui.button('Open', () => {
            grok.shell.addTableView(DG.DataFrame.fromObjects(res)!);
          });
          infoDiv.append(openResultsButton);
        }
        grok.shell.info(infoDiv);
      });
    }, 'Download patents');
    downloadPatentsButton.classList.add(className);
    return downloadPatentsButton;
  };

  const numPatentsByMol: {[key: string]: PatentInfo[]} = {};
  const similarities: {[key: string]: number} = {};
  const isSimilarity = table.col('similarity');
  const patentsIdsByMol: {[key: string]: string[]} = {};

  for (let i = 0; i < table.rowCount; i++) {
    const smiles: string = table.col('smiles')?.get(i);
    if (!numPatentsByMol[smiles])
      numPatentsByMol[smiles] = [];
    if (!patentsIdsByMol[smiles])
      patentsIdsByMol[smiles] = [];
    const patentId = table.col('doc_id')?.get(i);
    const patentInfo: PatentInfo = {
      smiles: smiles,
      title: table.col('title')?.get(i),
      id: `[${patentId}](${`https://www.surechembl.org/patent/${patentId}`})`,
      language: table.col('language')?.get(i),
      assign_applic: table.col('assign_applic')?.get(i),
      published: table.col('published')?.get(i).toString(),
    };
    numPatentsByMol[smiles].push(patentInfo);
    patentsIdsByMol[smiles].push(patentId);
    if (isSimilarity) {
      if (!similarities[smiles])
        similarities[smiles] = table.col('similarity')?.get(i);
    }
  }

  const addAllPatentsToWorkspace = createAddToWorkspaceButton(Object.values(numPatentsByMol).reduce((acc, val) => acc.concat(val), []),
    `All patents`, 'Add all found patents to workspace', 'surechembl-add-all-patents-to-workspace-button');
  const downloadAllPatentsButton = downloadPatents(Object.values(patentsIdsByMol).reduce((acc, val) => acc.concat(val), []),
    'surechembl-download-all-patents-to-workspace-button');
  compsHost.append(ui.divH([addAllPatentsToWorkspace, downloadAllPatentsButton]));

  Object.keys(numPatentsByMol).forEach((key: string) => {
    const molHost = ui.div();
    const res = grok.chem.drawMolecule(key, WIDTH, HEIGHT, false);
    res.style.float = 'left';
    molHost.append(res);
    const acc = ui.accordion();
    acc.root.style.paddingLeft = '25px';
    const accPane = acc.addPane(`Patents found: ${numPatentsByMol[key].length}`, () => {
      const df = DG.DataFrame.fromObjects(numPatentsByMol[key])!;
      df.columns.remove('smiles');
      return df.plot.grid().root;
    });

    const addToWorkspaceButton = createAddToWorkspaceButton(numPatentsByMol[key], `Patents for ${key}`,
      'Add table to workspace', 'surechembl-add-patents-to-workspace-button');
    const downloadPatentsButton = downloadPatents(patentsIdsByMol[key], 'surechembl-add-patents-to-workspace-button');

    const accPaneHeader = accPane.root.getElementsByClassName('d4-accordion-pane-header');
    if (accPaneHeader.length) {
      accPaneHeader[0].append(addToWorkspaceButton);
      accPaneHeader[0].append(downloadPatentsButton);
    }
    const molDiv = ui.divV([molHost], {style: {paddingBottom: '15px'}});
    if (isSimilarity)
      molDiv.prepend(ui.divText(similarities[key].toFixed(4).toString()));
    molDiv.append(acc.root);
    compsHost.append(molDiv);
  });
}


//name: sureChemblSubstructureSearch
//input: string molecule {semType: Molecule}
//input: int limit
//output: dataframe df
export async function sureChemblSubstructureSearch(molecule: string, limit: number): Promise<DG.DataFrame | null> {
  try {
    const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
    const smarts = mol.get_smarts();
    mol?.delete();
    const df: DG.DataFrame | null =
      await grok.data.query(`${_package.name}:searchPatentBySubstructure`, {'pattern': smarts, 'maxMols': limit});
    return df;
  } catch (e: any) {
    console.error('In SubstructureSearch: ' + e.toString());
    throw e;
  }
}

//name: sureChemblSimilaritySearch
//input: string molecule {semType: Molecule}
//input: int limit
//input: double similarityThreshold
//output: dataframe df
export async function sureChemblSimilaritySearch(molecule: string, limit: number, similarityThreshold?: number): Promise<DG.DataFrame | null> {
  try {
    const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
    const smiles = mol.get_smiles();
    mol?.delete();
    const df: DG.DataFrame | null =
      await grok.data.query(`${_package.name}:searchPatentBySimilarity`, {'pattern': smiles, 'maxMols': limit, 'threshold': similarityThreshold ?? defaultSimilarityThreshold});
    return df;
  } catch (e: any) {
    console.error('In SimilaritySearch: ' + e.toString());
    throw e;
  }
}

