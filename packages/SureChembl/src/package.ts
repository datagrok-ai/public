/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/surechembl.css';
import {SearchType} from './constants';
export * from './package.g';
export const _package = new DG.Package();

const WIDTH = 200;
const HEIGHT = 100;

const defaultMolLimit = 10;
const defaultSimilarityThreshold = 0.6;
const MAX_SMILES_LENGTH = 5000;

export function init() {
}


type PatentInfo = {
  title: string,
  id: string,
  language: string,
  assign_applic: string,
  published: string,
  smiles: string,
  fields: string
}


let _rdkitPromise: Promise<any> | null = null;
function getRdkitModule(): Promise<any> {
  if (!_rdkitPromise)
    _rdkitPromise = grok.functions.call('Chem:getRdKitModule');
  return _rdkitPromise;
}


function patentSearch(molecule: string, searchType: SearchType): DG.Widget {
  const isSimilarity = searchType === SearchType.similarity;
  const fnName = isSimilarity ?
    PackageFunctions.sureChemblSimilaritySearch.name :
    PackageFunctions.sureChemblSubstructureSearch.name;
  let currentSearch = '';

  const root = ui.divV([]);
  const widget = new DG.Widget(root);

  const compsHost = ui.div();
  compsHost.style.overflowY = 'scroll';

  let similarityThreshold: DG.InputBase | null = null;

  const runSearch = (limit: number, threshold: number) => {
    const token = `${molecule}|${limit}|${threshold}`;
    if (currentSearch === token)
      return;
    currentSearch = token;

    ui.empty(compsHost);
    compsHost.append(ui.loader());

    const params: {[key: string]: any} = {molecule: molecule, limit: limit};
    if (isSimilarity)
      params.similarityThreshold = threshold;

    grok.functions.call(`${_package.name}:${fnName}`, params)
      .then((table: DG.DataFrame | null) => {
        if (token === currentSearch)
          updateSearchPanel(table, compsHost);
      })
      .catch((err: any) => {
        if (token !== currentSearch)
          return;
        if (compsHost.children.length > 0)
          compsHost.removeChild(compsHost.firstChild!);
        const div = ui.divText('Error');
        grok.shell.error(err);
        ui.tooltip.bind(div, `${err}`);
        compsHost.appendChild(div);
      });
  };

  const molLimit = ui.input.int('Molecules limit', {value: defaultMolLimit});
  widget.subs.push(DG.debounce(molLimit.onChanged, 1000).subscribe(() => {
    runSearch(molLimit.value ?? defaultMolLimit,
      similarityThreshold?.value ?? defaultSimilarityThreshold);
  }));

  root.append(molLimit.root);

  if (isSimilarity) {
    similarityThreshold = ui.input.float('Similarity cutoff',
      {value: defaultSimilarityThreshold, min: 0, max: 1, step: 0.05, showSlider: true});
    widget.subs.push(DG.debounce(similarityThreshold.onChanged, 1000).subscribe(() => {
      runSearch(molLimit.value ?? defaultMolLimit,
        similarityThreshold!.value ?? defaultSimilarityThreshold);
    }));
    root.append(similarityThreshold.root);
  }

  root.append(compsHost);

  runSearch(molLimit.value ?? defaultMolLimit,
    similarityThreshold?.value ?? defaultSimilarityThreshold);
  return widget;
}

function updateSearchPanel(table: DG.DataFrame | null, compsHost: HTMLDivElement) {
  ui.empty(compsHost);
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

  const smilesCol = table.col('smiles')!;
  const titleCol = table.col('title')!;
  const docIdCol = table.col('doc_id')!;
  const langCol = table.col('language')!;
  const assignCol = table.col('assign_applic')!;
  const publishedCol = table.col('published')!;
  const fieldsCol = table.col('fields')!;
  const similarityCol = table.col('similarity');

  const patentsByMol: {[key: string]: PatentInfo[]} = {};
  const similarities: {[key: string]: number} = {};

  for (let i = 0; i < table.rowCount; i++) {
    const smiles: string = smilesCol.get(i);
    if (!patentsByMol[smiles])
      patentsByMol[smiles] = [];
    const patentId = docIdCol.get(i);
    const patentInfo: PatentInfo = {
      smiles: smiles,
      title: titleCol.get(i),
      id: `[${patentId}](${`https://www.surechembl.org/patent/${patentId}`})`,
      language: langCol.get(i),
      assign_applic: assignCol.get(i),
      published: publishedCol.get(i).toString(),
      fields: fieldsCol.get(i),
    };
    patentsByMol[smiles].push(patentInfo);
    if (similarityCol && !similarities[smiles])
      similarities[smiles] = similarityCol.get(i);
  }

  const addAllPatentsToWorkspace = createAddToWorkspaceButton(Object.values(patentsByMol).flat(),
    `All patents`, 'Add all found patents to workspace', 'surechembl-add-all-patents-to-workspace-button');
  compsHost.append(ui.divH([addAllPatentsToWorkspace]));

  Object.keys(patentsByMol).forEach((key: string) => {
    const molHost = ui.div();
    const res = grok.chem.drawMolecule(key, WIDTH, HEIGHT, false);
    res.style.float = 'left';
    molHost.append(res);
    const acc = ui.accordion();
    acc.root.style.paddingLeft = '25px';
    const accPane = acc.addPane(`Patents found: ${patentsByMol[key].length}`, () => {
      const df = DG.DataFrame.fromObjects(patentsByMol[key])!;
      df.columns.remove('smiles');
      return df.plot.grid().root;
    });

    const addToWorkspaceButton = createAddToWorkspaceButton(patentsByMol[key], `Patents for ${key}`,
      'Add table to workspace', 'surechembl-add-patents-to-workspace-button');

    const accPaneHeader = accPane.root.getElementsByClassName('d4-accordion-pane-header');
    if (accPaneHeader.length)
      accPaneHeader[0].append(addToWorkspaceButton);
    const molDiv = ui.divV([molHost], {style: {paddingBottom: '15px'}});
    if (similarityCol)
      molDiv.prepend(ui.divText(similarities[key].toFixed(4).toString()));
    molDiv.append(acc.root);
    compsHost.append(molDiv);
  });
}

async function patentDataSearch(molecule: string, limit: number, searchType: SearchType,
  threshold?: number): Promise<DG.DataFrame | null> {
  try {
    if (!grok.chem.isMolBlock(molecule) && molecule?.length > MAX_SMILES_LENGTH)
      throw new Error(`SMILES string longer than ${MAX_SMILES_LENGTH} characters not supported`);
    const rdkit = await getRdkitModule();
    const mol = rdkit.get_mol(molecule);
    let pattern: string;
    try {
      pattern = searchType === SearchType.substructure ? mol.get_smarts() : mol.get_smiles();
    } finally {
      mol?.delete();
    }
    const queryName = searchType === SearchType.substructure ? 'searchPatentBySubstructure' : 'searchPatentBySimilarity';
    const params: {[key: string]: any} = {pattern: pattern, maxMols: limit};
    if (searchType === SearchType.similarity)
      params.threshold = threshold ?? defaultSimilarityThreshold;
    return await grok.data.query(`${_package.name}:${queryName}`, params);
  } catch (e: any) {
    console.error(`In ${searchType} search: ${e.toString()}`);
    throw e;
  }
}


function buildSearchWidget(molecule: string, searchType: SearchType): DG.Widget {
  if (!molecule)
    return new DG.Widget(ui.divText('SMILES is empty'));
  return patentSearch(molecule, searchType);
}


export class PackageFunctions {
  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 * * *',
    },
  })
  static async sureChemblSubstructureSearch(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) molecule: string,
    @grok.decorators.param({'type': 'int'}) limit: number): Promise<DG.DataFrame | null> {
    return patentDataSearch(molecule, limit, SearchType.substructure);
  }


  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 * * *',
    },
  })
  static async sureChemblSimilaritySearch(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) molecule: string,
    @grok.decorators.param({'type': 'int'}) limit: number,
      similarityThreshold?: number): Promise<DG.DataFrame | null> {
    return patentDataSearch(molecule, limit, SearchType.similarity, similarityThreshold);
  }


  @grok.decorators.func({
    'name': 'Databases | SureChEMBL | Substructure Search',
    'condition': 'true',
    'meta': {role: 'panel,widgets'},
  })
  static sureChemblSubstructureSearchWidget(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) molecule: string): DG.Widget {
    return buildSearchWidget(molecule, SearchType.substructure);
  }


  @grok.decorators.func({
    'name': 'Databases | SureChEMBL | Similarity Search',
    'condition': 'true',
    'meta': {role: 'panel,widgets'},
  })
  static sureChemblSimilaritySearchWidget(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) molecule: string): DG.Widget {
    return buildSearchWidget(molecule, SearchType.similarity);
  }
}
