import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const WIDTH = 200;
const HEIGHT = 100;
const BASE_URL = 'https://www.ebi.ac.uk/chembl/api/data';

export enum SEARCH_TYPE {
  SUBSTRUCTURE = 'substructure',
  SIMILARITY = 'similarity',
}

export enum ELEMENTS {
  CHEMBL_ID = 'molecule_chembl_id',
  SMILES = 'canonical_smiles',
  PROPERTIES = 'molecule_properties',
  SIMILARITY = 'similarity',
}

export async function getData(searchType: SEARCH_TYPE, smiles: string, score: number | null = null):
  Promise<DG.DataFrame | null> {
  const response = await grok.dapi.fetchProxy(
    `${BASE_URL}/${searchType}/${smiles}${searchType === SEARCH_TYPE.SUBSTRUCTURE ? '' : `/${score}`}`);

  const responseText = await response.text();
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(responseText, 'text/xml');
  const molecules = xmlDoc.getElementsByTagName('molecule');
  const rowCount = Math.min(molecules.length, 20);

  const df = DG.DataFrame.create(rowCount);
  let validCount = rowCount;
  for (let i = 0; i < rowCount; i++) {
    const molecule = molecules[i];

    const chemblId = molecule.getElementsByTagName(ELEMENTS.CHEMBL_ID)[0];
    if (!chemblId) {
      validCount = i;
      break;
    }
    let col = df.columns.getOrCreate(ELEMENTS.CHEMBL_ID, DG.TYPE.STRING);
    col.set(i, chemblId.textContent);

    const smiles = molecule.getElementsByTagName(ELEMENTS.SMILES)[0];
    col = df.columns.getOrCreate(ELEMENTS.SMILES, DG.TYPE.STRING);
    col.set(i, smiles.textContent);

    if (searchType === SEARCH_TYPE.SIMILARITY) {
      const similarity = molecule.getElementsByTagName(ELEMENTS.SIMILARITY)[0];
      col = df.columns.getOrCreate(ELEMENTS.SIMILARITY, DG.TYPE.FLOAT);
      col.set(i, parseInt(similarity.textContent ?? '0') / 100);
    }
  }
  for (let i = rowCount - 1; i >= validCount; i--)
    df.rows.removeAt(i);
  return df;
}

export async function chemblSubstructureSearch(mol: string): Promise<DG.DataFrame | null> {
  return await getData(SEARCH_TYPE.SUBSTRUCTURE, mol);
}

export async function chemblSimilaritySearch(molecule: string): Promise<DG.DataFrame | null> {
  return await getData(SEARCH_TYPE.SIMILARITY, molecule, 40);
}

export async function getSmiles(molString: string): Promise<string> {
  return await grok.functions.call('Chem:convertMolNotation',
    {molecule: molString, sourceNotation: 'unknown', targetNotation: 'smiles'});
}

export async function chemblSearchWidget(mol: string, substructure: boolean = false): Promise<DG.Widget> {
  try {
    mol = await getSmiles(mol);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule string is malformed'));
  }
  const headerHost = ui.divH([], {style: {position: 'absolute'}});
  const compsHost = ui.div([ui.loader(), headerHost], 'd4-flex-wrap chem-viewer-grid chem-search-panel-wrapper');
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

    const moleculeCol = table.getCol(ELEMENTS.SMILES);
    const chemblIdCol = table.getCol(ELEMENTS.CHEMBL_ID);
    const molCount = Math.min(table.rowCount, 20);

    for (let i = 0; i < molCount; i++) {
      const molHost = ui.divV([]);
      const res = grok.chem.drawMolecule(moleculeCol.get(i), WIDTH, HEIGHT, true);
      molHost.append(res);
      if (!substructure)
        molHost.append(ui.divText(`Score: ${table.getCol(ELEMENTS.SIMILARITY).get(i)?.toFixed(2)}`));

      ui.tooltip.bind(molHost,
        () => ui.divText(`ChEMBL ID: ${chemblIdCol.get(i)}\nClick to open in ChEMBL Database`));
      molHost.addEventListener('click',
        () => window.open(`https://www.ebi.ac.uk/chembl/compound_report_card/${chemblIdCol.get(i)}`, '_blank'));
      compsHost.appendChild(molHost);
    }

    headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
      table.name = substructure ? 'ChEMBL Substructure Search' : 'ChEMBL Similarity Search';
      grok.shell.addTableView(table);
    }, 'Open compounds as table'));
    compsHost.style.overflowY = 'auto';
  })
    .catch((err: any) => {
      if (compsHost.children.length > 0)
        compsHost.removeChild(compsHost.firstChild!);

      const div = ui.divText('No matches');
      ui.tooltip.bind(div, `${err}`);
      compsHost.appendChild(div);
    });
  return new DG.Widget(panel);
}
