/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {registerChemblIdHandler} from './handlers';
import {chemblBioactivityForTargetsSearch, chemblPKForDrugSearch} from './search-scripts';

export const _package = new DG.Package();

const WIDTH = 200;
const HEIGHT = 100;

//tags: autostart
export function init() {
  //Register handlers
  // DG.ObjectHandler.register(new ChemblIdHandler());
  registerChemblIdHandler(_package);
}

export async function chemblSubstructureSearch(molecule: string): Promise<DG.DataFrame | null> {
  try {
    const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
    const smarts = mol.get_smarts();
    mol?.delete();
    const df: DG.DataFrame | null =
      await grok.data.query(`${_package.name}:patternSubstructureSearch`, {'pattern': smarts, 'maxRows': 100});
    return df;
  } catch (e: any) {
    console.error('In SubstructureSearch: ' + e.toString());
    throw e;
  }
}

export async function chemblSimilaritySearch(molecule: string): Promise<DG.DataFrame | null> {
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

//name: chemblSearchWidgetLocalDb
//tags: widgets
//input: string mol {semType: Molecule}
//input: bool substructure
//output: widget result
export function chemblSearchWidgetLocalDb(mol: string, substructure: boolean = false): DG.Widget {
  const headerHost = ui.div([]);
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid');
  const panel = ui.divV([headerHost, compsHost]);
  const searchFunc = substructure ? async () => chemblSubstructureSearch(mol) : async () => chemblSimilaritySearch(mol);

  searchFunc().then((table: DG.DataFrame | null) => {
    compsHost.removeChild(compsHost.firstChild!);
    if (table === null || table.rowCount === 0) {
      compsHost.appendChild(ui.divText('No matches'));
      return;
    }

    const moleculeCol = table.getCol('smiles');
    const chemblId = table.getCol('chembl_id');
    const molCount = Math.min(table.rowCount, 20);

    if (!substructure) {
      const similarityCol = table.getCol('similarity');
      const order = similarityCol.getSortedOrder();
      const descendingOrder = order.slice().sort((a, b) => similarityCol.get(b) - similarityCol.get(a));
      
      const reorderedMols: string[] = new Array<string>(descendingOrder.length);
      const reorderedScores: number[] = new Array<number>(descendingOrder.length);
      
      for (let i = 0; i < descendingOrder.length; ++i) {
        const index = descendingOrder[i];
        reorderedMols[i] = moleculeCol.get(index);
        reorderedScores[i] = similarityCol.get(index);
      }
      
      moleculeCol.init((i) => reorderedMols[i]);
      similarityCol.init((i) => reorderedScores[i]);
    }

    for (let i = 0; i < molCount; i++) {
      const molHost = ui.divV([]);
      grok.functions.call('Chem:drawMolecule', {'molStr': moleculeCol.get(i), 'w': WIDTH, 'h': HEIGHT, 'popupMenu': true})
        .then((res: HTMLElement) => {
          molHost.append(res);
          if (!substructure)
            molHost.append(ui.divText(`Score: ${table.getCol('similarity').get(i).toFixed(2)}`));
        });

      ui.tooltip.bind(molHost,
        () => ui.divText(`ChEMBL ID: ${chemblId.get(i)}\nClick to open in ChEMBL Database`));
      molHost.addEventListener('click',
        () => window.open(`https://www.ebi.ac.uk/chembl/compound_report_card/${chemblId.get(i)}`, '_blank'));
      compsHost.appendChild(molHost);
    }

    headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
      table.name = `"ChEMBL Similarity Search"`;
      grok.shell.addTableView(table);
    }, 'Open compounds as table'));
    compsHost.style.overflowY = 'auto';
  }).catch((err: any) => {
    if (compsHost.children.length > 0)
      compsHost.removeChild(compsHost.firstChild!);

    const div = ui.divText('No matches');
    ui.tooltip.bind(div, `${err}`);
    compsHost.appendChild(div);
  });
  return new DG.Widget(panel);
}

//name: Databases | ChEMBL | Substructure Search (Internal)
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSubstructureSearchPanel(mol: string): DG.Widget {
  return mol ? chemblSearchWidgetLocalDb(mol, true) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Databases | ChEMBL | Similarity Search (Internal)
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSimilaritySearchPanel(mol: string): DG.Widget {
  return mol ? chemblSearchWidgetLocalDb(mol) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Chembl targets by organism
//tags: HitTriageDataSource
//input: int maxNumberOfMolecules = 1000 [Maximum number of rows to return]
//input: string organism = "Shigella" [Organism name]
//output: dataframe compounds
export async function getChemblCompoundsByOrganism(maxNumberOfMolecules: number, organism: string): Promise<DG.DataFrame> {
  const df = await grok.data.query('Chembl:StructuresByOrganism', {maxNumberOfMolecules: maxNumberOfMolecules, organism: organism});
  return df;
}

//name: Chembl Compounds
//tags: HitTriageDataSource
//input: int maxNumberOfMolecules = 1000 [Maximum number of rows to return]
//output: dataframe compounds
export async function getChemblCompounds(maxNumberOfMolecules: number): Promise<DG.DataFrame> {
  const df = await grok.data.query('Chembl:ChemblNumberOfStructures', {maxNumberOfMolecules: maxNumberOfMolecules});
  return df;
}

//name: Chembl molregno
//tags: HitTriageFunction
//input: dataframe table [Input data table] {caption: Table}
//input: column molecules {caption: Molecules; semType: Molecule}
//output: dataframe result
export async function chemblMolregno(table: DG.DataFrame, molecules: DG.Column): Promise<DG.DataFrame> {
  const name = table.columns.getUnusedName('CHEMBL molregno');
  table.columns.addNewInt(name);
  for (let i = 0; i < molecules.length; i++) {
    const smile = molecules.get(i);
    if (!smile) {
      table.set(name, i, null);
      continue;
    }
    const canonical = grok.chem.convert(smile, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
    const resDf: DG.DataFrame = await grok.data.query('Chembl:ChemblMolregNoBySmiles', {smiles: canonical});
    const res: number = resDf.getCol('molregno').toList()[0];
    table.set(name, i, res);
  }
  return table;
}


//name: chemblIdToSmilesTs
//meta.role: converter
//meta.inputRegexp: (CHEMBL[0-9]+)
//input: string id = "CHEMBL1185" { semType: CHEMBL_ID }
//output: string smiles { semType: Molecule }
export async function chemblIdToSmilesTs(id: string): Promise<string> {
  // this function won't be needed once we include queries to functions in the startupData
  // damn, it returns null instead of the actual molecule
  return await grok.functions.call('Chembl:chemblIdToSmiles', {id: id});
  //return 'CN(C)CCc1c[nH]c2ccc(C[C@H]3COC(=O)N3)cc12';
}

//name: chemblBioactivitySearchWidget
//tags: search
//input: string s
//output: widget w
export async function chemblBioactivitySearchWidget(s: string) {
  return await chemblBioactivityForTargetsSearch(s);
}

//name: chemblPKForDrugSearchWidget
//tags: search
//input: string s
//output: widget w
export async function chemblPKForDrugSearchWidget(s: string) {
  return await chemblPKForDrugSearch(s);
}

