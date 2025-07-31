/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {registerChemblIdHandler} from './handlers';
import {_demoDatabasesChembl} from './demo';

export const _package = new DG.Package();

const WIDTH = 200;
const HEIGHT = 100;

export * from './package.g';

export class PackageFunctions {
  @grok.decorators.autostart()
  static init() {
    registerChemblIdHandler(_package);
  }


  @grok.decorators.func({
    'tags': [
      'widgets',
    ],
    'name': 'chemblSearchWidgetLocalDb',
  })
  static async chemblSearchWidgetLocalDb(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) mol: string,
      substructure: boolean = false): Promise<DG.Widget> {
    const headerHost = ui.div([]);
    const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid chem-search-panel-wrapper');
    const panel = ui.divV([headerHost, compsHost]);
    const searchFunc = substructure ? async () => chemblSubstructureSearch(mol) : async () => chemblSimilaritySearch(mol);

    try {
      const table = await searchFunc();
      compsHost.removeChild(compsHost.firstChild!);
      if (table === null || table.rowCount === 0) {
        compsHost.appendChild(ui.divText('No matches'));
        return new DG.Widget(panel);
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
      return new DG.Widget(panel);
    } catch (err: any) {
      if (compsHost.children.length > 0)
        compsHost.removeChild(compsHost.firstChild!);

      const div = ui.divText('No matches');
      ui.tooltip.bind(div, `${err}`);
      compsHost.appendChild(div);
      return new DG.Widget(panel);
    }
  }


  @grok.decorators.panel({
    'tags': [
      'widgets',
    ],
    'name': 'Databases | ChEMBL | Substructure Search (Internal)',
    'condition': 'true',
  })
  static async chemblSubstructureSearchPanel(
  @grok.decorators.param({'options': {'semType': 'Molecule'}}) mol: string): Promise<DG.Widget> {
    return mol ? await PackageFunctions.chemblSearchWidgetLocalDb(mol, true) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.func({
    'tags': [
      'widgets',
    ],
    'name': 'Databases | ChEMBL | Similarity Search (Internal)',
    'condition': 'true',
  })
  static async chemblSimilaritySearchPanel(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) mol: string): Promise<DG.Widget> {
    return mol ? await PackageFunctions.chemblSearchWidgetLocalDb(mol) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.func({
    'tags': [
      'HitTriageDataSource',
    ],
    'name': 'Chembl targets by organism',
  })
  static async getChemblCompoundsByOrganism(
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '1000'}}) maxNumberOfMolecules: number,
    @grok.decorators.param({'options': {'initialValue': '\'Shigella\''}}) organism: string): Promise<DG.DataFrame> {
    const df = await grok.data.query('Chembl:StructuresByOrganism', {maxNumberOfMolecules: maxNumberOfMolecules, organism: organism});
    return df;
  }


  @grok.decorators.func({
    'tags': [
      'HitTriageDataSource',
    ],
    'name': 'Chembl Compounds',
  })
  static async getChemblCompounds(
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '1000'}}) maxNumberOfMolecules: number): Promise<DG.DataFrame> {
    const df = await grok.data.query('Chembl:ChemblNumberOfStructures', {maxNumberOfMolecules: maxNumberOfMolecules});
    return df;
  }


  @grok.decorators.func({
    'tags': [
      'HitTriageFunction',
    ],
    'name': 'Chembl molregno',
  })
  static async chemblMolregno(
    @grok.decorators.param({'options': {'caption': 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({'options': {'caption': 'Molecules', 'semType': 'Molecule'}}) molecules: DG.Column): Promise<DG.DataFrame> {
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


  @grok.decorators.func({
    meta: {
      'role': 'converter',
      'inputRegexp': '(CHEMBL[0-9]+)',
    },
  })
  static async chemblIdToSmilesTs(
  @grok.decorators.param({'options': {'semType': 'CHEMBL_ID', 'initialValue': '\'CHEMBL1185\''}}) id: string): Promise<string> {
    return await grok.functions.call('Chembl:chemblIdToSmiles', {id: id});
  }


  @grok.decorators.func({
    'meta': {
      'demoPath': 'Cheminformatics | Database Queries',
    },
    'name': 'Database Queries',
    'description': 'Running various queries to chemical databases using convenient input forms',
  })
  static async demoDatabasesChembl(): Promise<void> {
    await _demoDatabasesChembl();
    grok.shell.windows.showContextPanel = true;
    grok.shell.windows.showHelp = true;
    setTimeout(() => grok.shell.windows.help.showHelp('/help/access/databases/databases#parameterized-queries'), 1000);
  }
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
