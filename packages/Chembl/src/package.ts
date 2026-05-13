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
    registerChemblIdHandler();
  }


  @grok.decorators.func({
    name: 'chemblSearchWidgetLocalDb',
    meta: {role: 'widgets'},
  })
  static async chemblSearchWidgetLocalDb(
    @grok.decorators.param({options: {semType: 'Molecule'}}) mol: string,
      substructure: boolean = false): Promise<DG.Widget> {
    const headerHost = ui.div([]);
    const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid chem-search-panel-wrapper');
    const panel = ui.divV([headerHost, compsHost]);

    try {
      const table = await (substructure ? chemblSubstructureSearch(mol) : chemblSimilaritySearch(mol));
      ui.empty(compsHost);
      if (table === null || table.rowCount === 0) {
        compsHost.appendChild(ui.divText('No matches'));
        return new DG.Widget(panel);
      }

      const moleculeCol = table.getCol('smiles');
      const chemblId = table.getCol('chembl_id');
      const similarityCol = substructure ? null : table.getCol('similarity');
      const molCount = Math.min(table.rowCount, 20);

      const displayOrder = similarityCol ? similarityCol.getSortedOrder().slice().reverse() : null;

      for (let i = 0; i < molCount; i++) {
        const rowIdx = displayOrder ? displayOrder[i] : i;
        const molHost = ui.divV([]);
        grok.functions.call('Chem:drawMolecule', {'molStr': moleculeCol.get(rowIdx), 'w': WIDTH, 'h': HEIGHT, 'popupMenu': true})
          .then((res: HTMLElement) => {
            molHost.append(res);
            if (similarityCol)
              molHost.append(ui.divText(`Score: ${similarityCol.get(rowIdx).toFixed(2)}`));
          });

        ui.tooltip.bind(molHost,
          () => ui.divText(`ChEMBL ID: ${chemblId.get(rowIdx)}\nClick to open in ChEMBL Database`));
        molHost.addEventListener('click',
          () => window.open(`https://www.ebi.ac.uk/chembl/compound_report_card/${chemblId.get(rowIdx)}`, '_blank'));
        compsHost.appendChild(molHost);
      }

      headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
        table.name = substructure ? 'ChEMBL Substructure Search' : 'ChEMBL Similarity Search';
        grok.shell.addTableView(table);
      }, 'Open compounds as table'));
      compsHost.style.overflowY = 'auto';
      return new DG.Widget(panel);
    } catch (err: any) {
      ui.empty(compsHost);
      const div = ui.divText('No matches');
      ui.tooltip.bind(div, `${err}`);
      compsHost.appendChild(div);
      return new DG.Widget(panel);
    }
  }


  @grok.decorators.panel({
    name: 'Databases | ChEMBL | Substructure Search (Internal)',
    meta: {role: 'widgets'},
  })
  static async chemblSubstructureSearchPanel(
  @grok.decorators.param({options: {semType: 'Molecule'}}) mol: string): Promise<DG.Widget> {
    return mol ? await PackageFunctions.chemblSearchWidgetLocalDb(mol, true) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.panel({
    name: 'Databases | ChEMBL | Similarity Search (Internal)',
    meta: {role: 'widgets'},
  })
  static async chemblSimilaritySearchPanel(
    @grok.decorators.param({options: {semType: 'Molecule'}}) mol: string): Promise<DG.Widget> {
    return mol ? await PackageFunctions.chemblSearchWidgetLocalDb(mol) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.func({
    name: 'Chembl targets by organism',
    meta: {role: 'hitTriageDataSource'},
  })
  static async getChemblCompoundsByOrganism(
    @grok.decorators.param({type: 'int', options: {initialValue: '1000', description: 'Maximum number of rows to return'}}) maxNumberOfMolecules: number,
    @grok.decorators.param({options: {initialValue: '\'Shigella\'', description: 'Organism name'}}) organism: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:StructuresByOrganism', {maxNumberOfMolecules: maxNumberOfMolecules, organism: organism});
  }


  @grok.decorators.func({
    name: 'Chembl Compounds',
    meta: {role: 'hitTriageDataSource'},
  })
  static async getChemblCompounds(
    @grok.decorators.param({type: 'int', options: {initialValue: '1000', description: 'Maximum number of rows to return'}}) maxNumberOfMolecules: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ChemblNumberOfStructures', {maxNumberOfMolecules: maxNumberOfMolecules});
  }


  @grok.decorators.func({
    name: 'Chembl molregno',
    meta: {role: 'hitTriageFunction'},
  })
  static async chemblMolregno(
    @grok.decorators.param({options: {caption: 'Table', description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({options: {caption: 'Molecules', semType: 'Molecule'}}) molecules: DG.Column): Promise<DG.DataFrame> {
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
      role: 'converter',
      inputRegexp: '(CHEMBL[0-9]+)',
    },
    outputs: [{name: 'result', type: 'string', options: {semType: 'Molecule'}}],
  })
  static async chemblIdToSmilesTs(
  @grok.decorators.param({options: {semType: 'CHEMBL_ID', initialValue: '\'CHEMBL1185\''}}) id: string): Promise<string> {
    return await grok.functions.call('Chembl:chemblIdToSmiles', {id: id});
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Cheminformatics | Database Queries',
    },
    name: 'Database Queries',
    description: 'Running various queries to chemical databases using convenient input forms',
  })
  static async demoDatabasesChembl(): Promise<void> {
    await _demoDatabasesChembl();
    grok.shell.windows.showContextPanel = true;
    grok.shell.windows.showHelp = true;
    setTimeout(() => grok.shell.windows.help.showHelp('/help/access/databases/databases#parameterized-queries'), 1000);
  }
}


async function chemblPatternSearch(molecule: string, mode: 'substructure' | 'similarity'): Promise<DG.DataFrame | null> {
  try {
    if (!grok.chem.isMolBlock(molecule) && molecule?.length > 5000)
      throw new Error('SMILES string longer than 5000 characters not supported');
    const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
    let pattern: string;
    try {
      pattern = mode === 'substructure' ? mol.get_smarts() : mol.get_smiles();
    } finally {
      mol?.delete();
    }
    const query = mode === 'substructure' ? 'patternSubstructureSearch' : 'patternSimilaritySearch';
    return await grok.data.query(`${_package.name}:${query}`, {'pattern': pattern, 'maxRows': 100});
  } catch (e: any) {
    console.error(`In ${mode === 'substructure' ? 'Substructure' : 'Similarity'}Search: ` + e.toString());
    throw e;
  }
}


export function chemblSubstructureSearch(molecule: string): Promise<DG.DataFrame | null> {
  return chemblPatternSearch(molecule, 'substructure');
}


export function chemblSimilaritySearch(molecule: string): Promise<DG.DataFrame | null> {
  return chemblPatternSearch(molecule, 'similarity');
}
