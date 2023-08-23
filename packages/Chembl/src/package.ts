import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

export const _package = new DG.Package();

const WIDTH = 200;
const HEIGHT = 100;

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
//input: string searchType
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
    const molregnoCol = table.getCol('molregno');
    const molCount = Math.min(table.rowCount, 20);

    for (let i = 0; i < molCount; i++) {
      const molHost = ui.divV([]);
      grok.functions.call('Chem:drawMolecule', {'molStr': moleculeCol.get(i), 'w': WIDTH, 'h': HEIGHT, 'popupMenu': true})
        .then((res: HTMLElement) => {
          molHost.append(res);
          if (!substructure)
            molHost.append(ui.divText(`Score: ${table.getCol('similarity').get(i).toFixed(2)}`));
        });

      ui.tooltip.bind(molHost,
        () => ui.divText(`ChEMBL ID: ${molregnoCol.get(i)}\nClick to open in ChEMBL Database`));
      molHost.addEventListener('click',
        () => window.open(`https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL${molregnoCol.get(i)}`, '_blank'));
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


/*
//name_: Chembl Browser
//tags_: app
export async function Browser() {
  // Filter inputs
  const molecule = ui.moleculeInput('Substructure', 'C1CCCCC1');
  const subName = ui.stringInput('Subname', '');
  const molregno = ui.stringInput('Molregno', '');
  const ro5Violation = ui.choiceInput('RO5 Violations', 'All', ['All', '0', '1', '2', '3', '4']);
  const maxPhase = ui.choiceInput('Max Phase', 'All', ['All', '0', '1', '2', '3', '4']);
  const moleculeType = ui.choiceInput('Molecule type', 'All',
    ['Protein', 'Oligonucleotide', 'Unknown', 'Antibody', 'Oligosaccharide', 'Unclassified', 'Enzyme', 'Cell', 'All']);
  // Never used. TODO: remove?
  // const clear = ui.button([ui.iconFA('trash-alt'), 'Clear filters'], () => clearFilters());

  const controlPanel = ui.form([molecule, subName, molregno, ro5Violation, maxPhase, moleculeType]);
  $(controlPanel).addClass('ui-form-condensed');

  // Filter handlers
  molecule.onChanged(() => update());
  subName.onChanged(() => update());
  molregno.onChanged(() => findByMolregno());
  ro5Violation.onChanged(() => update());
  maxPhase.onChanged(() => update());
  moleculeType.onChanged(() => update());

  let v: DG.View;
  let r: DG.Viewer;

  async function initView() {
    const parser = document.createElement('a');
    parser.href = window.location.toString();
    const pathSegments = parser.href.split('?');

    // if we came to app by link with molregno in URL
    if (pathSegments.length === 2 && pathSegments[1].includes('molregno')) {
      const parsedMolregno = parseInt(pathSegments[1].split('=')[1]);
      molregno.value = parsedMolregno.toString();
      const data = await grok.data.query(`${_package.name}:cbFindByMolregno`, {'molregno': parsedMolregno});
      data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
      v = grok.shell.newView('Chembl Browser');
      r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
      v.append(ui.divV([r]));
      v.toolbox = ui.div(controlPanel);
      v.box = true;
    } else if (pathSegments.length > 2 && pathSegments[1].includes('substructure')) {
      // if we came to app by link with set of parameters in URL
      const parsedSubstructure = pathSegments[1].split('=')[1];
      const parsedSubname = pathSegments[2].split('=')[1];
      const parsedRo5 = parseInt(pathSegments[3].split('=')[1]);
      const parsedMaxPhase = parseInt(pathSegments[4].split('=')[1]);
      const parsedMoleculeType = pathSegments[5].split('=')[1];

      const queryParameters = {
        'substructure': parsedSubstructure,
        'subname': parsedSubname,
        'num_ro5_violations': parsedRo5,
        'max_phase': parsedMaxPhase,
        'molecule_type': parsedMoleculeType,
      };

      molecule.value = parsedSubstructure;
      subName.value = parsedSubname;
      ro5Violation.value = parsedRo5.toString();
      maxPhase.value = parsedMaxPhase.toString();
      moleculeType.value = parsedMoleculeType;

      const data = await grok.data.query(`${_package.name}:cbChemblBrowserQuery`, queryParameters);
      data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
      v = grok.shell.newView('Chembl Browser');
      r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
      v.append(ui.divV([r]));
      v.toolbox = ui.div(controlPanel);
      v.box = true;
    } else {
      const data = await grok.data.query(`${_package.name}:cbAllChemblStructures`, {});
      data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
      v = grok.shell.newView('Chembl Browser');
      r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
      v.append(ui.divV([r]));
      v.toolbox = ui.div(controlPanel);
      v.box = true;
    }
  }

  async function update() {
    const ro5 = ro5Violation.value === 'All' ? -1 : parseInt(ro5Violation.value!);
    const qMaxPhase = maxPhase.value === 'All' ? -1 : parseInt(maxPhase.value!);
    const queryParameters = {
      'substructure': molecule.value,
      'subname': subName.value,
      'num_ro5_violations': ro5,
      'max_phase': qMaxPhase,
      'molecule_type': moleculeType.value,
    };
    const query = await grok.data.query(`${_package.name}:cbChemblBrowserQuery`, queryParameters);
    if (query.rowCount === 0)
      grok.shell.info('No results for this filter were found');
    else {
      query.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
      v.root.children[0].remove();
      r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
      r.setOptions({look: look});
      console.log(r.getOptions());
      v.toolbox = ui.div(controlPanel);
      v.append(ui.divV([r]));
      v.box = true;
      v.path = '';
      v.path = `/apps/Chemblbrowser?substructure=${molecule.value}?subname=${subName.value}?` +
        `num_ro5_violations=${ro5}?max_phase=${qMaxPhase}?molecule_type=${moleculeType.value}`;
    }
  }

  async function findByMolregno() {
    const query = await grok.data.query(`${_package.name}:cbFindByMolregno`, {'molregno': parseInt(molregno.value)});
    if (query.rowCount === 0)
      grok.shell.info('No results for this filter were found');
    else {
      query.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
      v.root.children[0].remove();
      r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
      v.toolbox = ui.div(controlPanel);
      v.append(ui.divV([r]));
      v.box = true;
      v.path = '/apps/Chemblbrowser?molregno=' + molregno.value;
    }
  }

  async function clearFilters() {
    molecule.value = '';
    subName.value = '';
    ro5Violation.value = 'All';
    maxPhase.value = 'All';
    moleculeType.value = 'All';
    molregno.value = '';
    v.path = '';
  }

  await initView();

  const look = {
    'sketchState': {
      '#type': 'SketchState',
      'elementStates': [
        {
          'left': 142,
          'top': 94,
          'width': 100,
          'height': 20,
          'type': 'field',
          'viewerSettings': {
            'table': 'FindByMolregno',
            'column': 'molregno',
            'format': null,
          },
        },
        {
          'left': 142,
          'top': 118,
          'width': 100,
          'height': 20,
          'type': 'field',
          'viewerSettings': {
            'table': 'FindByMolregno',
            'column': 'compound_name',
            'format': null,
          },
        },
        {
          'left': 10,
          'top': 190,
          'width': 112,
          'height': 20,
          'type': 'html',
          'viewerSettings': {
            'markup': '<label class="d4-sketch-column-name">max_phase</label>',
          },
        },
        {
          'left': 142,
          'top': 214,
          'width': 100,
          'height': 20,
          'type': 'field',
          'viewerSettings': {
            'table': 'FindByMolregno',
            'column': 'chembl_id',
            'format': null,
          },
        },
        {
          'left': 10,
          'top': 118,
          'width': 112,
          'height': 20,
          'type': 'html',
          'viewerSettings': {
            'markup': '<label class="d4-sketch-column-name">compound_name</label>',
          },
        },
        {
          'left': 142,
          'top': 142,
          'width': 100,
          'height': 20,
          'type': 'field',
          'viewerSettings': {
            'table': 'FindByMolregno',
            'column': 'num_ro5_violations',
            'format': null,
          },
        },
        {
          'left': 10,
          'top': 142,
          'width': 112,
          'height': 20,
          'type': 'html',
          'viewerSettings': {
            'markup': '<label class="d4-sketch-column-name">num_ro5_violations</label>',
          },
        },
        {
          'left': 10,
          'top': 94,
          'width': 112,
          'height': 20,
          'type': 'html',
          'viewerSettings': {
            'markup': '<label class="d4-sketch-column-name">molregno</label>',
          },
        },
        {
          'left': 142,
          'top': 166,
          'width': 100,
          'height': 20,
          'type': 'field',
          'viewerSettings': {
            'table': 'FindByMolregno',
            'column': 'synonyms',
            'format': null,
          },
        },
        {
          'left': 10,
          'top': 166,
          'width': 112,
          'height': 20,
          'type': 'html',
          'viewerSettings': {
            'markup': '<label class="d4-sketch-column-name">synonyms</label>',
          },
        },
        {
          'left': 10,
          'top': 214,
          'width': 112,
          'height': 20,
          'type': 'html',
          'viewerSettings': {
            'markup': '<label class="d4-sketch-column-name">chembl_id</label>',
          },
        },
        {
          'left': 142,
          'top': 190,
          'width': 100,
          'height': 20,
          'type': 'field',
          'viewerSettings': {
            'table': 'FindByMolregno',
            'column': 'max_phase',
            'format': null,
          },
        },
        {
          'left': 7,
          'top': 20,
          'width': 235,
          'height': 70,
          'type': 'field',
          'viewerSettings': {
            'table': 'FindByMolregno',
            'column': 'canonical_smiles',
            'format': null,
          },
        },
      ],
    },
  };
}
*/
