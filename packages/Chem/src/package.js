import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createRDKit} from './RDKit_minimal_2021.03_17.js';
import {getMolColumnPropertyPanel} from './chem_column_property_panel.js';
import * as chemSearches from './chem_searches.js';
import {setSearchesRdKitModule, moleculesToFingerprints} from './chem_searches.js';
import {RDKitCellRenderer} from './rdkit_cell_renderer.js';
import {SubstructureFilter} from './chem_substructure_filter.js';

let rdKitModule = null;
let rdKitWorkerWebRoot = null;
let initialized = false;
let structure = {};
export let _package = new DG.Package();


//name: initChem
export async function initChem() {
  if (!initialized) {
    // structure.name = "Chem";
    rdKitWorkerWebRoot = _package.webRoot;
    rdKitModule = await createRDKit(rdKitWorkerWebRoot);
    setSearchesRdKitModule(rdKitModule);
    console.log('RDKit (package) initialized');
    rdKitModule.prefer_coordgen(false);
    _package.STORAGE_NAME = 'rdkit_descriptors';
    _package.KEY = 'selected';
    _package.rdKitRendererCache = new DG.LruCache();
    _package.rdKitRendererCache.onItemEvicted = (mol) => {
      mol.delete();
    };
    initialized = true;
  }
}

//tags: init
export async function init() {
  return initChem();
}

//name: SubstructureFilter
//description: RDKit-based substructure filter
//tags: filter
//output: filter result
export function substructureFilter() {
  return new SubstructureFilter();
}

function _svgDiv(mol) {
  let root = ui.div();
  root.innerHTML = mol.get_svg();
  return root;
}

//name: getCLogP
//input: string smiles {semType: Molecule}
//output: double cLogP
export function getCLogP(smiles) {
  let mol = rdKitModule.get_mol(smiles);
  return JSON.parse(mol.get_descriptors()).CrippenClogP;
}

//name: RDKit Info
//tags: panel, widgets
//input: string smiles {semType: Molecule}
//output: widget result
export function rdkitInfoPanel(smiles) {
  let mol = rdKitModule.get_mol(smiles);
  return new DG.Widget(ui.divV([
    this._svgDiv(mol),
    ui.divText(`${this.getCLogP(smiles)}`)
  ]));
}

//name: RDKit Settings
//input: column molColumn {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molColumn) {
  return getMolColumnPropertyPanel(molColumn);
}

//name: rdkitCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//meta-cell-renderer-sem-type: Molecule
//output: grid_cell_renderer result
export async function rdkitCellRenderer() {
  //let props = DG.toJs(await this.getProperties());
  // if (props?.Renderer && props.Renderer === 'RDKit') {
  return new RDKitCellRenderer(rdKitModule);
  //}
}

//name: similarityScoring
//input: column molStringsColumn
//input: string molString
//input: bool sorted
//output: dataframe result
// deprecated
export async function similarityScoring(molStringsColumn, molString, sorted) {
  try {
    if (molStringsColumn === null || molString === null) throw "An input was null";
    let result = await _chemSimilarityScoring(molStringsColumn, molString, {'sorted': sorted});
    if (result == null) {
      return DG.DataFrame.create();
    }
    return (sorted ? result : DG.DataFrame.fromColumns([result]));
  } catch (e) {
    console.error("In similarityScoring: " + e.toString());
    throw e;
  }
}

//name: getSimilarities
//input: column molStringsColumn
//input: string molString
//output: dataframe result
export async function getSimilarities(molStringsColumn, molString) {
  try {
    if (molStringsColumn === null || molString === null) throw "An input was null";
    let result = await chemSearches.chemGetSimilarities(molStringsColumn, molString);
    // TODO: get rid of a wrapping DataFrame and be able to return Columns
    return result ? DG.DataFrame.fromColumns([result]) : DG.DataFrame.create();
  } catch (e) {
    console.error("In getSimilarities: " + e.toString());
    throw e;
  }
}

//name: getMorganFingerprints
//input: column molColumn {semType: Molecule}
//output: column result [fingerprints]
export function getMorganFingerprints(molColumn) {
  return moleculesToFingerprints(molColumn);
}

//name: findSimilar
//input: column molStringsColumn
//input: string molString
//input: int limit
//input: int cutoff
//output: dataframe result
export async function findSimilar(molStringsColumn, molString, aLimit, aCutoff) {
  try {
    if (molStringsColumn === null || molString === null || aLimit === null || aCutoff === null) throw "An input was null";
    let result = await chemSearches.chemFindSimilar(molStringsColumn, molString, {limit: aLimit, cutoff: aCutoff});
    return result ? result : DG.DataFrame.create();
  } catch (e) {
    console.error("In getSimilarities: " + e.toString());
    throw e;
  }
}

//name: searchSubstructure
//input: column molStringsColumn
//input: string molString
//input: bool substructLibrary
//input: string molStringSmarts
//output: column result
export async function searchSubstructure(molStringsColumn, molString, substructLibrary, molStringSmarts) {
  try {
    if (molStringsColumn === null || molString === null || substructLibrary === null || molStringSmarts === null)
      throw "An input was null";
    let result =
      substructLibrary ?
        await chemSearches.chemSubstructureSearchLibrary(molStringsColumn, molString, molStringSmarts, rdKitWorkerWebRoot) :
        chemSearches.chemSubstructureSearchGraph(molStringsColumn, molString);
    return DG.Column.fromList('object', 'bitset', [result]);
  } catch (e) {
    console.error("In substructureSearch: " + e.toString());
    throw e;
  }
}

//name: substructureSearch
//input: column molStringsColumn
//input: string molString
//input: bool substructLibrary
//output: column result
// deprecated
export async function substructureSearch(molStringsColumn, molString, substructLibrary) {
  return this.searchSubstructure(molStringsColumn, molString, substructLibrary);
}

//tags: app
function descriptorsApp(context) {
  let defaultSmiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
  let sketcherValue = defaultSmiles;

  let windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  let table = DG.DataFrame.create();
  table.name = 'Descriptors';
  let view = grok.shell.addTableView(table);

  let dsDiv = ui.divV([], 'grok-prop-panel');
  dsDiv.appendChild(this.descriptorsWidget(defaultSmiles).root);

  let sketcher = grok.chem.sketcher((smiles, molfile) => {
    sketcherValue = smiles;
    removeChildren(dsDiv);
    dsDiv.appendChild(this.descriptorsWidget(smiles).root);
  }, defaultSmiles);
  let addButton = ui.bigButton('ADD', async () => {
    this.getSelected().then(selected => {
      grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'smiles', selected).then(t => {
        let columnNames = table.columns.names();
        if ((table.columns.length !== selected.length + 1) || selected.some(s => !columnNames.includes(s))) {
          table = DG.DataFrame.create();
          table.name = 'Descriptors';
          view.dataFrame = table;
          for (let col of t.columns.toList())
            table.columns.addNew(col.name, col.type);
        }
        table.rows.addNew(t.columns.toList().map(c => c.get(0)));
      });
    });
  });
  addButton.style.marginTop = '12px';
  let skDiv = ui.divV([sketcher, addButton], 'grok-prop-panel,dlg-sketcher,pure-form');

  let skNode = view.dockManager.dock(skDiv, DG.DOCK_TYPE.RIGHT, null, 'Sketcher', 0.25);
  view.dockManager.dock(dsDiv, DG.DOCK_TYPE.DOWN, skNode, 'Descriptors', 0.5);

  grok.events.onViewRemoved.subscribe((v) => {
    if (v.name === view.name) {
      windows.showToolbox = true;
      windows.showHelp = true;
      windows.showProperties = true;
    }
  });
}

//name: Chem Descriptors
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function descriptorsWidget(smiles) {
  let widget = new DG.Widget(ui.div());
  let result = ui.div();
  let selectButton = ui.bigButton('SELECT', async () => {
    openDescriptorsDialog(await this.getSelected(), async (selected) => {
      await grok.dapi.userDataStorage.postValue(this.STORAGE_NAME, this.KEY, JSON.stringify(selected));
      update();
    });
  });
  selectButton.style.marginTop = '20px';

  let update = () => {
    removeChildren(result);
    result.appendChild(ui.loader());
    this.getSelected().then(selected => {
      grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${smiles}`), 'smiles', selected).then(table => {
        removeChildren(result);
        let map = {};
        for (let descriptor of selected)
          map[descriptor] = table.col(descriptor).get(0);
        result.appendChild(ui.tableFromMap(map));
      });
    });
  }

  widget.root.appendChild(result);
  widget.root.appendChild(selectButton);

  update();

  return widget;
}

//description: Get selected descriptors
export async function getSelected() {
  let str = await grok.dapi.userDataStorage.getValue(this.STORAGE_NAME, this.KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  if (selected.length === 0) {
    selected = (await grok.chem.descriptorsTree())['Lipinski']['descriptors'].slice(0, 3).map(p => p['name']);
    await grok.dapi.userDataStorage.postValue(this.STORAGE_NAME, this.KEY, JSON.stringify(selected));
  }
  return selected;
}

//description: Open descriptors selection dialog
function openDescriptorsDialog(selected, onOK) {
  grok.chem.descriptorsTree().then(descriptors => {
    let tree = ui.tree();
    tree.root.style.maxHeight = '400px';

    let groups = {};
    let items = [];

    for (let groupName in descriptors) {
      let group = tree.group(groupName, null, false);
      group.enableCheckBox();
      groups[groupName] = group;

      for (let descriptor of descriptors[groupName]['descriptors']) {
        let item = group.item(descriptor['name'], descriptor);
        item.enableCheckBox(selected.includes(descriptor['name']));
        items.push(item);
      }
    }

    let clear = ui.button('NONE', () => {
      for (let g in groups) groups[g].checked = false;
      for (let i of items) i.checked = false;
    });

    ui.dialog('Chem Descriptors')
      .add(clear)
      .add(tree.root)
      .onOK(() => onOK(items.filter(i => i.checked).map(i => i.value['name'])))
      .show();
  });
}

function renderMolRdKitCanvasCache(molString, canvas, x, y, w, h) {

  let mol = this.rdKitRendererCache.getOrCreate(molString, (s) => {
    try {
      return rdKitModule.get_mol(s);
    } catch (e) {
      return rdKitModule.get_mol("");
    }
  });

  const opts = {
    "clearBackground": false,
    "offsetx": Math.floor(x),
    "offsety": -Math.floor(y),
    "width": Math.floor(w),
    "height": Math.floor(h)
  };
  mol.draw_to_canvas_with_highlights(canvas, JSON.stringify(opts));

}

//description: Removes all children from node
function removeChildren(node) {
  while (node.firstChild)
    node.removeChild(node.firstChild);
}

//name: saveAsSdf
//description: Save as SDF
//tags: fileExporter
export function saveAsSdf() {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: open dialog
  //todo: UI for choosing structure column if necessary
  //todo: UI for choosing columns with properties

  let table = grok.shell.t;
  let structureColumn = table.columns.bySemType('Molecule');
  if (structureColumn == null)
    return;

  let result = '';

  for (let i = 0; i < table.rowCount; i++) {
    try {
      let mol = new OCL.Molecule.fromSmiles(structureColumn.get(i));
      result += `\n${mol.toMolfile()}\n`;

      // properties
      for (let col of table.columns)
        if (col !== structureColumn) {
          result += `>  <${col.name}>\n${col.get(i)}\n\n`;
        }

      result += '$$$$'
    }
    catch (error) {
      console.error(error);
    }
  }

  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}

//name: Structure WIP
//description: 2D molecule representation
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function structureWIP(smiles) {
  const mol = rdKitModule.get_mol(smiles);
  return new DG.Widget(this._svgDiv(mol));
}

//name: Properties WIP
//description: Basic molecule properties
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function propertiesWIP(smiles) {
  const mol = new OCL.Molecule.fromSmiles(smiles);
  const formula = mol.getMolecularFormula();
  const molProps = new OCL.MoleculeProperties(mol);

  //TODO: name, need PubChem
  const map = {
    'SMILES': smiles,
    'Formula': formula.formula,
    'MW': formula.absoluteWeight,
    'Number of HBA': molProps.acceptorCount,
    'Number of HBD': molProps.donorCount,
    'LogP': molProps.logP,
    'LogS': molProps.logS,
    'Polar Surface Area': molProps.polarSurfaceArea,
    'Number of rotatabe bonds': molProps.rotatableBondCount,
    'Number of stereo centers': molProps.stereoCenterCount,
  };

  return new DG.Widget(ui.tableFromMap(map));
}

//name: SDF WIP
//description: Molecule as SDF
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function sdfWIP(smiles) {
  const mol = new OCL.Molecule.fromSmiles(smiles);
  return new DG.Widget(ui.textInput('', mol.toMolfile()).root);
}

//name: 3D WIP
//description: 3D molecule representation
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function get3dWIP(smiles) {
  //TODO: implement
  //what is dml?
  return new DG.Widget();
}

//name: Toxicity WIP
//description: Toxicity prediction. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function toxicityWIP(smiles) {
  const mol = new OCL.Molecule.fromSmiles(smiles);
  const riskTypes = {
    0: 'Mutagenicity',
    1: 'Tumorigenicity',
    2: 'Irritating effects',
    3: 'Reproductive effects'
  };
  const riskLevels = {
    0: 'Unknown',
    1: 'None',
    2: 'Low',
    3: 'High'
  };

  const risks = {};
  Object.keys(riskTypes).forEach((typeId) => {
    risks[riskTypes[typeId]] = riskLevels[new OCL.ToxicityPredictor().assessRisk(mol, typeId)];
  });

  //FIXME: no such settings as in Dart: processValueElement, Aydar
  return new DG.Widget(ui.tableFromMap(risks));
}

//name: Drug Likeness WIP
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/drug-likeness.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function drugLikenessWIP(smiles) {
  //TODO: implement
  //what is chem?
  return new DG.Widget();
}

//name: Structural Alerts WIP
//description: Screening drug candidates against structural alerts, i.e. chemical fragments associated to a toxicological response
//help-url: /help/domains/chem/info-panels/structural-alerts.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function structuralAlertsWIP(smiles) {
  //TODO: implement
  //what is chem?
  return new DG.Widget();
}

//tags: unitTest
export async function _testSubstructureSearch() {
  await testSubstructureSearch();
}