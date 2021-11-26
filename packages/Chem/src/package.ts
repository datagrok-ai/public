import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// @ts-ignore
import {getMolColumnPropertyPanel} from './chem_column_property_panel';
import * as chemSearches from './chem_searches';
import {SubstructureFilter} from './chem_substructure_filter';
import {RDKitCellRenderer} from './rdkit_cell_renderer';
import * as OCL from 'openchemlib/full.js';
import {drugLikenessWidget} from './widgets/drug-likeness';
import {molfileWidget} from './widgets/molfile';
import {propertiesWidget} from './widgets/properties';
import {initStructuralAlertsContext, structuralAlertsWidget} from './widgets/structural-alerts-widget';
import {structure2dWidget} from './widgets/structure2d';
import {structure3dWidget} from './widgets/structure3d';
import {toxicityWidget} from './widgets/toxicity';
import {OCLCellRenderer} from './ocl_cell_renderer';
import {chemSpace} from './analysis/chem_space';
import {getDescriptors} from './descriptors/descriptors_calculation';
import * as chemCommonRdKit from './chem_common_rdkit';
import {rGroupAnalysis} from './analysis/r_group';
import {chemLock, chemUnlock} from './chem_common';
import {MoleculeViewer} from './chem_similarity_search';

let structure = {};
const _STORAGE_NAME = 'rdkit_descriptors';
const _KEY = 'selected';

const getRdKitModuleLocal = chemCommonRdKit.getRdKitModule;
const initRdKitService = chemCommonRdKit.initRdKitService;
const getRdKitService = chemCommonRdKit.getRdKitService;
const getRdKitWebRoot = chemCommonRdKit.getRdKitWebRoot;
const drawMoleculeToCanvas = chemCommonRdKit.drawMoleculeToCanvas;

/**
* Usage:
* let a = await grok.functions.call('Chem:getRdKitModule');
* let b = a.get_mol('C1=CC=CC=C1');
* alert(b.get_pattern_fp());
**/

//name: getRdKitModule
//output: object module
export function getRdKitModule() {
  return getRdKitModuleLocal();
}

export let _package: any = new DG.Package();

//name: initChem
export async function initChem() {
  chemLock('initChem');
  await initRdKitService(_package.webRoot);
  const path = getRdKitWebRoot() + 'data-samples/alert_collection.csv';
  const table = await grok.data.loadTable(path);
  const alertsSmartsList = table.columns['smarts'].toList();
  const alertsDescriptionsList = table.columns['description'].toList();
  await initStructuralAlertsContext(alertsSmartsList, alertsDescriptionsList);
  chemUnlock('initChem');
}

//tags: init
export async function init() {
  await initChem();
}

//name: initChemAutostart
//tags: autostart
export async function initChemAutostart() {
  await initChem();
}

//name: SubstructureFilter
//description: RDKit-based substructure filter
//tags: filter
//output: filter result
export function substructureFilter() {
  return new SubstructureFilter();
}

export function _svgDiv(mol: any) {
  let root = ui.div();
  root.innerHTML = mol.get_svg();
  return root;
}

//name: canvasMol
//input: int x
//input: int y
//input: int w
//input: int h
//input: object canvas
//input: string molString
//input: string scaffoldMolString
export function canvasMol(
  x: number, y: number, w: number, h: number, canvas: HTMLCanvasElement,
  molString: string, scaffoldMolString: string | null = null) {
  drawMoleculeToCanvas(x, y, w, h, canvas, molString, scaffoldMolString == '' ? null : scaffoldMolString);
}

//name: getCLogP
//input: string smiles {semType: Molecule}
//output: double cLogP
export function getCLogP(smiles: string) {
  let mol = getRdKitModuleLocal().get_mol(smiles);
  return JSON.parse(mol.get_descriptors()).CrippenClogP;
}

//name: RDKit Info
//tags: panel, widgets
//input: string smiles {semType: Molecule}
//output: widget result
export function rdkitInfoPanel(smiles: string) {
  let mol = getRdKitModuleLocal().get_mol(smiles);
  return new DG.Widget(ui.divV([
    _svgDiv(mol),
    ui.divText(`${getCLogP(smiles)}`)
  ]));
}

//name: RDKit Settings
//input: column molColumn {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molColumn: DG.Column) {
  return getMolColumnPropertyPanel(molColumn);
}

//name: rdkitCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//meta-cell-renderer-sem-type: Molecule
//output: grid_cell_renderer result
export async function rdkitCellRenderer() {
  //let props = DG.toJs(await this.getProperties());
  // if (props?.Renderer && props.Renderer === 'RDKit') {
  return new RDKitCellRenderer(getRdKitModuleLocal());
  //}
}

//name: oclCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//meta-cell-renderer-sem-type: Molecule
//output: grid_cell_renderer result
export async function oclCellRenderer() {
  return new OCLCellRenderer();
}

//name: getSimilarities
//input: column molStringsColumn
//input: string molString
//output: dataframe result
export async function getSimilarities(molStringsColumn: DG.Column, molString: string) {
  chemLock('getSimilarities');
  try {
    if (molStringsColumn === null || molString === null) throw "An input was null";
    // TODO: Make in the future so that the return type is always one
    const result = (await chemSearches.chemGetSimilarities(molStringsColumn, molString)) as unknown; // TODO: !
    // TODO: get rid of a wrapping DataFrame and be able to return Columns
    chemUnlock('getSimilarities');
    return result ? DG.DataFrame.fromColumns([result as DG.Column]) : DG.DataFrame.create();
  } catch (e: any) {
    console.error("In getSimilarities: " + e.toString());
    throw e;
  }
  chemUnlock('getSimilarities');
}

//name: getMorganFingerprints
//input: column molColumn {semType: Molecule}
//output: column result [fingerprints]
export async function getMorganFingerprints(molColumn: DG.Column) {
  chemLock('getMorganFingerprints');
  const fingerprints = await chemSearches.chemGetMorganFingerprints(molColumn);
  const fingerprintsBitsets: DG.BitSet[] = [];
  for (let i = 0; i < fingerprints.length; ++i) {
    const fingerprint = DG.BitSet.fromBytes(fingerprints[i].getRawData().buffer, fingerprints[i].length);
    fingerprintsBitsets.push(fingerprint);
  }
  chemUnlock('getMorganFingerprints');
  return DG.Column.fromList('object', 'fingerprints', fingerprintsBitsets);
}

//name: getMorganFingerprint
//input: string molString {semType: Molecule}
//output: object fingerprintBitset [Fingerprints]
export function getMorganFingerprint(molString: string) {
  const bitArray = chemSearches.chemGetMorganFingerprint(molString);
  return DG.BitSet.fromBytes(bitArray.getRawData(), bitArray.length);
}

//name: findSimilar
//input: column molStringsColumn
//input: string molString
//input: int limit
//input: int cutoff
//output: dataframe result
export async function findSimilar(molStringsColumn: DG.Column, molString: string, aLimit: number, aCutoff: number) {
  chemLock('chemFindSimilar');
  try {
    if (molStringsColumn === null || molString === null || aLimit === null || aCutoff === null) throw "An input was null";
    let result = await chemSearches.chemFindSimilar(molStringsColumn, molString, {limit: aLimit, cutoff: aCutoff});
    chemUnlock('chemFindSimilar');
    return result ? result : DG.DataFrame.create();
  } catch (e: any) {
    console.error("In getSimilarities: " + e.toString());
    throw e;
  }
  chemUnlock('chemFindSimilar');
}

//name: searchSubstructure
//input: column molStringsColumn
//input: string molString
//input: bool substructLibrary
//input: string molStringSmarts
//output: column result
export async function searchSubstructure(molStringsColumn: DG.Column, molString: string, substructLibrary: boolean, molStringSmarts: string) {
  chemLock('searchSubstructure');
  try {
    if (molStringsColumn === null || molString === null || substructLibrary === null || molStringSmarts === null)
      throw "An input was null";
    let result =
      substructLibrary ?
        await chemSearches.chemSubstructureSearchLibrary(molStringsColumn, molString, molStringSmarts/*, webRoot*/) :
        chemSearches.chemSubstructureSearchGraph(molStringsColumn, molString);
    chemUnlock('searchSubstructure');
    return DG.Column.fromList('object', 'bitset', [result]);
  } catch (e: any) {
    console.error("In substructureSearch: " + e.toString());
    throw e;
  }
  chemUnlock('searchSubstructure');
}

//name: Descriptors App
//tags: app
export function descriptorsApp(context: any) {
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
  dsDiv.appendChild(descriptorsWidget(defaultSmiles).root);

  let sketcher = grok.chem.sketcher((smiles: string, molfile: string) => {
    sketcherValue = smiles;
    removeChildren(dsDiv);
    dsDiv.appendChild(descriptorsWidget(smiles).root);
  }, defaultSmiles);
  let addButton = ui.bigButton('ADD', async () => {
    getSelected().then(selected => {
      grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'smiles', selected).then(t => {
        let columnNames = table.columns.names();
        if ((table.columns.length !== selected.length + 1) || selected.some((s: any) => !columnNames.includes(s))) {
          table = DG.DataFrame.create();
          table.name = 'Descriptors';
          view.dataFrame = table;
          for (let col of t.columns.toList())
            table.columns.addNew(col.name, col.type);
        }
        table.rows.addNew(t.columns.toList().map((c: any) => c.get(0)));
      });
    });
  });
  addButton.style.marginTop = '12px';
  let skDiv = ui.divV([sketcher, addButton], 'grok-prop-panel,dlg-sketcher,pure-form');

  let skNode = view.dockManager.dock(skDiv, DG.DOCK_TYPE.RIGHT, null, 'Sketcher', 0.25);
  view.dockManager.dock(dsDiv, DG.DOCK_TYPE.DOWN, skNode, 'Descriptors', 0.5);

  grok.events.onViewRemoved.subscribe((v: any) => {
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
export function descriptorsWidget(smiles: string) {
  let widget = new DG.Widget(ui.div());
  let result = ui.div();
  let selectButton = ui.bigButton('SELECT', async () => {
    openDescriptorsDialog(await getSelected(), async (selected: any) => {
      await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
      update();
    });
  });
  selectButton.style.marginTop = '20px';

  let update = () => {
    removeChildren(result);
    result.appendChild(ui.loader());
    getSelected().then(selected => {
      grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${smiles}`), 'smiles', selected).then((table: any) => {
        removeChildren(result);
        let map: { [_: string]: any } = {};
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
  let str = await grok.dapi.userDataStorage.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  if (selected.length === 0) {
    selected = (await grok.chem.descriptorsTree() as any)['Lipinski']['descriptors'].slice(0, 3).map((p: any) => p['name']);
    await grok.dapi.userDataStorage.postValue(_STORAGE_NAME, _KEY, JSON.stringify(selected));
  }
  return selected;
}

//description: Open descriptors selection dialog
function openDescriptorsDialog(selected: any, onOK: any) {
  grok.chem.descriptorsTree().then((descriptors: { [_: string]: any }) => {
    let tree = ui.tree();
    tree.root.style.maxHeight = '400px';

    let groups: { [_: string]: any } = {};
    let items: DG.TreeViewNode[] = [];

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
      .onOK(() => onOK(items.filter(i => i.checked).map((i: any) => i.value['name'])))
      .show();
  });
}

//description: Removes all children from node
function removeChildren(node: any) {
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
      let mol = OCL.Molecule.fromSmiles(structureColumn.get(i));
      result += `\n${mol.toMolfile()}\n`;

      // properties
      for (let col of table.columns)
        if (col !== structureColumn) {
          result += `>  <${col.name}>\n${col.get(i)}\n\n`;
        }

      result += '$$$$'
    } catch (error) {
      console.error(error);
    }
  }

  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}

//name: Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/drug-likeness.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function drugLikeness(smiles: string) {
  return smiles ? new DG.Widget(ui.divText('SMILES is empty')): drugLikenessWidget(smiles);
}

//name: Molfile
//description: Molecule as Molfile
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function molfile(smiles: string) {
  return smiles ? new DG.Widget(ui.divText('SMILES is empty')): molfileWidget(smiles);
}

//name: Properties
//description: Basic molecule properties
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function properties(smiles: string) {
  return smiles ? new DG.Widget(ui.divText('SMILES is empty')): propertiesWidget(smiles);
}

//name: Structural Alerts
//description: Screening drug candidates against structural alerts, i.e. chemical fragments associated to a toxicological response
//help-url: /help/domains/chem/info-panels/structural-alerts.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structuralAlerts(smiles: string) {
  return smiles ? new DG.Widget(ui.divText('SMILES is empty')): structuralAlertsWidget(smiles);
}

//name: Structure 2D
//description: 2D molecule representation
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function structure2d(smiles: string) {
  return smiles ? new DG.Widget(ui.divText('SMILES is empty')): structure2dWidget(smiles);
}

//name: Structure 3D
//description: 3D molecule representation
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structure3d(smiles: string) {
  return smiles ? new DG.Widget(ui.divText('SMILES is empty')): structure3dWidget(smiles);
}

//name: Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function toxicity(smiles: string) {
  return smiles ? new DG.Widget(ui.divText('SMILES is empty')): toxicityWidget(smiles);
}

//name: R-Groups Analysis
//top-menu: Chem | R-Groups Analysis
export function rGroupsAnalysisMenu() {
  const col = grok.shell.t.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col === null) {
    grok.shell.error('Current table does not contain molecules')
    return;
  }
  rGroupAnalysis(col);
}

//name: Chem | R-Groups Analysis
//friendly-name: Chem | R-Groups Analysis
//tags: panel, chem
//input: column col {semType: Molecule}
export function rGroupsAnalysisPanel(col: DG.Column) {
  rGroupAnalysis(col);
}

//top-menu: Chem | Chemical Space...
//name: Chem Space
//input: dataframe table
//input: column smiles { semType: Molecule }
//output: viewer result
export async function chemSpaceTopMenu(table: DG.DataFrame, smiles: DG.Column) {
  return new Promise<void>(async (resolve, reject) => {
    await chemSpace(table, smiles);
    resolve();
  });
}

//name: Chem | DescriptorsPort...
//tags: panel
//input: column smiles { semType: Molecule }
//output: string result
export async function descriptors(table: DG.DataFrame, smiles: DG.Column) {
  return(getDescriptors(smiles));
}

/*
//top-menu: Chem | Similarity Search...
export async function chemSimilaritySearch() {
  // Shouldn't name it similaritySearch
  grok.shell.addTableView(grok.data.demo.molecules(100)).addViewer(new MoleculeViewer());
}
 */