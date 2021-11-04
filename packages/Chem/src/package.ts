import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {createRDKit} from './RDKit_minimal_2021.03_17.js';
import {getMolColumnPropertyPanel} from './chem_column_property_panel';
import * as chemSearches from './chem_searches';
import {setSearchesRdKitModule, moleculesToFingerprints} from './chem_searches';
import {setCommonRdKitModule, drawMoleculeToCanvas} from './chem_common';
import {SubstructureFilter} from './chem_substructure_filter';
import {RDKitCellRenderer} from './rdkit_cell_renderer';
import * as OCL from 'openchemlib/full.js';
import {drugLikenessWidget} from './widgets/drug-likeness';
import {molfileWidget} from './widgets/molfile';
import {propertiesWidget} from './widgets/properties';
import {loadAlertsCollection, setStructuralAlertsRdKitModule, structuralAlertsWidget} from './widgets/structural-alerts';
import {structure2dWidget} from './widgets/structure2d';
import {structure3dWidget} from './widgets/structure3d';
import {toxicityWidget} from './widgets/toxicity';
import {OCLCellRenderer} from './ocl_cell_renderer';
import {getRGroups, getMCS} from "./chem_rgroup_analysis";

export let rdKitModule: any = null;
let rdKitWorkerWebRoot: string | undefined;
let initialized = false;
let structure = {};
const _STORAGE_NAME = 'rdkit_descriptors';
const _KEY = 'selected';

export const _package = new DG.Package();

//name: initChem
export async function initChem() {
  if (!initialized) {
    // structure.name = "Chem";
    rdKitWorkerWebRoot = _package.webRoot;
    // @ts-ignore
    rdKitModule = await createRDKit(rdKitWorkerWebRoot);
    setSearchesRdKitModule(rdKitModule);
    setCommonRdKitModule(rdKitModule);
    setStructuralAlertsRdKitModule(rdKitModule);
    await loadAlertsCollection();
    console.log('RDKit (package) initialized');
    rdKitModule.prefer_coordgen(false);
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
  let mol = rdKitModule.get_mol(smiles);
  return JSON.parse(mol.get_descriptors()).CrippenClogP;
}

//name: RDKit Info
//tags: panel, widgets
//input: string smiles {semType: Molecule}
//output: widget result
export function rdkitInfoPanel(smiles: string) {
  let mol = rdKitModule.get_mol(smiles);
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
  return new RDKitCellRenderer(rdKitModule);
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
  try {
    if (molStringsColumn === null || molString === null) throw "An input was null";
    // TODO: Make in the future so that the return type is always one
    let result = (await chemSearches.chemGetSimilarities(molStringsColumn, molString)) as DG.Column;
    // TODO: get rid of a wrapping DataFrame and be able to return Columns
    return result ? DG.DataFrame.fromColumns([result]) : DG.DataFrame.create();
  } catch (e: any) {
    console.error("In getSimilarities: " + e.toString());
    throw e;
  }
}

//name: getMorganFingerprints
//input: column molColumn {semType: Molecule}
//output: column result [fingerprints]
export function getMorganFingerprints(molColumn: DG.Column) {
  return moleculesToFingerprints(molColumn);
}

//name: findSimilar
//input: column molStringsColumn
//input: string molString
//input: int limit
//input: int cutoff
//output: dataframe result
export async function findSimilar(molStringsColumn: DG.Column, molString: string, aLimit: number, aCutoff: number) {
  try {
    if (molStringsColumn === null || molString === null || aLimit === null || aCutoff === null) throw "An input was null";
    let result = await chemSearches.chemFindSimilar(molStringsColumn, molString, {limit: aLimit, cutoff: aCutoff});
    return result ? result : DG.DataFrame.create();
  } catch (e: any) {
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
export async function searchSubstructure(molStringsColumn: DG.Column, molString: string, substructLibrary: boolean, molStringSmarts: string) {
  try {
    if (molStringsColumn === null || molString === null || substructLibrary === null || molStringSmarts === null)
      throw "An input was null";
    let result =
      substructLibrary ?
        await chemSearches.chemSubstructureSearchLibrary(molStringsColumn, molString, molStringSmarts, rdKitWorkerWebRoot) :
        chemSearches.chemSubstructureSearchGraph(molStringsColumn, molString);
    return DG.Column.fromList('object', 'bitset', [result]);
  } catch (e: any) {
    console.error("In substructureSearch: " + e.toString());
    throw e;
  }
}

//tags: app
function descriptorsApp(context: any) {
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
  return drugLikenessWidget(smiles);
}

//name: Molfile
//description: Molecule as Molfile
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function molfile(smiles: string) {
  return molfileWidget(smiles);
}

//name: Properties
//description: Basic molecule properties
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function properties(smiles: string) {
  return propertiesWidget(smiles);
}

//name: Structural Alerts
//description: Screening drug candidates against structural alerts, i.e. chemical fragments associated to a toxicological response
//help-url: /help/domains/chem/info-panels/structural-alerts.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structuralAlerts(smiles: string) {
  return structuralAlertsWidget(smiles);
}

//name: Structure 2D
//description: 2D molecule representation
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function structure2d(smiles: string) {
  return structure2dWidget(smiles);
}

//name: Structure 3D
//description: 3D molecule representation
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structure3d(smiles: string) {
  return structure3dWidget(smiles);
}

//name: Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function toxicity(smiles: string) {
  return toxicityWidget(smiles);
}

//name: rGroupsAnalytics
//input: dataframe df
//input: column col {semType: Molecule}
export function rGroupsAnalytics(df: DG.DataFrame, col: DG.Column) {
  let sketcherSmile = '';
  function onChanged(smiles: string) {
    sketcherSmile = smiles;
  }

  let sketcher = grok.chem.sketcher(onChanged, sketcherSmile);
  let columnPrefixInput = ui.stringInput('Column prefix', 'R');
  let visualAnalysisCheck = ui.boolInput('Visual analysis', true);

  let mcsButton = ui.button('MCS', async () => {
    let smiles = await getMCS(col);
    // sketcher.setSmiles(smiles);
    // sketcher.remove();
    // sketcherSmile = smiles;
    // sketcher = grok.chem.sketcher(onChanged, sketcherSmile);
    // mcsButton.insertAdjacentElement('beforebegin', sketcher);
  });

  let dlg = ui.dialog({
    title: 'R-Group Analysis',
    helpUrl: '/help/domains/chem/cheminformatics.md#r-group-analysis'
    })
    .add(ui.div([
      sketcher,
      ui.tooltip.bind(mcsButton, "Most Common Substructure"),
      columnPrefixInput,
      visualAnalysisCheck
    ]))
    .onOK(async () => {
      let res = await getRGroups(col, sketcherSmile, columnPrefixInput.value);
      for (let resCol of res.columns) {
        resCol.semType = DG.SEMTYPE.MOLECULE;
        col.dataFrame.columns.add(resCol);
      }
      if (res.columns.length == 0)
        grok.shell.error("None R-Groups were found");
      let view = null;
      for (let v of grok.shell.tableViews) {
        view = v;
        break;
      }
      if (visualAnalysisCheck.value && view) {
        let plot = view.trellisPlot({
          xColumnNames: [res.columns[0].name],
          yColumnNames: [res.columns[1].name]});
      }
    });
  dlg.show();
  dlg.initDefaultHistory();
}

//name: ChemSimilaritySpace
//input: dataframe table
//input: column molColumn {semType: Molecule}
//input: int cycleNum = 100
//input: bool allowLongParameters = false
//output: graphics
export async function chemSimilaritySpace(table: DG.DataFrame, molColumn: DG.Column, cycleNum: number, allowLongParameters: number) {
  const fpColumn = getMorganFingerprints(molColumn);
  if (fpColumn.stats.missingValueCount > 0) {
    throw new Error('Molecule column has a null entry');
  }

  if (cycleNum * fpColumn.length * 100 >= (1e9) && !allowLongParameters) {
    throw new Error('The given cycle and step numbers are too high to be runned. \
    If you want to run it anyway, please check the parameter allowLongParameters');
  }

  if (window.Worker) {
    const myWorker = new Worker(rdKitWorkerWebRoot + 'src/chem_stochastic_proximity_embedding.js');
    const fpBuffers = new Array(fpColumn.length);

    for (let i = 0; i < fpColumn.length; ++i) {
      const buffer = fpColumn.get(i).getBuffer();
      fpBuffers[i] = buffer;
    }

    myWorker.postMessage([fpColumn.length, fpBuffers,
      2, null, null, 1.0, 2.0, 0.01, fpColumn.length * 100, cycleNum]);

    return new Promise<void>((resolve, reject) => {
      myWorker.onmessage = function(event) {
        const coordinates = event.data;
        const coords = [
          DG.Column.fromFloat32Array('SPE_X', coordinates[0]),
          DG.Column.fromFloat32Array('SPE_Y', coordinates[1]),
        ];
        table = DG.DataFrame.fromColumns(table.columns.toList().concat(coords));
        const view = grok.shell.addTableView(table);
        view.scatterPlot({
          x: 'SPE_X',
          y: 'SPE_Y',
        });
        resolve();
      };
      myWorker.onerror = function(error) {
        reject(error.message);
      };
    });
  } else {
    throw new Error('Your browser doesn\'t support web workers.');
  }
}