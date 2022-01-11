import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// @ts-ignore
import {getMolColumnPropertyPanel} from './chem-column-property-panel';
import * as chemSearches from './chem-searches';
import {SubstructureFilter} from './chem-substructure-filter';
import {RDKitCellRenderer} from './rdkit-cell-renderer';
import * as OCL from 'openchemlib/full.js';
import {drugLikenessWidget} from './widgets/drug-likeness';
import {molfileWidget} from './widgets/molfile';
import {propertiesWidget} from './widgets/properties';
import {structuralAlertsWidget} from './widgets/structural-alerts';
import {structure2dWidget} from './widgets/structure2d';
import {structure3dWidget} from './widgets/structure3d';
import {toxicityWidget} from './widgets/toxicity';
import {chemSpace} from './analysis/chem-space';
import {getActivityCliffs} from './analysis/activity-cliffs';
import {getDescriptorsSingle} from './descriptors/descriptors-calculation';
import {addDescriptors} from './descriptors/descriptors-calculation';
import {getDescriptorsApp} from './descriptors/descriptors-calculation';
import {addMcs} from './panels/find-mcs';
import {addInchis} from './panels/inchi';
import {addInchiKeys} from './panels/inchi';
import * as chemCommonRdKit from './chem-common-rdkit';
import {convertToRDKit, rGroupAnalysis} from './analysis/r-group-analysis';
import {chemLock, chemUnlock} from './chem-common';
import {identifiersWidget} from './widgets/identifiers';
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;
import {oclMol} from './chem-common-ocl';
import $ from 'cash-dom';
import '../css/chem.css';
import {RDMol} from './rdkit-api';
import {isMolBlock} from './chem-utils';

const getRdKitModuleLocal = chemCommonRdKit.getRdKitModule;
const drawMoleculeToCanvas = chemCommonRdKit.drawMoleculeToCanvas;
let initialized: boolean = false;

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

export const _package: DG.Package = new DG.Package();
let _rdRenderer: RDKitCellRenderer;
let _renderer: GridCellRendererProxy;
let _renderers: Map<string, DG.GridCellRenderer>;
let _properties: any;

//tags: init
export async function initChem() {
  chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
  await chemCommonRdKit.initRdKitModuleLocal();
  _properties = await _package.getProperties();
  _rdRenderer = new RDKitCellRenderer(getRdKitModuleLocal());
  _renderer = new GridCellRendererProxy(_rdRenderer, 'Molecule');
  _renderers = new Map();
  _properties = {};
}

//tags: autostart
export async function initChemAutostart() { }

//name: SubstructureFilter
//description: RDKit-based substructure filter
//tags: filter
//output: filter result
//meta.semType: Molecule
export function substructureFilter() {
  return new SubstructureFilter();
}

export function _svgDiv(mol: any) {
  const root = ui.div();
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

export function renderMolecule(
  molStr: string,
  options: {renderer?: 'RDKit' | 'OpenChemLib', width?: number, height?: number},
) {
  options.renderer ??= _properties.Renderer as 'RDKit' | 'OpenChemLib';
  options.width ??= 200;
  options.height ??= 150;

  let mol: OCL.Molecule | RDMol;
  let molFile: string;
  let smiles: string;
  isMolBlock(molStr) ? molFile = molStr : smiles = molStr;

  const moleculeHost = ui.canvas(options.width, options.height);
  $(moleculeHost).addClass('chem-canvas');

  switch (options.renderer) {
  case 'RDKit':
    mol = getRdKitModuleLocal().get_mol(convertToRDKit(molStr));
    (mol as RDMol).draw_to_canvas(moleculeHost, options.width, options.height);
    molFile ??= (mol as RDMol).get_molblock();
    smiles ??= (mol as RDMol).get_smiles();
    break;
  case 'OpenChemLib':
    mol = oclMol(molStr);
    OCL.StructureView.drawMolecule(moleculeHost, mol as OCL.Molecule);
    molFile ??= mol.toMolfile();
    smiles ??= mol.toSmiles();
    break;
  default:
    throw new Error(`Renderer '${options.renderer}' is not supported.`);
  }

  const moreBtn = ui.iconFA(
    'ellipsis-v',
    () => {
      const menu = DG.Menu.popup();
      menu.item('Copy SMILES', () => {
        navigator.clipboard.writeText(smiles);
        grok.shell.info('SMILES copied!');
      });
      menu.item('Copy Molfile', () => {
        navigator.clipboard.writeText(molFile);
        grok.shell.info('Molfile copied!');
      });
      menu.item('Sketch', () => {
        const sketcher = new Sketcher();
        isMolBlock(molStr) ? sketcher.setMolFile(molStr) : sketcher.setSmiles(molStr);
        ui.dialog()
          .add(sketcher)
          .show();
      });
      menu.show();
    },
    'More',
  );
  $(moreBtn).addClass('chem-mol-view-icon pep-more-icon');

  return ui.divV([moreBtn, moleculeHost], 'chem-mol-box');
}

//name: getCLogP
//input: string smiles {semType: Molecule}
//output: double cLogP
export function getCLogP(smiles: string) {
  const mol = getRdKitModuleLocal().get_mol(smiles);
  return JSON.parse(mol.get_descriptors()).CrippenClogP;
}

export class GridCellRendererProxy extends DG.GridCellRenderer {
  renderer: DG.GridCellRenderer;
  _cellType: string;

  constructor(renderer: DG.GridCellRenderer, cellType: string) {
    super();
    this.renderer = renderer;
    this._cellType = cellType;
  }

  get defaultWidth(): number | null { return this.renderer.defaultWidth;  }
  get defaultHeight(): number | null { return this.renderer.defaultHeight; }

  get name(): string { return this.renderer.name; }
  get cellType(): string { return this._cellType; }

  renderInternal(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    this.renderer.renderInternal(g, x, y, w, h, gridCell, cellStyle);
  }
}

//name: rdkitCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: RDKit
export async function rdkitCellRenderer() {
  return new RDKitCellRenderer(getRdKitModuleLocal());
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//meta-cell-renderer-sem-type: Molecule
//output: grid_cell_renderer result
export async function chemCellRenderer() {
  const renderer = _properties.Renderer ?? 'RDKit';
  if (!_renderers.has(renderer)) {
    const renderFunctions = DG.Func.find({meta: {chemRendererName: renderer}});
    if (renderFunctions.length > 0) {
      const r = await renderFunctions[0].apply();
      _renderers.set(_properties.Renderer, r);
      return r;
    }
  }

  _renderer.renderer = renderer ?? _renderer.renderer;
  return _renderer;

  //return _renderers.get(renderer) ?? new RDKitCellRenderer(getRdKitModuleLocal());
}

//name: getSimilarities
//input: column molStringsColumn
//input: string molString
//output: dataframe result
export async function getSimilarities(molStringsColumn: DG.Column, molString: string) {
  chemLock('getSimilarities');
  try {
    if (molStringsColumn === null || molString === null) throw 'An input was null';
    // TODO: Make in the future so that the return type is always one
    const result = (await chemSearches.chemGetSimilarities(molStringsColumn, molString)) as unknown; // TODO: !
    // TODO: get rid of a wrapping DataFrame and be able to return Columns
    chemUnlock('getSimilarities');
    return result ? DG.DataFrame.fromColumns([result as DG.Column]) : DG.DataFrame.create();
  } catch (e: any) {
    console.error('In getSimilarities: ' + e.toString());
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
    //@ts-ignore
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
  //@ts-ignore
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
    if (molStringsColumn === null || molString === null || aLimit === null || aCutoff === null) throw 'An input was null';
    const result = await chemSearches.chemFindSimilar(molStringsColumn, molString, {limit: aLimit, cutoff: aCutoff});
    chemUnlock('chemFindSimilar');
    return result ? result : DG.DataFrame.create();
  } catch (e: any) {
    console.error('In getSimilarities: ' + e.toString());
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
      throw 'An input was null';

    const result =
      substructLibrary ?
        await chemSearches.chemSubstructureSearchLibrary(molStringsColumn, molString, molStringSmarts/*, webRoot*/) :
        chemSearches.chemSubstructureSearchGraph(molStringsColumn, molString);
    chemUnlock('searchSubstructure');
    return DG.Column.fromList('object', 'bitset', [result]);
  } catch (e: any) {
    console.error('In substructureSearch: ' + e.toString());
    throw e;
  }
  chemUnlock('searchSubstructure');
}

//name: Descriptors App
//tags: app
export function descriptorsApp(context: any) {
  getDescriptorsApp();
}


//name: saveAsSdf
//description: Save as SDF
//tags: fileExporter
export function saveAsSdf() {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: open dialog
  //todo: UI for choosing structure column if necessary
  //todo: UI for choosing columns with properties

  const table = grok.shell.t;
  const structureColumn = table.columns.bySemType('Molecule');
  if (structureColumn == null)
    return;


  let result = '';

  for (let i = 0; i < table.rowCount; i++) {
    try {
      const mol = OCL.Molecule.fromSmiles(structureColumn.get(i));
      result += `\n${mol.toMolfile()}\n`;

      // properties
      for (const col of table.columns) {
        if (col !== structureColumn)
          result += `>  <${col.name}>\n${col.get(i)}\n\n`;
      }

      result += '$$$$';
    } catch (error) {
      console.error(error);
    }
  }

  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}


//#region Top menu

//top-menu: Chem | Activity Cliffs...
//name: Activity Cliffs
//description: detect activity cliffs
//input: dataframe df [Input data table]
//input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
export async function activityCliffs(dataframe: DG.DataFrame, smiles: DG.Column, activities: DG.Column, similarity: number) {
  await getActivityCliffs(dataframe, smiles, activities, similarity);
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

//name: R-Groups Analysis
//top-menu: Chem | R-Groups Analysis...
export function rGroupsAnalysisMenu() {
  const col = grok.shell.t.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col === null) {
    grok.shell.error('Current table does not contain molecules');
    return;
  }
  rGroupAnalysis(col);
}

/*
//top-menu: Chem | Similarity Search...
export async function chemSimilaritySearch() {
  // Shouldn't name it similaritySearch
  grok.shell.addTableView(grok.data.demo.molecules(100)).addViewer(new MoleculeViewer());
}
 */
//#endregion

//#region Molecule column property panel

//name: Chem | Substructure search Port...
//friendly-name: Chem | Substructure search Port...
//tags: panel, chem
//input: column smiles { semType: Molecule }
export async function substructurePanel(smiles: DG.Column) {
  grok.shell.warning('ss');
}

//name: Chem | Descriptors Port...
//friendly-name: Chem | Descriptors Port...
//tags: panel, chem
//input: column smiles { semType: Molecule }
//output: string result
export async function descriptorsPanel(smiles: DG.Column) {
  const table: DG.DataFrame = grok.shell.t;
  addDescriptors(smiles, table);
}

//name: Chem | Find MCS Port
//friendly-name: Chem | Find MCS Port
//tags: panel, chem
//input: column col {semType: Molecule}
export function addMcsPanel(col: DG.Column) {
  addMcs(col);
}

//name: Chem | To InchI Port
//friendly-name: Chem | To InchI Port
//tags: panel, chem
//input: column col {semType: Molecule}
export function addInchisPanel(col: DG.Column) {
  addInchis(col);
}

//name: Chem | To InchI Keys Port
//friendly-name: Chem | To InchI Keys Port
//tags: panel, chem
//input: column col {semType: Molecule}
export function addInchisKeysPanel(col: DG.Column) {
  addInchiKeys(col);
}

//name: RDKit Settings
//input: column molColumn {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molColumn: DG.Column) {
  return getMolColumnPropertyPanel(molColumn);
}

//#endregion

//#region Single molecule property panel

//name: Chem Descriptors
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function descriptorsWidget(smiles: string) {
  return getDescriptorsSingle(smiles);
}

//name: Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/drug-likeness.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function drugLikeness(smiles: string) {
  return smiles ? drugLikenessWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Molfile
//description: Molecule as Molfile
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function molfile(smiles: string) {
  return smiles ? molfileWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Properties
//description: Basic molecule properties
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function properties(smiles: string) {
  return smiles ? propertiesWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Structural Alerts
//description: Screening drug candidates against structural alerts, i.e. chemical fragments associated to a toxicological response
//help-url: /help/domains/chem/info-panels/structural-alerts.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structuralAlerts(smiles: string) {
  return smiles ? structuralAlertsWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Structure 2D
//description: 2D molecule representation
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function structure2d(smiles: string) {
  return smiles ? structure2dWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Structure 3D
//description: 3D molecule representation
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structure3d(smiles: string) {
  return smiles ? structure3dWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function toxicity(smiles: string) {
  return smiles ? toxicityWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Identifiers
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function identifiers(smiles: string) {
  return smiles ? await identifiersWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: convertMolecule
//tags: unitConverter
//input: string molecule {semType: Molecule}
//input: string from {choices:["smiles", "molblock", "inchi", "v3Kmolblock"]}
//input: string to {choices:["smiles", "molblock", "inchi", "v3Kmolblock"]}
//output: string result {semType: Molecule}
export function convertMolecule(molecule: string, from: string, to: string): string {
  let mol;
  try {
    mol = getRdKitModule().get_mol(molecule);
    if (to === 'molblock')
      return mol.get_molblock();
    if (to === 'smiles')
      return mol.get_smiles();
    if (to === 'v3Kmolblock')
      return mol.get_v3Kmolblock();
    if (to == 'inchi')
      return mol.get_inchi();
    throw `Failed to convert molecule: unknown target unit: "${to}"`;
  } finally {
    mol?.delete();
  }
}

/*//tags: cellEditor
//description: Molecule
//input: grid_cell cell
export function editMoleculeCell(cell: DG.GridCell) {
  const sketcher = new Sketcher();
  sketcher.setMolecule(cell.cell.value);

  ui.dialog()
    .add(sketcher)
    .onOK(() => cell.cell.value = sketcher.getMolFile())
    .show();
}
*/
