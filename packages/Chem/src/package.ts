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
import {initStructuralAlertsContext, structuralAlertsWidget} from './widgets/structural-alerts';
import {structure2dWidget} from './widgets/structure2d';
import {structure3dWidget} from './widgets/structure3d';
import {toxicityWidget} from './widgets/toxicity';
import {OCLCellRenderer} from './ocl_cell_renderer';
import {chemSpace} from './analysis/chem_space';
import {getActivityCliffs} from './analysis/activity-cliffs';
import {getDescriptorsSingle} from './descriptors/descriptors_calculation';
import {addDescriptors} from './descriptors/descriptors_calculation';
import {getDescriptorsApp} from './descriptors/descriptors_calculation';
import {addMcs} from './panels/find-mcs';
import {addInchis} from './panels/inchi';
import {addInchiKeys} from './panels/inchi';
import * as chemCommonRdKit from './chem_common_rdkit';
import {rGroupAnalysis} from './analysis/r_group';
import {chemLock, chemUnlock} from './chem_common';
import {MoleculeViewer} from './chem_similarity_search';
import { identifiersWidget } from './widgets/identifiers';

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
//meta.semType: Molecule
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
  getDescriptorsApp();
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
//top-menu: Chem | R-Groups Analysis
export function rGroupsAnalysisMenu() {
  const col = grok.shell.t.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col === null) {
    grok.shell.error('Current table does not contain molecules')
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

//name: Chem | Descriptors Port...
//friendly-name: Chem | Descriptors Port...
//tags: panel, chem
//input: column smiles { semType: Molecule }
//output: string result
export async function descriptorsPanel(smiles: DG.Column) {
  let table: DG.DataFrame = grok.shell.t;
  addDescriptors(smiles, table);
}

//name: Chem | R-Groups Analysis Port
//friendly-name: Chem | R-Groups Analysis Port
//tags: panel, chem
//input: column col {semType: Molecule}
export function rGroupsAnalysisPanel(col: DG.Column) {
  rGroupAnalysis(col);
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

//#endregion