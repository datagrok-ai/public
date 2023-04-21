import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import '../css/chem.css';
import * as chemSearches from './chem-searches';
import {GridCellRendererProxy, RDKitCellRenderer} from './rendering/rdkit-cell-renderer';
import {getDescriptorsApp, getDescriptorsSingle} from './descriptors/descriptors-calculation';
import {assure, delay} from '@datagrok-libraries/utils/src/test';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {OpenChemLibSketcher} from './open-chem/ocl-sketcher';
import {_importSdf} from './open-chem/sdf-importer';
import {OCLCellRenderer} from './open-chem/ocl-cell-renderer';
import Sketcher = DG.chem.Sketcher;
import {activityCliffsIdx, CLIFFS_DF_NAME, getActivityCliffs, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {IUMAPOptions, ITSNEOptions} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {SequenceSpaceFunctionEditor} from '@datagrok-libraries/ml/src/functionEditors/seq-space-editor';
import {ActivityCliffsFunctionEditor} from '@datagrok-libraries/ml/src/functionEditors/activity-cliffs-editor';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils';
import {EMPTY_MOLECULE_MESSAGE, elementsTable} from './constants';
import {similarityMetric} from '@datagrok-libraries/ml/src/distance-metrics-methods';

//widget imports
import {SubstructureFilter} from './widgets/chem-substructure-filter';
import {drugLikenessWidget} from './widgets/drug-likeness';
import {propertiesWidget} from './widgets/properties';
import {structuralAlertsWidget} from './widgets/structural-alerts';
import {structure2dWidget} from './widgets/structure2d';
import {toxicityWidget} from './widgets/toxicity';

//panels imports
import {addInchiKeys, addInchis} from './panels/inchi';
import {getMolColumnPropertyPanel} from './panels/chem-column-property-panel';
import {checkForStructuralAlerts} from './panels/structural-alerts';
import { addDescriptors } from './descriptors/descriptors-calculation';

//utils imports
import { ScaffoldTreeViewer} from "./widgets/scaffold-tree";
import {Fingerprint} from './utils/chem-common';
import * as chemCommonRdKit from './utils/chem-common-rdkit';
import {IMolContext, getMolSafe, isSmarts} from './utils/mol-creation_rdkit';
import {checkMoleculeValid, checkMolEqualSmiles, _rdKitModule} from './utils/chem-common-rdkit';
import {_convertMolNotation} from './utils/convert-notation-utils';
import {molToMolblock} from './utils/convert-notation-utils';
import {getAtomsColumn, checkPackage} from './utils/elemental-analysis-utils';
import {saveAsSdfDialog} from './utils/sdf-utils';
import {getSimilaritiesMarix} from './utils/similarity-utils';

//analytical imports
import {createPropPanelElement, createTooltipElement} from './analysis/activity-cliffs';
import {chemDiversitySearch, ChemDiversityViewer} from './analysis/chem-diversity-viewer';
import {chemSimilaritySearch, ChemSimilarityViewer} from './analysis/chem-similarity-viewer';
import {chemSpace, getEmbeddingColsNames} from './analysis/chem-space';
import {rGroupAnalysis} from './analysis/r-group-analysis';
import {ChemSearchBaseViewer} from './analysis/chem-search-base-viewer';

//file importers
import {_importTripos} from './file-importers/mol2-importer';
import {_importSmi} from './file-importers/smi-importer';

//script api
import {generateScaffoldTree} from "./scripts-api";
import { renderMolecule } from './rendering/render-molecule';
import { RDKitReactionRenderer } from './rendering/rdkit-reaction-renderer';
import { structure3dWidget } from './widgets/structure3d';
import { identifiersWidget } from './widgets/identifiers';
import { closeAllAccordionPanes, demoScaffold, getAccordionPane, openSketcher, scrollTable } from './utils/demo-utils';

const drawMoleculeToCanvas = chemCommonRdKit.drawMoleculeToCanvas;
const SKETCHER_FUNCS_FRIENDLY_NAMES: {[key: string]: string} = {
  OpenChemLib: 'OpenChemLib',
  Ketcher: 'Ketcher',
  Marvin: 'Marvin',
  ChemDraw: 'ChemDraw'
}

const PREVIOUS_SKETCHER_NAMES: {[key: string]: string} = {
  'Open Chem Sketcher': 'OpenChemLib',
  'ketcherSketcher': 'Ketcher',
  'Marvin JS': 'Marvin',
  'Chem Draw': 'ChemDraw'
}

/**
 * Usage:
 * let a = await grok.functions.call('Chem:getRdKitModule');
 * let b = a.get_mol('C1=CC=CC=C1');
 * alert(b.get_pattern_fp());
 **/

//name: getRdKitModule
//output: object module
export function getRdKitModule() {
  return chemCommonRdKit.getRdKitModule();
}

export const _package: DG.Package = new DG.Package();
export let _properties: any;

let _rdRenderer: RDKitCellRenderer;
export let renderer: GridCellRendererProxy;
let _renderers: Map<string, DG.GridCellRenderer>;

//tags: init
export async function initChem(): Promise<void> {
  chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
  await chemCommonRdKit.initRdKitModuleLocal();
  _properties = await _package.getProperties();
  _rdRenderer = new RDKitCellRenderer(getRdKitModule());
  renderer = new GridCellRendererProxy(_rdRenderer, 'Molecule');
  let storedSketcherType = await grok.dapi.userDataStorage.getValue(DG.chem.STORAGE_NAME, DG.chem.KEY, true);
  if (PREVIOUS_SKETCHER_NAMES[storedSketcherType])
    storedSketcherType = PREVIOUS_SKETCHER_NAMES[storedSketcherType];
  if (!storedSketcherType && _properties.Sketcher)
    storedSketcherType = SKETCHER_FUNCS_FRIENDLY_NAMES[_properties.Sketcher]
  const sketcherFunc = DG.Func.find({tags: ['moleculeSketcher']}).find(e => e.name === storedSketcherType || e.friendlyName === storedSketcherType);
  if (sketcherFunc)
    DG.chem.currentSketcherType = sketcherFunc.friendlyName;
  else {
    if(!!storedSketcherType) {
      grok.shell.warning(`Package with ${storedSketcherType} function is not installed. Switching to ${DG.DEFAULT_SKETCHER}.`);
    }
    DG.chem.currentSketcherType = DG.DEFAULT_SKETCHER;
  }
  _renderers = new Map();
}

//tags: autostart
export async function initChemAutostart(): Promise<void> { }

//name: Chemistry | Most Diverse Structures
//tags: tooltip, panel
//input: column col {semType: Molecule}
//output: widget
export async function chemTooltip(col: DG.Column): Promise<DG.Widget | undefined> {
  const version = col.version;

  for (let i = 0; i < col.length; ++i) {
    if (!col.isNone(i) && isSmarts(col.get(i)))
      return;
  }

  const divMain = ui.div();
  divMain.append(ui.divText('Most diverse structures'));
  const divStructures = ui.div();
  divStructures.classList.add('d4-flex-wrap');

  if (col.temp['version'] !== version || col.temp['molIds'].length === 0) {
    const molIds = await chemDiversitySearch(
      col, similarityMetric['Tanimoto'], 9, 'Morgan' as Fingerprint, true);

    Object.assign(col.temp, {
      'version': version,
      'molIds': molIds
    });
  }

  const molIdsCached = col.temp['molIds'];
  for (let i = 0; i < molIdsCached.length; ++i) {
    divStructures.append(renderMolecule(col.get(molIdsCached[i]), {width: 150, height: 75}));
  }

  divMain.append(divStructures);
  return new DG.Widget(divMain);
}

//name: Scaffold Tree
//tags: viewer
//meta.trellisable: true
//meta.icon: files/icons/scaffold-tree-icon.svg
//output: viewer result
export function scaffoldTreeViewer() : ScaffoldTreeViewer {
 return new ScaffoldTreeViewer();
}

//name: SubstructureFilter
//description: RDKit-based substructure filter
//tags: filter
//output: filter result
//meta.semType: Molecule
export function substructureFilter(): SubstructureFilter {
  return new SubstructureFilter();
}

//name: canvasMol
//input: int x
//input: int y
//input: int w
//input: int h
//input: object canvas
//input: string molString
//input: string scaffoldMolString
//input: object options {optional: true}
export function canvasMol(
  x: number, y: number, w: number, h: number, canvas: HTMLCanvasElement,
  molString: string, scaffoldMolString: string | null = null,
  options = {normalizeDepiction: true, straightenDepiction: true}
): void {
  drawMoleculeToCanvas(x, y, w, h, canvas,
    molString, scaffoldMolString == '' ? null : scaffoldMolString,
    options);
}


//name: drawMolecule
//input: string molStr
//input: int w {optional: true}
//input: int h {optional: true}
//input: bool popupMenu {optional: true}
//output: object canvas
export function drawMolecule(molStr: string, w?: number, h?: number, popupMenu?: boolean): HTMLElement {
  return renderMolecule(molStr, {width: w, height: h, popupMenu: popupMenu});
}


//name: getCLogP
//input: string smiles {semType: Molecule}
//output: double cLogP
export function getCLogP(smiles: string): number {
  const mol = getRdKitModule().get_mol(smiles);
  const res = JSON.parse(mol.get_descriptors()).CrippenClogP;
  mol?.delete();
  return res;
}

//name: rdKitCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: RDKit
export async function rdKitCellRenderer(): Promise<RDKitCellRenderer> {
  return new RDKitCellRenderer(getRdKitModule());
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-ChemicalReaction
//meta.cellType: ChemicalReaction
//meta-cell-renderer-sem-type: ChemicalReaction
//output: grid_cell_renderer result
export async function rdKitReactionRenderer(): Promise<RDKitReactionRenderer> {
  return new RDKitReactionRenderer(getRdKitModule());
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//meta.cellType: Molecule
//meta-cell-renderer-sem-type: Molecule
//output: grid_cell_renderer result
export async function chemCellRenderer(): Promise<DG.GridCellRenderer> {
  const propertiesRenderer: string = _properties.Renderer ?? 'RDKit';
  if (!_renderers.has(propertiesRenderer)) {
    const renderFunctions = DG.Func.find({meta: {chemRendererName: propertiesRenderer}});
    if (renderFunctions.length > 0) {
      const r = await renderFunctions[0].apply();
      _renderers.set(_properties.Renderer, r);
      return r;
    }
  }

  renderer.renderer = _renderers.get(propertiesRenderer)!;
  return renderer;
}

export async function getMorganFingerprints(molColumn: DG.Column): Promise<DG.Column> {
  assure.notNull(molColumn, 'molColumn');

  try {
    const fingerprints = await chemSearches.chemGetFingerprints(molColumn, Fingerprint.Morgan);
    const fingerprintsBitsets: (DG.BitSet | null)[] = [];
    for (let i = 0; i < fingerprints.length; ++i) {
      const fingerprint = fingerprints[i] ? DG.BitSet.fromBytes(fingerprints[i]!.getRawData().buffer, fingerprints[i]!.length) : null;
      fingerprintsBitsets.push(fingerprint);
    }
    return DG.Column.fromList('object', 'fingerprints', fingerprintsBitsets);
  } catch (e: any) {
    console.error('Chem | Catch in getMorganFingerprints: ' + e.toString());
    throw e;
  }
}

//name: getMorganFingerprint
//input: string molString {semType: Molecule}
//output: object fingerprintBitset [Fingerprints]
export function getMorganFingerprint(molString: string): DG.BitSet {
  const bitArray = chemSearches.chemGetFingerprint(molString, Fingerprint.Morgan);
  return DG.BitSet.fromBytes(bitArray.getRawData().buffer, bitArray.length);
}

//name: getSimilarities
//input: column molStringsColumn
//input: string molString
//output: dataframe result
export async function getSimilarities(molStringsColumn: DG.Column, molString: string): Promise<DG.DataFrame> {
  try {
    const result = await chemSearches.chemGetSimilarities(molStringsColumn, molString);
    return result ? DG.DataFrame.fromColumns([result]) : DG.DataFrame.create();
  } catch (e: any) {
    console.error('Chem | Catch in getSimilarities: ' + e.toString());
    throw e;
  }
}

//name: getDiversities
//input: column molStringsColumn
//input: int limit
//output: dataframe result
export async function getDiversities(molStringsColumn: DG.Column, limit: number = Number.MAX_VALUE): Promise<DG.DataFrame> {
  try {
    const result = await chemSearches.chemGetDiversities(molStringsColumn, limit);
    return result ? DG.DataFrame.fromColumns([result]) : DG.DataFrame.create();
  } catch (e: any) {
    console.error('Chem | Catch in getDiversities: ' + e.toString());
    throw e;
  }
}

//name: findSimilar
//input: column molStringsColumn
//input: string molString
//input: int limit
//input: int cutoff
//output: dataframe result
export async function findSimilar(molStringsColumn: DG.Column, molString: string, limit: number = Number.MAX_VALUE, cutoff: number = 0.0)
  : Promise<DG.DataFrame> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(molString, 'molString');
  assure.notNull(limit, 'limit');
  assure.notNull(cutoff, 'cutoff');

  try {
    const result = await chemSearches.chemFindSimilar(molStringsColumn, molString, {limit: limit, cutoff: cutoff});
    return result ? result : DG.DataFrame.create();
  } catch (e: any) {
    console.error('Chem | In findSimilar: ' + e.toString());
    throw e;
  }
}

//name: searchSubstructure
//input: column molStringsColumn
//input: string molString
//input: string molBlockFailover
//output: column result
export async function searchSubstructure(
  molStringsColumn: DG.Column, molString: string,
  molBlockFailover: string): Promise<DG.Column<any>> {

  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(molString, 'molString');
  assure.notNull(molBlockFailover, 'molBlockFailover');

  try {
    const result = await chemSearches.chemSubstructureSearchLibrary(molStringsColumn, molString, molBlockFailover);
    return DG.Column.fromList('object', 'bitset', [result]); // TODO: should return a bitset itself
  } catch (e: any) {
    console.error('Chem | In substructureSearch: ' + e.toString());
    throw e;
  }
}

//name: Molecular Descriptors
//tags: app
export function descriptorsApp(): void {
  getDescriptorsApp();
}

//name: saveAsSdf
//description: As SDF
//tags: fileExporter
export async function saveAsSdf(): Promise<void> {
  const progressIndicator = DG.TaskBarProgressIndicator.create('Saving as SDF...');
  saveAsSdfDialog();
  progressIndicator.close();
}

//#region Top menu

//name: Chem Similarity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/chem-similarity-search-viewer.svg
export function similaritySearchViewer(): ChemSimilarityViewer {
  return new ChemSimilarityViewer();
}

//top-menu: Chem | Search | Similarity Search...
//name: Similarity Search
//description: finds the most similar molecule
export function similaritySearchTopMenu(): void {
  (grok.shell.v as DG.TableView).addViewer('Chem Similarity Search');
}

//name: Chem Diversity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/chem-diversity-search-viewer.svg
export function diversitySearchViewer(): ChemDiversityViewer {
  return new ChemDiversityViewer();
}

//top-menu: Chem | Search | Diversity Search...
//name: Diversity Search
//description: finds the most diverse molecules
export function diversitySearchTopMenu(): void {
  (grok.shell.v as DG.TableView).addViewer('Chem Diversity Search');
}


//name: ChemSpaceEditor
//tags: editor
//input: funccall call
export function ChemSpaceEditor(call: DG.FuncCall) {
  const funcEditor = new SequenceSpaceFunctionEditor(DG.SEMTYPE.MOLECULE);
  ui.dialog({title: 'Chemical Space'})
    .add(funcEditor.paramsUI)
    .onOK(async () => {
      call.func.prepare(funcEditor.funcParams).call(true);
    })
    .show();
}


//top-menu: Chem | Analyze Structure | Chemical Space...
//name: Chem Space
//input: dataframe table
//input: column molecules { semType: Molecule }
//input: string methodName { choices:["UMAP", "t-SNE"] }
//input: string similarityMetric { choices:["Tanimoto", "Asymmetric", "Cosine", "Sokal"] }
//input: bool plotEmbeddings = true
//input: object options {optional: true}
//editor: Chem:ChemSpaceEditor
export async function chemSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: string,
  similarityMetric: string = 'Tanimoto', plotEmbeddings: boolean, options?: IUMAPOptions | ITSNEOptions): Promise<DG.Viewer | undefined> {
  if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
    grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
    return;
  }

  const embedColsNames = getEmbeddingColsNames(table);

  const chemSpaceParams = {
    seqCol: molecules,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: [embedColsNames[0], embedColsNames[1]],
    options: options
  };
  const chemSpaceRes = await chemSpace(chemSpaceParams);
  const embeddings = chemSpaceRes.coordinates;

  for (const col of embeddings) {
    table.columns.add(col);
  }
  if (plotEmbeddings)
    return grok.shell
      .tableView(table.name)
      .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Chem space'});
}


//name: Chem Space Embeddings
//input: string col
//input: string methodName
//input: string similarityMetric
//input: string xAxis
//input: string yAxis
//input: object options {optional: true}
//output: object result
export async function getChemSpaceEmbeddings(col: DG.Column, methodName: string,
  similarityMetric: string = 'Tanimoto', xAxis: string, yAxis: string, options?: any): Promise<ISequenceSpaceResult> {
  //need to create dataframe to add fingerprints column
  if (!col.dataFrame) {
    const dfForFp = DG.DataFrame.create(col.length);
    dfForFp.columns.add(col);
  }
  const chemSpaceParams = {
    seqCol: col,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: [xAxis, yAxis],
    options: options
  };
  const chemSpaceRes = await chemSpace(chemSpaceParams);
  return chemSpaceRes;
}

//name: Chem Similarities Matrix
//input: int dim
//input: column col
//input: dataframe df
//input: string colName
//input: object simArr
//output: object res
export async function getChemSimilaritiesMatrix(dim: number, col: DG.Column,
  df: DG.DataFrame, colName: string, simArr: DG.Column[]): Promise<DG.Column[]> {
  //need to create dataframe to add fingerprints column
  if (!col.dataFrame) {
    const dfForFp = DG.DataFrame.create(col.length);
    dfForFp.columns.add(col);
  }
  return await getSimilaritiesMarix(dim, col, df, colName, simArr);
}

//top-menu: Chem | Analyze Structure | Elemental Analysis...
//name: Elemental Analysis
//description: function that implements elemental analysis
//input: dataframe table
//input: column molCol { semType: Molecule }
//input: bool radarView = false
//input: bool radarGrid = false
export function elementalAnalysis(table: DG.DataFrame, molCol: DG.Column, radarView: boolean, radarGrid: boolean): void {
  if (molCol.semType !== DG.SEMTYPE.MOLECULE) {
    grok.shell.info(`The column ${molCol.name} doesn't contain molecules`);
    return;
  }

  const [elements, invalid]: [Map<string, Int32Array>, number[]] = getAtomsColumn(molCol);
  let columnNames: string[] = [];

  if (invalid.filter((el) => el !== null).length > 0) {
    console.log(`Invalid rows ${invalid.map((i) => i.toString()).join(', ')}`);
    grok.shell.warning('Dataset contains malformed data!');
  }

  for (let elName of elementsTable) {
    const value = elements.get(elName);
    if (value) {
      let column = DG.Column.fromInt32Array(elName, value);
      column.name = table.columns.getUnusedName(column.name);
      table.columns.add(column);
      columnNames.push(column.name);
    }
  }

  let view = grok.shell.getTableView(table.name);

  if (radarView) {
    const packageExists = checkPackage('Charts', '_radarViewerDemo');
    if (packageExists) {
      let radarViewer = DG.Viewer.fromType('Radar', table, {
        valuesColumnNames: columnNames,
      });
      view.addViewer(radarViewer);
    } else {
      grok.shell.warning('Charts package is not installed');
    }
  }

  if (radarGrid) {
    const packageExists = checkPackage('PowerGrid', 'radarCellRenderer');
    if (packageExists) {
      let gc = view.grid.columns.add({gridColumnName: 'elementsRadar', cellType: 'radar'});
      gc.settings = {columnNames: Array.from(elements.keys())};
      gc.width = 300;
    } else {
      grok.shell.warning('PowerGrid package is not installed');
    }
  }
}

//name: R-Groups Analysis
//top-menu: Chem | Analyze SAR | R-Groups Analysis...

export function rGroupsAnalysisMenu(): void {
  const col = grok.shell.t.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col === null) {
    grok.shell.error('Current table does not contain molecules');
    return;
  }
  rGroupAnalysis(col);
}

//name: ActivityCliffsEditor
//tags: editor
//input: funccall call
export function ActivityCliffsEditor(call: DG.FuncCall) {
  const funcEditor = new ActivityCliffsFunctionEditor(DG.SEMTYPE.MOLECULE);
  ui.dialog({title: 'Activity Cliffs'})
    .add(funcEditor.paramsUI)
    .onOK(async () => {
      call.func.prepare(funcEditor.funcParams).call(true);
    })
    .show();
}

//top-menu: Chem | Analyze SAR | Activity Cliffs...
//name: Activity Cliffs
//description: detect activity cliffs
//input: dataframe table [Input data table]
//input: column molecules {type:categorical; semType: Molecule}
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE"] }
//input: object options {optional: true}
//editor: Chem:ActivityCliffsEditor
export async function activityCliffs(df: DG.DataFrame, molecules: DG.Column, activities: DG.Column,
  similarity: number, methodName: string, options?: IUMAPOptions | ITSNEOptions) {
  if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
    grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
    return;
  }
  if (activities.type !== DG.TYPE.INT && activities.type !== DG.TYPE.BIG_INT && activities.type !== DG.TYPE.FLOAT) {
    grok.shell.error(`Column ${activities.name} is not numeric`);
    return;
  }
  const axesNames = getEmbeddingColsNames(df);
  if (df.rowCount > 500) {
    ui.dialog().add(ui.divText(`Activity cliffs analysis might take several minutes.
    Do you want to continue?`))
      .onOK(async () => {
        const progressBar = DG.TaskBarProgressIndicator.create(`Activity cliffs running...`);
        await getActivityCliffs(df, molecules, null as any, axesNames, 'Activity cliffs', activities, similarity, 'Tanimoto',
          methodName, DG.SEMTYPE.MOLECULE, {'units': molecules.tags['units']}, chemSpace, getSimilaritiesMarix,
          createTooltipElement, createPropPanelElement, undefined, options);
        progressBar.close();
      })
      .show();
  } else {
    await getActivityCliffs(df, molecules, null as any, axesNames, 'Activity cliffs', activities, similarity, 'Tanimoto',
      methodName, DG.SEMTYPE.MOLECULE, {'units': molecules.tags['units']}, chemSpace, getSimilaritiesMarix,
      createTooltipElement, createPropPanelElement, undefined, options);
  }
}

//top-menu: Chem | Calculate | To InchI...
//name: To InchI
//input: dataframe table [Input data table]
//input: column molecules {type:categorical; semType: Molecule}
export function addInchisTopMenu(table: DG.DataFrame, col: DG.Column): void {
  addInchis(table, col);
}

//top-menu: Chem | Calculate | To InchI Keys...
//name: To InchI Keys
//input: dataframe table [Input data table]
//input: column molecules {type:categorical; semType: Molecule}
export function addInchisKeysTopMenu(table: DG.DataFrame, col: DG.Column): void {
  addInchiKeys(table, col);
}

//#endregion

//#region Molecule column property panel


//name: Chemistry | Rendering
//input: column molColumn {semType: Molecule}
//tags: panel, exclude-actions-panel
//output: widget result
export function molColumnPropertyPanel(molColumn: DG.Column): DG.Widget {
  return getMolColumnPropertyPanel(molColumn);
}

//name: Chemistry | Descriptors
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function descriptorsWidget(smiles: string): DG.Widget {
  return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles) ? getDescriptorsSingle(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}

//name: Biology | Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.
//help-url: /help/domains/chem/info-panels/drug-likeness.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function drugLikeness(smiles: string): DG.Widget {
  return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles) ? drugLikenessWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}


//name: Chemistry | Properties
//description: Basic molecule properties
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function properties(smiles: DG.SemanticValue): Promise<DG.Widget> {
  return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles.value) ? propertiesWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}

//name: Biology | Structural Alerts
//description: Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response
//help-url: /help/domains/chem/info-panels/structural-alerts.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structuralAlerts(smiles: string): Promise<DG.Widget> {
  return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles) ? structuralAlertsWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}


//name: Structure | Identifiers
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function identifiers(smiles: string): Promise<DG.Widget> {
  return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles) ? await identifiersWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}


//name: Structure | 3D Structure
//description: 3D molecule representation
//tags: panel, chem, widgets
//input: string molecule { semType: Molecule }
//output: widget result
export async function structure3D(molecule: string): Promise<DG.Widget> {
  return molecule && !DG.chem.Sketcher.isEmptyMolfile(molecule) ? structure3dWidget(molecule) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}


//name: Structure | 2D Structure
//description: 2D molecule representation
//tags: panel, chem, widgets
//input: string molecule { semType: Molecule }
//output: widget result
export function structure2d(molecule: string): DG.Widget {
  return molecule && !DG.chem.Sketcher.isEmptyMolfile(molecule) ? structure2dWidget(molecule) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}


//name: Biology | Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function toxicity(smiles: string): DG.Widget {
  return smiles  && !DG.chem.Sketcher.isEmptyMolfile(smiles) ? toxicityWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
}


//name: convertMolNotation
//description: RDKit-based conversion for SMILES, SMARTS, InChi, Molfile V2000 and Molfile V3000
//tags: unitConverter
//input: string molecule {semType: Molecule}
//input: string sourceNotation {choices:["smiles", "smarts", "molblock", "inchi", "v3Kmolblock"]}
//input: string targetNotation {choices:["smiles", "smarts", "molblock", "inchi", "v3Kmolblock"]}
//output: string result {semType: Molecule}
export function convertMolNotation(molecule: string, sourceNotation: string, targetNotation: string): string {
  return _convertMolNotation(molecule, sourceNotation, targetNotation, getRdKitModule());
}

//tags: cellEditor
//description: Molecule
//input: grid_cell cell
export async function editMoleculeCell(cell: DG.GridCell): Promise<void> {
  const sketcher = new Sketcher();
  const unit = cell.cell.column.tags[DG.TAGS.UNITS];
  let molecule = cell.cell.value;
  if (unit === DG.chem.Notation.Smiles) {
    //convert to molFile to draw in coordinates similar to dataframe cell
    molecule = convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
  }
  sketcher.setMolecule(molecule);
  const dlg = ui.dialog()
    .add(sketcher)
    .onOK(() => {
      if (unit === DG.chem.Notation.Smiles) {
        //set new cell value only in case smiles has been edited (to avoid undesired molecule orientation change)
        const newValue = sketcher.getSmiles();
        const mol = checkMoleculeValid(cell.cell.value);
        if (!checkMolEqualSmiles(mol, newValue))
          cell.cell.value = newValue;
        mol?.delete();
      } else
        cell.cell.value = sketcher.getMolFile();
      Sketcher.addToCollection(Sketcher.RECENT_KEY, sketcher.getMolFile());
    })
  .show({resizable: true});
  ui.onSizeChanged(dlg.root).subscribe((_) => {
    if(!sketcher.sketcher?.isInitialized)
      return;
    sketcher._autoResized ? sketcher._autoResized = false : sketcher.resize();
  });
}

//name: OpenChemLib
//tags: moleculeSketcher
//output: widget sketcher
export function openChemLibSketcher(): OpenChemLibSketcher {
  return new OpenChemLibSketcher();
}

//name: importSdfs
//description: Opens SDF file
//tags: file-handler
//meta.ext: sdf,mol
//input: list bytes
//output: list tables
export function importSdf(bytes: Uint8Array): DG.DataFrame[] | void {
  try {
    return _importSdf(Uint8Array.from(bytes));
  } catch(e:any){
    grok.shell.warning('file is not supported or malformed');
    grok.shell.error(e);
  }
}

//name: importSmi
//description: Opens smi file
//tags: file-handler
//meta.ext: smi
//input: list bytes
//output: list tables
export function importSmi(bytes: Uint8Array): DG.DataFrame[] | void {
  try {
    return _importSmi(Uint8Array.from(bytes));
  } catch(e:any){
    grok.shell.warning('file is not supported or malformed');
    grok.shell.error(e);
  }
}

//name: importMol2
//description: Opens smi file
//tags: file-handler
//meta.ext: mol2
//input: list bytes
//output: list tables
export function importMol2(bytes: Uint8Array): DG.DataFrame[] | void {
  try {
    return _importTripos(Uint8Array.from(bytes));
  } catch(e:any){
    grok.shell.warning('file is not supported or malformed');
    grok.shell.error(e);
  }
}

//name: importMol
//description: Opens MOL file
//tags: file-handler
//meta.ext: mol
//input: string content
//output: list tables
export function importMol(content: string): DG.DataFrame[] | void {
  try {
    const molCol = DG.Column.string('molecule', 1).init((_) => content);
    return [DG.DataFrame.fromColumns([molCol])];
  } catch(e:any){
    grok.shell.warning('file is not supported or malformed');
    grok.shell.error(e);
  }
}

//name: oclCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: OpenChemLib
export async function oclCellRenderer(): Promise<OCLCellRenderer> {
  return new OCLCellRenderer();
}

//name: Sort by similarity
//description: Sorts a molecular column by similarity
//tags: exclude-actions-panel
//meta.action: Sort by similarity
//input: semantic_value value { semType: Molecule }
export async function sortBySimilarity(value: DG.SemanticValue): Promise<void> {
  const molCol = value.cell.column;
  const tableRowIdx = value.cell.rowIndex;
  const dframe = molCol.dataFrame;
  const smiles = molCol.get(tableRowIdx);

  const grid = value.viewer as DG.Grid;
  ui.setUpdateIndicator(grid.root, true);
  const progressBar = DG.TaskBarProgressIndicator.create('Sorting Structures...');
  progressBar.update(0, 'Installing ScaffoldGraph..: 0% completed');
  const fingerprints : DG.DataFrame = await callChemSimilaritySearch(dframe, molCol, smiles, 'Tanimoto', 1000000, 0.0, Fingerprint.Morgan);
  ui.setUpdateIndicator(grid.root, false);
  progressBar.update(100, 'Sort completed');
  progressBar.close();

  const idxCol = fingerprints.columns.byName('indexes');
  grid.sort([], []);
  grid.setRowOrder(idxCol.toList());
  grid.props.pinnedRows = [tableRowIdx];
  grid.scrollToPixels(0,0); //to address the bug in the core
}

//name: Use as filter
//description: Adds this structure as a substructure filter
//tags: exclude-actions-panel
//meta.action: Use as filter
//input: semantic_value value { semType: Molecule }
export function useAsSubstructureFilter(value: DG.SemanticValue): void {
  const tv = grok.shell.tv;
  if (tv == null)
    throw 'Requires an open table view.';

  const molCol = value.cell.column;
  const molecule = value.value;
  if (molCol == null)
    throw 'Molecule column not found.';

  let molblock;

  //in case molecule is smiles setting correct coordinates to save molecule orientation in filter
  if (value.cell.column.tags[DG.TAGS.UNITS] == DG.chem.Notation.Smiles)
    molblock = convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
  else
    molblock = molToMolblock(molecule, getRdKitModule());

  tv.getFiltersGroup({createDefaultFilters: false}).add({
    type: DG.FILTER_TYPE.SUBSTRUCTURE,
    column: molCol.name,
    columnName: molCol.name,
    molBlock: molblock,
  });
}

//name: Copy as SMILES
//description: Copies structure as smiles
//tags: exclude-actions-panel
//meta.action: Copy as SMILES
//input: semantic_value value { semType: Molecule }
export function copyAsSmiles(value: DG.SemanticValue): void {
  const smiles = !DG.chem.isMolBlock(value.value) && !isSmarts(value.value) ? value.value :
    _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, getRdKitModule());
  navigator.clipboard.writeText(smiles);
  grok.shell.info('Smiles copied to clipboard');
}

//name: Copy as MOLFILE V2000
//description: Copies structure as molfile V2000
//tags: exclude-actions-panel
//meta.action: Copy as MOLFILE V2000
//input: semantic_value value { semType: Molecule }
export function copyAsMolfileV2000(value: DG.SemanticValue): void {
  const molfileV2000 = DG.chem.isMolBlock(value.value) && !value.value.includes('V3000') ? value.value :
    _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, getRdKitModule());
  navigator.clipboard.writeText(molfileV2000);
  grok.shell.info('Molfile V2000 copied to clipboard');
}


//name: Copy as MOLFILE V3000
//description: Copies structure as molfile V3000
//tags: exclude-actions-panel
//meta.action: Copy as MOLFILE V3000
//input: semantic_value value { semType: Molecule }
export function copyAsMolfileV3000(value: DG.SemanticValue): void {
  const molfileV3000 = DG.chem.isMolBlock(value.value) && value.value.includes('V3000') ? value.value :
    _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.V3KMolBlock, getRdKitModule());
  navigator.clipboard.writeText(molfileV3000);
  grok.shell.info('Molfile V3000 copied to clipboard');
}

//name: Copy as SMARTS
//description: Copies structure as smarts
//tags: exclude-actions-panel
//meta.action: Copy as SMARTS
//input: semantic_value value { semType: Molecule }
export function copyAsSmarts(value: DG.SemanticValue): void {
  const smarts = !DG.chem.isMolBlock(value.value) && isSmarts(value.value) ? value.value :
    _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.Smarts, getRdKitModule());
  navigator.clipboard.writeText(smarts);
  grok.shell.info('Smarts copied to clipboard');
}


//name: isSmiles
//input: string s
//output: bool res
export function isSmiles(s: string) : boolean {
  const ctx: IMolContext = getMolSafe(s, {}, _rdKitModule, true);
  if (ctx.mol !== null) {
    ctx.mol.delete();
    return true;
  }
 return false;
}

//name: detectSmiles
//input: column col
//input: int min
export function detectSmiles(col: DG.Column, min: number) : void {
  if (DG.Detector.sampleCategories(col, isSmiles, min, 10, 0.8)) {
    col.tags[DG.TAGS.UNITS] = DG.UNITS.Molecule.SMILES;
    col.semType = DG.SEMTYPE.MOLECULE;
  }
}

//name: chemSimilaritySearch
//input: dataframe df
//input: column col
//input: string molecule
//input: string metricName
//input: int limit
//input: double minScore
//input: string fingerprint
//output: dataframe result
export async function callChemSimilaritySearch(
  df: DG.DataFrame,
  col: DG.Column,
  molecule: string,
  metricName: string,
  limit: number,
  minScore: number,
  fingerprint: string): Promise<DG.DataFrame> {
  return await chemSimilaritySearch(df, col, molecule, metricName, limit, minScore, fingerprint as Fingerprint);
}


//name: chemDiversitySearch
//input: column col
//input: string metricName
//input: int limit
//input: string fingerprint
//output: dataframe result
export async function callChemDiversitySearch(
  col: DG.Column,
  metricName: string,
  limit: number,
  fingerprint: string): Promise<number[]> {
  return await chemDiversitySearch(col, similarityMetric[metricName], limit, fingerprint as Fingerprint);
}


//top-menu: Chem | Analyze Structure | Scaffold Tree
//name: addScaffoldTree
export function addScaffoldTree(): void {
  grok.shell.tv.addViewer(ScaffoldTreeViewer.TYPE);
}


//name: getScaffoldTree
//input: dataframe data
//input: int ringCutoff = 10 [Ignore molecules with # rings > N]
//input: bool dischargeAndDeradicalize = false [Remove charges and radicals from scaffolds]
//output: string result
export async function getScaffoldTree(data: DG.DataFrame,
                                      ringCutoff: number = 0,
                                      dischargeAndDeradicalize: boolean = false
                                      ): Promise<string> {
  const molColumn = data.columns.bySemType(DG.SEMTYPE.MOLECULE);
  const invalid: number[] = new Array<number>(data.columns.length);
  const smiles = molColumn?.getTag(DG.TAGS.UNITS) === DG.UNITS.Molecule.SMILES;
  const smilesList: string[] = new Array<string>(data.columns.length);
  for (let rowI = 0; rowI < molColumn!.length; rowI++) {
    let el: string = molColumn?.get(rowI);
    if (!smiles)
      try {
        el = convertMolNotation(el, DG.UNITS.Molecule.MOLBLOCK, DG.UNITS.Molecule.SMILES);
      }
      catch {
        invalid[rowI] = rowI;
      }

    smilesList[rowI] = el;
  }
  const smilesColumn: DG.Column = DG.Column.fromStrings('smiles', smilesList);
  smilesColumn.name = data.columns.getUnusedName(smilesColumn.name);
  data.columns.add(smilesColumn);
  const scriptRes = await generateScaffoldTree(data, smilesColumn!.name, ringCutoff, dischargeAndDeradicalize);
  return scriptRes;
}


//name: filterMoleculeDuplicates
//input: list molecules
//input: string molecule
//output: list result
export function removeDuplicates(molecules: string[], molecule: string): string[] {
  const mol1 = checkMoleculeValid(molecule);
  if (!mol1) throw (`Molecule is possibly malformed`);
  const filteredMolecules = molecules.filter((smiles) => !checkMolEqualSmiles(mol1, smiles));
  mol1.delete();
  return filteredMolecules;
}

//name: demoCh00
//meta.demoPath: Cheminformatics | Overview
export async function demoCh00(): Promise<void> {
  const table = DG.DataFrame.fromCsv(await _package.files.readAsText('sar-small.csv'));
  grok.shell.windows.showProperties = true;
  const sketcherType = DG.chem.currentSketcherType;
  DG.chem.currentSketcherType = 'OpenChemLib';

  //1. Open table
  const tv = grok.shell.addTableView(table);
  await delay(2000);
  const propPanel = document.getElementsByClassName('grok-entity-prop-panel')[0];
  closeAllAccordionPanes(propPanel!);

   //2. Show some molecule properties on context pane
   grok.shell.windows.showHelp = false;
   const structurePaneContent = getAccordionPane('Structure', propPanel!);
   getAccordionPane('3D Structure', structurePaneContent!);
   const biologyPaneContent = getAccordionPane('Biology', propPanel!);
   getAccordionPane('Toxicity', biologyPaneContent!);
   getAccordionPane('Drug Likeness', biologyPaneContent!);
   await delay(3000);
   table.currentRowIdx = 3;
   await delay(3000);
   table.currentRowIdx = 5;
   await delay(3000);
   grok.shell.windows.showHelp = true;

   //3. Scroll a little bit
  const canvas = tv.grid.root.getElementsByTagName('canvas')[2];
  await scrollTable(canvas, 300, 15, 100);
  await delay(2000);

  //4. Filtering
  const filters = tv.getFiltersGroup();
  await delay(1000);
  const sketcherDlg = await openSketcher(filters.root);
  const sketcherInput = sketcherDlg!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;
  sketcherInput.value = 'C1CCCCC1';
  await delay(1000);
  sketcherInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
  await delay(1000);
  await scrollTable(canvas, 200, 10, 200);
  await delay(1000);
  Array.from(sketcherDlg!.getElementsByTagName('span')).find(el => el.textContent === 'CANCEL')?.click();
  delay(500);
  filters.close();
  delay(1000);

  //5.Align to scaffold
  grok.shell.o = tv.dataFrame.col('smiles');
  await delay(2000);
  closeAllAccordionPanes(propPanel!);
  const chemistryPaneContent = getAccordionPane('Chemistry', propPanel!);
  const renderingPaneContent = getAccordionPane('Rendering', chemistryPaneContent!) as HTMLElement;
  await delay(1000);
  const scaffoldSketcher = await openSketcher(renderingPaneContent);
  const scaffoldSketcherInput = scaffoldSketcher!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;

  let dT = null;
  try {dT = new DataTransfer();} catch (e) { }
  const evt = new ClipboardEvent('paste', {clipboardData: dT});
  evt.clipboardData!.setData('text/plain', demoScaffold);
  scaffoldSketcherInput.value = demoScaffold;
  await delay(100);
  scaffoldSketcherInput.dispatchEvent(evt);
  await delay(1000);
  await scrollTable(canvas, 200, 10, 200);
  await delay(1000);
  Array.from(scaffoldSketcher!.getElementsByTagName('span')).find(el => el.textContent === 'CANCEL')?.click();

  DG.chem.currentSketcherType = sketcherType;

}


//name: demoCh01
//meta.demoPath: Cheminformatics | Similarity Search
export async function demoCh01(): Promise<void> {
  const table = DG.DataFrame.fromCsv(await _package.files.readAsText('sar-small.csv'));
  grok.shell.windows.showProperties = true;
  grok.shell.windows.showHelp = true;

  //1. Open table
  const tv = grok.shell.addTableView(table);
  await delay(2000);

  //2. Add similarity search viewer
  const similarityViewer = tv.addViewer('Chem Similarity Search');
  grok.shell.o = similarityViewer;
  await delay(2000);

  //3. Change target molecule
  table.currentRowIdx = 1;
  await delay(1000);
  table.currentRowIdx = 2;
  await delay(1000);
  table.currentRowIdx = 3;

}


//name: demoCh02
//meta.demoPath: Cheminformatics | R Group Analysis
export async function demoCh02(): Promise<void> {
  const table = DG.DataFrame.fromCsv(await _package.files.readAsText('sar-small.csv'));
  grok.shell.windows.showProperties = true;

  //1. Open table
  const tv = grok.shell.addTableView(table);
  await delay(2000);

  //2. Open RGroup analysis sketcher and run analysis
  rGroupAnalysis(table.col('smiles')!);
  await delay(2000);
  const sketcher = document.getElementsByClassName('d4-dialog')[0];
  const sketcherInput = sketcher!.getElementsByClassName('grok-sketcher-input')[0]?.children[0] as HTMLInputElement;
  sketcherInput.value = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
  await delay(1000);
  sketcherInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
  await delay(1000);
  Array.from(sketcher!.getElementsByTagName('span')).find(el => el.textContent === 'OK')?.click();
}


//name: demoCh03
//meta.demoPath: Cheminformatics | Activity Cliffs
export async function demoCh03(): Promise<void> {
  const table = DG.DataFrame.fromCsv(await _package.files.readAsText('activity_cliffs.csv'));
  grok.shell.windows.showProperties = true;

  //1. Open table
  const tv = grok.shell.addTableView(table);
  await delay(1000);

  //2. Run Activity cliffs
  const molecules = table.col('smiles')!
  const progressBar = DG.TaskBarProgressIndicator.create(`Activity cliffs running...`);
  const axesNames = getEmbeddingColsNames(table); //@ts-ignore
  const scatterPlot = await getActivityCliffs(table, molecules, null as any, axesNames, 'Activity cliffs', table.col('Activity')!, 80, 'Tanimoto',
    't-SNE', DG.SEMTYPE.MOLECULE, {'units': molecules.tags['units']}, chemSpace, getSimilaritiesMarix,
    createTooltipElement, createPropPanelElement, undefined, undefined, 0.5);
  progressBar.close();
  await delay(1000);

  //3. Open cliffs table
  (Array.from(scatterPlot!.root.children)
    .filter((it) => it.className === 'ui-btn ui-btn-ok scatter_plot_link cliffs_grid')[0] as HTMLElement).click();
  await delay(1000);

  //4. Select some random cliffs
  let cliffsGrid: DG.Viewer | null = null;
  for (const i of tv.viewers) {
    if (i.dataFrame.name === `${CLIFFS_DF_NAME}${activityCliffsIdx}`)
      cliffsGrid = i;
  }
  cliffsGrid!.dataFrame.currentRowIdx = 0;
  await delay(2000);
  cliffsGrid!.dataFrame.currentRowIdx = 1;
  await delay(2000);
  cliffsGrid!.dataFrame.currentRowIdx = 2;
}


//name: demoCh04
export async function demoCh04(): Promise<void> { //Databases integration in property panel
  const table = _importSdf(await _package.files.readAsBytes('mol1K.sdf'))[0];
  grok.shell.windows.showProperties = true;
  grok.shell.windows.showHelp = false;

  //1. Open table
  const tv = grok.shell.addTableView(table);
  await delay(2000);

  //2. Open tabs on property panel
  const propPanel = document.getElementsByClassName('grok-entity-prop-panel')[0];
  closeAllAccordionPanes(propPanel!);
  const databasesPaneContent = getAccordionPane('Databases', propPanel!);
  getAccordionPane('ChEMBL (Internal) Substructure Search', databasesPaneContent!) as HTMLElement;
  getAccordionPane('ChEMBL (Internal) Similarity Search', databasesPaneContent!) as HTMLElement;
  await delay(3000);
  table.currentRowIdx = 2;
  await delay(3000);
  table.currentRowIdx = 5;
}


//name: demoCh05
//meta.demoPath: Cheminformatics | Databases
export async function demoCh05(): Promise<void> {
  const ids = ['CHEMBL1827', 'CHEMBL1829', 'CHEMBL1830'];
  const query = `SELECT m.chembl_id AS compound_chembl_id, s.canonical_smiles, act.standard_type, act.standard_value
  FROM compound_structures s, molecule_dictionary m, compound_records r, docs d, activities act, assays a, target_dictionary t
  WHERE s.molregno     = m.molregno
  AND m.molregno       = r.molregno
  AND r.record_id      = act.record_id
  AND r.doc_id         = d.doc_id
  AND act.assay_id     = a.assay_id
  AND a.tid            = t.tid
  AND act.standard_type = 'IC50'
  AND t.chembl_id      = '~id~';`

  const connection = await grok.functions.eval('Chembl:Chembl');
  const queryPanel = ui.box();
  const gridDiv = ui.div();
  const scatterPlot = ui.div();
  const barchart = ui.div();

  const totalDiv = ui.splitV([
    queryPanel,
    gridDiv,
    ui.splitH([
      scatterPlot,
      barchart
    ])
  ], {style: {height: '100%', width: '100%'}});

  const loading = (isLoading: boolean) => {
    ui.setUpdateIndicator(gridDiv, isLoading);
    ui.setUpdateIndicator(scatterPlot, isLoading);
    ui.setUpdateIndicator(barchart, isLoading);
  }

  const loadNewQuery = (id: string) => {
    ui.empty(queryPanel);
    const queryDiv = ui.textInput(`Compound activity details for target = ${id}`, query.replace(`~id~`, id));
    queryDiv.input.style.height = '100%';
    queryPanel.append(queryDiv.root);
    const dBQuery = connection!.query('', query.replace(`~id~`, id));
    dBQuery.adHoc = true;
    ui.empty(gridDiv);
    ui.empty(scatterPlot);
    ui.empty(barchart);
    loading(true);
    return dBQuery;
  }

  const loadQueryResults = async (t: DG.DataFrame) => {
    await grok.data.detectSemanticTypes(t);
    const grid = t.plot.grid().root;
    grid.style.width = '100%';
    loading(false);
    gridDiv.append(grid);
    barchart.append(t.plot.bar().root);
    scatterPlot.append(t.plot.scatter().root);
    await delay(1500);
  }


  grok.shell.newView('Databases', [totalDiv]);

  setTimeout(async () => {
    for (const id of ids) {
      const t = await loadNewQuery(id).executeTable();
      await loadQueryResults(t);
    }
  }, 500);

}
