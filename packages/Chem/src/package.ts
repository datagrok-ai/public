import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getMolColumnPropertyPanel} from './panels/chem-column-property-panel';
import * as chemSearches from './chem-searches';
import {SubstructureFilter} from './widgets/chem-substructure-filter';
import {GridCellRendererProxy, RDKitCellRenderer} from './rendering/rdkit-cell-renderer';
import {drugLikenessWidget} from './widgets/drug-likeness';
import {molfileWidget} from './widgets/molfile';
import {propertiesWidget} from './widgets/properties';
import {structuralAlertsWidget} from './widgets/structural-alerts';
import {structure2dWidget} from './widgets/structure2d';
import {structure3dWidget} from './widgets/structure3d';
import {toxicityWidget} from './widgets/toxicity';
import {chemSpace, getEmbeddingColsNames} from './analysis/chem-space';
import {getDescriptorsApp, getDescriptorsSingle} from './descriptors/descriptors-calculation';
import {addInchiKeys, addInchis} from './panels/inchi';
import {addMcs} from './panels/find-mcs';
import * as chemCommonRdKit from './utils/chem-common-rdkit';
import {_rdKitModule} from './utils/chem-common-rdkit';
import {rGroupAnalysis} from './analysis/r-group-analysis';
import {identifiersWidget} from './widgets/identifiers';
import {_convertMolNotation, isMolBlock, MolNotation} from './utils/convert-notation-utils';
import '../css/chem.css';
import {chemSimilaritySearch, ChemSimilarityViewer} from './analysis/chem-similarity-viewer';
import {chemDiversitySearch, ChemDiversityViewer} from './analysis/chem-diversity-viewer';
import {saveAsSdfDialog} from './utils/sdf-utils';
import {Fingerprint} from './utils/chem-common';
import {assure} from '@datagrok-libraries/utils/src/test';
import {OpenChemLibSketcher} from './open-chem/ocl-sketcher';
import {_importSdf} from './open-chem/sdf-importer';
import {OCLCellRenderer} from './open-chem/ocl-cell-renderer';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import Sketcher = DG.chem.Sketcher;
import {getActivityCliffs, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils';
import {checkForStructuralAlerts} from './panels/structural-alerts';
import {createPropPanelElement, createTooltipElement} from './analysis/activity-cliffs';
import {getAtomsColumn, checkPackage} from './utils/elemental-analysis-utils';
import {elementsTable} from './constants';
import {getSimilaritiesMarix} from './utils/similarity-utils';
import {molToMolblock} from './utils/convert-notation-utils';
import {similarityMetric} from '@datagrok-libraries/utils/src/similarity-metrics';
import {_importSmi} from './file-importers/smi-importer';
import {scaffoldTreeGeneration} from './scripts-api';
import { gasteigerChargesWidget } from './widgets/gasteiger-charges';

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
  _renderers = new Map();
}

//tags: autostart
export async function initChemAutostart(): Promise<void> { }

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
    const fingerprintsBitsets: DG.BitSet[] = [];
    for (let i = 0; i < fingerprints.length; ++i) {
      const fingerprint = DG.BitSet.fromBytes(fingerprints[i].getRawData().buffer, fingerprints[i].length);
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

//name: Descriptors App
//tags: app
export function descriptorsApp(): void {
  getDescriptorsApp();
}

//name: saveAsSdf
//description: Save as SDF
//tags: fileExporter
export function saveAsSdf(): void {
  saveAsSdfDialog();
}

//#region Top menu

//top-menu: Chem | Activity Cliffs...
//name: Activity Cliffs
//description: detect activity cliffs
//input: dataframe table [Input data table]
//input: column molecules {type:categorical; semType: Molecule}
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
export async function activityCliffs(df: DG.DataFrame, molecules: DG.Column, activities: DG.Column,
  similarity: number, methodName: string): Promise<DG.Viewer | undefined> {

  if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
    grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
    return;
  }

  const axesNames = getEmbeddingColsNames(df);
  const options: { [key: string]: any } = {
    'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005},
  };
  return await getActivityCliffs(
    df,
    molecules,
    null as any,
    axesNames,
    'Activity cliffs',
    activities,
    similarity,
    'Tanimoto',
    methodName,
    DG.SEMTYPE.MOLECULE,
    {'units': molecules.tags['units']},
    chemSpace,
    getSimilaritiesMarix,
    createTooltipElement,
    createPropPanelElement,
    options[methodName]);
}

//top-menu: Chem | Chemical Space...
//name: Chem Space
//input: dataframe table
//input: column molecules { semType: Molecule }
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
//input: string similarityMetric { choices:["Tanimoto", "Asymmetric", "Cosine", "Sokal"] }
//input: bool plotEmbeddings = true
export async function chemSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: string,
  similarityMetric: string = 'Tanimoto', plotEmbeddings: boolean): Promise<DG.Viewer | undefined> {

  if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
    grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
    return;
  }

  const embedColsNames = getEmbeddingColsNames(table);

  // dimensionality reducing algorithm doesn't handle empty values correctly so remove empty values at this step
  const withoutEmptyValues = DG.DataFrame.fromColumns([molecules]).clone();
  const emptyValsIdxs = removeEmptyStringRows(withoutEmptyValues, molecules);
  const chemSpaceParams = {
    seqCol: withoutEmptyValues.col(molecules.name)!,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: [embedColsNames[0], embedColsNames[1]],
  };

  const chemSpaceRes = await chemSpace(chemSpaceParams);
  const embeddings = chemSpaceRes.coordinates;

  //inserting empty values back into results
  for (const col of embeddings) {
    const listValues = col.toList();
    emptyValsIdxs.forEach((ind: number) => listValues.splice(ind, 0, null));
    table.columns.add(DG.Column.float(col.name, table.rowCount).init((i)=> listValues[i]));
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
//output: object result
export async function getChemSpaceEmbeddings(col: DG.Column, methodName: string,
  similarityMetric: string = 'Tanimoto', xAxis: string, yAxis: string): Promise<ISequenceSpaceResult> {
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

//name: R-Groups Analysis
//top-menu: Chem | R-Groups Analysis...

export function rGroupsAnalysisMenu(): void {
  const col = grok.shell.t.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col === null) {
    grok.shell.error('Current table does not contain molecules');
    return;
  }
  rGroupAnalysis(col);
}

//top-menu: Chem | Substituent Analysis...
//name: substituentAnalysis
//input: dataframe table
//output: viewer result
export function substituentAnalysisMenu(table: DG.DataFrame): void {
  const packageExists = checkPackage('Charts', '_SubstituentAnalysisViewer');
  if (packageExists) {
    const substituentAnalysisViewer = DG.Viewer.fromType('SubstituentAnalysisViewer', table, {});
    grok.shell.tv.addViewer(substituentAnalysisViewer);
  } else {
    grok.shell.warning('Charts package is not installed');
  }
}


//#endregion

//#region Molecule column property panel

//name: Chem | Find MCS
//friendly-name: Chem | Find MCS
//tags: panel, chem
//input: column col {semType: Molecule}
export function addMcsPanel(col: DG.Column): void {
  addMcs(col);
}

//name: Chem | To InchI
//friendly-name: Chem | To InchI
//tags: panel, chem
//input: column col {semType: Molecule}
export function addInchisPanel(col: DG.Column): void {
  addInchis(col);
}

//top-menu: Chem | Elemental Analysis...
//name: Elemental analysis
//description: function that implements elemental analysis
//input: dataframe table
//input: column molCol { semType: Molecule }
//input: bool radarView = false
//input: bool radarGrid = false
export function elementalAnalysis(table: DG.DataFrame, molCol: DG.Column, radarView: boolean, radarGrid: boolean): void {
  const [elements, invalid]: [Map<string, Int32Array>, number[]] = getAtomsColumn(molCol);
  let columnNames: string[] = [];

  if (invalid.length > 0)
    console.log(`Invalid rows ${invalid.map((i) => i.toString()).join(', ')}`);

  for (let elName of elementsTable) {
    const value = elements.get(elName);
    if (value) {
      if (table.columns.contains(elName)) {
        break;
      } else {
        let column = DG.Column.fromInt32Array(elName, value);
        table.columns.add(column);
        columnNames.push(elName);
      }
    }
  }

  let view = grok.shell.getTableView(table.name);

  if (radarView) {
    const packageExists = checkPackage('Charts', 'radarViewerDemo');
    if (packageExists) {
      let radarViewer = DG.Viewer.fromType('RadarViewer', table, {
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


//name: Chem | To InchI Keys
//friendly-name: Chem | To InchI Keys
//tags: panel, chem
//input: column col {semType: Molecule}
export function addInchisKeysPanel(col: DG.Column): void {
  addInchiKeys(col);
}


//name: Chem
//input: column molColumn {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molColumn: DG.Column): DG.Widget {
  return getMolColumnPropertyPanel(molColumn);
}

//name: Chem Descriptors
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function descriptorsWidget(smiles: string): DG.Widget {
  return getDescriptorsSingle(smiles);
}

//name: Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.
//help-url: /help/domains/chem/info-panels/drug-likeness.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function drugLikeness(smiles: string): DG.Widget {
  return smiles ? drugLikenessWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Molfile
//description: Molecule as Molfile
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function molfile(smiles: string): DG.Widget {
  return smiles ? molfileWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Properties
//description: Basic molecule properties
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function properties(smiles: DG.SemanticValue): Promise<DG.Widget> {
  return smiles ? propertiesWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Structural Alerts
//description: Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response
//help-url: /help/domains/chem/info-panels/structural-alerts.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structuralAlerts(smiles: string): Promise<DG.Widget> {
  return smiles ? structuralAlertsWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Structure 2D
//description: 2D molecule representation
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function structure2d(smiles: string): DG.Widget {
  return smiles ? structure2dWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Structure 3D
//description: 3D molecule representation
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function structure3d(smiles: string): Promise<DG.Widget> {
  if (isMolBlock(smiles)) {
    const mol = getRdKitModule().get_mol(smiles);
    smiles = mol.get_smiles();
    mol?.delete();
  }

  return smiles ? structure3dWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function toxicity(smiles: string): DG.Widget {
  return smiles ? toxicityWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Identifiers
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function identifiers(smiles: string): Promise<DG.Widget> {
  return smiles ? await identifiersWidget(smiles) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Gasteiger Partial Charges
//tags: demo, chem, rdkit, panel, widgets
//description: The Gasteiger partial charges visualization, RDKit based
//input: string molString {semType: Molecule}
//output: widget result
export async function gasteigerPartialCharges(molString: string): Promise<DG.Widget> {
  return gasteigerChargesWidget(molString);
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
  sketcher.setMolecule(cell.cell.value);
  ui.dialog()
    .add(sketcher)
    .onOK(() => {
      cell.cell.value = unit == 'molblock' ? sketcher.getMolFile() : sketcher.getSmiles();
      Sketcher.addRecent(sketcher.getMolFile());
    })
    .show();
}

//name: SimilaritySearchViewer
//tags: viewer
//output: viewer result
export function similaritySearchViewer(): ChemSimilarityViewer {
  return new ChemSimilarityViewer();
}

//top-menu: Chem | Similarity Search...
//name: similaritySearch
//description: finds the most similar molecule
export function similaritySearchTopMenu(): void {
  (grok.shell.v as DG.TableView).addViewer('SimilaritySearchViewer');
}

//name: DiversitySearchViewer
//tags: viewer
//output: viewer result
export function diversitySearchViewer(): ChemDiversityViewer {
  return new ChemDiversityViewer();
}

//top-menu: Chem | Diversity Search...
//name: diversitySearch
//description: finds the most diverse molecules
export function diversitySearchTopMenu() {
  (grok.shell.v as DG.TableView).addViewer('DiversitySearchViewer');
}

//name: inchiToSmiles
//input: string id
//output: string smiles {semType: Molecule}
//meta.role: converter
//meta.inputRegexp: (InChI\=.+)
export async function inchiToSmiles(id: string) {
  const mol = getRdKitModule().get_mol(id);
  return mol.get_smiles();
}

//name: Open Chem Sketcher
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

//name: Use as filter
//description: Adds this structure as a substructure filter
//meta.action: Use as filter
//input: string mol { semType: Molecule }
export function useAsSubstructureFilter(mol: string): void {
  const tv = grok.shell.tv;
  if (tv == null)
    throw 'Requires an open table view.';

  const molCol = tv.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (molCol == null)
    throw 'Molecule column not found.';

  tv.getFiltersGroup({createDefaultFilters: false}).add({
    type: DG.FILTER_TYPE.SUBSTRUCTURE,
    column: molCol.name,
    columnName: molCol.name,
    molBlock: molToMolblock(mol, getRdKitModule())
  });
}

//name: detectSmiles
//input: column col
//input: int min
export function detectSmiles(col: DG.Column, min: number) {
  function isSmiles(s: string) {
    let d: RDMol | null = null;
    try {
      d = _rdKitModule.get_mol(s);
      return true;
    } catch {
      return false;
    } finally {
      d?.delete();
    }
  }

  if (DG.Detector.sampleCategories(col, isSmiles, min, 10, 0.8)) {
    col.tags[DG.TAGS.UNITS] = DG.UNITS.Molecule.SMILES;
    col.semType = DG.SEMTYPE.MOLECULE;
  }
}

//name: Chem | Structural Alerts...
//tags: panel, chem
//input: column col { semType: Molecule }
export async function getStructuralAlerts(col: DG.Column<string>): Promise<void> {
  await checkForStructuralAlerts(col);
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

//name: scaffoldTree
//input: dataframe data
//output: string result
export async function scaffoldTree(data: DG.DataFrame) {
  const smilesColumn = data.columns.bySemType(DG.SEMTYPE.MOLECULE);
  const scriptRes = await scaffoldTreeGeneration(data, smilesColumn!.name, smilesColumn!.name);
  return scriptRes;
}
