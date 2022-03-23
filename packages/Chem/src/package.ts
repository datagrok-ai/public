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
import {chemSpace} from './analysis/chem-space';
import {getActivityCliffs} from './analysis/activity-cliffs';
import {addDescriptors, getDescriptorsApp, getDescriptorsSingle} from './descriptors/descriptors-calculation';
import {addInchiKeys, addInchis} from './panels/inchi';
import {addMcs} from './panels/find-mcs';
import * as chemCommonRdKit from './utils/chem-common-rdkit';
import {rGroupAnalysis} from './analysis/r-group-analysis';
import {identifiersWidget} from './widgets/identifiers';
import {convertMoleculeImpl, isMolBlock} from './utils/chem-utils';
import '../css/chem.css';
import {ChemSimilarityViewer} from './analysis/chem-similarity-viewer';
import {ChemDiversityViewer} from './analysis/chem-diversity-viewer';
import {_saveAsSdf} from './utils/sdf-utils';
import {Fingerprint} from './utils/chem-common';
import {assure} from '@datagrok-libraries/utils/src/test';
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;

import {OpenChemLibSketcher} from './open-chem/ocl-sketcher';
import {_importSdf} from './open-chem/sdf-importer';
import {OCLCellRenderer} from './open-chem/ocl-cell-renderer';

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
export async function initChem() {
  chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
  await chemCommonRdKit.initRdKitModuleLocal();
  _properties = await _package.getProperties();
  _rdRenderer = new RDKitCellRenderer(getRdKitModule());
  renderer = new GridCellRendererProxy(_rdRenderer, 'Molecule');
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
  const mol = getRdKitModule().get_mol(smiles);
  const res = JSON.parse(mol.get_descriptors()).CrippenClogP;
  mol?.delete();
  return res;
}

//name: rdKitCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: RDKit
export async function rdKitCellRenderer() {
  return new RDKitCellRenderer(getRdKitModule());
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

  renderer.renderer = renderer ?? renderer.renderer;
  return renderer;
}

//name: getMorganFingerprints
//input: column molColumn {semType: Molecule}
//output: column result [fingerprints]
export async function getMorganFingerprints(molColumn: DG.Column) {
  assure.notNull(molColumn, 'molColumn');

  try {
    const fingerprints = await chemSearches.chemGetFingerprints(molColumn, Fingerprint.Morgan);
    const fingerprintsBitsets: DG.BitSet[] = [];
    for (let i = 0; i < fingerprints.length; ++i) {
      //@ts-ignore
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
export function getMorganFingerprint(molString: string) {
  const bitArray = chemSearches.chemGetFingerprint(molString, Fingerprint.Morgan);
  //@ts-ignore
  return DG.BitSet.fromBytes(bitArray.getRawData(), bitArray.length);
}

//name: getSimilarities
//input: column molStringsColumn
//input: string molString
//output: dataframe result
export async function getSimilarities(molStringsColumn: DG.Column, molString: string) {
  try {
    const result = await chemSearches.chemGetSimilarities(molStringsColumn, molString);
    return result ? DG.DataFrame.fromColumns([result as DG.Column]) : DG.DataFrame.create();
  } catch (e: any) {
    console.error('Chem | Catch in getSimilarities: ' + e.toString());
    throw e;
  }
}

//name: findSimilar
//input: column molStringsColumn
//input: string molString
//input: int limit
//input: int cutoff
//output: dataframe result
export async function findSimilar(molStringsColumn: DG.Column, molString: string, limit: number, cutoff: number) {
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
//input: bool substructLibrary
//input: string molStringSmarts
//output: column result
export async function searchSubstructure(
  molStringsColumn: DG.Column, molString: string,
  substructLibrary: boolean, molStringSmarts: string) {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(molString, 'molString');
  assure.notNull(substructLibrary, 'substructLibrary');
  assure.notNull(molStringSmarts, 'molStringSmarts');

  try {
    const result =
      substructLibrary ?
        await chemSearches.chemSubstructureSearchLibrary(molStringsColumn, molString, molStringSmarts) :
        chemSearches.chemSubstructureSearchGraph(molStringsColumn, molString);
    return DG.Column.fromList('object', 'bitset', [result]); // TODO: should return a bitset itself
  } catch (e: any) {
    console.error('Chem | In substructureSearch: ' + e.toString());
    throw e;
  }
}

//name: Descriptors App
//tags: app
export function descriptorsApp(context: any) {
  getDescriptorsApp();
}

//name: saveAsSdf
//description: Save as SDF
//tags: fileExporter
export function saveAsSdf() {_saveAsSdf();}

//#region Top menu

//top-menu: Chem | Activity Cliffs...
//name: Activity Cliffs
//description: detect activity cliffs
//input: dataframe df [Input data table]
//input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
export async function activityCliffs(df: DG.DataFrame,
  smiles: DG.Column,
  activities: DG.Column,
  similarity: number,
  methodName: string) {
  await getActivityCliffs(df, smiles, activities, similarity, methodName);
}

//top-menu: Chem | Chemical Space...
//name: Chem Space
//input: dataframe table
//input: column smiles { semType: Molecule }
//input: string methodName { choices:["UMAP", "t-SNE", "SPE", "pSPE", "OriginalSPE"] }
//input: string similarityMetric { choices:["Tanimoto", "Asymmetric", "Cosine", "Sokal"] }
//input: bool plotEmbeddings = true
//output: viewer result
export async function chemSpaceTopMenu(table: DG.DataFrame,
  smiles: DG.Column,
  methodName: string,
  similarityMetric: string = 'Tanimoto',
  plotEmbeddings: boolean) {
  return new Promise<void>(async (resolve, reject) => {
    const embeddings = await chemSpace(smiles, methodName, similarityMetric);
    const cols = table.columns as DG.ColumnList;
    for (const col of embeddings)
      cols.add(col);

    const view = grok.shell.addTableView(table);
    if (plotEmbeddings)
      view.scatterPlot({x: 'Embed_X', y: 'Embed_Y'});
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

//#endregion

//#region Molecule column property panel

//name: Chem | Find MCS
//friendly-name: Chem | Find MCS
//tags: panel, chem
//input: column col {semType: Molecule}
export function addMcsPanel(col: DG.Column) {
  addMcs(col);
}

//name: Chem | To InchI
//friendly-name: Chem | To InchI
//tags: panel, chem
//input: column col {semType: Molecule}
export function addInchisPanel(col: DG.Column) {
  addInchis(col);
}

//name: Chem | To InchI Keys
//friendly-name: Chem | To InchI Keys
//tags: panel, chem
//input: column col {semType: Molecule}
export function addInchisKeysPanel(col: DG.Column) {
  addInchiKeys(col);
}

//name: Chem
//input: column molColumn {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molColumn: DG.Column) {
  return getMolColumnPropertyPanel(molColumn);
}

//name: Chem Descriptors
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function descriptorsWidget(smiles: string) {
  return getDescriptorsSingle(smiles);
}

//name: Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.
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
//description: Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response
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
//input: string molecule { semType: Molecule }
//output: widget result
export async function structure3d(molecule: string) {
  if (isMolBlock(molecule)) {
    const mol = getRdKitModule().get_mol(molecule);
    molecule = mol.get_smiles();
    mol?.delete();
  }

  return molecule ? structure3dWidget(molecule) : new DG.Widget(ui.divText('SMILES is empty'));
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
  return convertMoleculeImpl(molecule, from, to, getRdKitModule());
}


//tags: cellEditor
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

//name: SimilaritySearchViewer
//tags: viewer
//output: viewer result
export function similaritySearchViewer() {
  return new ChemSimilarityViewer();
}

//top-menu: Chem | Similarity Search...
//name: similaritySearch
//description: finds the most similar molecule
//output: viewer result
export function similaritySearchTopMenu() {
  (grok.shell.v as DG.TableView).addViewer('SimilaritySearchViewer');
}

//name: DiversitySearchViewer
//tags: viewer
//output: viewer result
export function diversitySearchViewer() {
  return new ChemDiversityViewer();
}

//top-menu: Chem | Diversity Search...
//name: diversitySearch
//description: finds the most diverse molecules
//output: viewer result
export function diversitySearchTopMenu() {
  (grok.shell.v as DG.TableView).addViewer('DiversitySearchViewer');
}

//input: string id
//output: string smiles { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: InChI\=.+
export function inchiToSmiles(id: string) {
  const mol = chemCommonRdKit.getRdKitModule().get_mol(id);
  const smiles = mol.get_smiles();
  mol.delete();
  return smiles;
}

//name: openChemLibSketch
//description: Sketches a molecule
//top-menu: Chem | OpenChemLib Sketch
export function openChemLibSketch() {
  ui.dialog()
    .add(openChemLibSketcher().root)
    .showModal(true);
}

//name: openChemLibSketcher
//tags: moleculeSketcher
//output: widget sketcher
export function openChemLibSketcher() {
  return new OpenChemLibSketcher();
}

//name: importSdfs
//description: Opens SDF file
//tags: file-handler
//meta.ext: sdf
//input: list bytes
//output: list tables
export function importSdf(bytes: Uint8Array) {
  return _importSdf(Uint8Array.from(bytes));
}

//name: oclCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: OpenChemLib
export async function oclCellRenderer() {
  return new OCLCellRenderer();
}
