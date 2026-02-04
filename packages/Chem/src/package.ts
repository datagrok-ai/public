/* eslint-disable camelcase */
/* eslint-disable guard-for-in */
/* eslint-disable max-params */
/* eslint-disable max-len */
/* eslint-disable max-lines */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as chemSearches from './chem-searches';
import {GridCellRendererProxy, RDKitCellRenderer} from './rendering/rdkit-cell-renderer';
import {assure} from '@datagrok-libraries/test/src/test';
import {OpenChemLibSketcher} from './open-chem/ocl-sketcher';
import {_importSdf} from './open-chem/sdf-importer';
import {OCLCellRenderer} from './open-chem/ocl-cell-renderer';
import Sketcher = DG.chem.Sketcher;
import {runActivityCliffs, getActivityCliffsEmbeddings, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {ActivityCliffsEditor as ActivityCliffsFunctionEditor}
  from '@datagrok-libraries/ml/src/functionEditors/activity-cliffs-function-editor';
import {MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT, EMPTY_MOLECULE_MESSAGE,
  SMARTS_MOLECULE_MESSAGE, elementsTable} from './constants';
import {similarityMetric} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {calculateDescriptors, getDescriptorsTree} from './docker/api';
import {addDescriptorsColsToDf, getDescriptorsSingle, getSelected, openDescriptorsDialogDocker} from './descriptors/descriptors-calculation';
import {identifiersWidget, openMapIdentifiersDialog, textToSmiles} from './widgets/identifiers';

//widget imports
import {SubstructureFilter} from './widgets/chem-substructure-filter';
import {drugLikenessWidget} from './widgets/drug-likeness';
import {addPropertiesAsColumns, getChemPropertyFunc, getPropertiesAsColumns, propertiesWidget} from './widgets/properties';
import {structuralAlertsWidget} from './widgets/structural-alerts';
import {structure2dWidget} from './widgets/structure2d';
import {getToxicityRisksColumns, toxicityWidget} from './widgets/toxicity';

//panels imports
import {getInchiKeysImpl, getInchisImpl} from './panels/inchi';
import {getMolColumnPropertyPanel} from './panels/chem-column-property-panel';
import {ScaffoldTreeViewer} from './widgets/scaffold-tree';
import {ScaffoldTreeFilter} from './widgets/scaffold-tree-filter';
import {Fingerprint} from './utils/chem-common';
import * as chemCommonRdKit from './utils/chem-common-rdkit';
import {IMolContext, getMolSafe, isFragment, _isSmarts} from './utils/mol-creation_rdkit';
import {checkMoleculeValid, checkMolEqualSmiles, _rdKitModule} from './utils/chem-common-rdkit';
import {_convertMolNotation, convertNotationForColumn} from './utils/convert-notation-utils';
import {molToMolblock} from './utils/convert-notation-utils';
import {getAtomsColumn, checkPackage} from './utils/elemental-analysis-utils';
import {saveAsSdfDialog} from './utils/sdf-utils';
import {getSimilaritiesMarix} from './utils/similarity-utils';

//analytical imports
import {ActivityCliffsParams, createPropPanelElement, createTooltipElement} from './analysis/activity-cliffs';
import {chemDiversitySearch, ChemDiversityViewer} from './analysis/chem-diversity-viewer';
import {chemSimilaritySearch, ChemSimilarityViewer} from './analysis/chem-similarity-viewer';
import {chemSpace, runChemSpace} from './analysis/chem-space';
import {RGroupDecompRes, RGroupParams, rGroupAnalysis, rGroupDecomp} from './analysis/r-group-analysis';
import {MatchedMolecularPairsViewer} from './analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer';

//file importers
import {_importTripos} from './file-importers/mol2-importer';
import {_importSmi} from './file-importers/smi-importer';

import {Chem} from './scripts-api';
import {renderMolecule} from './rendering/render-molecule';
import {RDKitReactionRenderer} from './rendering/rdkit-reaction-renderer';
import {structure3dWidget} from './widgets/structure3d';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {_demoActivityCliffs, _demoActivityCliffsLayout, _demoChemicalSpace, _demoChemOverview, _demoMMPA,
  _demoRgroupAnalysis, _demoRGroups, _demoScaffoldTree, _demoSimilarityDiversitySearch} from './demo/demo';
import {getStructuralAlertsByRules, RuleSet, STRUCT_ALERTS_RULES_NAMES} from './panels/structural-alerts';
import {getmolColumnHighlights} from './widgets/col-highlights';
import {RDLog, RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {malformedDataWarning} from './utils/malformed-data-utils';
import {DimReductionBaseEditor} from '@datagrok-libraries/ml/src/functionEditors/dimensionality-reduction-editor';
import {getEmbeddingColsNames}
  from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import {ITSNEOptions, IUMAPOptions} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/multi-column-dim-reducer';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {getMCS} from './utils/most-common-subs';
import JSZip from 'jszip';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {fetchWrapper} from '@datagrok-libraries/utils/src/fetch-utils';
import {CHEM_PROP_MAP} from './open-chem/ocl-service/calculations';
import {oclMol} from './utils/chem-common-ocl';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';
import {OCLService} from './open-chem/ocl-service';
import {MpoProfileDialog} from './analysis/mpo';
import {MmmpFunctionEditor, MmpDiffTypes} from './analysis/molecular-matched-pairs/mmp-function-editor';
import {SCALING_METHODS} from './analysis/molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from './analysis/molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {MixtureCellRenderer} from './rendering/mixture-cell-renderer';
import {createComponentPane, createMixtureWidget, Mixfile} from './utils/mixfile';
import {biochemicalPropertiesDialog} from './widgets/biochem-properties-widget';
import {checkCurrentView} from './utils/ui-utils';
import {DesirabilityProfile, mpo, PropertyDesirability, WEIGHTED_AGGREGATIONS_LIST, WeightedAggregation} from '@datagrok-libraries/statistics/src/mpo/mpo';
//@ts-ignore
import '../css/chem.css';
import {addDeprotectedColumn, DeprotectEditor} from './analysis/deprotect';
import {MpoProfilesView} from './mpo/mpo-profiles-view';

import $ from 'cash-dom';
import {MpoProfileCreateView} from './mpo/mpo-create-profile';
import {calculateMpoCore, findSuitableProfiles, loadMpoProfiles, MPO_PROFILE_CHANGED_EVENT, MPO_TEMPLATE_PATH} from './mpo/utils';

export {getMCS};
export * from './package.g';

/** Temporary polyfill */

function getDecoratorFunc() {
  return function(_args: any) {
    return function(
      _target: any,
      _propertyKey: string,
      _descriptor: PropertyDescriptor,
    ) { };
  };
}

// Ensure decorators object exists and polyfill missing decorators
if (!grok.decorators)
  (grok as any).decorators = {};

const decorators = [
  'func', 'init', 'param', 'panel', 'editor', 'demo', 'app',
  'appTreeBrowser', 'fileHandler', 'fileExporter', 'model', 'viewer', 'filter', 'cellRenderer', 'autostart',
  'dashboard', 'folderViewer', 'semTypeDetector', 'packageSettingsEditor', 'functionAnalysis', 'converter',
  'fileViewer', 'model', 'treeBrowser', 'polyfill',
];

decorators.forEach((decorator) => {
  if (!(grok.decorators as any)[decorator])
    (grok.decorators as any)[decorator] = getDecoratorFunc();
});

/** End temporary polyfill */

const drawMoleculeToCanvas = chemCommonRdKit.drawMoleculeToCanvas;
const SKETCHER_FUNCS_FRIENDLY_NAMES: {[key: string]: string} = {
  OpenChemLib: 'OpenChemLib',
  Ketcher: 'Ketcher',
  Marvin: 'Marvin',
  ChemDraw: 'ChemDraw',
};

const PREVIOUS_SKETCHER_NAMES: {[key: string]: string} = {
  'Open Chem Sketcher': 'OpenChemLib',
  'ketcherSketcher': 'Ketcher',
  'Marvin JS': 'Marvin',
  'Chem Draw': 'ChemDraw',
};
let container: DG.DockerContainer;

export const _package: DG.Package = new DG.Package();
export let _properties: any;

let _rdRenderer: RDKitCellRenderer;
export let renderer: GridCellRendererProxy;
let _renderers: Map<string, DG.GridCellRenderer>;
let _initChemPromise: Promise<void> | null = null;

async function initChemInt(): Promise<void> {
  chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
  await chemCommonRdKit.initRdKitModuleLocal();
  _properties = await _package.getProperties();
  _rdRenderer = new RDKitCellRenderer(PackageFunctions.getRdKitModule());
  renderer = new GridCellRendererProxy(_rdRenderer, 'Molecule');
  let storedSketcherType = grok.userSettings.getValue(DG.chem.STORAGE_NAME, DG.chem.KEY) ?? '';
  if (PREVIOUS_SKETCHER_NAMES[storedSketcherType])
    storedSketcherType = PREVIOUS_SKETCHER_NAMES[storedSketcherType];
  if (!storedSketcherType && _properties.Sketcher)
    storedSketcherType = SKETCHER_FUNCS_FRIENDLY_NAMES[_properties.Sketcher];

  const sketcherFunctions = DG.Func.find({meta: {role: DG.FUNC_TYPES.MOLECULE_SKETCHER}});
  const sketcherFunc = sketcherFunctions.find((e) => e.name === storedSketcherType || e.friendlyName === storedSketcherType);
  if (sketcherFunc)
    DG.chem.currentSketcherType = sketcherFunc.friendlyName;
  else {
    if (!!storedSketcherType) {
      grok.shell.warning(
        `Package with ${storedSketcherType} function is not installed.Switching to ${DG.DEFAULT_SKETCHER}.`);
    }

    DG.chem.currentSketcherType = DG.DEFAULT_SKETCHER;
  }
  _renderers = new Map();
}

export class PackageFunctions {
  /**
   * Usage:
   * let a = await grok.functions.call('Chem:getRdKitModule');
   * let b = a.get_mol('C1=CC=CC=C1');
   * alert(b.get_pattern_fp());
   **/

  @grok.decorators.func({outputs: [{type: 'object', name: 'module'}]})
  static getRdKitModule(): RDModule {
    return chemCommonRdKit.getRdKitModule();
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'handler'}]})
  static getMolFileHandler(molString: string): MolfileHandlerBase {
    return MolfileHandler.getInstance(molString);
  }

  @grok.decorators.init()
  static async init(): Promise<void> {
    if (!_initChemPromise)
      _initChemPromise = initChemInt();
    await _initChemPromise;
  }

  @grok.decorators.autostart()
  static async initChemAutostart(): Promise<void> { }

  @grok.decorators.func({'top-menu': 'Chem | Transform | Recalculate Coordinates...', 'name': 'Recalculate Coordinates',
    'description': 'Recalculates 2D coordinates for molecules in the column using Open Chem Lib',
    'meta': {'role': 'transform'},
  })
  static async recalculateCoordsViaOCL(@grok.decorators.param({}) table: DG.DataFrame,
  @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
  @grok.decorators.param({options: {initialValue: 'true'}}) join: boolean = true,
  ): Promise<DG.Column | void> {
    const oclService = new OCLService();
    const newMols = await oclService.recalculateCoordinates(molecules.toList());
    oclService.terminate();
    const newCol = DG.Column.fromStrings(table.columns.getUnusedName(`${molecules.name}_recalcCoords`), newMols);
    newCol.semType = 'Molecule';
    if (!join)
      return newCol;
    table.columns.add(newCol);
  }

  @grok.decorators.func({
    name: 'Chemistry | Most Diverse Structures',
    meta: {role: 'tooltip'},
  })
  static async chemTooltip(
    @grok.decorators.param({options: {semType: 'Molecule'}}) col: DG.Column): Promise<DG.Widget | undefined> {
    const initialWidth = 255;
    const initialHeight = 90;
    const tooltipMaxWidth = 500;
    const version = col.version;
    const colCategories = col.categories;
    for (let i = 0; i < Math.min(colCategories.length, 100); ++i) {
      if (!!colCategories[i] && _isSmarts(colCategories[i]))
        return;
    }

    const divMain = ui.div();
    divMain.append(ui.divText('Most diverse structures', 'chem-tooltip-text'));
    const divStructures = ui.div([ui.loader()]);
    divStructures.classList.add('chem-tooltip-structure-div');
    const getDiverseStructures = async (): Promise<void> => {
      if (col.temp['version'] !== version || col.temp['molIds'].length === 0) {
        const molIds = await chemDiversitySearch(
          col, similarityMetric[BitArrayMetricsNames.Tanimoto], Math.min(6, colCategories.length), Fingerprint.Morgan, DG.BitSet.create(col.length).setAll(true), true);

        Object.assign(col.temp, {
          'version': version,
          'molIds': molIds,
        });
      }
      ui.empty(divStructures);
      const molIdsCached = col.temp['molIds'];
      for (let i = 0; i < molIdsCached.length; ++i)
        divStructures.append(renderMolecule(col.categories[molIdsCached[i]], {width: 75, height: 32}));
    };

    divMain.append(divStructures);
    const widget = new DG.Widget(divMain);

    Object.assign(widget.root.style, {
      position: 'relative',
      width: `${initialWidth}px`,
      height: `${initialHeight}px`,
    });

    const tooltip = document.querySelector('.d4-tooltip');
    if (tooltip) {
      const {width, height} = tooltip.getBoundingClientRect();
      const isWideTooltip = width + initialWidth > tooltipMaxWidth;

      Object.assign(widget.root.style, {
        left: isWideTooltip ? '0' : `${width - widget.root.offsetWidth}px`,
        top: isWideTooltip ? '0' : `${30 - height}px`,
        width: `${isWideTooltip ? initialWidth : initialWidth + width}px`,
        height: isWideTooltip ? `${initialHeight}px` : '0',
      });
    }

    setTimeout(() => getDiverseStructures(), 10);
    return widget;
  }

  @grok.decorators.func({
    name: 'Scaffold Tree',
    meta: {icon: 'files/icons/scaffold-tree-icon.svg', role: 'viewer'},
    outputs: [{name: 'result', type: 'viewer'}],
  })
  static scaffoldTreeViewer() : ScaffoldTreeViewer {
    return new ScaffoldTreeViewer();
  }

  @grok.decorators.func({
    name: 'Substructure Filter',
    description: 'RDKit-based substructure filter',
    meta: {semType: 'Molecule', primaryFilter: 'true', allowMultipleFiltersForColumn: 'false', role: 'filter'},
    outputs: [{name: 'result', type: 'filter'}],
  })
  static substructureFilter(): SubstructureFilter {
    return new SubstructureFilter();
  }

  @grok.decorators.func()
  static canvasMol(
    @grok.decorators.param({type: 'int'}) x: number,
    @grok.decorators.param({type: 'int'}) y: number,
    @grok.decorators.param({type: 'int'}) w: number,
    @grok.decorators.param({type: 'int'}) h: number,
    @grok.decorators.param({type: 'object'}) canvas: HTMLCanvasElement,
    @grok.decorators.param({type: 'string'}) molString: string,
    @grok.decorators.param({type: 'string'}) scaffoldMolString: string | null = null,
    @grok.decorators.param({type: 'object', options: {optional: true}}) options : any = {normalizeDepiction: true, straightenDepiction: true},
  ): void {
    drawMoleculeToCanvas(x, y, w, h, canvas,
      molString, scaffoldMolString == '' ? null : scaffoldMolString,
      options);
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'canvas'}]})
  static drawMolecule(
    @grok.decorators.param({type: 'string'}) molStr: string,
    @grok.decorators.param({type: 'int', options: {optional: true}}) w?: number,
    @grok.decorators.param({type: 'int', options: {optional: true}}) h?: number,
    @grok.decorators.param({type: 'bool', options: {optional: true}}) popupMenu?: boolean,
  ): HTMLElement {
    return renderMolecule(molStr, {width: w, height: h, popupMenu: popupMenu});
  }

  @grok.decorators.func()
  static getCLogP(
    @grok.decorators.param({type: 'string', options: {semType: 'Molecule'}}) smiles: string): number {
    const mol = PackageFunctions.getRdKitModule().get_mol(smiles);
    const res = JSON.parse(mol.get_descriptors()).CrippenClogP;
    mol?.delete();
    return res;
  }

  @grok.decorators.func({
    outputs: [{name: 'result', type: 'grid_cell_renderer'}],
    meta: {chemRendererName: 'RDKit'},
  })
  static async rdKitCellRenderer(): Promise<RDKitCellRenderer> {
    return new RDKitCellRenderer(PackageFunctions.getRdKitModule());
  }

  @grok.decorators.func({
    name: 'chemCellRenderer',
    meta: {'cellType': 'ChemicalReaction', 'role': 'cellRenderer'},
    outputs: [{name: 'result', type: 'grid_cell_renderer'}],
  })
  static async rdKitReactionRenderer(): Promise<RDKitReactionRenderer> {
    return new RDKitReactionRenderer(PackageFunctions.getRdKitModule());
  }

  @grok.decorators.func({
    name: 'chemMixtureRenderer',
    meta: {'cellType': 'ChemicalMixture', 'role': 'cellRenderer'},
    outputs: [{name: 'result', type: 'grid_cell_renderer'}],
  })
  static async rdKitMixtureRenderer(): Promise<MixtureCellRenderer> {
    return new MixtureCellRenderer(PackageFunctions.getRdKitModule());
  }

  @grok.decorators.func({
    name: 'chemCellRenderer',
    meta: {'cellType': 'Molecule', 'role': 'cellRenderer'},
    outputs: [{name: 'result', type: 'grid_cell_renderer'}],
  })
  static async chemCellRenderer(): Promise<DG.GridCellRenderer> {
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

  @grok.decorators.func({
    name: 'getMorganFingerprints',
    meta: {vectorFunc: 'true'},
  })
  static async getMorganFingerprints(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molColumn: DG.Column): Promise<DG.Column> {
    assure.notNull(molColumn, 'molColumn');
    try {
      const fingerprints = await chemSearches.chemGetFingerprints(molColumn, Fingerprint.Morgan, false);
      const fingerprintsBitsets: (DG.BitSet | null)[] = [];
      for (let i = 0; i < fingerprints.length; ++i) {
        const fingerprint = fingerprints[i] ?
          DG.BitSet.fromBytes(fingerprints[i]!.getRawData().buffer as ArrayBuffer, fingerprints[i]!.length) : null;
        fingerprintsBitsets.push(fingerprint);
      }
      return DG.Column.fromList('object', 'fingerprints', fingerprintsBitsets);
    } catch (e: any) {
      console.error('Chem | Catch in getMorganFingerprints: ' + e.toString());
      throw e;
    }
  }

  @grok.decorators.func({
    outputs: [{type: 'object', name: 'fingerprintBitset', options: {description: 'Fingerprints'}}],
  })
  static getMorganFingerprint(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molString: string): DG.BitSet {
    const bitArray = chemSearches.chemGetFingerprint(molString, Fingerprint.Morgan);
    return DG.BitSet.fromBytes(bitArray.getRawData().buffer as ArrayBuffer, bitArray.length);
  }

  @grok.decorators.func()
  static async getSimilarities(
    molStringsColumn: DG.Column,
    molString: string): Promise<DG.DataFrame> {
    try {
      const result = await chemSearches.chemGetSimilarities(molStringsColumn, molString);
      return result ? DG.DataFrame.fromColumns([result]) : DG.DataFrame.create();
    } catch (e: any) {
      console.error('Chem | Catch in getSimilarities: ' + e.toString());
      throw e;
    }
  }

  @grok.decorators.func()
  static async getDiversities(
    molStringsColumn: DG.Column,
    @grok.decorators.param({type: 'int'}) limit: number = Number.MAX_VALUE): Promise<DG.DataFrame> {
    try {
      const result = await chemSearches.chemGetDiversities(molStringsColumn, limit);
      return result ? DG.DataFrame.fromColumns([result]) : DG.DataFrame.create();
    } catch (e: any) {
      console.error('Chem | Catch in getDiversities: ' + e.toString());
      throw e;
    }
  }

  @grok.decorators.func()
  static async findSimilar(
    molStringsColumn: DG.Column,
    molString: string,
    @grok.decorators.param({type: 'int'}) limit: number = Number.MAX_VALUE,
    @grok.decorators.param({type: 'int'}) cutoff: number = 0.0): Promise<DG.DataFrame> {
    assure.notNull(molStringsColumn, 'molStringsColumn');
    assure.notNull(molString, 'molString');
    assure.notNull(limit, 'limit');
    assure.notNull(cutoff, 'cutoff');
    try {
      const result = await chemSearches.chemFindSimilar(molStringsColumn, molString, {limit: limit, minScore: cutoff});
      return result ? result : DG.DataFrame.create();
    } catch (e: any) {
      console.error('Chem | In findSimilar: ' + e.toString());
      throw e;
    }
  }

  @grok.decorators.func()
  static async searchSubstructure(
    molStringsColumn: DG.Column,
    molString: string,
    molBlockFailover: string): Promise<DG.Column<any>> {
    assure.notNull(molStringsColumn, 'molStringsColumn');
    assure.notNull(molString, 'molString');
    assure.notNull(molBlockFailover, 'molBlockFailover');
    try {
      const result = await chemSearches.chemSubstructureSearchLibrary(molStringsColumn, molString, molBlockFailover);
      const resBitset = DG.BitSet.fromBytes(result.buffer.buffer as ArrayBuffer, molStringsColumn.length);
      return DG.Column.fromList('object', 'bitset', [resBitset]); // TODO: should return a bitset itself
    } catch (e: any) {
      console.error('Chem | In substructureSearch: ' + e.toString());
      throw e;
    }
  }

  @grok.decorators.fileExporter({
    description: 'As SDF...',
  })
  static async saveAsSdf(): Promise<void> {
    saveAsSdfDialog();
  }

  //#region Top menu

  @grok.decorators.func({
    name: 'Chem Similarity Search',
    outputs: [{name: 'result', type: 'viewer'}],
    meta: {icon: 'files/icons/chem-similarity-search-viewer.svg', role: 'viewer'},
  })
  static similaritySearchViewer(): ChemSimilarityViewer {
    return new ChemSimilarityViewer();
  }

  @grok.decorators.func({
    'name': 'Similarity Search',
    'top-menu': 'Chem | Search | Similarity Search...',
  })
  static similaritySearchTopMenu(): void {
    (grok.shell.v as DG.TableView).addViewer('Chem Similarity Search');
  }

  @grok.decorators.func({
    name: 'Chem Diversity Search',
    outputs: [{name: 'result', type: 'viewer'}],
    meta: {icon: 'files/icons/chem-diversity-search-viewer.svg', role: 'viewer'},
  })
  static diversitySearchViewer(): ChemDiversityViewer {
    return new ChemDiversityViewer();
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Search | Diversity Search...',
    'name': 'Diversity Search',
  })
  static diversitySearchTopMenu(): void {
    (grok.shell.v as DG.TableView).addViewer('Chem Diversity Search');
  }


  @grok.decorators.func({
    'top-menu': 'Chem | Calculate | Descriptors...',
  })
  static async descriptorsDocker(): Promise<void> {
    await openDescriptorsDialogDocker();
  }

  //function with tranfrom tag to be able to run within data sync projects, adds column to dataframe
  @grok.decorators.func({
    name: 'calculateDescriptorsTransform',
    meta: {role: 'transform'},
  })
  static async calculateDescriptorsTransform(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'list<string>'}) selected: string[],
  ): Promise<void> {
    const cols = await calculateDescriptors(molecules, selected);
    addDescriptorsColsToDf(table, cols);
  }

  //vector function to run in add new column dialog
  @grok.decorators.func({
    name: 'getDescriptors',
    meta: {vectorFunc: 'true'},
  })
  static async getDescriptors(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'list<string>', options: {optional: true}}) selected?: string[],
  ): Promise<DG.DataFrame> {
    if (!selected || selected.length === 0)
      selected = await getSelected();
    const cols = (await calculateDescriptors(molecules, selected)).filter((col) => col != null);
    return DG.DataFrame.fromColumns(cols);
  }

  @grok.decorators.func({outputs: [{name: 'descriptors', type: 'object'}]})
  static async chemDescriptorsTree(): Promise<object> {
    return await fetchWrapper(() => getDescriptorsTree());
  }

  @grok.decorators.func({
    'name': 'Map Identifiers',
    'top-menu': 'Chem | Calculate | Map Identifiers...',
  })
  static async getMapIdentifiers() {
    await openMapIdentifiersDialog();
  }

  @grok.decorators.func()
  static async freeTextToSmiles(molfile: string): Promise<string | null> {
    return await textToSmiles(molfile);
  }

  @grok.decorators.func()
  static async chemDescriptors(table: DG.DataFrame, molecules: DG.Column, descriptors: string[]): Promise<void> {
    const descCols = await fetchWrapper(() => calculateDescriptors(molecules, descriptors));
    addDescriptorsColsToDf(table, descCols);
  }

  @grok.decorators.editor({name: 'SearchSubstructureEditor'})
  static searchSubstructureEditor(call: DG.FuncCall): void {
    if (grok.shell.tv.dataFrame.rowCount > MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT) {
      grok.shell.warning(`Too many rows, maximum for substructure search is ${MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT}`);
      return;
    }
    const molColumns = grok.shell.tv.dataFrame.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE);
    if (!molColumns.length) {
      grok.shell.warning(`Data doesn't contain molecule columns`);
      return;
    } else if (molColumns.length === 1)
      call.func.prepare({molecules: molColumns[0]}).call(true);
    else {
      const colInput = ui.input.column('Molecules', {table: grok.shell.tv.dataFrame, value: molColumns[0]});
      ui.dialog({title: 'Substructure search'});
      ui.dialog({title: 'Substructure search'})
        .add(colInput)
        .onOK(async () => {
          call.func.prepare({molecules: colInput.value}).call(true);
        })
        .show();
    }
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Search | Substructure Search...',
    'name': 'Substructure Search',
    'editor': 'Chem:SearchSubstructureEditor',
  })
  static SubstructureSearchTopMenu(
    @grok.decorators.param({type: 'column', options: {semType: 'Molecule'}}) molecules: DG.Column): void {
    const fg = grok.shell.tv.getFiltersGroup({createDefaultFilters: false});
    grok.shell.tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: molecules.name,
      columnName: molecules.name,
      molBlock: DG.WHITE_MOLBLOCK,
    });
    grok.shell.tv.grid.scrollToCell(molecules, 0);
    const filterHeader = Array.from(fg.root!.getElementsByClassName('d4-filter-header'))
      .find((el) => Array.from(el!.getElementsByTagName('label')).find((it) => it.textContent === molecules.name));
    if (filterHeader) {
      setTimeout(() => {
        const sketchLink = (filterHeader.parentElement as HTMLElement).getElementsByClassName('sketch-link')[0];
        const element = sketchLink ?? (filterHeader.parentElement as HTMLElement)
          .getElementsByClassName('chem-canvas')[0];
        (element as HTMLElement).click();
      }, 500);
    }
  }

  @grok.decorators.func({
    'name': 'Cluster MCS',
    'top-menu': 'Chem | Calculate | Cluster MCS...',
    'description': 'Calculates most common substructures for each cluster',
  })
  static async clusterMCSTopMenu(
    table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Molecule'}}) molCol: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'categorical'}}) clusterCol: DG.Column): Promise<void> {
    const c = await PackageFunctions.performClusterMCS(molCol, clusterCol);
    c.name = table.columns.getUnusedName(c.name);
    table.columns.add(c);
  }

  @grok.decorators.func({
    name: 'clusterMCS',
    friendlyName: 'Cluster MCS',
    description: 'Calculates most common substructures for each cluster',
    meta: {vectorFunc: 'true'},
    outputs: [{name: 'result', type: 'column', options: {semType: 'Molecule'}}],
  })
  static async performClusterMCS(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molCol: DG.Column,
    @grok.decorators.param({options: {type: 'string'}}) clusterCol: DG.Column): Promise<DG.Column> {
    const PG = DG.TaskBarProgressIndicator.create('Most common substructures...');
    const mcsCol = DG.Column.string('Cluster MCS', molCol.length);
    try {
      const clusteredMols = new Array(clusterCol.categories.length).fill(null).map(() => [] as string[]);
      const indexes = clusterCol.getRawData();
      const mols = molCol.toList();
      //optimization to avoid pushing to arrays
      for (let i = 0; i < indexes.length; i++)
        clusteredMols[indexes[i]].push(mols[i]);

      const mcsResult = await (await chemCommonRdKit.getRdKitService()).clusterMCS(clusteredMols, true, true);
      mcsCol.semType = DG.SEMTYPE.MOLECULE;
      mcsCol.init((i) => mcsResult[indexes[i]]);
      mcsCol.setTag('.structure-filter-type', 'Categorical');
      mcsCol.setTag('.ignore-custom-filter', 'true');
    } catch (e) {
      grok.shell.error('Cluster MCS Error');
      console.error(e);
    }
    PG.close();
    return mcsCol;
  }

  @grok.decorators.editor()
  static ChemSpaceEditor(
    @grok.decorators.param({type: 'funccall'}) call: DG.FuncCall): void {
    const funcEditor = new DimReductionBaseEditor({semtype: DG.SEMTYPE.MOLECULE});
    const clusterMCS = ui.input.bool('Cluster MCS', {value: false, tooltipText: 'Perform MCS on clustered data'});
    const editor = funcEditor.getEditor();
    editor.appendChild(clusterMCS.root);
    const dialog = ui.dialog({title: 'Chemical space'})
      .add(editor)
      .onOK(async () => {
        const params = funcEditor.getParams();
        return call.func.prepare({
          molecules: params.col,
          table: params.table,
          methodName: params.methodName,
          similarityMetric: params.similarityMetric,
          plotEmbeddings: params.plotEmbeddings,
          options: params.options,
          preprocessingFunction: params.preprocessingFunction,
          clusterEmbeddings: params.clusterEmbeddings,
          clusterMCS: !!clusterMCS.value,
        }).call();
      });
    dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
    dialog.show();
  }

  @grok.decorators.func({
    name: 'Fingerprints',
    meta: {supportedSemTypes: 'Molecule', supportedDistanceFunctions: 'Tanimoto,Asymmetric,Cosine,Sokal', role: 'dimRedPreprocessingFunction'},
    outputs: [{type: 'object', name: 'result'}],
  })
  static async getFingerprints(
    @grok.decorators.param({options: {semType: 'Molecule'}}) col: DG.Column,
    @grok.decorators.param({options: {optional: true}}) _metric: string | undefined = undefined,
    @grok.decorators.param({type: 'string', options: {caption: 'Fingerprint type', optional: true,
      choices: ['Morgan', 'RDKit', 'Pattern', 'AtomPair', 'MACCS', 'TopologicalTorsion'], initialValue: 'Morgan'}}) fingerprintType: Fingerprint = Fingerprint.Morgan,
  ) {
    //TODO: get rid of fallback
    let fingerprintTypeStr = fingerprintType as string;
    if ((fingerprintTypeStr.startsWith('\'') || fingerprintTypeStr.startsWith('"')) &&
      fingerprintTypeStr.endsWith('\'') || fingerprintTypeStr.endsWith('"'))
      fingerprintTypeStr = fingerprintTypeStr.slice(1, -1);

    const fpColumn = await chemSearches.chemGetFingerprints(col, fingerprintTypeStr as Fingerprint, false);
    malformedDataWarning(fpColumn, col);
    return {entries: fpColumn, options: {}};
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Analyze | Chemical Space...',
    'name': 'Chem Space',
    'description': 'Maps the dataset to 2D plot based on similarity',
    'editor': 'Chem:ChemSpaceEditor',
    'outputs': [{type: 'viewer', name: 'result'}],
  })
  static async chemSpaceTopMenu(
    @grok.decorators.param({type: 'dataframe'}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'string', options: {choices: ['UMAP', 't-SNE']}}) methodName: DimReductionMethods,
    @grok.decorators.param({type: 'string', options: {choices: ['Tanimoto', 'Asymmetric', 'Cosine', 'Sokal']}}) similarityMetric: BitArrayMetrics = BitArrayMetricsNames.Tanimoto,
    @grok.decorators.param({type: 'bool', options: {initialValue: 'true'}}) plotEmbeddings: boolean,
    @grok.decorators.param({type: 'object', options: {optional: true}}) options?: (IUMAPOptions | ITSNEOptions) & Options,
    @grok.decorators.param({type: 'func', options: {optional: true}}) preprocessingFunction?: DG.Func,
    @grok.decorators.param({type: 'bool', options: {optional: true}}) clusterEmbeddings?: boolean,
    @grok.decorators.param({type: 'bool', options: {optional: true}}) clusterMCS?: boolean): Promise<DG.Viewer | undefined> {
    //workaround for functions which add viewers to tableView (can be run only on active table view)
    checkCurrentView(table);
    // after this check we are sure that tableView is active
    const tv = grok.shell.tv;
    if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
      grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
      return;
    }
    const clusterColName = table.columns.getUnusedName('Cluster (DBSCAN)');
    const embedColsNames: string[] = getEmbeddingColsNames(table);
    const funcCall = await DG.Func.find({name: 'chemSpaceTransform'})[0].prepare({
      table: table,
      molecules: molecules,
      methodName: methodName,
      similarityMetric: similarityMetric,
      plotEmbeddings: false,
      options: JSON.stringify(options),
      preprocessingFunction: preprocessingFunction,
      clusterEmbeddings: clusterEmbeddings,
    }).call(undefined, undefined, {processed: false});
    let res: DG.ScatterPlotViewer = funcCall.getOutputParamValue();

    if (plotEmbeddings) {
      res = tv.scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Chemical space'});
      const description = `Molecules column: ${molecules.name}, method: ${methodName}, ${options ? `options: ${JSON.stringify(options)},` : ``} similarity: ${similarityMetric}`;
      res.setOptions({description: description, descriptionVisibilityMode: 'Never'});
      //temporary fix (to save backward compatibility) since labels option type has been changed from string to array in 1.23 platform version
      if (Object.keys(res.props).includes('labelColumnNames')) { //@ts-ignore
        if (res.props['labelColumnNames'].constructor.name == 'Array')
          res.setOptions({labelColumnNames: [molecules.name]});
      }
      if (clusterEmbeddings) {
        res.props.colorColumnName = clusterColName;
        if (clusterMCS) {
          const clusterCol = table.col(clusterColName)!;
          const mcsCol = await PackageFunctions.performClusterMCS(molecules, clusterCol);
          mcsCol.name = table.columns.getUnusedName(mcsCol.name);
          table.columns.add(mcsCol);
        }
      }
    }
    return res;
  }

  @grok.decorators.func({
    outputs: [{type: 'viewer', name: 'result'}],
    meta: {role: 'transform'},
  })
  static async chemSpaceTransform(
    table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Molecule'}})molecules: DG.Column,
    @grok.decorators.param({type: 'string'}) methodName: DimReductionMethods,
    @grok.decorators.param({type: 'string'}) similarityMetric: BitArrayMetrics = BitArrayMetricsNames.Tanimoto,
    @grok.decorators.param({options: {initialValue: 'true'}}) plotEmbeddings: boolean,
    @grok.decorators.param({options: {optional: true}}) options?: string,
    @grok.decorators.param({options: {optional: true}}) clusterEmbeddings?: boolean,
  ): Promise<DG.Viewer | undefined> {
    const res = await runChemSpace(table, molecules, methodName, similarityMetric, plotEmbeddings, JSON.parse(options ?? '{}'),
      undefined, clusterEmbeddings);
    console.log(`returned from runChemSpace`);
    return res;
  }

  @grok.decorators.func({
    name: 'Chem Space Embeddings',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async getChemSpaceEmbeddings(
    @grok.decorators.param({type: 'string'}) col: DG.Column,
    @grok.decorators.param({type: 'string'}) methodName: DimReductionMethods,
    @grok.decorators.param({type: 'string'}) similarityMetric: BitArrayMetrics = BitArrayMetricsNames.Tanimoto,
    xAxis: string,
    yAxis: string,
    @grok.decorators.param({type: 'object', options: {optional: true}}) options?: any): Promise<ISequenceSpaceResult> {
    //need to create dataframe to add fingerprints column
    if (!col.dataFrame) {
      const dfForFp = DG.DataFrame.create(col.length);
      dfForFp.columns.add(col);
    }
    const chemSpaceParams = {
      seqCol: col,
      methodName: methodName,
      similarityMetric: similarityMetric as BitArrayMetrics,
      embedAxesNames: [xAxis, yAxis],
      options: options,
    };
    const chemSpaceRes = await chemSpace(chemSpaceParams);
    return chemSpaceRes;
  }

  @grok.decorators.func({
    name: 'Chem Similarities Matrix',
    outputs: [{name: 'result', type: 'object'}],
  })
  static async getChemSimilaritiesMatrix(
    @grok.decorators.param({type: 'int'}) dim: number,
      col: DG.Column,
      df: DG.DataFrame,
      colName: string,
    @grok.decorators.param({type: 'object'}) simArr: DG.Column[]): Promise<(DG.Column | null)[]> {
    //need to create dataframe to add fingerprints column
    if (!col.dataFrame) {
      const dfForFp = DG.DataFrame.create(col.length);
      dfForFp.columns.add(col);
    }
    return await getSimilaritiesMarix(dim, col, df, colName, simArr);
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Analyze | Elemental Analysis...',
    'name': 'Elemental Analysis',
  })
  static async elementalAnalysis(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({options: {description: 'Add a standalone radar viewer', initialValue: 'false'}}) radarViewer: boolean,
    @grok.decorators.param({options: {description: 'Show radar in grid cells', initialValue: 'false'}}) radarGrid: boolean): Promise<void> {
    if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
      grok.shell.info(`The column ${molecules.name} doesn't contain molecules`);
      return;
    }
    //workaround for functions which add viewers to tableView (can be run only on active table view)
    checkCurrentView(table);
    const view = grok.shell.tv;
    const funcCall = await DG.Func.find({name: 'runElementalAnalysis'})[0].prepare({
      table: table,
      molecules: molecules,
    }).call(undefined, undefined, {processed: false});
    const columnNames: string[] = funcCall.getOutputParamValue();

    if (radarViewer) {
      let packageExists = checkPackage('Charts', 'radarViewerDemo');
      if (!packageExists) //for compatibility with previous versions of Charts
        packageExists = checkPackage('Charts', '_radarViewerDemo');
      if (packageExists) {
        const radarViewer = DG.Viewer.fromType('Radar', table, {
          valuesColumnNames: columnNames,
        });
        view.addViewer(radarViewer);
      } else
        grok.shell.warning('Charts package is not installed');
    }
    if (radarGrid) {
      const packageExists = checkPackage('PowerGrid', 'radarCellRenderer');
      if (packageExists) {
        const gc = view.grid.columns.add({gridColumnName: `elements (${molecules.name})`, cellType: 'radar'});
        gc.settings = {columnNames: columnNames};
        gc.width = 300;
      } else
        grok.shell.warning('PowerGrid package is not installed');
    }
  }

  @grok.decorators.func({
    name: 'runElementalAnalysis',
    outputs: [{name: 'res', type: 'list'}],
    meta: {role: 'transform'},
  })
  static runElementalAnalysis(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column): string[] {
    const [elements, invalid]: [Map<string, Int32Array>, number[]] = getAtomsColumn(molecules);
    const columnNames: string[] = [];
    if (invalid.filter((el) => el !== null).length > 0) {
      console.log(`Invalid rows ${invalid.map((i) => i.toString()).join(', ')}`);
      grok.shell.warning('Dataset contains malformed data!');
    }
    const extendedElementsTable = ['R'].concat(elementsTable).concat(['Molecule Charge']);
    for (const elName of extendedElementsTable) {
      const value = elements.get(elName);
      if (value) {
        const column = DG.Column.fromInt32Array(elName, value);
        column.name = table.columns.getUnusedName(column.name);
        invalid.map((i) => {
          column.set(i, null);
        });
        table.columns.add(column);
        columnNames.push(column.name);
      }
    }
    return columnNames;
  }

  @grok.decorators.func({
    'name': 'R-Groups Analysis',
    'top-menu': 'Chem | Analyze | R-Groups Analysis...',
  })
  static rGroupsAnalysisMenu(): void {
    const col = grok.shell.t.columns.bySemType(DG.SEMTYPE.MOLECULE);
    if (col === null) {
      grok.shell.error('Current table does not contain molecules');
      return;
    }
    rGroupAnalysis(col);
  }


  @grok.decorators.func({
    outputs: [{name: 'result', type: 'object'}],
    meta: {role: 'transform'},
  })
  static async rGroupDecomposition(
    df: DG.DataFrame,
    molColName: string,
    core: string,
    rGroupName: string,
    rGroupMatchingStrategy: string,
    @grok.decorators.param({options: {optional: true, initialValue: 'false'}}) onlyMatchAtRGroups: boolean): Promise<RGroupDecompRes | undefined> {
    const params: RGroupParams = {
      molColName: molColName,
      core: core,
      rGroupName: rGroupName,
      rGroupMatchingStrategy: rGroupMatchingStrategy,
      onlyMatchAtRGroups: onlyMatchAtRGroups,
    };
    const col = df.col(molColName);
    if (col === null)
      throw new Error(`Current table does not contain ${params.molColName} column`); //exception
    return await rGroupDecomp(col, params);
  }

  @grok.decorators.editor()
  static ActivityCliffsEditor(
    call: DG.FuncCall): void {
    const funcEditor = new ActivityCliffsFunctionEditor({semtype: DG.SEMTYPE.MOLECULE});
    const dialog = ui.dialog({title: 'Activity Cliffs'})
      .add(funcEditor.getEditor())
      .onOK(async () => {
        const params = funcEditor.getParams();
        if (params.activities) {
          call.func.prepare({
            table: params.table,
            molecules: params.col,
            activities: params.activities,
            similarity: params.similarityThreshold,
            methodName: params.methodName,
            similarityMetric: params.similarityMetric,
            preprocessingFunction: params.preprocessingFunction,
            options: params.options,
          }).call(true);
        } else
          grok.shell.error(`Column with activities has not been selected. Table contains no numeric columns.`);
      });
    dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
    dialog.show();
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Analyze | Activity Cliffs...',
    'name': 'Activity Cliffs',
    'description': 'Detects pairs of molecules with similar structure and significant difference in any given property',
    'editor': 'Chem:ActivityCliffsEditor',
  })
  static async activityCliffs(
    @grok.decorators.param({options: {description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {type: 'categorical', semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'numerical'}}) activities: DG.Column,
    @grok.decorators.param({options: {description: 'Similarity cutoff', initialValue: '80'}}) similarity: number,
    @grok.decorators.param({type: 'string', options: {choices: ['UMAP', 't-SNE']}}) methodName: DimReductionMethods,
    @grok.decorators.param({type: 'string', options: {choices: ['Tanimoto', 'Asymmetric', 'Cosine', 'Sokal']}}) similarityMetric: BitArrayMetrics,
    @grok.decorators.param({type: 'func', options: {optional: true}}) preprocessingFunction: DG.Func,
    @grok.decorators.param({type: 'object', options: {optional: true}}) options?: (IUMAPOptions | ITSNEOptions) & Options,
    @grok.decorators.param({options: {optional: true}}) isDemo?: boolean,
    @grok.decorators.param({options: {optional: true}}) isTest?: boolean): Promise<void> {
    if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
      grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
      return;
    }
    if (activities.type !== DG.TYPE.INT && activities.type !== DG.TYPE.BIG_INT && activities.type !== DG.TYPE.FLOAT) {
      grok.shell.error(`Column ${activities.name} is not numeric`);
      return;
    }
    //workaround for functions which add viewers to tableView (can be run only on active table view)
    checkCurrentView(table);

    const allowedRowCount = 10000;
    const fastRowCount = methodName === DimReductionMethods.UMAP ? 5000 : 2000;
    if (table.rowCount > allowedRowCount) {
      grok.shell.warning(`Too many rows, maximum for activity cliffs is ${allowedRowCount}`);
      return;
    }

    const runActCliffs = async (): Promise<void> => {
      await DG.Func.find({name: 'activityCliffsTransform'})[0].prepare({
        table: table,
        molecules: molecules,
        activities: activities,
        similarity: similarity,
        methodName: methodName,
        similarityMetric: similarityMetric,
        options: JSON.stringify(options),
        isDemo: isDemo,
      }).call(undefined, undefined, {processed: false});

      const view = grok.shell.tv;

      const description = `Molecules: ${molecules.name}, activities: ${activities.name}, method: ${methodName}, ${options ? `options: ${JSON.stringify(options)},` : ``} similarity: ${similarityMetric}, similarity cutoff: ${similarity}`;
      view.addViewer(DG.VIEWER.SCATTER_PLOT, {
        xColumnName: axesNames[0],
        yColumnName: axesNames[1],
        color: activities.name,
        showXSelector: false,
        showYSelector: false,
        showSizeSelector: false,
        showColorSelector: false,
        markerMinSize: 5,
        markerMaxSize: 25,
        title: 'Activity cliffs',
        initializationFunction: 'activityCliffsInitFunction',
        description: description,
        descriptionVisibilityMode: 'Never',
      }) as DG.ScatterPlotViewer;
    };

    const axesNames = getEmbeddingColsNames(table);
    if (table.rowCount > fastRowCount && !isTest) {
      ui.dialog().add(ui.divText(`Activity cliffs analysis might take several minutes.
      Do you want to continue?`))
        .onOK(async () => {
          const progressBar = DG.TaskBarProgressIndicator.create(`Activity cliffs running...`);
          await runActCliffs();
          progressBar.close();
        })
        .show();
    } else
      await runActCliffs();
  }

  @grok.decorators.func({
    name: 'activityCliffsInitFunction',
  })
  static async activityCliffsInitFunction(
    @grok.decorators.param({type: 'viewer'}) sp: DG.ScatterPlotViewer): Promise<void> {
    const tag = sp.dataFrame.getTag('activityCliffsParams');
    if (!tag) {
      grok.shell.error(`Activity cliffs parameters not found in table tags`);
      return;
    }
    const actCliffsParams: ActivityCliffsParams = JSON.parse(tag);
    const molCol = sp.dataFrame.col(actCliffsParams.molColName)!;
    const actCol = sp.dataFrame.col(actCliffsParams.activityColName)!;
    const encodedColWithOptions = await PackageFunctions.getFingerprints(molCol);

    const axesNames = [sp.getOptions().look['xColumnName'], sp.getOptions().look['yColumnName']];

    await runActivityCliffs(sp, sp.dataFrame, molCol, encodedColWithOptions, actCol, axesNames,
      actCliffsParams.similarity, actCliffsParams.similarityMetric, actCliffsParams.options, DG.SEMTYPE.MOLECULE,
      {'units': molCol.meta.units!}, createTooltipElement, createPropPanelElement, undefined, undefined, actCliffsParams.isDemo);
    //temporary fix (to save backward compatibility) since labels option type has been changed from string to array in 1.23 platform version
    if (Object.keys(sp.props).includes('labelColumnNames')) { //@ts-ignore
      if (sp.props['labelColumnNames'].constructor.name == 'Array')
        sp.setOptions({labelColumnNames: [molCol.name]});
    }
    sp.render(sp.getInfo()['canvas'].getContext('2d'));
  }

  @grok.decorators.func({
    meta: {role: 'transform'},
  })
  static async activityCliffsTransform(
    @grok.decorators.param({options: {description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {type: 'categorical', semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'column', options: {type: 'numerical'}}) activities: DG.Column,
    @grok.decorators.param({options: {description: 'Similarity cutoff', initialValue: '80'}}) similarity: number,
    @grok.decorators.param({type: 'string', options: {choices: ['UMAP', 't-SNE']}}) methodName: DimReductionMethods,
    @grok.decorators.param({type: 'string', options: {choices: ['Tanimoto', 'Asymmetric', 'Cosine', 'Sokal']}}) similarityMetric: BitArrayMetrics,
    @grok.decorators.param({options: {optional: true}}) options?: string,
    @grok.decorators.param({options: {optional: true}}) isDemo?: boolean): Promise<void> {
    const preprocessingFunction = DG.Func.find({name: 'getFingerprints', package: 'Chem'})[0];
    const axesNames = getEmbeddingColsNames(table);
    await getActivityCliffsEmbeddings(table, molecules, axesNames, similarity,
      similarityMetric, methodName, JSON.parse(options ?? '{}'), preprocessingFunction);
    const tagContent: ActivityCliffsParams = {
      molColName: molecules.name,
      activityColName: activities.name,
      similarityMetric: similarityMetric,
      similarity: similarity,
      options: options ?? {},
      isDemo: isDemo,
    };
    table.setTag('activityCliffsParams', JSON.stringify(tagContent));
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Calculate | To InchI...',
    'name': 'To InchI',
    'meta': {'role': 'transform'},
  })
  static addInchisTopMenu(
    table: DG.DataFrame,
    @grok.decorators.param({name: 'molecules', options: {semType: 'Molecule'}}) col: DG.Column): void {
    const inchiCol = getInchisImpl(col);
    inchiCol.name = table.columns.getUnusedName(inchiCol.name);
    table.columns.add(inchiCol);
  }

  @grok.decorators.func({
    name: 'getInchis',
    meta: {vectorFunc: 'true'},
  })
  static getInchis(
    @grok.decorators.param({type: 'column<string>', options: {semType: 'Molecule'}}) molecules: DG.Column): DG.Column {
    return getInchisImpl(molecules);
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Calculate | To InchI Keys...',
    'name': 'To InchI Keys',
    'meta': {'role': 'transform'},
  })
  static addInchisKeysTopMenu(
    @grok.decorators.param({options: {description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({name: 'molecules', options: {semType: 'Molecule'}}) col: DG.Column): void {
    const inchiKeyCol = getInchiKeysImpl(col);
    inchiKeyCol.name = table.columns.getUnusedName(inchiKeyCol.name);
    table.columns.add(inchiKeyCol);
  }

  @grok.decorators.func({
    name: 'getInchiKeys',
    meta: {vectorFunc: 'true'},
  })
  static getInchiKeys(
    @grok.decorators.param({type: 'column<string>', options: {semType: 'Molecule'}}) molecules: DG.Column): DG.Column {
    return getInchiKeysImpl(molecules);
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Analyze | Structural Alerts...',
    'name': 'Structural Alerts',
    'description': 'Highlights the fragments that could lead to potential chemical hazards',
    'meta': {'role': 'hitTriageFunction'},
  })
  static async structuralAlertsTopMenu(
    @grok.decorators.param({options: {description: 'Input data table', caption: 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {caption: 'Molecules', semType: 'Molecule', type: 'categorical'}}) molecules: DG.Column,
    @grok.decorators.param({options: {caption: 'PAINS', initialValue: 'true', description: '"Pan Assay Interference Compounds filters"'}}) pains: boolean,
    @grok.decorators.param({options: {caption: 'BMS', initialValue: 'false', description: '"Bristol-Myers Squibb HTS Deck filters"'}}) bms: boolean,
    @grok.decorators.param({options: {caption: 'SureChEMBL', initialValue: 'false', description: '"MedChem unfriendly compounds from SureChEMBL"'}}) sureChembl: boolean,
    @grok.decorators.param({options: {caption: 'MLSMR', initialValue: 'false', description: '"NIH MLSMR Excluded Functionality filters"'}}) mlsmr: boolean,
    @grok.decorators.param({options: {caption: 'Dundee', initialValue: 'false', description: '"University of Dundee NTD Screening Library filters"'}}) dundee: boolean,
    @grok.decorators.param({options: {caption: 'Inpharmatica', initialValue: 'false', description: '"Inpharmatica filters"'}}) inpharmatica: boolean,
    @grok.decorators.param({options: {caption: 'LINT', initialValue: 'false', description: '"Pfizer LINT filters"'}}) lint: boolean,
    @grok.decorators.param({options: {caption: 'Glaxo', initialValue: 'false', description: '"Glaxo Wellcome Hard filters"'}}) glaxo: boolean): Promise<DG.DataFrame | void> {
    if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
      grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
      return;
    }

    await DG.Func.find({name: 'runStructuralAlerts'})[0].prepare({
      table: table,
      molecules: molecules,
      pains: pains,
      bms: bms,
      sureChembl: sureChembl,
      mlsmr: mlsmr,
      dundee: dundee,
      inpharmatica: inpharmatica,
      lint: lint,
      glaxo: glaxo,
    }).call(undefined, undefined, {processed: false});

    return table;
  }

  @grok.decorators.func({
    name: 'runStructuralAlerts',
    outputs: [{name: 'result', type: 'dataframe'}],
    meta: {role: 'transform'},
  })
  static async runStructuralAlerts(
    @grok.decorators.param({options: {caption: 'Table', description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {caption: 'Molecules', type: 'categorical', semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({options: {caption: 'PAINS', initialValue: 'true', description: '"Pan Assay Interference Compounds filters"'}}) pains: boolean,
    @grok.decorators.param({options: {caption: 'BMS', initialValue: 'false', description: '"Bristol-Myers Squibb HTS Deck filters"'}}) bms: boolean,
    @grok.decorators.param({options: {caption: 'SureChEMBL', initialValue: 'false', description: '"MedChem unfriendly compounds from SureChEMBL"'}}) sureChembl: boolean,
    @grok.decorators.param({options: {caption: 'MLSMR', initialValue: 'false', description: '"NIH MLSMR Excluded Functionality filters"'}}) mlsmr: boolean,
    @grok.decorators.param({options: {caption: 'Dundee', initialValue: 'false', description: '"University of Dundee NTD Screening Library filters"'}}) dundee: boolean,
    @grok.decorators.param({options: {caption: 'Inpharmatica', initialValue: 'false', description: '"Inpharmatica filters"'}}) inpharmatica: boolean,
    @grok.decorators.param({options: {caption: 'LINT', initialValue: 'false', description: '"Pfizer LINT filters"'}}) lint: boolean,
    @grok.decorators.param({options: {caption: 'Glaxo', initialValue: 'false', description: '"Glaxo Wellcome Hard filters"'}}) glaxo: boolean): Promise<void> {
    if (table.rowCount > 5000)
      grok.shell.info('Structural Alerts detection will take a while to run');

    const ruleSet: RuleSet = {'PAINS': pains, 'BMS': bms, 'SureChEMBL': sureChembl, 'MLSMR': mlsmr,
      'Dundee': dundee, 'Inpharmatica': inpharmatica, 'LINT': lint, 'Glaxo': glaxo};
    const resultDf = await getStructuralAlertsByRules(molecules, ruleSet);

    if (resultDf) {
      for (const resultCol of resultDf.columns) {
        resultCol.name = table.columns.getUnusedName(`${resultCol.name} (${molecules.name})`);
        table.columns.add(resultCol.clone());
      }
    }
  }

  @grok.decorators.func({
    name: 'getStructuralAlerts',
    meta: {vectorFunc: 'true'},
  })
  static async getStructuralAlerts(
    @grok.decorators.param({type: 'column<string>', options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'list<string>', options: {optional: true}}) alerts?: string[]): Promise<DG.DataFrame> {
    const lowerCaseAlerts = alerts?.map((it) => it.toLowerCase());
    const ruleSet: {[key: string]: boolean} = {};
    for (const rule of STRUCT_ALERTS_RULES_NAMES)
      ruleSet[rule] = !lowerCaseAlerts || lowerCaseAlerts.length === 0 ? true : lowerCaseAlerts.includes(rule.toLowerCase());

    const resultDf = await getStructuralAlertsByRules(molecules, ruleSet as RuleSet);
    if (!resultDf)
      throw Error(`Structural alerts haven't been calculated`);
    return resultDf;
  }

  //#endregion

  //#region Molecule column property panel


  @grok.decorators.panel({
    name: 'Chemistry | Rendering',
    meta: {'exclude-actions-panel': 'true'},
  })
  static molColumnPropertyPanel(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molColumn: DG.Column): DG.Widget {
    return getMolColumnPropertyPanel(molColumn);
  }

  @grok.decorators.panel({
    name: 'Chemistry | Highlight',
    meta: {'exclude-actions-panel': 'true'},
  })
  static molColumnHighlights(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molColumn: DG.Column): DG.Widget {
    return getmolColumnHighlights(molColumn);
  }

  @grok.decorators.panel({
    name: 'Chemistry | Descriptors',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static descriptorsWidget(
    @grok.decorators.param({options: {semType: 'Molecule'}}) smiles: string): DG.Widget {
    if (!smiles || DG.chem.Sketcher.isEmptyMolfile(smiles))
      return new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
    return _isSmarts(smiles) || isFragment(smiles) ?
      new DG.Widget(ui.divText(SMARTS_MOLECULE_MESSAGE)) :
      getDescriptorsSingle(smiles);
  }

  @grok.decorators.panel({
    'name': 'Biology | Drug Likeness',
    'description': 'Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.',
    'help-url': '/help/domains/chem/info-panels/drug-likeness.md',
    'meta': {'role': 'widgets', 'domain': 'chem'},
  })
  static drugLikeness(
    @grok.decorators.param({options: {semType: 'Molecule'}}) smiles: DG.SemanticValue): DG.Widget {
    return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles.value) ?
      _isSmarts(smiles.value) || isFragment(smiles.value) ? new DG.Widget(ui.divText(SMARTS_MOLECULE_MESSAGE)) :
        drugLikenessWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
  }


  @grok.decorators.panel({
    name: 'Chemistry | Properties',
    description: 'Basic molecule properties',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static properties(
    @grok.decorators.param({options: {semType: 'Molecule'}}) smiles: DG.SemanticValue): DG.Widget {
    return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles.value) ?
      _isSmarts(smiles.value) || isFragment(smiles.value)? new DG.Widget(ui.divText(SMARTS_MOLECULE_MESSAGE)) :
        propertiesWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
  }

  @grok.decorators.func({
    name: 'getChemPropertyFunction',
    description: 'Return chem property function',
    outputs: [{name: 'result', type: 'object'}],
  })
  static getChemPropertyFunction(
    @grok.decorators.param({type: 'string'}) name: string): null | ((smiles: string) => any) {
    return getChemPropertyFunc(name);
  }

  @grok.decorators.panel({
    'name': 'Biology | Structural Alerts',
    'description': 'Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response',
    'help-url': '/help/domains/chem/info-panels/structural-alerts.md',
    'meta': {'role': 'widgets', 'domain': 'chem'},
  })
  static async structuralAlerts(
    @grok.decorators.param({options: {semType: 'Molecule'}}) smiles: string): Promise<DG.Widget> {
    return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles) ?
      _isSmarts(smiles) || isFragment(smiles) ? new DG.Widget(ui.divText(SMARTS_MOLECULE_MESSAGE)) :
        structuralAlertsWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
  }

  @grok.decorators.panel({
    name: 'Structure | Identifiers',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async identifiers(
    @grok.decorators.param({options: {semType: 'Molecule'}}) smiles: string): Promise<DG.Widget> {
    return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles) ?
      _isSmarts(smiles) || isFragment(smiles) ? new DG.Widget(ui.divText(SMARTS_MOLECULE_MESSAGE)) :
        await identifiersWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
  }

  @grok.decorators.panel({
    name: 'Structure | 3D Structure',
    description: '3D molecule representation',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async structure3D(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecule: string): Promise<DG.Widget> {
    return molecule && !DG.chem.Sketcher.isEmptyMolfile(molecule) ?
      _isSmarts(molecule) || isFragment(molecule) ? new DG.Widget(ui.divText(SMARTS_MOLECULE_MESSAGE)) :
        structure3dWidget(molecule) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
  }

  @grok.decorators.panel({
    name: 'Structure | 2D Structure',
    description: '2D molecule representation',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static structure2d(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecule: string): DG.Widget {
    return molecule && !DG.chem.Sketcher.isEmptyMolfile(molecule) ?
      structure2dWidget(molecule) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
  }

  @grok.decorators.panel({
    'name': 'Biology | Toxicity',
    'description': 'Toxicity prediction. Calculated by openchemlib',
    'help-url': '/help/domains/chem/info-panels/toxicity-risks.md',
    'meta': {'role': 'widgets', 'domain': 'chem'},
  })
  static toxicity(
    @grok.decorators.param({options: {semType: 'Molecule'}}) smiles: DG.SemanticValue): DG.Widget {
    return smiles && !DG.chem.Sketcher.isEmptyMolfile(smiles.value) ?
      _isSmarts(smiles.value) || isFragment(smiles.value) ? new DG.Widget(ui.divText(SMARTS_MOLECULE_MESSAGE)) :
        toxicityWidget(smiles) : new DG.Widget(ui.divText(EMPTY_MOLECULE_MESSAGE));
  }

  @grok.decorators.func({
    name: 'convertMoleculeNotation',
    meta: {vectorFunc: 'true'},
  })
  static async convertMoleculeNotation(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecule: DG.Column,
    @grok.decorators.param({type: 'string'}) targetNotation: DG.chem.Notation,
    @grok.decorators.param({options: {initialValue: 'false', optional: true, nullable: true}}) kekulize?: boolean,
  ): Promise<DG.Column> {
    let col: DG.Column;
    let newColName = `${molecule.name}_${targetNotation}`;
    try {
      if (!!molecule.dataFrame?.columns)
        newColName = molecule.dataFrame.columns.getUnusedName(newColName);
    } catch (e) {}

    try {
      const res = await convertNotationForColumn(molecule, targetNotation, kekulize ?? false);
      col = DG.Column.fromStrings(newColName, res);
      col.semType = DG.SEMTYPE.MOLECULE;
    } catch (e: any) {
      col = DG.Column.string(newColName, molecule.length).init((_) => e?.message);
    }
    return col;
  }

  @grok.decorators.func({
    name: 'convertMolNotation',
    description: 'RDKit-based conversion for SMILES, SMARTS, InChi, Molfile V2000 and Molfile V3000',
    outputs: [{name: 'result', type: 'string', options: {semType: 'Molecule'}}],
    meta: {role: 'unitConverter'},
  })
  static convertMolNotation(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecule: string,
    @grok.decorators.param({type: 'string', options: {choices: ['smiles', 'cxsmiles', 'smarts', 'cxsmarts', 'molblock', 'v3Kmolblock']}}) sourceNotation: DG.chem.Notation,
    @grok.decorators.param({type: 'string', options: {choices: ['smiles', 'cxsmiles', 'smarts', 'cxsmarts', 'molblock', 'v3Kmolblock']}}) targetNotation: DG.chem.Notation ): string {
    return _convertMolNotation(molecule, sourceNotation, targetNotation, PackageFunctions.getRdKitModule());
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Transform | Convert Notation...',
    'name': 'Convert Notation',
    'meta': {'role': 'transform'},
  })
  static async convertNotation(
    data: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Molecule'}}) molecules: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {choices: ['smiles', 'cxsmiles', 'smarts', 'cxsmarts', 'molblock', 'v3Kmolblock'], initialValue: 'smiles'}}) targetNotation: DG.chem.Notation,
    @grok.decorators.param({options: {initialValue: 'false'}}) overwrite: boolean = false,
    @grok.decorators.param({options: {initialValue: 'true'}}) join: boolean = true,
    @grok.decorators.param({options: {initialValue: 'false', optional: true, nullable: true}}) kekulize?: boolean): Promise<DG.Column<string> | void> {
    const res = await convertNotationForColumn(molecules, targetNotation, kekulize ?? false);
    const units = targetNotation === DG.chem.Notation.MolBlock ? DG.UNITS.Molecule.MOLBLOCK :
      targetNotation === DG.chem.Notation.V3KMolBlock ? DG.UNITS.Molecule.V3K_MOLBLOCK : DG.UNITS.Molecule.SMILES;
    if (overwrite) {
      for (let i = 0; i < molecules.length; i++)
        molecules.set(i, res[i], false);
      molecules.meta.units = units;
    } else {
      const colName = data.columns.getUnusedName(`${molecules.name}_${targetNotation}`);
      const col = DG.Column.fromStrings(colName, res);
      col.meta.units = units;
      col.semType = DG.SEMTYPE.MOLECULE;
      if (!join)
        return col;
      data.columns.add(col);
    }
  }

  @grok.decorators.func({
    name: 'Convert Notation...',
    meta: {action: 'Convert Notation...'},
  })
  static convertMolNotationAction(
    @grok.decorators.param({options: {semType: 'Molecule'}}) col: DG.Column) {
    const func = DG.Func.find({name: 'convertNotation', package: 'Chem'})[0];
    if (!func || !col?.dataFrame)
      return;
    func.prepare({data: col.dataFrame, molecules: col}).edit();
  }

  @grok.decorators.func({
    name: 'Convert Mixture To Smiles...',
    meta: {action: 'Convert mixture to smiles...'},
  })
  static convertMixtureToSmiles(
    @grok.decorators.param({options: {semType: 'ChemicalMixture'}}) col: DG.Column): void {
    // Each cell is a Mixfile JSON string
    const smilesArr: string[] = [];

    function collectSmiles(component: any, smilesList: string[]) {
      if (component.smiles)
        smilesList.push(component.smiles);
      else if (component.molfile) {
        try {
          const smiles = PackageFunctions.convertMolNotation(component.molfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
          smilesList.push(smiles);
        } catch {}
      }
      if (Array.isArray(component.contents)) {
        for (const child of component.contents)
          collectSmiles(child, smilesList);
      }
    }

    for (let i = 0; i < col.length; ++i) {
      let mixfile: Mixfile;
      try {
        mixfile = JSON.parse(col.get(i));
      } catch {
        smilesArr.push('');
        continue;
      }
      const smilesList: string[] = [];
      collectSmiles(mixfile, smilesList);
      smilesArr.push(smilesList.join('.'));
    }

    if (col.dataFrame) {
      const name = col.dataFrame.columns.getUnusedName(`${col.name}_smiles`);
      const smilesCol = DG.Column.fromStrings(name, smilesArr);
      smilesCol.meta.units = DG.UNITS.Molecule.SMILES;
      smilesCol.semType = DG.SEMTYPE.MOLECULE;
      col.dataFrame.columns.add(smilesCol);
    }
  }

  @grok.decorators.func({
    description: 'Molecule',
    meta: {role: 'cellEditor'},
  })
  static async editMoleculeCell(
    @grok.decorators.param({type: 'grid_cell'}) cell: DG.GridCell): Promise<void> {
    const sketcher = new Sketcher(undefined, PackageFunctions.validateMolecule);
    const unit = cell.cell.column.meta.units;
    let molecule = cell.cell.value;
    if (unit === DG.chem.Notation.Smiles) {
      //convert to molFile to draw in coordinates similar to dataframe cell
      molecule = PackageFunctions.convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    }
    sketcher.setMolecule(molecule);
    const dlg = ui.dialog()
      .add(sketcher)
      .onOK(() => {
        if (unit === DG.chem.Notation.Smiles) {
          //set new cell value only in case smiles has been edited (to avoid undesired molecule orientation change)
          const newValue = sketcher.getSmiles();
          const mol = checkMoleculeValid(cell.cell.value);
          if (!checkMolEqualSmiles(mol, newValue)) {
            try {
              //@ts-ignore TODO Remove on js-api update
              cell.setValue(newValue, true);
            } catch {
              cell.cell.value = newValue;
            }
          }
          mol?.delete();
        } else {
          try {
            //@ts-ignore TODO Remove on js-api update
            cell.setValue(sketcher.getMolFile(), true);
          } catch {
            cell.cell.value = sketcher.getMolFile();
          }
        }
        Sketcher.addToCollection(Sketcher.RECENT_KEY, sketcher.getMolFile());
      })
      .show({resizable: true});
    ui.onSizeChanged(dlg.root).subscribe((_) => {
      if (!sketcher.sketcher?.isInitialized)
        return;
      sketcher._autoResized ? sketcher._autoResized = false : sketcher.resize();
    });
  }

  @grok.decorators.func({
    name: 'OpenChemLib',
    outputs: [{name: 'sketcher', type: 'widget'}],
    meta: {role: 'moleculeSketcher'},
  })
  static openChemLibSketcher(): OpenChemLibSketcher {
    return new OpenChemLibSketcher();
  }

  @grok.decorators.fileHandler({
    description: 'Opens SDF file',
    ext: 'sdf,mol',
  })
  static importSdf(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array): DG.DataFrame[] | void {
    try {
      return _importSdf(Uint8Array.from(bytes));
    } catch (e:any) {
      grok.shell.warning('file is not supported or malformed');
      grok.shell.error(e);
    }
  }

  @grok.decorators.fileHandler({
    description: 'Opens smi file',
    ext: 'smi',
  })
  static importSmi(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array): DG.DataFrame[] | void {
    try {
      return _importSmi(Uint8Array.from(bytes));
    } catch (e:any) {
      grok.shell.warning('file is not supported or malformed');
      grok.shell.error(e);
    }
  }

  @grok.decorators.fileHandler({
    description: 'Opens smi file',
    ext: 'mol2',
  })
  static importMol2(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array): DG.DataFrame[] | void {
    try {
      return _importTripos(Uint8Array.from(bytes));
    } catch (e:any) {
      grok.shell.warning('file is not supported or malformed');
      grok.shell.error(e);
    }
  }

  @grok.decorators.fileHandler({
    description: 'Opens MOL file',
    ext: 'mol',
  })
  static importMol(content: string): DG.DataFrame[] | void {
    try {
      const molCol = DG.Column.string('molecule', 1).init((_) => content);
      return [DG.DataFrame.fromColumns([molCol])];
    } catch (e:any) {
      grok.shell.warning('file is not supported or malformed');
      grok.shell.error(e);
    }
  }

  @grok.decorators.func({
    outputs: [{name: 'result', type: 'grid_cell_renderer'}],
    meta: {chemRendererName: 'OpenChemLib'},
  })
  static async oclCellRenderer(): Promise<OCLCellRenderer> {
    return new OCLCellRenderer();
  }

  @grok.decorators.func({
    name: 'Sort by similarity',
    description: 'Sorts a molecular column by similarity',
    meta: {action: 'Sort by similarity'},
  })
  static async sortBySimilarity(
    @grok.decorators.param({options: {semType: 'Molecule'}}) value: DG.SemanticValue): Promise<void> {
    if (!value || !value.cell || !value.cell.dart ||!value.cell.column) {
      grok.shell.error('Sorting by similarity requires a valid table cell');
      return;
    }
    const molCol = value.cell.column;
    const tableRowIdx = value.cell.rowIndex;
    const dframe = molCol.dataFrame;
    const smiles = molCol.get(tableRowIdx);
    checkCurrentView(dframe);
    const grid = value.viewer as DG.Grid ?? grok.shell.tv?.grid;
    if (!grid)
      throw new Error('Cannnot access value grid');
    ui.setUpdateIndicator(grid.root, true);
    const progressBar = DG.TaskBarProgressIndicator.create('Sorting Structures...');
    progressBar.update(0, 'Installing ScaffoldGraph..: 0% completed');
    const fingerprints : DG.DataFrame = await PackageFunctions.callChemSimilaritySearch(dframe, molCol, smiles,
      BitArrayMetricsNames.Tanimoto, Fingerprint.Morgan, 1000000, 0.0);
    ui.setUpdateIndicator(grid.root, false);
    progressBar.update(100, 'Sort completed');
    progressBar.close();

    const idxCol = fingerprints.columns.byName('indexes');
    grid.sort([], []);
    grid.setRowOrder(idxCol.toList());
    //grid.props.pinnedRows = [tableRowIdx];
    //next two rows can be added after Chem is updated to version 1.7 API
    grid.props.pinnedRowColumnNames = [molCol.name];
    grid.props.pinnedRowValues = [value.value];
    grid.scrollToPixels(0, 0); //to address the bug in the core
  }

  @grok.decorators.func({
    name: 'Use as filter',
    description: 'Adds this structure as a substructure filter',
    meta: {action: 'Use as filter'},
  })
  static useAsSubstructureFilter(
    @grok.decorators.param({options: {semType: 'Molecule'}}) value: DG.SemanticValue): void {
    const tv = grok.shell.tv;
    if (tv == null)
      throw new Error('Requires an open table view.');

    const molCol = value.cell?.column;
    const molecule = value.value;
    if (molCol == null)
      throw new Error('Molecule column not found.');

    let molblock;

    //in case molecule is smiles setting correct coordinates to save molecule orientation in filter
    if (value.cell.column.meta.units == DG.chem.Notation.Smiles)
      molblock = PackageFunctions.convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    else
      molblock = molToMolblock(molecule, PackageFunctions.getRdKitModule());

    tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: molCol.name,
      columnName: molCol.name,
      molBlock: molblock,
    }, false);
  }

  @grok.decorators.func({
    name: 'Copy as...',
    description: 'Copies structure in different formats',
    meta: {'action': 'Copy as...', 'exclude-current-value-menu': 'true'},
  })
  static copyAsAction(
    @grok.decorators.param({options: {semType: 'Molecule'}}) value: DG.SemanticValue) {
    const formats = ['Smiles', 'MolfileV2000', 'MolfileV3000', 'Smarts'];
    const menu = DG.Menu.popup();

    formats.forEach((format) => {
      const func = DG.Func.find({package: 'Chem', name: `copyAs${format}`})[0];
      if (func) {
        menu.item(format, () => {
          func.apply({value: value});
        });
      }
    });
    menu.show();
  }

  @grok.decorators.func({
    name: 'Copy as SMILES',
    description: 'Copies structure as smiles',
    meta: {'action': 'Copy as SMILES', 'exclude-actions-panel': 'true'},
  })
  static copyAsSmiles(
    @grok.decorators.param({options: {semType: 'Molecule'}}) value: DG.SemanticValue): void {
    const smiles = !DG.chem.isMolBlock(value.value) && !_isSmarts(value.value) ? value.value :
      _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, PackageFunctions.getRdKitModule());
    navigator.clipboard.writeText(smiles);
    grok.shell.info('Smiles copied to clipboard');
  }

  @grok.decorators.func({
    name: 'Copy as MOLFILE V2000',
    description: 'Copies structure as molfile V2000',
    meta: {'action': 'Copy as MOLFILE V2000', 'exclude-actions-panel': 'true'},
  })
  static copyAsMolfileV2000(
    @grok.decorators.param({options: {semType: 'Molecule'}}) value: DG.SemanticValue): void {
    const molfileV2000 = DG.chem.isMolBlock(value.value) && !value.value.includes('V3000') ? value.value :
      _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, PackageFunctions.getRdKitModule());
    navigator.clipboard.writeText(molfileV2000);
    grok.shell.info('Molfile V2000 copied to clipboard');
  }

  @grok.decorators.func({
    name: 'Copy as MOLFILE V3000',
    description: 'Copies structure as molfile V3000',
    meta: {'action': 'Copy as MOLFILE V3000', 'exclude-actions-panel': 'true'},
  })
  static copyAsMolfileV3000(
    @grok.decorators.param({options: {semType: 'Molecule'}}) value: DG.SemanticValue): void {
    const molfileV3000 = DG.chem.isMolBlock(value.value) && value.value.includes('V3000') ? value.value :
      _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.V3KMolBlock, PackageFunctions.getRdKitModule());
    navigator.clipboard.writeText(molfileV3000);
    grok.shell.info('Molfile V3000 copied to clipboard');
  }

  @grok.decorators.func({
    name: 'Copy as SMARTS',
    description: 'Copies structure as smarts',
    meta: {'action': 'Copy as SMARTS', 'exclude-actions-panel': 'true'},
  })
  static copyAsSmarts(
    @grok.decorators.param({options: {semType: 'Molecule'}})value: DG.SemanticValue): void {
    const smarts = !DG.chem.isMolBlock(value.value) && _isSmarts(value.value) ? value.value :
      _convertMolNotation(value.value, DG.chem.Notation.Unknown, DG.chem.Notation.Smarts, PackageFunctions.getRdKitModule());
    navigator.clipboard.writeText(smarts);
    grok.shell.info('Smarts copied to clipboard');
  }

  @grok.decorators.func({
    name: 'Copy as IMAGE',
    description: 'Copies structure as Image',
    meta: {'action': 'Copy as Image', 'exclude-actions-panel': 'true'},
  })
  static copyAsImage(
    @grok.decorators.param({options: {semType: 'Molecule'}}) value: DG.SemanticValue): void {
    if (!value?.value || !value?.gridCell?.bounds || !value?.gridCell?.renderer)
      return;
    const gridCellBounds = value.gridCell.bounds;
    const w = 600;
    const heightMultiplier = w / gridCellBounds.width;
    const h = gridCellBounds.height * heightMultiplier;

    const renderer = value.gridCell.renderer;

    const canvas = ui.canvas(w * window.devicePixelRatio, h * window.devicePixelRatio);
    renderer.render(canvas.getContext('2d')!, 0, 0, w, h, value.gridCell, value.gridCell.style);
    canvas.toBlob((blob) => {
      if (!blob)
        return;
      navigator.clipboard.write([new ClipboardItem({'image/png': blob})]).then(() => {
        grok.shell.info('Image copied to clipboard');
      });
    });
  }

  @grok.decorators.func()
  static isSmiles(
    @grok.decorators.param({type: 'string'}) s: string) : boolean {
    const ctx: IMolContext = getMolSafe(s, {}, _rdKitModule, true);
    if (ctx.mol !== null) {
      ctx.mol.delete();
      return true;
    }
    return false;
  }

  @grok.decorators.func()
  static isSmarts(s: string): boolean {
    return !!s.match(/\[.?#\d|\$|&|;|,|!.?]/g);
  }

  @grok.decorators.func()
  static detectSmiles(col: DG.Column,
    @grok.decorators.param({type: 'int'}) min: number) : void {
    if (DG.Detector.sampleCategories(col, PackageFunctions.isSmiles, min, 10, 0.8)) {
      col.meta.units = DG.UNITS.Molecule.SMILES;
      col.semType = DG.SEMTYPE.MOLECULE;
    }
  }

  @grok.decorators.func({
    name: 'chemSimilaritySearch',
    outputs: [{name: 'result', type: 'dataframe'}],
  })
  static async callChemSimilaritySearch(
    df: DG.DataFrame,
    col: DG.Column,
    molecule: string,
    @grok.decorators.param({type: 'string'}) metricName: BitArrayMetrics,
    fingerprint: string,
    @grok.decorators.param({type: 'int'}) limit: number,
    minScore: number): Promise<DG.DataFrame> {
    const res = await chemSimilaritySearch(df, col, molecule, metricName, limit, minScore,
      fingerprint as Fingerprint, DG.BitSet.create(col.length).setAll(true));
    return res ?? DG.DataFrame.create();
  }


  @grok.decorators.func({
    name: 'chemDiversitySearch',
    outputs: [{name: 'result', type: 'dataframe'}],
  })
  static async callChemDiversitySearch(
    col: DG.Column,
    @grok.decorators.param({type: 'string'}) metricName: BitArrayMetrics,
    fingerprint: string, @grok.decorators.param({type: 'int'}) limit: number): Promise<number[]> {
    return await chemDiversitySearch(col, similarityMetric[metricName], limit,
      fingerprint as Fingerprint, DG.BitSet.create(col.length).setAll(true));
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Calculate | Chemical Properties...',
    'name': 'Chemical Properties',
    'description': 'Calculates chemical properties and adds them as columns to the input table. properties include Molecular Weight (MW), Hydrogen Bond Acceptors (HBA), Hydrogen Bond Donors (HBD), LogP (Partition), LogS (Solubility), Polar Surface Area (PSA), Rotatable Bonds, Stereo Centers, Molecule Charge.',
    'meta': {'function_family': 'biochem-calculator', 'method_info.author': 'Open Chem Lib Team', 'method_info.year': '2024', 'method_info.github': 'https://github.com/actelion/openchemlib', 'role': 'hitTriageFunction,transform'}})
  static async addChemPropertiesColumns(
    @grok.decorators.param({type: 'dataframe', options: {description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({options: {initialValue: 'true'}}) MW?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) HBA?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) HBD?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) logP?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) logS?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) PSA?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) rotatableBonds?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) stereoCenters?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) moleculeCharge?: boolean,
  ): Promise<void> {
    const propArgs: string[] = ([] as string[]).concat(MW ? ['MW'] : [], HBA ? ['HBA'] : [],
      HBD ? ['HBD'] : [], logP ? ['LogP'] : [], logS ? ['LogS'] : [], PSA ? ['PSA'] : [],
      rotatableBonds ? ['Rotatable bonds'] : [], stereoCenters ? ['Stereo centers'] : [],
      moleculeCharge ? ['Molecule charge'] : []);
    const pb = DG.TaskBarProgressIndicator.create('Chemical properties ...');
    try {
      await addPropertiesAsColumns(table, molecules, propArgs);
    } finally {
      pb.close();
    }
  }

  @grok.decorators.func({
    name: 'getProperties',
    meta: {vectorFunc: 'true'},
  })
  static async getProperties(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'list<string>', options: {optional: true}}) selected?: string[]): Promise<DG.DataFrame> {
    const propNames = Object.keys(CHEM_PROP_MAP);
    let props: string[] = [];
    if (!selected || selected.length === 0)
      props = propNames;
    else {
      for (const propName of propNames)
        props = props.concat(selected.includes(propName) ? [propName] : []);
    }

    const cols = await getPropertiesAsColumns(molecules, props);

    return DG.DataFrame.fromColumns(cols);
  }


  @grok.decorators.func({
    'top-menu': 'Chem | Calculate | Toxicity Risks...',
    'name': 'Toxicity Risks',
    'meta': {'role': 'hitTriageFunction,transform'},
  })
  static async addChemRisksColumns(
    @grok.decorators.param({options: {description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({options: {initialValue: 'true'}}) mutagenicity?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) tumorigenicity?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) irritatingEffects?: boolean,
    @grok.decorators.param({options: {initialValue: 'false'}}) reproductiveEffects?: boolean,
  ): Promise<void> {
    const pb = DG.TaskBarProgressIndicator.create('Toxicity risks ...');
    try {
      const toxCols = await getToxicityRisksColumns(molecules, {mutagenicity, tumorigenicity, irritatingEffects, reproductiveEffects});
      toxCols.forEach((col) => table.columns.add(col));
    } finally {
      pb.close();
    }
  }


  @grok.decorators.func({
    name: 'getToxicityRisks',
    meta: {vectorFunc: 'true'},
  })
  static async getToxicityRisks(
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'list<string>', options: {optional: true}}) risks?: string[]): Promise<DG.DataFrame> {
    const toxCols = await getToxicityRisksColumns(molecules, {
      mutagenicity: !risks || !risks.length ? true : risks.includes('mutagenicity'),
      tumorigenicity: !risks || !risks.length ? true : risks.includes('tumorigenicity'),
      irritatingEffects: !risks || !risks.length ? true : risks.includes('irritatingEffects'),
      reproductiveEffects: !risks || !risks.length ? true : risks.includes('reproductiveEffects'),
    });
    return DG.DataFrame.fromColumns(toxCols);
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Analyze | Scaffold Tree',
    'name': 'addScaffoldTree',
    'description': 'Generates a hierarchical tree based on the scaffolds presented in dataset'})
  static addScaffoldTree(): void {
    DG.ObjectPropertyBag.setDefaultProperty('Scaffold Tree', 'allowGenerate', true);
    grok.shell.tv.addViewer(ScaffoldTreeViewer.TYPE);
  }


  @grok.decorators.func({
    name: 'Matched Molecular Pairs Analysis',
    meta: {showInGallery: 'false', role: 'viewer'},
    outputs: [{name: 'result', type: 'viewer'}],
  })
  static mmpViewer(): MatchedMolecularPairsViewer {
    return new MatchedMolecularPairsViewer();
  }

  @grok.decorators.editor({
    name: 'MMPEditor',
  })
  static MMPEditor(call: DG.FuncCall): void {
    const funcEditor = new MmmpFunctionEditor();
    const editor = funcEditor.getEditor();
    const dialog = ui.dialog({title: 'Matched Molecular Pairs'})
      .add(editor)
      .onOK(async () => {
        const params = funcEditor.getParams();
        return call.func.prepare(params).call();
      });
    // dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
    dialog.show();
  }

  @grok.decorators.func({
    'name': 'Matched Molecular Pairs',
    'editor': 'Chem:MMPEditor',
    'top-menu': 'Chem | Analyze | Matched Molecular Pairs...'})
  static async mmpAnalysis(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'column_list', options: {type: 'numerical'}}) activities: DG.Column[],
    @grok.decorators.param({type: 'string_list'}) diffTypes: MmpDiffTypes[],
    @grok.decorators.param({type: 'string_list'}) scalings: SCALING_METHODS[],
    @grok.decorators.param({type: 'double',
      options: {description: 'Maximum fragment size relative to core', initialValue: '0.4'}}) fragmentCutoff: number = 0.4): Promise<void> {
    if (activities.length < 1) {
      grok.shell.warning('MMP analysis requires at least one activity');
      return;
    }

    //workaround for functions which add viewers to tableView (can be run only on active table view)
    checkCurrentView(table);

    const view = grok.shell.tv as DG.TableView;

    const activityColsNames = [];
    for (let i = 0; i < scalings.length; i++) {
      if (scalings[i] === SCALING_METHODS.NONE)
        activityColsNames.push(activities[i].name);
      else {
        const scaledCol = scaleActivity(activities[i], scalings[i]);
        const name = grok.shell.tv.dataFrame.columns.getUnusedName(scaledCol.name);
        scaledCol.name = name;
        grok.shell.tv.dataFrame.columns.add(scaledCol);
        activityColsNames.push(name);
      }
    }

    const viewer = view.addViewer('Matched Molecular Pairs Analysis');
    viewer.setOptions({moleculesColumnName: molecules.name, activities: activityColsNames, diffTypes: diffTypes,
      scalings: scalings, fragmentCutoff});
    viewer.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/refs/heads/master/help/datagrok/solutions/domains/chem/chem.md#matched-molecular-pairs';
  }

  @grok.decorators.func({
    name: 'Scaffold Tree Filter',
    description: 'Scaffold Tree filter',
    outputs: [{name: 'result', type: 'filter'}],
    meta: {semType: 'Molecule', allowMultipleFiltersForColumn: 'false', role: 'filter'},
  })
  static scaffoldTreeFilter(): ScaffoldTreeFilter {
    return new ScaffoldTreeFilter();
  }

  @grok.decorators.func()
  static async getScaffoldTree(data: DG.DataFrame,
    @grok.decorators.param({type: 'int', options: {description: 'Ignore molecules with # rings > N', initialValue: '10'}}) ringCutoff: number = 0,
    @grok.decorators.param({options: {description: 'Remove charges and radicals from scaffolds', initialValue: 'false'}}) dischargeAndDeradicalize: boolean = false,
  ): Promise<string> {
    const molColumn = data.columns.bySemType(DG.SEMTYPE.MOLECULE);

    const categories = molColumn?.categories;
    if (categories?.length === 1 && !categories[0]?.trim())
      throw new Error('Molecule column is empty.');

    const smiles = molColumn?.meta.units === DG.UNITS.Molecule.SMILES;
    const invalid: number[] = [];
    const smilesList: string[] = [];
    for (let rowI = 0; rowI < molColumn!.length; rowI++) {
      const mol = molColumn?.get(rowI);
      try {
        smilesList[smilesList.length] = smiles ?
          mol :
          PackageFunctions.convertMolNotation(mol, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
      } catch {
        invalid[invalid.length] = rowI;
      }
    }
    const smilesColumn: DG.Column = DG.Column.fromStrings('smiles', smilesList);
    const dataFrame: DG.DataFrame = DG.DataFrame.fromColumns([smilesColumn]);
    const scriptBlob = await Chem.generateScaffoldTree(dataFrame, smilesColumn!.name, ringCutoff, dischargeAndDeradicalize);
    const scriptRes = new TextDecoder().decode(scriptBlob.data);
    return scriptRes;
  }


  @grok.decorators.func({
    name: 'filterMoleculeDuplicates',
  })
  static removeDuplicates(molecules: string[], molecule: string): string[] {
    const mol1 = checkMoleculeValid(molecule);
    if (!mol1)
      throw new Error(`Molecule is possibly malformed`);
    const filteredMolecules = molecules.filter((smiles) => !checkMolEqualSmiles(mol1, smiles));
    mol1.delete();
    return filteredMolecules;
  }


  @grok.decorators.func({
    name: 'Demo Chem Overview',
    description: 'Overview of Cheminformatics functionality',
    meta: {isDemoScript: 'true', demoSkip: 'GROK-14320', demoPath: 'Cheminformatics | Overview'}})
  static async demoChemOverview(): Promise<void> {
    await _demoChemOverview();
  }


  @grok.decorators.func({
    name: 'Demo Similarity Search',
    description: 'Searching for most similar or diverse molecules in dataset',
    meta: {demoPath: 'Cheminformatics | Similarity & Diversity Search'},
  })
  static async demoSimilarityDiversitySearch(): Promise<void> {
    await _demoSimilarityDiversitySearch();
  }

  @grok.decorators.func({
    name: 'Demo Matched Molecular Pairs',
    description: 'Detect matched molecule pairs calculate the difference in activity values between them',
    meta: {demoPath: 'Cheminformatics | Matched Molecular Pairs'},
  })
  static async demoMMPA(): Promise<void> {
    await _demoMMPA();
  }

  @grok.decorators.func({
    name: 'Demo R Group Analysis',
    description: 'R Group Analysis including R-group decomposition and  visual analysis of the obtained R-groups',
    meta: {demoPath: 'Cheminformatics | R-Group Analysis', demoSkip: 'GROK-14320', isDemoDashboard: 'true'},
  })
  static async demoRgroupAnalysis(): Promise<void> {
    await _demoRGroups();
  }


  @grok.decorators.func({
    name: 'Demo Activity Cliffs',
    description: 'Searching similar structures with significant activity difference',
    meta: {demoPath: 'Cheminformatics | Molecule Activity Cliffs', demoSkip: 'GROK-14320', isDemoDashboard: 'true'},
  })
  static async demoMoleculeActivityCliffs(): Promise<void> {
    await _demoActivityCliffsLayout();
  }

  @grok.decorators.func({
    name: 'Demo Chemical Space',
    description: 'Maps the dataset to 2D plot based on similarity',
    meta: {demoPath: 'Cheminformatics | Chemical Space', demoSkip: 'GROK-14320'},
  })
  static async demoChemicalSpace(): Promise<void> {
    await _demoChemicalSpace();
  }


  @grok.decorators.func({
    name: 'Demo Scaffold Tree',
    description: 'Running scaffold analysis with hierarchical tree',
    meta: {demoPath: 'Cheminformatics | Scaffold Tree'},
  })
  static async demoScaffold(): Promise<void> {
    await _demoScaffoldTree();
  }

  @grok.decorators.func({
    'name': 'Names To Smiles',
    'top-menu': 'Chem | Transform | Names To Smiles...',
    'meta': {'role': 'transform'},
  })
  static async namesToSmiles(
    data: DG.DataFrame,
    @grok.decorators.param({type: 'column'}) names: DG.Column<string>): Promise<void> {
    const namesList = names.toList();
    const res = await grok.functions.call('Chembl:namesToSmiles', {names: namesList});
    const col = res.col('canonical_smiles');
    col.meta.units = DG.UNITS.Molecule.SMILES;
    col.semType = DG.SEMTYPE.MOLECULE;
    data.columns.add(col);
  }

  @grok.decorators.func({
    outputs: [{name: 'smiles', type: 'string', options: {semType: 'Molecule'}}],
    meta: {role: 'canonicalizer'},
  })
  static canonicalize(
    @grok.decorators.param({type: 'string', options: {semType: 'Molecule'}}) molecule: string): string {
    return PackageFunctions.convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
  }

  @grok.decorators.func({
    outputs: [{name: 'molecularFormula', type: 'string'}],
  })
  static getMolecularFormula(
    @grok.decorators.param({type: 'string', options: {semType: 'Molecule'}}) molecule: string): string {
    return oclMol(molecule).getMolecularFormula().formula;
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static validateMolecule(s: string): string | null {
    let logHandle: RDLog | null = null;
    let mol: RDMol | null = null;
    try {
      logHandle = _rdKitModule.set_log_capture('rdApp.error');
      mol = getMolSafe(s, {}, _rdKitModule, true).mol;
      let logBuffer = logHandle?.get_buffer();
      logHandle?.clear_buffer();
      if (!mol && !logBuffer)
        logBuffer = 'unknown error';
      return logBuffer ?? null;
    } finally {
      logHandle?.delete();
      mol?.delete();
    }
  }

  static async getContainer() {
    if (!container)
      container = await grok.dapi.docker.dockerContainers.filter('chemprop').first();
    return container;
  }

  static async getChempropError(response: Response): Promise<string> {
    const match = (await response.text()).match(/[\w.]+Error:\s*(.*)/);
    return match ? match[1] : response.statusText;
  }

  static async trainModelChemprop(table: string,
    predict: string,
    parameterValues: Record<string, any>): Promise<Uint8Array> {
    let chempropContainer = await PackageFunctions.getContainer();

    const body = {
      type: 'Chemprop',
      table: table,
      predict: predict,
      parameters: parameterValues,
    };

    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/modeling/train_chemprop', {
      method: 'POST',
      body: JSON.stringify(body),
      headers: {'Content-Type': 'application/json'},
    });

    if (response.status !== 201) {
      chempropContainer = await grok.dapi.docker.dockerContainers.find(chempropContainer.id);
      const started = chempropContainer.status.startsWith('started') || chempropContainer.status.startsWith('checking');
      if (!started)
        throw new Error(`Failed to start container: ${container.friendlyName}`);
      throw new Error(await PackageFunctions.getChempropError(response));
    }
    return new Uint8Array(await response.arrayBuffer());
  }

  static async applyModelChemprop(modelBlob: Uint8Array, table: string): Promise<DG.Column> {
    let chempropContainer = await PackageFunctions.getContainer();

    const body = {
      modelBlob: Array.from(modelBlob),
      table: table,
    };

    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/modeling/predict_chemprop', {
      method: 'POST',
      body: JSON.stringify(body),
      headers: {'Content-Type': 'application/json'},
    });

    if (response.status !== 201) {
      chempropContainer = await grok.dapi.docker.dockerContainers.find(chempropContainer.id);
      const started = chempropContainer.status.startsWith('started') || chempropContainer.status.startsWith('checking');
      if (!started)
        throw new Error(`Failed to start container: ${container.friendlyName}`);
      throw new Error(await PackageFunctions.getChempropError(response));
    }

    const data = await response.json();
    return DG.Column.fromStrings('outcome', data['outcome'].map((v: any) => v?.toString()));
  }

  @grok.decorators.func({
    name: 'trainChemprop',
    description: 'To be added',
    meta: {mlname: 'Chemprop', mlrole: 'train'},
    outputs: [{name: 'model', type: 'dynamic'}],
  })
  static async trainChemprop(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({options: {category: 'General', choices: ['regression', 'classification'], initialValue: 'regression', description: 'Type of dataset,e.g. classification or regression. This determines the loss function used during training.'}})
    dataset_type: string,
    @grok.decorators.param({options: {category: 'General', choices: ['mse', 'mae', 'rmse', 'bounded-mse', 'bounded-mae', 'bounded-rmse', 'r2', 'binary-mcc', 'multiclass-mcc', 'roc', 'prc', 'accuracy', 'f1'], initialValue: 'rmse', description: 'Metric to use during evaluation. Note:Does NOT affect loss function used during training (loss is determined by the `dataset_type` argument).'}})
    metric: string,
    @grok.decorators.param({type: 'int', options: {category: 'General', initialValue: '3', description: 'Number of classes when running multiclass classification'}})
    multiclass_num_classes: number,
    @grok.decorators.param({type: 'int', options: {category: 'General', initialValue: '1', description: 'Number of folds when performing cross validation'}})
    num_folds: number,
    @grok.decorators.param({type: 'int', options: {category: 'General', initialValue: '0', description: 'Random seed to use when splitting data into train/val/test sets. When `num_folds` > 1,the first fold uses this seed and all subsequent folds add 1 to the seed.'}})
    data_seed: number,
    @grok.decorators.param({type: 'list', options: {category: 'General', initialValue: '[0.8,0.1,0.1]', description: 'Split proportions for train/validation/test sets'}})
    split_sizes: any,
    @grok.decorators.param({options: {category: 'General', choices: ['random', 'scaffold_balanced', 'cv', 'cv_no_val', 'kennard_stone', 'kmeans', 'random_with_repeated_smiles'], initialValue: 'random', description: 'Method of splitting the data into train/val/test'}})
    split_type: string,
    @grok.decorators.param({options: {category: 'Model', choices: ['ReLU', 'LeakyReLU', 'PReLU', 'tanh', 'SELU', 'ELU'], initialValue: 'ReLU', description: 'Activation function'}})
    activation: string,
    @grok.decorators.param({options: {category: 'Model', initialValue: 'false', description: 'Use messages on atoms instead of messages on bonds'}})
    atom_messages: boolean,
    @grok.decorators.param({options: {category: 'Model', initialValue: 'false', description: 'Whether to add bias to linear layers'}})
    message_bias: boolean,
    @grok.decorators.param({type: 'int', options: {category: 'Model', initialValue: '1', description: 'Number of models in ensemble'}})
    ensemble_size: number,
    @grok.decorators.param({type: 'int', options: {category: 'Model', initialValue: '300', description: 'Dimensionality of hidden layers in MPN'}})
    message_hidden_dim: number,
    @grok.decorators.param({type: 'int', options: {category: 'Model', initialValue: '3', description: 'Number of message passing step'}})
    depth: number,
    @grok.decorators.param({options: {category: 'Model', initialValue: '0.0', description: 'Dropout probability'}})
    dropout: number,
    @grok.decorators.param({type: 'int', options: {category: 'Model', initialValue: '300', description: 'Hidden dim for higher-capacity FFN (defaults to hidden_size)'}})
    ffn_hidden_dim: number,
    @grok.decorators.param({type: 'int', options: {category: 'Model', initialValue: '2', description: 'Number of layers in FFN after MPN encoding'}})
    ffn_num_layers: number,
    @grok.decorators.param({type: 'int', options: {category: 'Training', initialValue: '50', description: 'Number of epochs to run'}})
    epochs: number,
    @grok.decorators.param({type: 'int', options: {category: 'Training', initialValue: '64', description: 'Batch size'}})
    batch_size: number,
    @grok.decorators.param({options: {category: 'Training', initialValue: '2.0', description: 'Number of epochs during which learning rate increases linearly from init_lr to max_lr. Afterwards,learning rate decreases exponentially from max_lr to final_lr.'}})
    warmup_epochs: number,
    @grok.decorators.param({options: {category: 'Training', initialValue: '0.001', description: 'Initial learning rate'}})
    init_lr: number,
    @grok.decorators.param({options: {category: 'Training', initialValue: '0.001', description: 'Maximum learning rate'}})
    max_lr: number,
    @grok.decorators.param({options: {category: 'Training', initialValue: '0.0001', description: 'Final learning rate'}})
    final_lr: number,
    @grok.decorators.param({options: {category: 'Training', initialValue: 'false', description: 'Turn off scaling of features'}})
    no_descriptor_scaling: boolean): Promise<Uint8Array | undefined> {
    const parameterValues = {
      'dataset_type': dataset_type,
      'metric': metric,
      'multiclass_num_classes': multiclass_num_classes,
      'activation': activation,
      'atom_messages': atom_messages,
      'batch_size': batch_size,
      'message_bias': message_bias,
      'depth': depth,
      'dropout': dropout,
      'ensemble_size': ensemble_size,
      'epochs': epochs,
      'ffn_hidden_dim': ffn_hidden_dim,
      'ffn_num_layers': ffn_num_layers,
      'final_lr': final_lr,
      'message_hidden_dim': message_hidden_dim,
      'init_lr': init_lr,
      'max_lr': max_lr,
      'no_descriptor_scaling': no_descriptor_scaling,
      'num_folds': num_folds,
      'data_seed': data_seed,
      'split_sizes': split_sizes,
      'split_type': split_type,
      'warmup_epochs': warmup_epochs,
    };
    predictColumn.name = df.columns.getUnusedName(predictColumn.name);
    df.columns.add(predictColumn);
    try {
      const modelBlob = await PackageFunctions.trainModelChemprop(df.toCsv(), predictColumn.name, parameterValues);
      const zip = new JSZip();
      const archive = await zip.loadAsync(modelBlob);
      const file = archive.file('blob.bin');
      const binBlob = await file?.async('uint8array')!;
      return binBlob;
    } catch (error: any) {
      grok.shell.error(error);
    }
  }

  @grok.decorators.func({
    meta: {mlname: 'Chemprop', mlrole: 'apply'},
    outputs: [{name: 'data_out', type: 'dataframe'}],
  })
  static async applyChemprop(df: DG.DataFrame,
    @grok.decorators.param({type: 'dynamic'}) model: Uint8Array) {
    try {
      const column = await PackageFunctions.applyModelChemprop(model, df.toCsv());
      return DG.DataFrame.fromColumns([column]);
    } catch (error: any) {
      grok.shell.error(error);
    }
  }

  @grok.decorators.func({
    meta: {mlname: 'Chemprop', mlrole: 'isApplicable'},
    outputs: [{name: 'result', type: 'bool'}],
  })
  static async isApplicableNN(df: DG.DataFrame, predictColumn: DG.Column) {
    if (df.columns.length > 1)
      return false;
    const featureColumn = df.columns.byIndex(0);
    if (featureColumn.semType != 'Molecule')
      return false;
    if (!predictColumn.matches('numerical'))
      return false;
    return true;
  }


  @grok.decorators.func({
    meta: {mlname: 'Chemprop', mlrole: 'isInteractive', mlupdate: 'false'},
    outputs: [{name: 'result', type: 'bool'}],
  })
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  static async isInteractiveNN(df: DG.DataFrame, predictColumn: DG.Column) {
    return true;
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Transform | Deprotect...',
    'name': 'Deprotect',
    'description': 'Removes drawn protecting groups / fragments from molecules',
    'editor': 'Chem:DeprotectEditor',
    'meta': {'role': 'transform'},
  })
  static async deprotect(
    @grok.decorators.param({options: {description: 'Input data table'}}) table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({options: {semType: 'Molecule', initialValue: 'O=C([N:1])OCC1c2ccccc2-c2ccccc21'}}) fragment: string): Promise<void> {
    const rdModule = PackageFunctions.getRdKitModule();
    addDeprotectedColumn(table, molecules, fragment, rdModule);
  }

  @grok.decorators.editor({
    name: 'Deprotect Editor',
    outputs: [{name: 'result', type: 'widget'}],
  })
  static deprotectEditor(call: DG.FuncCall): DG.Widget {
    return new DeprotectEditor(call);
  }

  @grok.decorators.func({
    name: 'beautifyMols',
    description: 'Beautifies the list of molecules and returns the list of beautified molecules',
  })
  static async beautifyMols(mols: string[]): Promise<string[]> {
    return await (await chemCommonRdKit.getRdKitService()).beautifyMolsV3K(mols);
  }

  @grok.decorators.func({
    description: 'Converts the list of molecules to V3K format using OCL',
  })
  static async convertToV3KViaOCL(mols: string[]): Promise<string[]> {
    const oc = new OCLService();
    const result = await oc.molfileToV3K(mols);
    oc.terminate();
    return result;
  }

  @grok.decorators.func({
    'name': 'mpo',
    'top-menu': 'Chem | Calculate | MPO Score...',
    'description': 'Calculates the MPO score for the column of molecules',
  })
  static async _mpo(): Promise<void> {
    const mpoDialog = new MpoProfileDialog();
    await mpoDialog.showDialog();
  }

  @grok.decorators.func({
    outputs: [{name: 'res', type: 'list'}],
    meta: {role: 'transform'},
  })
  static async mpoTransformFunction(
    df: DG.DataFrame,
    profileName: string,
    aggregation: WeightedAggregation,
    @grok.decorators.param({type: 'object'}) currentProperties: { [key: string]: PropertyDesirability },
  ): Promise<string[]> {
    const result = calculateMpoCore(df, profileName, currentProperties, aggregation);

    for (const warning of result.warnings)
      grok.shell.warning(warning);

    if (result.error)
      grok.shell.error(result.error);

    return result.columnNames;
  }

  @grok.decorators.fileViewer({
    fileViewer: 'json',
    fileViewerCheck: 'Chem:checkJsonMpoProfile',
  })
  static mpoProfileEditor(file: DG.FileInfo): DG.View {
    const view = DG.View.create();
    const saveButton = ui.bigButton('SAVE', () => {});
    saveButton.style.display = 'none';
    view.name = file.name;
    view.setRibbonPanels([[saveButton]]);

    file.readAsString().then((s) => {
      const mpoEditor = new MpoProfileEditor(undefined, true);
      mpoEditor.setProfile(JSON.parse(s));
      view.append(mpoEditor.root);

      mpoEditor.onChanged.subscribe((_) => saveButton.style.display = 'initial');
      saveButton.onclick = () => {
        grok.dapi.files.writeAsText(file, JSON.stringify(mpoEditor.getProfile()));
        saveButton.style.display = 'none';
      };
    });

    return view;
  }

  @grok.decorators.func({
    outputs: [{name: 'result', type: 'bool'}],
  })
  static checkJsonMpoProfile(content: string) {
    return JSON.parse(content)['type'] === 'MPO Desirability Profile';
  }

  @grok.decorators.panel({
    name: 'Chemistry | Mixture',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async mixtureWidget(
    @grok.decorators.param({type: 'string', options: {semType: 'ChemicalMixture'}}) mixture: string): Promise<DG.Widget> {
    return await createMixtureWidget(mixture);
  }

  @grok.decorators.panel({
    name: 'Chemistry | MixtureTree',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async mixtureTreeWidget(
    @grok.decorators.param({options: {semType: 'ChemicalMixture'}}) mixture: string): Promise<DG.Widget> {
    const mixtureObj = JSON.parse(mixture) as Mixfile;
    const resDiv = ui.divV([]);
    resDiv.append(ui.divText(`mixfileVersion: ${mixtureObj.mixfileVersion}`));
    if (mixtureObj.contents && mixtureObj.contents.length) {
      const contentsAcc = ui.accordion('contents');
      for (let i = 0; i < mixtureObj.contents.length; i++)
        contentsAcc.addPane(mixtureObj.contents![i].name ?? `component ${i + 1}`, () => createComponentPane(mixtureObj.contents![i]));
      resDiv.append(contentsAcc.root);
    }
    return new DG.Widget(resDiv);
  }

  @grok.decorators.func({
    'name': 'Biochemical Properties',
    'description': 'Dynamically discovers and executes tagged biochemical calculators',
    'top-menu': 'Chem | Calculate | Biochemical Properties',
  })
  static async biochemPropsWidget(): Promise<void> {
    await biochemicalPropertiesDialog();
  }

  @grok.decorators.app({
    'name': 'MPO profiles',
    'meta': {browsePath: 'Chem', icon: 'images/mpo.png'},
  })
  static async mpoProfilesApp(
    @grok.decorators.param({options: {metaUrl: true, optional: true}}) path?: string,
  ): Promise<DG.View> {
    const url = new URL(window.location.href);
    const params = url.searchParams;

    const hasPath = !!path;

    if (hasPath && url.pathname.endsWith('/Mpo/create-profile')) {
      const view = new MpoProfileCreateView();
      return view.tableView!;
    }

    const profileId = params.get('profileId');
    const profiles = await loadMpoProfiles();
    if (hasPath && profileId) {
      const profile = profiles.find((p) => p.name === decodeURIComponent(profileId));
      const view = new MpoProfileCreateView(profile, false, profile?.fileName);
      return view.view;
    }

    const infoView = new MpoProfilesView();
    await infoView.render();
    return infoView.view;
  }

  @grok.decorators.func()
  static async mpoProfilesAppTreeBrowser(
    @grok.decorators.param({type: 'dynamic'}) treeNode: DG.TreeViewGroup,
    @grok.decorators.param({type: 'view'}) browseView: any,
  ) {
    let openedView: DG.ViewBase | null = null;

    const refresh = async () => {
      treeNode.items.forEach((item) => item.remove());

      const profileFiles = await grok.dapi.files.list(MPO_TEMPLATE_PATH);
      for (const profile of profileFiles) {
        const content =
        JSON.parse(await grok.dapi.files.readAsText(
          `${MPO_TEMPLATE_PATH}/${profile.name}`,
        )) as DesirabilityProfile;

        treeNode.item(content.name).onSelected.subscribe(() => {
          openedView?.close();
          const editView = new MpoProfileCreateView(content, false, profile.name);
          openedView = editView.view;
          grok.shell.addPreview(editView.view);
        });
      }
    };

    await refresh();

    grok.events.onCustomEvent(MPO_PROFILE_CHANGED_EVENT).subscribe(async () => {
      await refresh();
    });
  }

  @grok.decorators.panel({
    name: 'Chemistry | MPO',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async mpoWidget(
  @grok.decorators.param({name: 'smiles', options: {semType: 'Molecule'}})
    semValue: DG.SemanticValue,
  ): Promise<DG.Widget<any>> {
    const container = ui.divV([]);
    const resultDiv = ui.div();
    const loader = ui.loader();
    resultDiv.appendChild(loader);

    const dataFrame = semValue.cell.dataFrame;
    const profiles = await loadMpoProfiles();
    const suitableProfiles = await findSuitableProfiles(dataFrame, profiles);

    if (suitableProfiles.length === 0) {
      container.appendChild(ui.divText('No suitable profiles available.'));
      const createButton = ui.button('Manage Profiles', async () => {
        const profilesView = new MpoProfilesView();
        await profilesView.show();
      });
      container.appendChild(createButton);

      return DG.Widget.fromRoot(container);
    }

    const profileInput = ui.input.choice('Profile', {
      items: suitableProfiles.map((p) => p.fileName),
      value: suitableProfiles[0].fileName,
      onValueChanged: () => calculateMpo(),
    });

    const aggregationInput = ui.input.choice('Aggregation', {
      items: WEIGHTED_AGGREGATIONS_LIST,
      nullable: false,
      onValueChanged: () => calculateMpo(),
    });

    async function calculateMpo() {
      if (!profileInput.value)
        return;

      const selected = suitableProfiles.find((p) => p.fileName === profileInput.value);
      if (!selected)
        return;

      ui.empty(resultDiv);
      resultDiv.appendChild(loader);

      const columns: DG.Column[] = [];

      for (const propertyName in selected.properties) {
        const column = dataFrame.columns.byName(propertyName);
        if (!column)
          continue;

        columns.push(column);

        const newTagValue = JSON.stringify(selected.properties[propertyName]);
        const existingTagValue = column.getTag('desirabilityTemplate');
        if (existingTagValue !== newTagValue)
          column.setTag('desirabilityTemplate', newTagValue);
      }

      const score = await mpo( //@ts-ignore
        dataFrame,
        columns,
        selected.name,
        aggregationInput.value as WeightedAggregation,
        false,
      );

      ui.empty(resultDiv);

      const mpoValue = score.get(semValue.cell.rowIndex);
      const valueHost = ui.divH([], {style: {position: 'relative', alignItems: 'center'}});

      const addColumnIcon = ui.iconFA(
        'plus',
        async () => {
          await mpo( //@ts-ignore
            dataFrame,
            columns,
            selected.name,
            aggregationInput.value as WeightedAggregation,
          );
        },
        'Add MPO score as a column',
      );

      ui.tools.setHoverVisibility(valueHost, [addColumnIcon]);
      $(addColumnIcon).css('color', '#2083d5').css('margin-right', '4px');

      valueHost.appendChild(addColumnIcon);
      valueHost.appendChild(ui.divText(mpoValue != null ? mpoValue.toFixed(2) : 'N/A'));

      const table = ui.tableFromMap({'MPO Score': valueHost});
      resultDiv.appendChild(table);
    }

    container.appendChild(ui.divV([profileInput.root, aggregationInput.root, resultDiv]));

    await calculateMpo();
    return DG.Widget.fromRoot(container);
  }
}
