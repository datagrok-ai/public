/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {BioTags, NOTATION, NOTATION_PROVIDER_CONSTRUCTOR_ROLE} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

import {OligoToolkitPackage} from './apps/common/model/oligo-toolkit-package';
import {APP_NAME} from './apps/common/view/const';
import {getSpecifiedAppUI} from './apps/common/view/utils';
import {CombinedAppUI} from './apps/common/view/combined-app-ui';
import {linkStrandsV3000} from './apps/structure/model/mol-transformations';
import {SequenceToMolfileConverter} from './apps/structure/model/sequence-to-molfile';
import {demoOligoPatternUI, demoOligoStructureUI, demoOligoTranslatorUI} from './demo/demo-st-ui';
import {getExternalAppViewFactories} from './plugins/mermade';
import {defaultErrorHandler} from './utils/err-info';

import {OligoNucleotideCellRenderer} from './oligo-renderer/cell-renderer';
import {buildOligoPanel} from './oligo-renderer/legend-panel';
import {buildOligoStructuresPanel} from './oligo-renderer/structures-panel';
import {combineSenseAntisenseToOligo, convertHelmColumnToOligo} from './oligo-renderer/converters';
import {
  openOligoCanvasDialog, openOligoHelmEditorDialog,
  copyHelmToClipboard, copyDuplexImageToClipboard,
} from './oligo-renderer/cell-actions';

import {polyToolConvert, polyToolConvertUI} from './polytool/pt-dialog';
import {polyToolEnumerateChemApp, polyToolEnumerateChemUI} from './polytool/pt-chem-enum-dialog';
import {polyToolEnumerateHelmUI, polyToolEnumerateSeq} from './polytool/pt-enumerate-seq-dialog';
import {_setPeptideColumn} from './polytool/utils';
import {PolyToolCsvLibHandler} from './polytool/csv-to-json-monomer-lib-converter';
import {ITranslationHelper} from './types';
import {PolyToolConvertFuncEditor} from './polytool/pt-convert-editor';
import {CyclizedNotationProvider} from './utils/cyclized';
import {getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {PolyToolDataRole, PolyToolTags} from './consts';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getPTCombineDialog} from './polytool/pt-combine-dialog';
import {PolyToolEnumeratorTypes} from './polytool/types';
import {splitterAsHelm} from '@datagrok-libraries/bio/src/utils/macromolecule';

export * from './package.g';

/** Temporary polyfill */

function getDecoratorFunc() {
  return function(args: any) {
    return function(
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
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
  'fileViewer', 'model', 'treeBrowser', 'polyfill'
];

decorators.forEach((decorator) => {
  if (!(grok.decorators as any)[decorator])
    (grok.decorators as any)[decorator] = getDecoratorFunc();
});

/** End temporary polyfill */

export const _package: OligoToolkitPackage = new OligoToolkitPackage({debug: true}/**/);

let initSequenceTranslatorPromise: Promise<void> | null = null;

async function initSequenceTranslatorInt(): Promise<void> {
  const [helmHelper, bioLibHelper] = await Promise.all([
    getHelmHelper(),
    getMonomerLibHelper(),
  ]);
  _package.completeInit(helmHelper, bioLibHelper.getMonomerLib());
}

export class PackageFunctions {
  @grok.decorators.app({
    icon: 'img/icons/toolkit.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Toolkit',
    tags: ['app']
  })
  static async oligoToolkitApp(): Promise<DG.ViewBase> {
    await _package.initLibData();
    const externalViewFactories = await getExternalAppViewFactories(_package);
    if (!externalViewFactories)
      throw new Error('External app view factories not loaded');
    const appUI = new CombinedAppUI(externalViewFactories!, _package);
    const view = await appUI.getAppView();
    return view;
  }


  @grok.decorators.init({tags: ['init']})
  static async init(): Promise<void> {
    if (initSequenceTranslatorPromise === null)
      _package.startInit(initSequenceTranslatorPromise = initSequenceTranslatorInt());

    return initSequenceTranslatorPromise;
  }

  @grok.decorators.app({
    icon: 'img/icons/translator.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Translator',
    tags: ['app']
  })
  static async oligoTranslatorApp(): Promise<DG.ViewBase> {
    const view = await getSpecifiedAppView(APP_NAME.TRANSLATOR);
    return view;
  }


  @grok.decorators.app({
    icon: 'img/icons/pattern.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Pattern',
    tags: ['app']
  })
  static async oligoPatternApp(): Promise<DG.ViewBase> {
    const view = await getSpecifiedAppView(APP_NAME.PATTERN);
    return view;
  }


  @grok.decorators.app({
    icon: 'img/icons/structure.png',
    browsePath: 'Peptides | Oligo Toolkit',
    name: 'Oligo Structure',
    tags: ['app']
  })
  static async oligoStructureApp(): Promise<DG.ViewBase> {
    const view = await getSpecifiedAppView(APP_NAME.STRUCTURE);
    return view;
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static async getTranslationHelper(): Promise<ITranslationHelper> {
    await _package.initLibData();
    return _package;
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static getCodeToWeightsMap(): Record<string, number> {
    const monomerLibWrapper = _package.monomerLibWrapper;
    const map = monomerLibWrapper.getCodesToWeightsMap();
    return Object.fromEntries(map);
  }


  @grok.decorators.func()
  static validateSequence(
    sequence: string): boolean {
    const validator = _package.createSequenceValidator(sequence);
    const format = _package.createFormatDetector(sequence).getFormat();
    return (format === null) ? false : validator.isValidSequence(format!);
  }


  @grok.decorators.func({
    name: 'validateSequence'
  })
  static getMolfileFromGcrsSequence(
    sequence: string,
    invert: boolean): string {
    return (new SequenceToMolfileConverter(sequence, invert, 'GCRS')).convert();
  }


  @grok.decorators.func()
  static linkStrands(
    @grok.decorators.param({type: 'object'}) strands: { senseStrands: string[], antiStrands: string[] }): string {
    return linkStrandsV3000(strands, true);
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Bioinformatics | Oligo Toolkit | Translator',
      path: '/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Translate',
      demoSkip: 'GROK-14320'
    },
    name: 'demoOligoTranslator',
    description: 'Translate oligonucleotide sequences across various formats accepted by different synthesizers'
  })
  static async demoTranslateSequence(): Promise<void> {
    await demoOligoTranslatorUI();
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Bioinformatics | Oligo Toolkit | Pattern',
      path: '%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex'
    },
    description: 'Design a modification pattern for an oligonucleotide sequence'
  })
  static async demoOligoPattern(): Promise<void> {
    await demoOligoPatternUI();
  }


  @grok.decorators.func({
    meta: {
      demoPath: 'Bioinformatics | Oligo Toolkit | Structure',
      path: '%20/apps/Tutorials/Demo/Bioinformatics/Oligonucleotide%20Sequence:%20Visualize%20duplex'
    },
    description: 'Visualize duplex and save SDF'
  })
  static async demoOligoStructure(): Promise<void> {
    await demoOligoStructureUI();
  }


  @grok.decorators.func()
  static async translateOligonucleotideSequence(
    sequence: string, sourceFormat: string, targetFormat: string
  ): Promise<string> {
    await _package.initLibData();
    return _package.createFormatConverter(sequence, sourceFormat).convertTo(targetFormat);
  }

  @grok.decorators.func({
    'top-menu': 'Bio | PolyTool | Convert...',
    'name': 'polyToolConvert',
    'description': 'editor for Performing conversion of sequences in custom notation to molfiles'
  })
  static async polyToolConvertTopMenu(): Promise<void> {
    await polyToolConvertUI();
  }

  @grok.decorators.editor({tags: ['editor']})
  static async getPolyToolConvertEditor(
    call: DG.FuncCall): Promise<DG.Column<string> | null> {
    const funcEditor = await PolyToolConvertFuncEditor.create(call);
    return await funcEditor.showDialog();
  }


  @grok.decorators.func({
    editor: 'SequenceTranslator:getPolyToolConvertEditor',
  })
  static async polyToolConvert2(
    table: DG.DataFrame,
    @grok.decorators.param({options: {caption: 'Sequence'}}) seqCol: DG.Column,
    @grok.decorators.param({options: {initialValue: 'true'}}) generateHelm: boolean,
    @grok.decorators.param({options: {initialValue: 'true'}}) chiralityEngine: boolean,
    @grok.decorators.param({type: 'object'}) rules: string[]
  ): Promise<DG.Column<string>> {
    const ptConvertRes = await polyToolConvert(seqCol, generateHelm, false, chiralityEngine, false, rules);
    return ptConvertRes[0];
  }


  @grok.decorators.func({
    'top-menu': 'Bio | PolyTool | Enumerate HELM...',
    'name': 'polyToolEnumerateHelm',
    'description': 'Dialog for configuring enumeration of a HELM sequence'
  })
  static async polyToolEnumerateHelmTopMenu(): Promise<void> {
    await polyToolEnumerateHelmUI(grok.shell.tv?.dataFrame.currentCell);
  }


  @grok.decorators.func({
    'top-menu': 'Chem | Transform | Markush Enumeration...',
    'name': 'Markush Enumerator',
    'description': 'Enumerate cores and R-group lists into a molecule table (Zip or Cartesian)'
  })
  static async chemEnumerateMarkushTopMenu(): Promise<void> {
    let cell: DG.Cell | undefined = undefined;
    if (grok.shell.tv?.dataFrame && grok.shell.tv.dataFrame.currentCell && grok.shell.tv.dataFrame.currentCell.column.semType === DG.SEMTYPE.MOLECULE)
      cell = grok.shell.tv.dataFrame.currentCell;
    polyToolEnumerateChemUI(cell);
  }


  @grok.decorators.func()
  static async polyToolColumnChoice(
    @grok.decorators.param({options: {description: 'Input data table'}}) df: DG.DataFrame,
      macroMolecule: DG.Column): Promise<void> {
    _setPeptideColumn(macroMolecule);
    await grok.data.detectSemanticTypes(df);
  }


  @grok.decorators.func()
  static async createMonomerLibraryForPolyTool(
    file: DG.FileInfo) {
    const fileContent = await file.readAsString();
    const libHandler = new PolyToolCsvLibHandler(file.fileName, fileContent);
    const libObject = await libHandler.getJson();
    const jsonFileName = file.fileName.replace(/\.csv$/, '.json');
    const jsonFileContent = JSON.stringify(libObject, null, 2);
    DG.Utils.download(jsonFileName, jsonFileContent);
  }


  @grok.decorators.func({
    meta: {
      icon: 'img/icons/structure.png',
      browsePath: 'Peptides | PolyTool',
      role: 'app'
    },
    name: 'HELM Enumerator',
    tags: ['app']
  })
  static async ptEnumeratorHelmApp(): Promise<void> {
    await polyToolEnumerateHelmUI();
  }


  @grok.decorators.func({
    meta: {
      icon: 'img/icons/structure.png',
      browsePath: 'Chem',
      role: 'app'
    },
    name: 'Markush Enumerator',
    tags: ['app'],
    outputs: [{type: 'view', name: 'result'}]
  })
  static async ptEnumeratorChemApp() {
    return await polyToolEnumerateChemApp();
  }


  @grok.decorators.func({
    name: 'Polytool Helm Enumerator dialog'
  })
  static async getPtHelmEnumeratorDialog(
    @grok.decorators.param({type: 'object', options: {nullable: true}}) cell?: DG.Cell) {
    return polyToolEnumerateHelmUI(cell);
  }


  @grok.decorators.func({
    name: 'Polytool Chem Enumerator dialog'
  })
  static async getPtChemEnumeratorDialog(
    @grok.decorators.param({type: 'object', options: {nullable: true}}) cell?: DG.Cell) {
    return polyToolEnumerateChemUI(cell);
  }


  /** Enumerator entry for OligoNucleotide cells.
   *
   * The cell value is HELM (under the hood). The enumerator dialog is built
   * around `Macromolecule` cells, so we wrap the oligo HELM in a temp
   * Macromolecule column and pass that cell in. The `outputAsOligo` flag
   * makes the dialog tag the enumerated result column as OligoNucleotide so
   * the duplex renderer picks it up automatically. */
  @grok.decorators.func({
    name: 'Polytool Oligo Enumerator dialog'
  })
  static async getPtOligoEnumeratorDialog(
    @grok.decorators.param({type: 'object', options: {nullable: true}}) cell?: DG.Cell) {
    if (!cell || cell.value == null)
      return polyToolEnumerateHelmUI(undefined, true);

    const helm = String(cell.value);
    const tempCol = DG.Column.fromStrings('helm', [helm]);
    tempCol.semType = DG.SEMTYPE.MACROMOLECULE;
    tempCol.meta.units = 'helm';
    tempCol.setTag('aligned', 'SEQ');
    tempCol.setTag('alphabet', 'RNA');
    tempCol.setTag('cell.renderer', 'helm');
    const tempDf = DG.DataFrame.fromColumns([tempCol]);
    const tempCell = tempDf.cell(0, 'helm');
    return polyToolEnumerateHelmUI(tempCell, true);
  }

  @grok.decorators.func({
    name: 'Enumerate Single HELM Sequence',
    description: 'Enumerate provided HELM sequence on provided positions with provided monomers and generates new table',
    outputs: [{type: 'dataframe', name: 'result'}]
  })
  static async enumerateSingleHelmSequence(
    helmSequence: string, positions: number[], monomerLists: string[][], toAtomicLevel: boolean = false
  ): Promise<DG.DataFrame> {
    return await polyToolEnumerateSeq(helmSequence, PolyToolDataRole.macromolecule, null, {
      type: PolyToolEnumeratorTypes.Single,
      placeholders: positions.map((pos, i) => ({position: pos, monomers: monomerLists[i]}))
    }, toAtomicLevel ? {
      generateHelm: true,
      chiralityEngine: true,
      highlightMonomers: false,
      rules: []
    } : false, _package.helmHelper);
  }

  @grok.decorators.func({
    name: 'Enumerate Single HELM Sequence with natural amino acids',
    description: 'Enumerate provided HELM sequence on all positions with natural amino acids and generates new table. Generated table has sequence column called "Enumerated", and molecule column called "Molfile(Enumerated) if toAtomicLevel is set to true. Keywords: Optimize, enumerate, HELM optimization, Maximize Minimize property. When you want to optimize certain peptide using for example logS, set toAtomicLevel to true and use generated molecule column to calculate given property using chem package functions.',
    outputs: [{type: 'dataframe', name: 'result'}]
  })
  static async enumerateSingleHelmSequenceWithNaturalAAs(
    helmSequence: string, toAtomicLevel: boolean = false
  ): Promise<DG.DataFrame> {
    const splitt = splitterAsHelm(helmSequence);
    const l = splitt.length;
    const positions = Array.from({length: l}, (_, i) => i);
    const monomerLists = positions.map((_part) => {
      return [...['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']];
    });
    return await PackageFunctions.enumerateSingleHelmSequence(helmSequence, positions, monomerLists, toAtomicLevel);
  }

  @grok.decorators.func({
    'name': 'Combine Sequences',
    'top-menu': 'Bio | PolyTool | Combine Sequences...'
  })
  static async getPolyToolCombineDialog() {
    getPTCombineDialog();
  }


  // ---- OligoNucleotide cell renderer + panel + converters ----

  @grok.decorators.func({
    name: 'oligoNucleotideCellRenderer',
    description: 'Renders OligoNucleotide (siRNA / ASO) duplex view in grid cells',
    tags: ['cellRenderer'],
    meta: {
      cellType: 'OligoNucleotide',
      columnTags: 'quality=OligoNucleotide',
      role: 'cellRenderer',
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}],
  })
  static oligoNucleotideCellRenderer(): DG.GridCellRenderer {
    return new OligoNucleotideCellRenderer();
  }

  /** Double-click cell editor: opens a full-screen modal with a nicely-rendered
   * canvas view of the duplex. Hover interactions (monomer/linkage tooltip
   * with cached RDKit structures) mirror what works in the grid cell. Editing
   * the HELM itself happens through the separate `Open HELM Editor` action. */
  @grok.decorators.func({
    name: 'editOligoNucleotideCell',
    description: 'OligoNucleotide',
    tags: ['cellEditor'],
    meta: {
      role: 'cellEditor',
    },
  })
  static editOligoNucleotideCell(
    @grok.decorators.param({type: 'grid_cell'}) cell: DG.GridCell,
  ): void {
    openOligoCanvasDialog(cell);
  }

  /** Cell context-menu action: open the HELM Web Editor for the cell's
   * sequence and write the edited HELM back on OK. Lives in the "Actions"
   * group on the cell's context menu (same surfacing convention as
   * `Copy as HELM`). */
  @grok.decorators.func({
    name: 'Open HELM Editor',
    description: 'Edit the oligonucleotide HELM in the HELM Web Editor',
    meta: {'action': 'Edit HELM'},
  })
  static openOligoHelmEditor(
    @grok.decorators.param({options: {semType: 'OligoNucleotide'}}) value: DG.SemanticValue,
  ): Promise<void> {
    return openOligoHelmEditorDialog(value);
  }

  @grok.decorators.func({
    name: 'Oligo-Nucleotide',
    description: 'Modifications, lengths, conjugates and color legend for an OligoNucleotide cell',
    tags: ['panel', 'widgets'],
    outputs: [{type: 'widget', name: 'result'}],
  })
  static oligoNucleotidePanel(
    @grok.decorators.param({type: 'semantic_value', options: {semType: 'OligoNucleotide'}}) value: DG.SemanticValue,
  ): DG.Widget {
    return buildOligoPanel(value);
  }

  @grok.decorators.func({
    name: 'Oligo Structures',
    description: 'Sense and antisense full molecular structures rendered separately',
    tags: ['panel', 'widgets'],
    outputs: [{type: 'widget', name: 'result'}],
  })
  static oligoNucleotideStructuresPanel(
    @grok.decorators.param({type: 'semantic_value', options: {semType: 'OligoNucleotide'}}) value: DG.SemanticValue,
  ): DG.Widget {
    return buildOligoStructuresPanel(value);
  }

  /** Cell context-menu action: copy the raw HELM string. Surfaced automatically
   * by the platform under the cell's "Copy" submenu because of `meta.action:
   * 'Copy as HELM'`; we set `exclude-actions-panel` so it doesn't also show up
   * in the right-side actions panel. */
  @grok.decorators.func({
    name: 'Copy as HELM',
    description: 'Copy the HELM string of an oligo cell to the clipboard',
    meta: {'action': 'Copy as HELM'},
  })
  static copyOligoAsHelm(
    @grok.decorators.param({options: {semType: 'OligoNucleotide'}}) value: DG.SemanticValue,
  ): void {
    copyHelmToClipboard(value);
  }

  /** Cell context-menu action: render the duplex to a high-resolution PNG with
   * transparent background and copy it to the system clipboard. Canvas pixel
   * dimensions are scaled up but the logical layout sees the original
   * gridCell bounds — so chip sizes match what's on-screen, just at higher
   * pixel density. drawDuplex itself never paints a backdrop, which keeps
   * the alpha channel clean. */
  @grok.decorators.func({
    name: 'Copy as Image',
    description: 'Copy a high-resolution image of the oligo duplex',
    meta: {'action': 'Copy as Image'},
  })
  static copyOligoAsImage(
    @grok.decorators.param({options: {semType: 'OligoNucleotide'}}) value: DG.SemanticValue,
  ): void {
    copyDuplexImageToClipboard(value);
  }

  // Invoked from the column / cell context menu via detectors.js (no top-menu).
  @grok.decorators.func({
    name: 'convertHelmToOligoNucleotide',
    description: 'Create a new column tagged as OligoNucleotide so HELM duplex cells render with the oligo view',
  })
  static async convertHelmToOligoNucleotide(
    table: DG.DataFrame,
    @grok.decorators.param({options: {caption: 'HELM column', semType: 'Macromolecule'}}) helmCol: DG.Column,
  ): Promise<DG.Column> {
    const out = convertHelmColumnToOligo(table, helmCol);
    await grok.data.detectSemanticTypes(table);
    return out;
  }

  // Invoked from the column / cell context menu via detectors.js — opens the
  // function editor (`.prepare().edit()`) so the user can pick antisense column.
  @grok.decorators.func({
    name: 'combineSenseAntisenseToOligoNucleotide',
    description: 'Combine separate sense + antisense HELM columns into one OligoNucleotide column',
  })
  static async combineSenseAntisenseToOligoNucleotide(
    table: DG.DataFrame,
    @grok.decorators.param({options: {caption: 'Sense', semType: 'Macromolecule'}}) senseCol: DG.Column,
    @grok.decorators.param({options: {caption: 'Antisense', semType: 'Macromolecule'}}) antiCol: DG.Column,
  ): Promise<DG.Column> {
    const out = combineSenseAntisenseToOligo(table, senseCol, antiCol);
    await grok.data.detectSemanticTypes(table);
    return out;
  }

  @grok.decorators.func({
    name: 'applyNotationProviderForHarmonizedSequence'
  })
  static applyNotationProviderForCyclized(
    @grok.decorators.param({type: 'column'}) col: DG.Column<string>,
      separator: string) {
    col.setTag('aligned', 'SEQ');
    col.setTag('alphabet', 'UN');
    col.setTag('.alphabetIsMultichar', 'true');
    col.meta.units = NOTATION.CUSTOM;
    if (separator && !col.tags[BioTags.separator])
      col.tags[BioTags.separator] = separator;
    col.tags[PolyToolTags.dataRole] = 'template';
    col.temp[SeqTemps.notationProvider] = new CyclizedNotationProvider(separator, _package.helmHelper);
  }

  @grok.decorators.func({
    name: 'harmonizedSequenceNotationProviderConstructor',
    meta: {role: 'notationProviderConstructor'}
  })
  static async harmonizedSequenceNotationProviderConstructor(): Promise<typeof CyclizedNotationProvider> {
    return CyclizedNotationProvider;
  }
}

//name: getSpecifiedAppView
async function getSpecifiedAppView(appName: string): Promise<DG.ViewBase> {
  await _package.initLibData();
  const appUI = getSpecifiedAppUI(appName, _package);
  const view = await appUI.getAppView();
  return view;
}
