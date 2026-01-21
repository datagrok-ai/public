/* eslint-disable max-lines-per-function */
/* eslint-disable rxjs/no-nested-subscribe */
/* eslint-disable max-params */
/* eslint-disable max-len */
/* eslint max-lines: 'off' */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import {DimReductionBaseEditor, PreprocessFunctionReturnType} from '@datagrok-libraries/ml/src/functionEditors/dimensionality-reduction-editor';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics, KnownMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {IMonomerLib, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {FastaFileHandler} from '@datagrok-libraries/bio/src/utils/fasta-handler';
import {SCORE} from '@datagrok-libraries/bio/src/utils/macromolecule/scoring';
import {createJsonMonomerLibFromSdf,} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ActivityCliffsEditor} from '@datagrok-libraries/ml/src/functionEditors/activity-cliffs-function-editor';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {getEmbeddingColsNames, multiColReduceDimensionality} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {ITSNEOptions, IUMAPOptions} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/multi-column-dim-reducer';
import {generateLongSequence, generateLongSequence2} from '@datagrok-libraries/bio/src/utils/generator';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {RDModule as _RDMoule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {ISeqHandler, SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';

import {getMacromoleculeColumns} from './utils/ui-utils';
import {MacromoleculeDifferenceCellRenderer, MacromoleculeSequenceCellRenderer,} from './utils/cell-renderer';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {SequenceAlignment} from './seq_align';
import {getEncodedSeqSpaceCol} from './analysis/sequence-space';
import {createLinesGrid, createPropPanelElement, createTooltipElement,} from './analysis/sequence-activity-cliffs';
import {SequenceSimilarityViewer} from './analysis/sequence-similarity-viewer';
import {SequenceDiversityViewer} from './analysis/sequence-diversity-viewer';
import {invalidateMols, MONOMERIC_COL_TAGS, SubstructureSearchDialog} from './substructure-search/substructure-search';
import {convert} from './utils/convert';
import {getMacromoleculeColumnPropertyPanel} from './widgets/representations';
import {saveAsFastaUI} from './utils/save-as-fasta';
import {BioSubstructureFilter} from './widgets/bio-substructure-filter';
import {WebLogoViewer} from './viewers/web-logo-viewer';
import {MonomerLibManager} from './utils/monomer-lib/lib-manager';
import {getMonomerLibraryManagerLink, showManageLibrariesDialog, showManageLibrariesView} from './utils/monomer-lib/library-file-manager/ui';
import {demoBioSimDiv} from './demo/bio01-similarity-diversity';
import {demoSeqSpace} from './demo/bio01a-hierarchical-clustering-and-sequence-space';
import {demoActivityCliffsCyclic} from './demo/bio01b-hierarchical-clustering-and-activity-cliffs';
import {demoToAtomicLevel} from './demo/bio03-atomic-level';
import {checkInputColumnUI} from './utils/check-input-column';
import {MsaWarning} from './utils/multiple-sequence-alignment';
import {multipleSequenceAlignmentUI} from './utils/multiple-sequence-alignment-ui';
import {WebLogoApp} from './apps/web-logo-app';
import {SplitToMonomersFunctionEditor} from './function-edtiors/split-to-monomers-editor';
import {splitToMonomersUI} from './utils/split-to-monomers';
import {MonomerCellRenderer} from './utils/monomer-cell-renderer';
import {BioPackage, BioPackageProperties} from './package-types';
import {getCompositionAnalysisWidget} from './widgets/composition-analysis-widget';
import {MacromoleculeColumnWidget} from './utils/macromolecule-column-widget';
import {addCopyMenuUI} from './utils/context-menu';
import {getRegionDo} from './utils/get-region';
import {GetRegionApp} from './apps/get-region-app';
import {GetRegionFuncEditor} from './utils/get-region-func-editor';
import {sequenceToMolfile} from './utils/sequence-to-mol';
import {detectMacromoleculeProbeDo} from './utils/detect-macromolecule-probe';
import {getMolColumnFromHelm} from './utils/helm-to-molfile/utils';
import {matchMoleculesWithMonomers, MonomerManager, standardizeMonomerLibrary} from './utils/monomer-lib/monomer-manager/monomer-manager';
import {calculateScoresWithEmptyValues} from './utils/calculate-scores';
import {SeqHelper} from './utils/seq-helper/seq-helper';
import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {molecular3DStructureWidget, toAtomicLevelWidget, toAtomicLevelSingle} from './widgets/to-atomic-level-widget';
import {handleSequenceHeaderRendering} from './widgets/sequence-scrolling-widget';
import {PolymerType} from '@datagrok-libraries/js-draw-lite/src/types/org';
import {BilnNotationProvider} from './utils/biln';

import * as api from './package-api';
export const _package = new BioPackage(/*{debug: true}/**/);
export * from './package.g';


// /** Avoid reassigning {@link monomerLib} because consumers subscribe to {@link IMonomerLib.onChanged} event */
// let monomerLib: MonomerLib | null = null;
let initBioPromise: Promise<void> | null = null;
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

export class PackageFunctions {
  @grok.decorators.func({description: 'Returns an instance of the monomer library helper', outputs: [{type: 'object', name: 'result'}]})
  static async getMonomerLibHelper(): Promise<IMonomerLibHelper> {
    return await MonomerLibManager.getInstance();
  }

  @grok.decorators.init({})
  static async initBio(): Promise<void> {
    if (initBioPromise == null)
      initBioPromise = initBioInt();

    await initBioPromise;
  }

  @grok.decorators.func({
    meta: {role: 'tooltip'},
  })
  static sequenceTooltip(
    @grok.decorators.param({options: {semType: 'Macromolecule'}}) col: DG.Column): DG.Widget<any> {
    const resWidget = new MacromoleculeColumnWidget(col, _package.seqHelper);
    const _resPromise = resWidget.init().then(() => { })
      .catch((err: any) => {
        const errMsg = err instanceof Error ? err.message : err.toString();
        grok.shell.error(errMsg);
      });
    return resWidget;
  }

  @grok.decorators.func({})
  static async standardiseMonomerLibrary(library: string): Promise<string> {
    return await standardizeMonomerLibrary(library);
  }

  @grok.decorators.func({'top-menu': 'Bio | Manage | Match with Monomer Library...', description: 'Matches molecules in a column with monomers from the selected library(s)',})
  static async matchWithMonomerLibrary(table: DG.DataFrame,
      @grok.decorators.param({type: 'column', options: {semType: 'Molecule'}})molecules: DG.Column,
      @grok.decorators.param({type: 'string', options: {choices: ['PEPTIDE', 'RNA', 'CHEM'], initialValue: 'PEPTIDE', caption: 'Polymer Type'}})polymerType: PolymerType = 'PEPTIDE') {
    const matchDF = await matchMoleculesWithMonomers(table, molecules.name, _package.monomerLib, polymerType);
    grok.shell.addTableView(matchDF);
  }

  // Keep for backward compatibility
  @grok.decorators.func({outputs: [{type: 'object', name: 'monomerLib'}]})
  static getBioLib(): IMonomerLib {
    return _package.monomerLib;
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static getSeqHandler(
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}}) sequence: DG.Column<string>): ISeqHandler {
    return _package.seqHelper.getSeqHandler(sequence);
  }

  // -- Panels --

  @grok.decorators.panel({name: 'Bioinformatics | Get Region', description: 'Creates a new column with sequences of the region between start and end'})
  static getRegionPanel(
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}}) seqCol: DG.Column<string>): DG.Widget {
    const funcName: string = 'getRegionTopMenu';
    const funcList = DG.Func.find({package: _package.name, name: funcName});
    if (funcList.length !== 1) throw new Error(`Package '${_package.name}' func '${funcName}' not found`);
    const func = funcList[0];
    const funcCall = func.prepare({table: seqCol.dataFrame, sequence: seqCol});
    const funcEditor = new GetRegionFuncEditor(funcCall, _package.seqHelper);
    return funcEditor.widget();
  }

  @grok.decorators.panel({
    name: 'Bioinformatics | Manage Monomer Libraries',
    meta: {'exclude-actions-panel': 'true'}
  })
  static async libraryPanel(
    @grok.decorators.param({name: 'seqColumn', options: {semType: 'Macromolecule'}}) _seqColumn: DG.Column): Promise<DG.Widget> {
    // return getLibraryPanelUI();
    return getMonomerLibraryManagerLink();
  }

  // -- Func Editors --

  @grok.decorators.editor({})
  static GetRegionEditor(call: DG.FuncCall): void {
    try {
      const funcEditor = new GetRegionFuncEditor(call, _package.seqHelper);
      funcEditor.dialog();
    } catch (err: any) {
      const errMsg = err instanceof Error ? err.message : err.toString();
      const errStack = err instanceof Error ? err.stack : undefined;
      grok.shell.error(`Get region editor error: ${errMsg}`);
      _package.logger.error(errMsg, undefined, errStack);
    }
  }

  @grok.decorators.editor({})
  static SplitToMonomersEditor(call: DG.FuncCall): void {
    const funcEditor = new SplitToMonomersFunctionEditor();
    ui.dialog({title: 'Split to Monomers'})
      .add(funcEditor.paramsUI)
      .onOK(async () => {
        return call.func.prepare(funcEditor.funcParams).call(true);
      })
      .show();
  }

  @grok.decorators.editor({})
  static SequenceSpaceEditor(call: DG.FuncCall) {
    const funcEditor = new DimReductionBaseEditor({semtype: DG.SEMTYPE.MACROMOLECULE});
    const dialog = ui.dialog({title: 'Sequence Space'})
      .add(funcEditor.getEditor())
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
        }).call();
      });
    dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
    dialog.show();
  }

  @grok.decorators.editor({})
  static SeqActivityCliffsEditor(call: DG.FuncCall) {
    const funcEditor = new ActivityCliffsEditor({semtype: DG.SEMTYPE.MACROMOLECULE});
    const dialog = ui.dialog({title: 'Activity Cliffs'})
      .add(funcEditor.getEditor())
      .onOK(async () => {
        const params = funcEditor.getParams();
        return call.func.prepare({
          table: params.table,
          molecules: params.col,
          activities: params.activities,
          similarity: params.similarityThreshold,
          methodName: params.methodName,
          similarityMetric: params.similarityMetric,
          preprocessingFunction: params.preprocessingFunction,
          options: params.options,
        }).call();
      });
    dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
    dialog.show();
  }


  // -- Cell renderers --

  @grok.decorators.func({
    name: 'customSequenceCellRenderer',
    meta: {
      cellType: 'sequence',
      columnTags: 'quality=Macromolecule, units=custom',
      role: 'cellRenderer'
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static customSequenceCellRenderer(): DG.GridCellRenderer {
    return new MacromoleculeSequenceCellRenderer();
  }

  @grok.decorators.func({
    name: 'fastaSequenceCellRenderer',
    meta: {
      cellType: 'sequence',
      columnTags: 'quality=Macromolecule, units=fasta',
      role: 'cellRenderer'
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static fastaSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
    return new MacromoleculeSequenceCellRenderer();
  }

  @grok.decorators.func({
    name: 'separatorSequenceCellRenderer',
    meta: {
      cellType: 'sequence',
      columnTags: 'quality=Macromolecule, units=separator',
      role: 'cellRenderer'
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static separatorSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
    return new MacromoleculeSequenceCellRenderer();
  }

  @grok.decorators.func({
    name: 'bilnSequenceCellRenderer',
    meta: {
      cellType: 'sequence',
      columnTags: 'quality=Macromolecule, units=biln',
      role: 'cellRenderer'
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static bilnSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
    return new MacromoleculeSequenceCellRenderer();
  }

  @grok.decorators.func({
    name: 'refineNotationProviderForBiln',
    outputs: [{type: 'bool', name: 'result'}],
    meta: {role: 'notationRefiner'},
  })
  static refineNotationProviderForBiln(
    @grok.decorators.param({type: 'column'}) col: DG.Column<string>,
    @grok.decorators.param({type: 'object'}) stats: {freq: { [key: string]: number; }, sameLength: boolean},
    @grok.decorators.param({type: 'string', options: {nullable: true, optional: true}}) separator: string | null
  ): boolean {
    if (separator !== '-')
      return false;// biln uses '-' as a separator
    const reCons = Object.keys(stats.freq).some((om) => om.match(/^.+\(\d{1,2},\d{1,2}\)$/));
    if (!reCons) {
      // biln might also encode monomers with hyphens in names encoded by []
      // here we know that there are no monomers with connections like (1,2) in names, so we can check for []
      const reBrackets = Object.keys(stats.freq).some((om) => om.includes('[') || om.includes(']'));
      if (!reBrackets)
        return false;
    }
    // refine the notation provider
    col.setTag('aligned', 'SEQ');
    col.setTag('alphabet', 'UN');
    col.setTag('.alphabetIsMultichar', 'true');
    col.meta.units = NOTATION.BILN;
    col.temp[SeqTemps.notationProvider] = new BilnNotationProvider(separator, _package.seqHelper, col);

    return true;
  }

  // // -- Property panels --

  @grok.decorators.panel({name: 'Bioinformatics | Sequence Renderer'})
  static macroMolColumnPropertyPanel(
    @grok.decorators.param({options: {semType: 'Macromolecule'}}) molColumn: DG.Column): DG.Widget {
    return getMacromoleculeColumnPropertyPanel(molColumn);
  }

  @grok.decorators.panel({name: 'Composition analysis',
    meta: {role: 'widgets', domain: 'bio'},
  })
  static compositionAnalysisWidget(
    @grok.decorators.param({options: {semType: 'Macromolecule'}}) sequence: DG.SemanticValue): DG.Widget {
    return getCompositionAnalysisWidget(sequence, _package.monomerLib, _package.seqHelper);
  }

  @grok.decorators.func({
    name: 'MacromoleculeDifferenceCellRenderer',
    meta: {
      cellType: 'MacromoleculeDifference',
      columnTags: 'quality=MacromoleculeDifference',
      role: 'cellRenderer'
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static macromoleculeDifferenceCellRenderer(): MacromoleculeDifferenceCellRenderer {
    return new MacromoleculeDifferenceCellRenderer();
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static sequenceAlignment(
    @grok.decorators.param({options: {choices: ['Local alignment', 'Global alignment']}}) alignType: string,
    @grok.decorators.param({options: {choices: ['AUTO', 'NUCLEOTIDES', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250', 'SCHNEIDER', 'TRANS']}}) alignTable: string,
      gap: number,
      seq1: string,
      seq2: string) {
    const toAlign = new SequenceAlignment(seq1, seq2, gap, alignTable);
    const res = alignType == 'Local alignment' ? toAlign.smithWaterman() : toAlign.needlemanWunsch();
    return res;
  }

  // -- Viewers --

  @grok.decorators.panel({
    name: 'WebLogo',
    description: 'WebLogo',
    meta: {icon: 'files/icons/weblogo-viewer.svg', role: 'viewer'},
    outputs: [{type: 'viewer', name: 'result'}]
  })
  static webLogoViewer() {
    return new WebLogoViewer();
  }

  @grok.decorators.panel({
    name: 'VdRegions',
    description: 'V-Domain regions viewer',
    meta: {icon: 'files/icons/vdregions-viewer.svg', role: 'viewer'},
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static vdRegionsViewer() {
    return new VdRegionsViewer();
  }

  // -- Top menu --

  @grok.decorators.func({name: 'getRegion', description: 'Gets a new column with sequences of the region between start and end'})
  static getRegion(
    @grok.decorators.param({type: 'column'})sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {optional: true}}) start: string | undefined,
    @grok.decorators.param({type: 'string', options: {optional: true}}) end: string | undefined,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'Name of the column to be created'}}) name: string | undefined): DG.Column<string> {
    return getRegionDo(sequence,
      start ?? null, end ?? null, name ?? null);
  }

  @grok.decorators.func({
    name: 'Get Region Top Menu',
    description: 'Get sequences for a region specified from a Macromolecule',
    'top-menu': 'Bio | Calculate | Get Region...',
    editor: 'Bio:GetRegionEditor'})
  static async getRegionTopMenu(
    @grok.decorators.param({options: {description: 'Input data table'}})table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Macromolecule', description: 'Sequence column'}}) sequence: DG.Column,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'Region start position name'}}) start: string | undefined,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'Region end position name'}}) end: string | undefined,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'Region column name'}}) name: string | undefined
  ): Promise<void> {
    const regCol = getRegionDo(sequence, start ?? null, end ?? null, name ?? null);
    sequence.dataFrame.columns.add(regCol);
    await grok.data.detectSemanticTypes(sequence.dataFrame); // to set renderer
  }

  @grok.decorators.func({
    name: 'Sequence Activity Cliffs',
    description: 'Detects pairs of molecules with similar structure and significant difference in any given property',
    'top-menu': 'Bio | Analyze | Activity Cliffs...',
    editor: 'Bio:SeqActivityCliffsEditor',
    outputs: []
  })
  static async activityCliffs(
    @grok.decorators.param({options: {description: 'Input data table'}})table: DG.DataFrame,
    @grok.decorators.param({type: 'string', options: {semType: 'Macromolecule', description: 'Input data table'}}) molecules: DG.Column<string>,
      activities: DG.Column,
    @grok.decorators.param({options: {initialValue: '80', description: 'Similarity cutoff'}}) similarity: number,
    @grok.decorators.param({type: 'string', options: {choices: ['UMAP', 't-SNE']}}) methodName: DimReductionMethods,
    @grok.decorators.param({type: 'string', options: {choices: ['Hamming', 'Levenshtein', 'Monomer chemical distance']}}) similarityMetric: MmDistanceFunctionsNames | BitArrayMetrics,
    @grok.decorators.param({type: 'func'}) preprocessingFunction: DG.Func,
    @grok.decorators.param({type: 'object', options: {optional: true}}) options?: (IUMAPOptions | ITSNEOptions) & Options,
    @grok.decorators.param({options: {optional: true}}) demo?: boolean): Promise<DG.Viewer | undefined> {
    //workaround for functions which add viewers to tableView (can be run only on active table view)
    if (table.name !== grok.shell.tv.dataFrame.name) {
      grok.shell.error(`Table ${table.name} is not a current table view`);
      return;
    }
    if (!checkInputColumnUI(molecules, 'Activity Cliffs'))
      return;
    const axesNames = getEmbeddingColsNames(table);
    const tags = {
      'units': molecules.meta.units!,
      'aligned': molecules.getTag(bioTAGS.aligned),
      'separator': molecules.getTag(bioTAGS.separator),
      'alphabet': molecules.getTag(bioTAGS.alphabet),
    };
    const columnDistanceMetric: MmDistanceFunctionsNames | BitArrayMetrics = similarityMetric;
    const seqCol = molecules;

    const runCliffs = async () => {
      const sp = await getActivityCliffs(
        table,
        seqCol,
        axesNames,
        'Activity cliffs', //scatterTitle
        activities,
        similarity,
        columnDistanceMetric, //similarityMetric
        methodName,
        {...(options ?? {})},
        DG.SEMTYPE.MACROMOLECULE,
        tags,
        preprocessingFunction,
        createTooltipElement,
        createPropPanelElement,
        createLinesGrid,
        undefined,
        demo
      );
      return sp;
    };

    const allowedRowCount = methodName === DimReductionMethods.UMAP ? 200_000 : 20_000;
    const fastRowCount = methodName === DimReductionMethods.UMAP ? 5_000 : 2_000;
    if (table.rowCount > allowedRowCount) {
      grok.shell.warning(`Too many rows, maximum for sequence activity cliffs is ${allowedRowCount}`);
      return;
    }

    const pi = DG.TaskBarProgressIndicator.create(`Running sequence activity cliffs ...`);
    const scRes = (await new Promise<DG.Viewer | undefined>((resolve, reject) => {
      if (table.rowCount > fastRowCount && !options?.[BYPASS_LARGE_DATA_WARNING]) {
        ui.dialog().add(ui.divText(`Activity cliffs analysis might take several minutes.
    Do you want to continue?`))
          .onOK(async () => {
            runCliffs().then((res) => resolve(res)).catch((err) => reject(err));
          })
          .onCancel(() => { resolve(undefined); })
          .show();
      } else
        runCliffs().then((res) => resolve(res)).catch((err) => reject(err));
    }).catch((err: any) => {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
      throw err;
    }).finally(() => { pi.close(); })) as DG.ScatterPlotViewer | undefined;
    if (scRes?.props?.xColumnName && scRes?.props?.yColumnName && table.col(scRes.props.xColumnName) && table.col(scRes.props.yColumnName)) {
    table.col(scRes.props.xColumnName)!.set(0, table.col(scRes.props.xColumnName)!.get(0)); // to trigger rendering
    table.col(scRes.props.yColumnName)!.set(0, table.col(scRes.props.yColumnName)!.get(0)); // to trigger rendering
    }

    return scRes;
  }

  @grok.decorators.func({
    name: 'Encode Sequences',
    meta: {
      supportedSemTypes: 'Macromolecule',
      supportedTypes: 'string',
      supportedDistanceFunctions: 'Hamming,Levenshtein,Monomer chemical distance,Needlemann-Wunsch',
      role: 'dimRedPreprocessingFunction'
    },
    outputs: [{type: 'object', name: 'result'}],
  })
  static async macromoleculePreprocessingFunction(
    @grok.decorators.param({options: {semType: 'Macromolecule'}})col: DG.Column,
    @grok.decorators.param({type: 'string'}) metric: MmDistanceFunctionsNames,
    @grok.decorators.param({options: {initialValue: '1', caption: 'Gap open penalty', optional: true}}) gapOpen: number = 1,
    @grok.decorators.param({options: {initialValue: '0.6', caption: 'Gap extension penalty', optional: true}}) gapExtend: number = 0.6,
    @grok.decorators.param({options: {caption: 'Fingerprint type', initialValue: 'Morgan', choices: ['Morgan', 'RDKit', 'Pattern', 'AtomPair', 'MACCS', 'TopologicalTorsion'], optional: true}}) fingerprintType : string = 'Morgan'): Promise<PreprocessFunctionReturnType> {
    if (col.semType !== DG.SEMTYPE.MACROMOLECULE)
      return {entries: col.toList(), options: {}};

    const {seqList, options} = await getEncodedSeqSpaceCol(col, metric, fingerprintType, gapOpen, gapExtend);
    return {entries: seqList, options};
  }

  @grok.decorators.func({name: 'Helm Fingerprints',
    meta: {
      supportedSemTypes: 'Macromolecule',
      supportedTypes: 'string',
      supportedUnits: 'helm',
      supportedDistanceFunctions: 'Tanimoto,Asymmetric,Cosine,Sokal'
    },
    outputs: [{type: 'object', name: 'result'}],
  })
  static async helmPreprocessingFunction(
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}}) col: DG.Column<string>,
    @grok.decorators.param({type: 'string'})_metric: BitArrayMetrics): Promise<PreprocessFunctionReturnType> {
    if (col.version !== col.temp[MONOMERIC_COL_TAGS.LAST_INVALIDATED_VERSION])
      await invalidateMols(col, _package.seqHelper, false);
    const molCol = col.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS];
    const fingerPrints: DG.Column<DG.BitSet | null> =
      await grok.functions.call('Chem:getMorganFingerprints', {molColumn: molCol});

    const entries: Array<BitArray | null> = new Array(fingerPrints.length).fill(null);
    for (let i = 0; i < fingerPrints.length; i++) {
      if (fingerPrints.isNone(i) || !fingerPrints.get(i))
        continue;
      const fp = fingerPrints.get(i)!;
      entries[i] = BitArray.fromUint32Array(fp.length, new Uint32Array(fp.getBuffer().buffer));
    }
    return {entries, options: {}};
  }

  @grok.decorators.func({
    name: 'Sequence Space',
    description: 'Creates 2D sequence space with projected sequences by pairwise distance',
    'top-menu': 'Bio | Analyze | Sequence Space...',
    editor: 'Bio:SequenceSpaceEditor',
    outputs: [],
  })
  static async sequenceSpaceTopMenu(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Macromolecule'}}) molecules: DG.Column,
    @grok.decorators.param({type: 'string', options: {choices: ['UMAP', 't-SNE']}}) methodName: DimReductionMethods,
    @grok.decorators.param({type: 'string', options: {choices: ['Hamming', 'Levenshtein', 'Monomer chemical distance']}}) similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames,
    @grok.decorators.param({options: {initialValue: 'true'}}) plotEmbeddings: boolean,
    @grok.decorators.param({type: 'func', options: {optional: true}}) preprocessingFunction?: DG.Func,
    @grok.decorators.param({type: 'object', options: {optional: true}}) options?: (IUMAPOptions | ITSNEOptions) & Options,
    @grok.decorators.param({options: {optional: true, initialValue: 'true'}}) clusterEmbeddings?: boolean,
    @grok.decorators.param({options: {optional: true}}) isDemo?: boolean
  ): Promise<DG.ScatterPlotViewer | undefined> {
    //workaround for functions which add viewers to tableView (can be run only on active table view)
    if (table.name !== grok.shell.tv.dataFrame.name) {
      grok.shell.error(`Table ${table.name} is not a current table view`);
      return;
    }
    const tableView =
      grok.shell.tv.dataFrame == table ? grok.shell.tv : undefined;
    if (!checkInputColumnUI(molecules, 'Sequence Space'))
      return;
    if (!preprocessingFunction)
      preprocessingFunction = DG.Func.find({name: 'macromoleculePreprocessingFunction', package: 'Bio'})[0];
    options ??= {};
    const res = await multiColReduceDimensionality(table, [molecules], methodName,
      [similarityMetric as KnownMetrics], [1], [preprocessingFunction], 'MANHATTAN',
      plotEmbeddings, clusterEmbeddings ?? false,
      /* dimRedOptions */ {...options, preprocessingFuncArgs: [options.preprocessingFuncArgs ?? {}]},
      /* uiOptions */{
        fastRowCount: 10000,
        scatterPlotName: 'Sequence space',
        bypassLargeDataWarning: options?.[BYPASS_LARGE_DATA_WARNING],
        tableView: tableView,
      });
    return res;
  }

  @grok.decorators.func({
    name: 'Molecules to HELM',
    'top-menu': 'Bio | Transform | Molecules to HELM...',
    description: 'Converts Peptide molecules to HELM notation by matching with monomer library',
  })
  static async moleculesToHelmTopMenu(
    @grok.decorators.param({name: 'table', options: {description: 'Input data table'}})table: DG.DataFrame,
    @grok.decorators.param({name: 'molecules', options: {semType: 'Molecule', description: 'Molecule column'}})molecules: DG.Column,
  ) {
    // collect current monomer library
    const monomerLib = _package.monomerLib;
    const libJSON = JSON.stringify(monomerLib.toJSON());
    await api.scripts.molToHelmConverterPy(table, molecules, libJSON);

    // semtype is not automatically set, so we set it manually
    const newCol = table.columns.toList().find((c) => c.name.toLowerCase().includes('regenerated sequence') && c.semType !== DG.SEMTYPE.MACROMOLECULE);
    if (newCol) {
      newCol.meta.units = NOTATION.HELM;
      newCol.semType = DG.SEMTYPE.MACROMOLECULE;
      newCol.setTag('cell.renderer', 'helm');
    }
  }

  @grok.decorators.func({name: 'Molecule to HELM Single', description: 'Converts a single molecule to HELM notation without requiring a table or column', outputs: [{name: 'result', type: 'string', options: {semType: 'Macromolecule', units: 'helm'}}]})
  static async moleculeToHelmSingle(
    @grok.decorators.param({name: 'molecule', options: {semType: 'Molecule', description: 'Input molecule'}})molecule: string,
  ): Promise<string> {
    // create temporary dataframe
    const tempCol = DG.Column.fromStrings('molecule', [molecule]);
    tempCol.semType = DG.SEMTYPE.MOLECULE;
    const tempDF = DG.DataFrame.fromColumns([tempCol]);
    // call converter
    await PackageFunctions.moleculesToHelmTopMenu(tempDF, tempCol);
    // get result
    const result = tempDF.columns.toList().find((c) => c.name.toLowerCase().includes('regenerated sequence'))?.get(0);
    return result ?? '';
  }

  @grok.decorators.func({
    name: 'To Atomic Level',
    description: 'Converts sequences to molblocks',
    'top-menu': 'Bio | Transform | To Atomic Level...',
  })
  static async toAtomicLevel(
    @grok.decorators.param({options: {description: 'Input data table'}})table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Macromolecule', caption: 'Sequence'}})seqCol: DG.Column,
    @grok.decorators.param({options: {initialValue: 'true', caption: 'Non-linear', description: 'Slower mode for cycling/branching HELM structures'}}) nonlinear: boolean = true,
    @grok.decorators.param({options: {initialValue: 'false', caption: 'Highlight monomers', description: 'Highlight monomers\' substructures of the molecule'}}) highlight: boolean = false
  ): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Converting to atomic level ...');
    try {
      await initBioPromise;
      const monomerLib = seqCol.temp[MmcrTemps.overriddenLibrary] ?? _package.monomerLib;
      const seqHelper = _package.seqHelper;
      const rdKitModule = _package.rdKitModule;
      await sequenceToMolfile(table, seqCol, nonlinear, highlight, monomerLib, seqHelper, rdKitModule);
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({
    name: 'To Atomic Level...',
    meta: {action: 'to atomic level'}
  })
  static async toAtomicLevelAction(
    @grok.decorators.param({options: {semType: 'Macromolecule'}}) seqCol: DG.Column) {
    if (!seqCol?.dataFrame)
      throw new Error('Sequence column is not found or its data frame is not empty');
    const func = DG.Func.find({name: 'toAtomicLevel', package: 'Bio'})[0];
    if (!func) throw new Error('To Atomic level Function not found');
    func.prepare({table: seqCol.dataFrame, seqCol: seqCol}).edit();
  }

  @grok.decorators.panel({
    name: 'Molecular Structure',
    meta: {role: 'widgets', domain: 'bio'},
  })
  static async toAtomicLevelPanel(
    @grok.decorators.param({name: 'sequence', type: 'semantic_value', options: {semType: 'Macromolecule'}})
      sequence: DG.SemanticValue
  ) : Promise<DG.Widget> {
    return toAtomicLevelWidget(sequence);
  }

  @grok.decorators.func({
    name: 'To Atomic Level Single sequence',
    description: 'Converts a single sequence to molblock',
    outputs: [{name: 'molfile', type: 'string', options: {semType: 'Molecule'}}]
  })
  static async toAtomicLevelSingleSeq(
    @grok.decorators.param({name: 'sequence', type: 'string', options: {semType: 'Macromolecule'}}) sequence: string,
  ) : Promise<string> {
    // create temporary column and table
    const isHelm = sequence.includes('$$');
    const isSeparator = sequence.split('').filter((c) => c == '/').length > 2;
    const isBiln = sequence.split('').filter((c) => c == '-').length > 2; // biln is super separator basically :D
    // const isFasta = !isHelm && !isSeparator && !isBiln;
    const tempCol = DG.Column.fromStrings('sequence', [sequence]);
    tempCol.semType = DG.SEMTYPE.MACROMOLECULE;
    tempCol.meta.units = isHelm ? NOTATION.HELM : isSeparator ? NOTATION.SEPARATOR : isBiln ? NOTATION.BILN : NOTATION.FASTA;
    tempCol.setTag(bioTAGS.aligned, 'SEQ');
    if (isSeparator)
      tempCol.setTag(bioTAGS.separator, '/');
    if (isBiln)
      tempCol.setTag(bioTAGS.separator, '-');
    // detect alphabet
    const dnaAlphabet = 'AGCT';
    const RNAAlphabet = 'AGCU';
    const isDNA = sequence.split('').every((c) => dnaAlphabet.includes(c) || c === '/' || c === '-');
    const isRNA = !isDNA && sequence.split('').every((c) => RNAAlphabet.includes(c) || c === '/' || c === '-');

    tempCol.setTag(bioTAGS.alphabet, isDNA ? ALPHABET.DNA : isRNA ? ALPHABET.RNA : 'UN');
    tempCol.setTag(bioTAGS.alphabetIsMultichar, 'true');

    const tempDF = DG.DataFrame.fromColumns([tempCol]);

    // get the sell as semantic value
    const cell = tempDF.cell(0, 'sequence');
    const semValue = DG.SemanticValue.fromTableCell(cell);

    const res = await toAtomicLevelSingle(semValue);
    if (res.errorText || !res.mol)
      throw new Error(res.errorText);
    return res.mol;
  }

  @grok.decorators.panel({
    name: 'Molecular 3D Structure',
    meta: {role: 'widgets', domain: 'bio'},
  })
  static async sequence3dStructureWidget(
    @grok.decorators.param({
      type: 'semantic_value',
      options: {semType: 'Macromolecule'}
    })
      sequence: DG.SemanticValue
  ): Promise<DG.Widget> {
    return molecular3DStructureWidget(sequence);
  }

  @grok.decorators.panel({
    name: 'MSA',
    description: 'Performs multiple sequence alignment',
    meta: {domain: 'bio'},
    'top-menu': 'Bio | Analyze | MSA...'
  })
  static multipleSequenceAlignmentDialog(): void {
    multipleSequenceAlignmentUI({}, _package.seqHelper)
      .catch((err: any) => {
        const [errMsg, errStack] = errInfo(err);
        if (err instanceof MsaWarning) {
          grok.shell.warning((err as MsaWarning).element);
          _package.logger.warning(errMsg);
          return;
        }
        grok.shell.error(errMsg);
        _package.logger.error(errMsg, undefined, errStack);
        // throw err; // This error throw is not handled
      });
  }

  @grok.decorators.func({
    name: 'Multiple Sequence Alignment',
    description: 'Multiple sequence alignment',
    meta: {domain: 'bio'}
  })
  static async alignSequences(
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}}) sequenceCol: DG.Column<string> | null = null,
    @grok.decorators.param({type: 'column'}) clustersCol: DG.Column | null = null,
    @grok.decorators.param({type: 'object', options: {optional: true}}) options?: any
  ): Promise<DG.Column<string>> {
    return multipleSequenceAlignmentUI({col: sequenceCol, clustersCol: clustersCol, ...options}, _package.seqHelper);
  }

  @grok.decorators.func({
    name: 'Composition Analysis',
    description: 'Visualizes sequence composition on a WebLogo plot',
    'top-menu': 'Bio | Analyze | Composition',
    meta: {
      icon: 'files/icons/composition-analysis.svg'
    },
    outputs: [{name: 'result', type: 'viewer'}]
  })
  static async compositionAnalysis(): Promise<void> {
    // Higher priority for columns with MSA data to show with WebLogo.
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    //@ts-ignore
    const colList: DG.Column[] = df.columns.toList().filter((col) => {
      if (col.semType != DG.SEMTYPE.MACROMOLECULE)
        return false;

      const _colSh = _package.seqHelper.getSeqHandler(col);
      // TODO: prevent for cyclic, branched or multiple chains in Helm
      return true;
    });

    const handler = async (col: DG.Column) => {
      if (!checkInputColumnUI(col, 'Composition'))
        return;

      const wlViewer = tv.addViewer('WebLogo', {sequenceColumnName: col.name});
      grok.shell.tv.dockManager.dock(wlViewer, DG.DOCK_TYPE.DOWN, null, 'Composition analysis', 0.25);
    };

    let col: DG.Column | null = null;
    if (colList.length == 0) {
      grok.shell.error('Current table does not contain sequences');
      return;
    } else if (colList.length > 1) {
      const colListNames: string [] = colList.map((col) => col.name);
      const selectedCol = colList.find((c) => { return _package.seqHelper.getSeqHandler(c).isMsa(); });
      const colInput: DG.InputBase = ui.input.choice(
        'Column', {value: selectedCol ? selectedCol.name : colListNames[0], items: colListNames});
      ui.dialog({
        title: 'Composition Analysis',
        helpUrl: 'https://datagrok.ai/help/datagrok/solutions/domains/bio/#sequence-composition',
      })
        .add(ui.div([
          colInput,
        ]))
        .onOK(async () => {
          const col: DG.Column | null = colList.find((col) => col.name == colInput.value) ?? null;

          if (col)
            await handler(col);
        })
        .show();
    } else
      col = colList[0];

    if (!col)
      return;

    await handler(col);
  }

  // -- Package settings editor --


  @grok.decorators.fileHandler({
    name: 'importFasta',
    description: 'Opens FASTA file',
    ext: 'fasta, fna, ffn, faa, frn, fa, fst',
  })
  static importFasta(
    fileContent: string
  ): DG.DataFrame[] {
    const ffh = new FastaFileHandler(fileContent);
    return ffh.importFasta();
  }

  @grok.decorators.fileHandler({
    name: 'importBam',
    description: 'Opens Bam file',
    ext: 'bam, bai',
  })
  static importBam(fileContent: string): DG.DataFrame[] {
    console.log(fileContent);
    return [];
  }

  @grok.decorators.func({
    name: 'convertDialog',
    'top-menu': 'Bio | Transform | Convert Sequence Notation...'
  })
  static convertDialog() {
    const col: DG.Column<string> | undefined = getMacromoleculeColumns()[0];
    convert(col, _package.seqHelper);
  }

  @grok.decorators.func({
    name: 'Convert Notation...',
    meta: {
      action: 'Convert Notation...'
    }
  })
  static convertColumnAction(
    @grok.decorators.param({
      options: {semType: 'Macromolecule'}
    })
      col: DG.Column
  ) {
    convert(col, _package.seqHelper);
  }

  @grok.decorators.func({
    name: 'monomerCellRenderer',
    meta: {
      cellType: 'Monomer',
      columnTags: 'quality=Monomer',
      role: 'cellRenderer'
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static monomerCellRenderer(): MonomerCellRenderer {
    return new MonomerCellRenderer();
  }

  @grok.decorators.func({
    name: 'testDetectMacromolecule',
  })
  static async testDetectMacromolecule(
    @grok.decorators.param({
      options: {choices: ['Demo:Files/', 'System:AppData/']}
    })
      path: string
  ): Promise<DG.DataFrame> {
    const pi = DG.TaskBarProgressIndicator.create('Test detectMacromolecule...');

    const fileList = await grok.dapi.files.list(path, true, '');
    //@ts-ignore
    const fileListToTest = fileList.filter((fi) => fi.fileName.endsWith('.csv'));

    let readyCount = 0;
    const res = [];

    for (const fileInfo of fileListToTest) {
      try {
        const csv = await grok.dapi.files.readAsText(path + fileInfo.fullPath);
        const df = DG.DataFrame.fromCsv(csv);

        for (const col of df.columns) {
          const semType = await grok.functions.call('Bio:detectMacromolecule', {col: col});
          if (semType === DG.SEMTYPE.MACROMOLECULE) {
            //console.warn(`file: ${fileInfo.path}, column: ${col.name}, ` +
            //  `semType: ${semType}, units: ${col.meta.units}`);
            // console.warn('file: '' + fileInfo.path + '', semType: '' + semType + '', ' +
            //   'units: '' + col.meta.units + ''');

            res.push({
              file: fileInfo.path, result: 'detected', column: col.name,
              message: `units: ${col.meta.units}`,
            });
          }
        }
      } catch (err: unknown) {
        // console.error('file: ' + fileInfo.path + ', error: ' + ex.toString());
        res.push({
          file: fileInfo.path, result: 'error', column: null,
          message: err instanceof Error ? err.message : (err as Object).toString(),
        });
      } finally {
        readyCount += 1;
        pi.update(100 * readyCount / fileListToTest.length, `Test ${fileInfo.fileName}`);
      }
    }

    grok.shell.info('Test Demo:Files for detectMacromolecule finished.');
    pi.close();
    const resDf = DG.DataFrame.fromObjects(res)!;
    resDf.name = `datasets_detectMacromolecule_${path}`;
    return resDf;
  }

  @grok.decorators.func({
    name: 'Split to Monomers',
    'top-menu': 'Bio | Transform | Split to Monomers...',
    editor: 'Bio:SplitToMonomersEditor',
  })
  static async splitToMonomersTopMenu(
    table: DG.DataFrame,
    @grok.decorators.param({
      options: {semType: 'Macromolecule'}
    })
    sequence: DG.Column
  ): Promise<DG.DataFrame> {
    return await splitToMonomersUI(table, sequence);
  }

  @grok.decorators.func({
    name: 'Bio: getHelmMonomers',
    outputs: [{name: 'result', type: 'object'}]
  })
  static getHelmMonomers(
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}})
      sequence: DG.Column<string>
  ): string[] {
    return _package.seqHelper.getSeqMonomers(sequence);
  }

  @grok.decorators.func({
    name: 'Sequence Similarity Search',
    meta: {icon: 'files/icons/sequence-similarity-viewer.svg', role: 'viewer'},
    outputs: [{name: 'result', type: 'viewer'}]
  })
  static similaritySearchViewer(): SequenceSimilarityViewer {
    return new SequenceSimilarityViewer(_package.seqHelper);
  }

  @grok.decorators.func({
    name: 'similaritySearch',
    description: 'Finds similar sequences',
    'top-menu': 'Bio | Search | Similarity Search',
    outputs: [{name: 'result', type: 'viewer'}]
  })
  static similaritySearchTopMenu(): void {
    const view = (grok.shell.v as DG.TableView);
    const viewer = view.addViewer('Sequence Similarity Search');
    view.dockManager.dock(viewer, 'down');
  }

  @grok.decorators.func({
    name: 'Sequence Diversity Search',
    meta: {icon: 'files/icons/sequence-diversity-viewer.svg', role: 'viewer'},
    outputs: [{name: 'result', type: 'viewer'}]
  })
  static diversitySearchViewer(): SequenceDiversityViewer {
    return new SequenceDiversityViewer(_package.seqHelper);
  }

  @grok.decorators.func({
    name: 'diversitySearch',
    description: 'Finds the most diverse sequences',
    'top-menu': 'Bio | Search | Diversity Search',
    outputs: [{name: 'result', type: 'viewer'}]
  })
  static diversitySearchTopMenu() {
    const view = (grok.shell.v as DG.TableView);
    const viewer = view.addViewer('Sequence Diversity Search');
    view.dockManager.dock(viewer, 'down');
  }

  @grok.decorators.editor({
    name: 'SearchSubsequenceEditor'
  })
  static searchSubsequenceEditor(call: DG.FuncCall) {
    const columns = getMacromoleculeColumns();
    if (columns.length === 1)
      call.func.prepare({macromolecules: columns[0]}).call(true);
    else
      new SubstructureSearchDialog(columns, _package.seqHelper);
  }

  @grok.decorators.func({
    name: 'Subsequence Search',
    'top-menu': 'Bio | Search | Subsequence Search ...',
    editor: 'Bio:SearchSubsequenceEditor'
  })
  static SubsequenceSearchTopMenu(macromolecules: DG.Column): void {
    grok.shell.tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: 'Bio:bioSubstructureFilter',
      column: macromolecules.name,
      columnName: macromolecules.name,
    });
    grok.shell.tv.grid.scrollToCell(macromolecules, 0);
  }

  @grok.decorators.func({
    name: 'Identity',
    description: 'Adds a column with fraction of matching monomers',
    'top-menu': 'Bio | Calculate | Identity...',
  })
  static async sequenceIdentityScoring(
    @grok.decorators.param({options: {description: 'Table containing Macromolecule column'}}) table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Macromolecule', description: 'Sequences to score'}}) macromolecule: DG.Column,
    @grok.decorators.param({options: {description: 'Sequence,matching column format'}}) reference: string
  ): Promise<DG.Column<number>> {
    const seqHelper = _package.seqHelper;
    const scores = calculateScoresWithEmptyValues(table, macromolecule, reference, SCORE.IDENTITY, seqHelper);
    return scores;
  }

  @grok.decorators.func({
    name: 'Similarity',
    description: 'Adds a column with similarity scores, calculated as sum of monomer fingerprint similarities',
    'top-menu': 'Bio | Calculate | Similarity...',
  })
  static async sequenceSimilarityScoring(
    @grok.decorators.param({options: {description: 'Table containing Macromolecule column'}}) table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Macromolecule', description: 'Sequences to score'}}) macromolecule: DG.Column,
    @grok.decorators.param({options: {description: 'Sequence,matching column format'}}) reference: string
  ): Promise<DG.Column<number>> {
    const seqHelper = _package.seqHelper;
    const scores = calculateScoresWithEmptyValues(table, macromolecule, reference, SCORE.SIMILARITY, seqHelper);
    return scores;
  }

  @grok.decorators.func({
    name: 'Manage Monomer Libraries',
    description: 'Manage HELM monomer libraries'
  })
  static async manageMonomerLibraries(): Promise<void> {
    showManageLibrariesDialog();
  }

  @grok.decorators.func({
    name: 'Manage Monomer Libraries View',
    'top-menu': 'Bio | Manage | Monomer Libraries'
  })
  static async manageLibrariesView(): Promise<void> {
    await showManageLibrariesView();
  }

  @grok.decorators.func({
    name: 'manageMonomersView',
    description: 'Edit and create monomers',
    'top-menu': 'Bio | Manage | Monomers'
  })
  static async manageMonomersView() {
    const monomerManager = await MonomerManager.getInstance();
    await monomerManager.getViewRoot();
  }

  @grok.decorators.app({
    name: 'Manage Monomer Libraries',
    browsePath: 'Peptides',
    icon: 'files/icons/monomers.png',
  })
  static async manageMonomerLibrariesView(): Promise<DG.View> {
    return await showManageLibrariesView(false);
  }

  // @grok.decorators.func({tags: ['monomer-lib-provider'], result: {type: 'object', name: 'result'}})
  // static async getMonomerLibFileProvider(): Promise<MonomerLibFromFilesProvider> {
  //   return
  // }

  @grok.decorators.func({name: 'Monomer Manager Tree Browser', meta: {role: 'appTreeBrowser'}})
  static async manageMonomerLibrariesViewTreeBrowser(treeNode: DG.TreeViewGroup) {
    const libraries = (await (await MonomerLibManager.getInstance()).getAvaliableLibraryNames());
    libraries.forEach((libName) => {
      const nodeName = libName.endsWith('.json') ? libName.substring(0, libName.length - 5) : libName;
      const libNode = treeNode.item(nodeName);
      // eslint-disable-next-line rxjs/no-ignored-subscription, rxjs/no-async-subscribe
      libNode.onSelected.subscribe(async () => {
        const monomerManager = await MonomerManager.getInstance();
        await monomerManager.getViewRoot(libName, true);
        monomerManager.resetCurrentRowFollowing();
      });
    });
  }

  @grok.decorators.fileExporter({description: 'As FASTA...'})
  static saveAsFasta() {
    saveAsFastaUI();
  }

  @grok.decorators.func({
    name: 'Bio Substructure Filter',
    description: 'Substructure filter for macromolecules',
    meta: {semType: 'Macromolecule', role: 'filter'},
    outputs: [{type: 'filter', name: 'result'}],
  })
  static bioSubstructureFilter(): BioSubstructureFilter {
    return new BioSubstructureFilter(_package.seqHelper, _package.logger);
  }

  @grok.decorators.func({
    name: 'Bio Substructure Filter Test',
    description: 'Substructure filter for Helm package tests',
    outputs: [{name: 'result', type: 'object'}]
  })
  static bioSubstructureFilterTest(): BioSubstructureFilter {
    return new BioSubstructureFilter(_package.seqHelper, _package.logger);
  }

  // -- Test apps --

  @grok.decorators.func()
  static async webLogoLargeApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('WebLogo');
    try {
      const urlParams = new URLSearchParams(window.location.search);
      const app = new WebLogoApp(urlParams, 'webLogoLargeApp');
      const df: DG.DataFrame = await _package.files.readCsv('data/sample_PT_100000x5.csv');
      await grok.data.detectSemanticTypes(df);
      await app.init(df);
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async webLogoAggApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('WebLogo ...');
    try {
      const urlParams = new URLSearchParams(window.location.search);
      const app = new WebLogoApp(urlParams, 'webLogoAggApp');
      const df: DG.DataFrame = await _package.files.readCsv('samples/FASTA_PT_activity.csv');
      await grok.data.detectSemanticTypes(df);
      await app.init(df);
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async getRegionApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('getRegion ...');
    try {
      const urlParams = new URLSearchParams(window.location.search);
      const app = new GetRegionApp(urlParams, 'getRegionApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async getRegionHelmApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('getRegion ...');
    try {
      const urlParams = new URLSearchParams(window.location.search);
      const df = await _package.files.readCsv('samples/HELM_empty_vals.csv');
      const app = new GetRegionApp(urlParams, 'getRegionHelmApp');
      await app.init({df: df, colName: 'HELM'});
    } finally {
      pi.close();
    }
  }

  // -- Tests long seq --

  //name: longSeqTableSeparator
  @grok.decorators.func()
  static longSeqTableSeparator(): void {
    const df = DG.DataFrame.fromColumns(generateLongSequence());
    grok.shell.addTableView(df);
  }

  //name: longSeqTableFasta
  @grok.decorators.func()
  static longSeqTableFasta(): void {
    const df = DG.DataFrame.fromColumns([generateLongSequence2(_package.seqHelper, NOTATION.FASTA)]);
    grok.shell.addTableView(df);
  }

  //name: longSeqTableHelm
  @grok.decorators.func()
  static longSeqTableHelm(): void {
    const df = DG.DataFrame.fromColumns([generateLongSequence2(_package.seqHelper, NOTATION.HELM)]);
    grok.shell.addTableView(df);
  }

  // -- Handle context menu --

  @grok.decorators.func()
  static addCopyMenu(
    @grok.decorators.param({type: 'object'})cell: DG.Cell,
    @grok.decorators.param({type: 'object'}) menu: DG.Menu): void {
    addCopyMenuUI(cell, menu, _package.seqHelper);
  }

  // -- Demo --
  // demoBio01

  @grok.decorators.demo({
    description: 'Sequence similarity tracking and evaluation dataset diversity',
    demoPath: 'Bioinformatics | Similarity, Diversity',
    path: '/apps/Tutorials/Demo/Bioinformatics/Similarity,%20Diversity',
  })
  static async demoBioSimilarityDiversity(): Promise<void> {
    await demoBioSimDiv();
  }

  @grok.decorators.demo({
    description: 'Exploring sequence space of Macromolecules, comparison with hierarchical clustering results',
    demoPath: 'Bioinformatics | Sequence Space',
    path: '/apps/Tutorials/Demo/Bioinformatics/Sequence%20Space',
    meta: {
      isDemoDashboard: 'true'
    }
  })
  static async demoBioSequenceSpace(): Promise<void> {
    await demoSeqSpace();
  }

  @grok.decorators.demo({
    description: 'Activity Cliffs analysis on Macromolecules data',
    demoPath: 'Bioinformatics | Activity Cliffs',
    path: '/apps/Tutorials/Demo/Bioinformatics/Activity%20Cliffs',
  })
  static async demoBioActivityCliffs(): Promise<void> {
    await demoActivityCliffsCyclic();
  }

  @grok.decorators.demo({
    description: 'Atomic level structure of Macromolecules',
    demoPath: 'Bioinformatics | Atomic Level',
    path: '/apps/Tutorials/Demo/Bioinformatics/Atomic%20Level',
  })
  static async demoBioAtomicLevel(): Promise<void> {
    await demoToAtomicLevel();
  }

  @grok.decorators.func({name: 'SDF to JSON Library'})
  static async sdfToJsonLib(table: DG.DataFrame) {
    const _jsonMonomerLibrary = createJsonMonomerLibFromSdf(table);
    const jsonMonomerLibrary = JSON.stringify(_jsonMonomerLibrary);
    DG.Utils.download(`${table.name}.json`, jsonMonomerLibrary);
  }

  // -- Utils --

  @grok.decorators.func({
    friendlyName: 'seq2atomic',
    description: 'Converts a `Macromolecule` sequence to its atomic level `Molecule` representation',
    outputs: [{name: 'molfile', type: 'string', options: {semType: 'Molecule'}}]
  })
  static async seq2atomic(
    @grok.decorators.param({options: {semType: 'Macromolecule'}})
      seq: string,
      nonlinear: boolean
  ): Promise<string | undefined> {
    if (!(seq.trim())) return '';
    return PackageFunctions.toAtomicLevelSingleSeq(seq);
  }

  // //description: Gets similarity to a reference sequence
  // //input: string seq { semType: Macromolecule }
  // //input: string ref { semType: Macromolecule }
  // //output: double result
  // export async function seqSimilarity(seq: string, ref: string): Promise<number> {
  //   // if (!(seq.trim())) return null;
  //   try {
  //     const seqCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, `seq`, [seq]);
  //     const df = DG.DataFrame.fromColumns([seqCol]);
  //     const semType = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
  //     if (semType) seqCol.semType = semType;
  //
  //     const resCol = await calculateScoresWithEmptyValues(df, seqCol, ref, SCORE.SIMILARITY);
  //     return resCol.get(0)!;
  //   } catch (err: any) {
  //     const [errMsg, errStack] = errInfo(err);
  //     _package.logger.error(errMsg, undefined, errStack);
  //     throw err;
  //   }
  // }

  @grok.decorators.func({
    name: 'seqIdentity',
    friendlyName: 'seqIdentity',
    description: 'Gets identity to a reference sequence',
  })
  static async seqIdentity(
    @grok.decorators.param({options: {semType: 'Macromolecule'}})
      seq: string,
    @grok.decorators.param({options: {semType: 'Macromolecule'}})
      ref: string
  ): Promise<number | null> {
    if (!(seq.trim())) return null;
    try {
      const seqCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, `seq`, [seq]);
      const df = DG.DataFrame.fromColumns([seqCol]);
      const semType = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
      if (!semType) throw new Error('Macromolecule required');

      const resCol = await calculateScoresWithEmptyValues(df, seqCol, ref, SCORE.IDENTITY, _package.seqHelper);
      return resCol.get(0);
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
      throw err;
    }
  }

  @grok.decorators.func()
  static async detectMacromoleculeProbe(
    file: DG.FileInfo,
    colName: string = '',
    @grok.decorators.param({options: {initialValue: '100'}})
    probeCount: number = 100
  ): Promise<void> {
    const csv: string = await file.readAsString();
    await detectMacromoleculeProbeDo(csv, colName, probeCount);
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static async getSeqHelper(): Promise<ISeqHelper> {
    await PackageFunctions.initBio();
    return _package.seqHelper;
  }

  @grok.decorators.func()
  static async getMolFromHelm(
    df: DG.DataFrame,
    @grok.decorators.param({type: 'column'})helmCol: DG.Column<string>,
    @grok.decorators.param({options: {initialValue: 'true'}})
    chiralityEngine: boolean = true
  ): Promise<DG.Column<string>> {
    return getMolColumnFromHelm(df, helmCol, chiralityEngine, _package.monomerLib);
  }
}


//export let hydrophobPalette: SeqPaletteCustom | null = null;

export class SeqPaletteCustom implements SeqPalette {
  private readonly _palette: { [m: string]: string };

  constructor(palette: { [m: string]: string }) {
    this._palette = palette;
  }

  public get(m: string): string {
    return this._palette[m];
  }
}

async function initBioInt() {
  const logPrefix = 'Bio: _package.initBio()';
  _package.logger.debug(`${logPrefix}, start`);
  const t1: number = window.performance.now();
  // very important that loading should happen in correct order!
  // first make sure chem and rdkit module are loaded
  const rdKitModule = await getRdKitModule();
  // then load package settings
  const pkgProps = await _package.getProperties();
  const bioPkgProps = new BioPackageProperties(pkgProps);
  _package.properties = bioPkgProps;
  // then load monomer lib
  const libHelper = await MonomerLibManager.getInstance();
  // Fix user lib settings for explicit stuck from a terminated test
  const libSettings = await getUserLibSettings();
  if (libSettings.explicit) {
    libSettings.explicit = [];
    await setUserLibSettings(libSettings);
  }
  await libHelper.awaitLoaded(Infinity);
  if (!libHelper.initialLoadCompleted)
    await libHelper.loadMonomerLib();
  // Do not wait for monomers and sets loaded
  libHelper.loadMonomerSets();
  const monomerLib = libHelper.getMonomerLib();
  const monomerSets = libHelper.getMonomerSets();
  // finally log
  const t2: number = window.performance.now();
  _package.logger.debug(`${logPrefix}, loading ET: ${t2 - t1} ms`);

  // const monomers: string[] = [];
  // const logPs: number[] = [];

  const seqHelper = new SeqHelper(libHelper, rdKitModule);
  _package.completeInit(seqHelper, monomerLib, monomerSets, rdKitModule);

  // NB! do not delete the code below. not used now but in future we might use hydrophobicity palette
  // const series = monomerLib!.getMonomerMolsByPolymerType('PEPTIDE')!;
  // Object.keys(series).forEach((symbol) => {
  //   monomers.push(symbol);
  //   const block = series[symbol].replaceAll('#R', 'O ');
  //   const mol = rdKitModule.get_mol(block);
  //   const logP = JSON.parse(mol.get_descriptors()).CrippenClogP;
  //   logPs.push(logP);
  //   mol?.delete();
  // });

  // const sum = logPs.reduce((a, b) => a + b, 0);
  // const avg = (sum / logPs.length) || 0;

  // const palette: { [monomer: string]: string } = {};
  // for (let i = 0; i < monomers.length; i++)
  //   palette[monomers[i]] = logPs[i] < avg ? '#4682B4' : '#DC143C';

  // hydrophobPalette = new SeqPaletteCustom(palette);

  _package.logger.debug(`${logPrefix}, end`);
  handleSequenceHeaderRendering();
}

// -- Package settings editor --

// //name: packageSettingsEditor
// //description: The database connection
// //tags: packageSettingsEditor
// //input: object propList
// //output: widget result
// export function packageSettingsEditor(propList: DG.Property[]): DG.Widget {
//   const widget = new PackageSettingsEditorWidget(propList);
//   widget.init().then(); // Ignore promise returned
//   return widget as DG.Widget;
// }

//name: test1
//output: object result
export function test1(): any {
  _package.logger.debug('Bio:test1() function');
  return {value: 'value1'};
}
