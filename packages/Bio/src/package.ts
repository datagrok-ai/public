/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  MacromoleculeDifferenceCellRenderer, MacromoleculeSequenceCellRenderer,
} from './utils/cell-renderer';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {SequenceAlignment} from './seq_align';
import {ISequenceSpaceResult, getEmbeddingColsNames, getSequenceSpace} from './analysis/sequence-space';
import {ISequenceSpaceParams, getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {
  createLinesGrid, createPropPanelElement, createTooltipElement, getChemSimilaritiesMatrix,
} from './analysis/sequence-activity-cliffs';
import {convert} from './utils/convert';
import {getMacromoleculeColumnPropertyPanel} from './widgets/representations';
import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {FastaFileHandler} from '@datagrok-libraries/bio/src/utils/fasta-handler';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils';

import {SequenceSimilarityViewer} from './analysis/sequence-similarity-viewer';
import {SequenceDiversityViewer} from './analysis/sequence-diversity-viewer';
import {SubstructureSearchDialog} from './substructure-search/substructure-search';
import {saveAsFastaUI} from './utils/save-as-fasta';
import {BioSubstructureFilter} from './widgets/bio-substructure-filter';
import {delay} from '@datagrok-libraries/utils/src/test';
import {
  TAGS as bioTAGS, ALPHABET, NOTATION,
} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {WebLogoViewer} from './viewers/web-logo-viewer';
import {createJsonMonomerLibFromSdf, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  MonomerLibHelper,
  getUserLibSettings,
  setUserLibSetting,
  getLibFileNameList,
  getLibraryPanelUI
} from './utils/monomer-lib';
import {getMacromoleculeColumns} from './utils/ui-utils';
import {DimReductionMethods, ITSNEOptions, IUMAPOptions} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {SequenceSpaceFunctionEditor} from '@datagrok-libraries/ml/src/functionEditors/seq-space-editor';
import {ActivityCliffsFunctionEditor} from '@datagrok-libraries/ml/src/functionEditors/activity-cliffs-editor';
import {SCORE, calculateScores} from '@datagrok-libraries/bio/src/utils/macromolecule/scoring';

import {demoBio01UI} from './demo/bio01-similarity-diversity';
import {demoBio01aUI} from './demo/bio01a-hierarchical-clustering-and-sequence-space';
import {demoBio01bUI} from './demo/bio01b-hierarchical-clustering-and-activity-cliffs';
import {demoBio03UI} from './demo/bio03-atomic-level';
import {demoBio05UI} from './demo/bio05-helm-msa-sequence-space';
import {checkInputColumnUI} from './utils/check-input-column';
import {multipleSequenceAlignmentUI} from './utils/multiple-sequence-alignment-ui';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {WebLogoApp} from './apps/web-logo-app';
import {SplitToMonomersFunctionEditor} from './function-edtiors/split-to-monomers-editor';
import {splitToMonomersUI} from './utils/split-to-monomers';
import {MonomerCellRenderer} from './utils/monomer-cell-renderer';
import {BioPackage, BioPackageProperties} from './package-types';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {PackageSettingsEditorWidget} from './widgets/package-settings-editor-widget';
import {getCompositionAnalysisWidget} from './widgets/composition-analysis-widget';
import {MacromoleculeColumnWidget} from './utils/macromolecule-column-widget';
import {addCopyMenuUI} from './utils/context-menu';
import {_getEnumeratorWidget, _setPeptideColumn} from './utils/enumerator-tools';
import {getRegionDo} from './utils/get-region';
import {GetRegionApp} from './apps/get-region-app';
import {GetRegionFuncEditor} from './utils/get-region-func-editor';
import {HelmToMolfileConverter} from './utils/helm-to-molfile';
import {DIMENSIONALITY_REDUCER_TERMINATE_EVENT}
  from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import { sequenceToMolfile } from './utils/sequence-to-mol';

export const _package = new BioPackage();

export const BYPASS_LARGE_DATA_WARNING = 'bypassLargeDataWarning';
// /** Avoid reassigning {@link monomerLib} because consumers subscribe to {@link IMonomerLib.onChanged} event */
// let monomerLib: MonomerLib | null = null;

//name: getMonomerLibHelper
//description:
//output: object result
export function getMonomerLibHelper(): IMonomerLibHelper {
  return MonomerLibHelper.instance;
}

export let hydrophobPalette: SeqPaletteCustom | null = null;

export class SeqPaletteCustom implements SeqPalette {
  private readonly _palette: { [m: string]: string };

  constructor(palette: { [m: string]: string }) {
    this._palette = palette;
  }

  public get(m: string): string {
    return this._palette[m];
  }
}

//tags: init
export async function initBio() {
  _package.logger.debug('Bio: initBio(), started');
  const module = await grok.functions.call('Chem:getRdKitModule');
  await Promise.all([
    (async () => { await MonomerLibHelper.instance.loadLibraries(); })(),
    (async () => {
      const pkgProps = await _package.getProperties();
      const bioPkgProps = new BioPackageProperties(pkgProps);
      _package.properties = bioPkgProps;
    })(),
  ]).finally(() => {
    _package.completeInit();
  });

  const monomerLib = MonomerLibHelper.instance.getBioLib();
  const monomers: string[] = [];
  const logPs: number[] = [];

  const series = monomerLib!.getMonomerMolsByPolymerType('PEPTIDE')!;
  Object.keys(series).forEach((symbol) => {
    monomers.push(symbol);
    const block = series[symbol].replaceAll('#R', 'O ');
    const mol = module.get_mol(block);
    const logP = JSON.parse(mol.get_descriptors()).CrippenClogP;
    logPs.push(logP);
    mol?.delete();
  });

  const sum = logPs.reduce((a, b) => a + b, 0);
  const avg = (sum / logPs.length) || 0;

  const palette: { [monomer: string]: string } = {};
  for (let i = 0; i < monomers.length; i++)
    palette[monomers[i]] = logPs[i] < avg ? '#4682B4' : '#DC143C';

  hydrophobPalette = new SeqPaletteCustom(palette);

  _package.logger.debug('Bio: initBio(), completed');
}

//name: sequenceTooltip
//tags: tooltip
//input: column col {semType: Macromolecule}
//output: widget result
export function sequenceTooltip(col: DG.Column): DG.Widget<any> {
  const resWidget = new MacromoleculeColumnWidget(col);
  const _resPromise = resWidget.init().then(() => { })
    .catch((err: any) => {
      const errMsg = err instanceof Error ? err.message : err.toString();
      grok.shell.error(errMsg);
    });
  return resWidget;
}

//name: getBioLib
//output: object monomerLib
export function getBioLib(): IMonomerLib {
  return MonomerLibHelper.instance.getBioLib();
}

// -- Panels --

//name: Get Region
//description: Creates a new column with sequences of the region between start and end
//tags: panel
//input: column seqCol {semType: Macromolecule}
//output: widget result
export function getRegionPanel(seqCol: DG.Column<string>): DG.Widget {
  const funcName: string = 'getRegionTopMenu';
  const funcList = DG.Func.find({package: _package.name, name: funcName});
  if (funcList.length !== 1) throw new Error(`Package '${_package.name}' func '${funcName}' not found`);
  const func = funcList[0];
  const funcCall = func.prepare({table: seqCol.dataFrame, sequence: seqCol});
  const funcEditor = new GetRegionFuncEditor(funcCall);
  return funcEditor.widget();
}

//name: Manage Libraries
//description:
//tags: panel, exclude-actions-panel
//input: column seqColumn {semType: Macromolecule}
//output: widget result
export async function libraryPanel(_seqColumn: DG.Column): Promise<DG.Widget> {
  return getLibraryPanelUI();
}

// -- Func Editors --

//name: GetRegionEditor
//tags: editor
//input: funccall call
export function GetRegionEditor(call: DG.FuncCall): void {
  try {
    const funcEditor = new GetRegionFuncEditor(call);
    funcEditor.dialog();
  } catch (err: any) {
    const errMsg = err instanceof Error ? err.message : err.toString();
    const errStack = err instanceof Error ? err.stack : undefined;
    grok.shell.error(`Get region editor error: ${errMsg}`);
    _package.logger.error(errMsg, undefined, errStack);
  }
}

//name: SplitToMonomersEditor
//tags: editor
//input: funccall call
export function SplitToMonomersEditor(call: DG.FuncCall): void {
  const funcEditor = new SplitToMonomersFunctionEditor();
  ui.dialog({title: 'Split to Monomers'})
    .add(funcEditor.paramsUI)
    .onOK(async () => {
      return call.func.prepare(funcEditor.funcParams).call(true);
    })
    .show();
}

//name: SequenceSpaceEditor
//tags: editor
//input: funccall call
export function SequenceSpaceEditor(call: DG.FuncCall) {
  const funcEditor = new SequenceSpaceFunctionEditor(DG.SEMTYPE.MACROMOLECULE);
  ui.dialog({title: 'Sequence Space'})
    .add(funcEditor.paramsUI)
    .onOK(async () => {
      return call.func.prepare(funcEditor.funcParams).call();
    })
    .show();
}

//name: SeqActivityCliffsEditor
//tags: editor
//input: funccall call
export function SeqActivityCliffsEditor(call: DG.FuncCall) {
  const funcEditor = new ActivityCliffsFunctionEditor(DG.SEMTYPE.MACROMOLECULE);
  ui.dialog({title: 'Activity Cliffs'})
    .add(funcEditor.paramsUI)
    .onOK(async () => {
      return call.func.prepare(funcEditor.funcParams).call(true);
    })
    .show();
}


// -- Package settings editor --

//name: packageSettingsEditor
//description: The database connection
//tags: packageSettingsEditor
//input: object propList
//output: widget result
export function packageSettingsEditor(propList: DG.Property[]): DG.Widget {
  const widget = new PackageSettingsEditorWidget(propList);
  widget.init().then(); // Ignore promise returned
  return widget as DG.Widget;
}

// -- Cell renderers --

//name: fastaSequenceCellRenderer
//tags: cellRenderer
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=fasta
//output: grid_cell_renderer result
export function fastaSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
  return new MacromoleculeSequenceCellRenderer();
}

// -- Property panels --

//name: Sequence Renderer
//input: column molColumn {semType: Macromolecule}
//tags: panel
//output: widget result
export function macroMolColumnPropertyPanel(molColumn: DG.Column): DG.Widget {
  return getMacromoleculeColumnPropertyPanel(molColumn);
}

//name: Composition analysis
//tags: panel, bio, widgets
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export function compositionAnalysisWidget(sequence: DG.SemanticValue): DG.Widget {
  return getCompositionAnalysisWidget(sequence);
}

//name: separatorSequenceCellRenderer
//tags: cellRenderer
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=separator
//output: grid_cell_renderer result
export function separatorSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
  return new MacromoleculeSequenceCellRenderer();
}

//name: MacromoleculeDifferenceCellRenderer
//tags: cellRenderer
//meta.cellType: MacromoleculeDifference
//meta.columnTags: quality=MacromoleculeDifference
//output: grid_cell_renderer result
export function macromoleculeDifferenceCellRenderer(): MacromoleculeDifferenceCellRenderer {
  return new MacromoleculeDifferenceCellRenderer();
}


//name: sequenceAlignment
//input: string alignType {choices: ['Local alignment', 'Global alignment']}
// eslint-disable-next-line max-len
//input: string alignTable {choices: ['AUTO', 'NUCLEOTIDES', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62','BLOSUM80','BLOSUM90','PAM30','PAM70','PAM250','SCHNEIDER','TRANS']}
//input: double gap
//input: string seq1
//input: string seq2
//output: object res
export function sequenceAlignment(alignType: string, alignTable: string, gap: number, seq1: string, seq2: string) {
  const toAlign = new SequenceAlignment(seq1, seq2, gap, alignTable);
  const res = alignType == 'Local alignment' ? toAlign.smithWaterman() : toAlign.needlemanWunch();
  return res;
}

// -- Viewers --

//name: WebLogo
//description: WebLogo
//tags: viewer, panel
//output: viewer result
//meta.icon: files/icons/weblogo-viewer.svg
export function webLogoViewer() {
  return new WebLogoViewer();
}

//name: VdRegions
//description: V-Domain regions viewer
//tags: viewer, panel
//meta.icon: files/icons/vdregions-viewer.svg
//output: viewer result
export function vdRegionsViewer() {
  return new VdRegionsViewer();
}


// -- Top menu --

//name: getRegion
//description: Gets a new column with sequences of the region between start and end
//input: column sequence
//input: string start  {optional: true}
//input: string end    {optional: true}
//input: string name   {optional: true} [Name of the column to be created]
//output: column result
export function getRegion(
  sequence: DG.Column<string>, start: string | undefined, end: string | undefined, name: string | undefined
): DG.Column<string> {
  return getRegionDo(sequence,
    start ?? null, end ?? null, name ?? null);
}

//top-menu: Bio | Convert | Get Region...
//name: Get Region Top Menu
//description: Get sequences for a region specified from a Macromolecule
//input: dataframe table                           [Input data table]
//input: column sequence  {semType: Macromolecule} [Sequence column]
//input: string start     {optional: true}         [Region start position name]
//input: string end       {optional: true}         [Region end position name]
//input: string name      {optional: true}         [Region column name]
//editor: Bio:GetRegionEditor
export async function getRegionTopMenu(
  table: DG.DataFrame, sequence: DG.Column,
  start: string | undefined, end: string | undefined, name: string | undefined
): Promise<void> {
  const regCol = getRegionDo(sequence, start ?? null, end ?? null, name ?? null);
  sequence.dataFrame.columns.add(regCol);
  await grok.data.detectSemanticTypes(sequence.dataFrame); // to set renderer
}

//top-menu: Bio | Analyze | Activity Cliffs...
//name: Sequence Activity Cliffs
//description: Detects pairs of molecules with similar structure and significant difference in any given property
//input: dataframe table [Input data table]
//input: column molecules {semType: Macromolecule}
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE"] }
//input: object options {optional: true}
//output: viewer result
//editor: Bio:SeqActivityCliffsEditor
export async function activityCliffs(df: DG.DataFrame, macroMolecule: DG.Column<string>, activities: DG.Column,
  similarity: number, methodName: DimReductionMethods, options?: (IUMAPOptions | ITSNEOptions) & Options,
): Promise<DG.Viewer | undefined> {
  if (!checkInputColumnUI(macroMolecule, 'Activity Cliffs'))
    return;
  const axesNames = getEmbeddingColsNames(df);
  const tags = {
    'units': macroMolecule.getTag(DG.TAGS.UNITS),
    'aligned': macroMolecule.getTag(bioTAGS.aligned),
    'separator': macroMolecule.getTag(bioTAGS.separator),
    'alphabet': macroMolecule.getTag(bioTAGS.alphabet),
  };
  const ncUH = UnitsHandler.getOrCreate(macroMolecule);
  let columnDistanceMetric: BitArrayMetrics | MmDistanceFunctionsNames = BitArrayMetricsNames.Tanimoto;
  let seqCol = macroMolecule;
  if (ncUH.isFasta() || (ncUH.isSeparator() && ncUH.alphabet && ncUH.alphabet !== ALPHABET.UN)) {
    if (ncUH.isFasta()) {
      columnDistanceMetric = ncUH.getDistanceFunctionName();
    } else {
      seqCol = ncUH.convert(NOTATION.FASTA);
      const uh = UnitsHandler.getOrCreate(seqCol);
      columnDistanceMetric = uh.getDistanceFunctionName();
      tags.units = NOTATION.FASTA;
    }
  }
  const runCliffs = async () => {
    const sp = await getActivityCliffs(
      df,
      seqCol,
      null,
      axesNames,
      'Activity cliffs', //scatterTitle
      activities,
      similarity,
      columnDistanceMetric, //similarityMetric
      methodName,
      DG.SEMTYPE.MACROMOLECULE,
      tags,
      getSequenceSpace,
      getChemSimilaritiesMatrix,
      createTooltipElement,
      createPropPanelElement,
      createLinesGrid,
      options);
    return sp;
  };

  const allowedRowCount = 20000;
  const fastRowCount = methodName === DimReductionMethods.UMAP ? 5000 : 2000;
  if (df.rowCount > allowedRowCount) {
    grok.shell.warning(`Too many rows, maximum for sequence activity cliffs is ${allowedRowCount}`);
    return;
  }

  if (df.rowCount > fastRowCount && !options?.[BYPASS_LARGE_DATA_WARNING]) {
    ui.dialog().add(ui.divText(`Activity cliffs analysis might take several minutes.
    Do you want to continue?`))
      .onOK(async () => {
        const progressBar = DG.TaskBarProgressIndicator.create(`Running sequence activity cliffs ...`);
        const res = await runCliffs();
        progressBar.close();
        return res;
      })
      .show();
  } else {
    const res = await runCliffs();
    return res;
  }
}

//top-menu: Bio | Analyze | Sequence Space...
//name: Sequence Space
//description: Creates 2D sequence space with projected sequences by pairwise distance
//input: dataframe table
//input: column molecules { semType: Macromolecule }
//input: string methodName { choices:["UMAP", "t-SNE"] }
//input: string similarityMetric { choices:["Hamming", "Levenshtein", "Monomer chemical distance"] }
//input: bool plotEmbeddings = true
//input: double sparseMatrixThreshold = 0.8 [Similarity Threshold for sparse matrix calculation]
//input: object options {optional: true}
//editor: Bio:SequenceSpaceEditor
export async function sequenceSpaceTopMenu(
  table: DG.DataFrame, macroMolecule: DG.Column, methodName: DimReductionMethods,
  similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames = MmDistanceFunctionsNames.LEVENSHTEIN,
  plotEmbeddings: boolean, sparseMatrixThreshold?: number, options?: (IUMAPOptions | ITSNEOptions) & Options,
): Promise<DG.Viewer | undefined> {
  // Delay is required for initial function dialog to close before starting invalidating of molfiles.
  // Otherwise, dialog is freezing
  await delay(10);
  if (!checkInputColumnUI(macroMolecule, 'Sequence space')) return;
  let scatterPlot: DG.ScatterPlotViewer | undefined = undefined;
  const pg = DG.TaskBarProgressIndicator.create('Initializing sequence space ...');
  // function for progress of umap
  try {
    function progressFunc(_nEpoch: number, epochsLength: number, embeddings: number[][]) {
      let embedXCol: DG.Column | null = null;
      let embedYCol: DG.Column | null = null;
      if (!table.columns.names().includes(embedColsNames[0])) {
        embedXCol = table.columns.add(DG.Column.float(embedColsNames[0], table.rowCount));
        embedYCol = table.columns.add(DG.Column.float(embedColsNames[1], table.rowCount));
        if (plotEmbeddings) {
          scatterPlot = grok.shell
            .tableView(table.name)
            .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Sequence space'});
        }
      } else {
        embedXCol = table.columns.byName(embedColsNames[0]);
        embedYCol = table.columns.byName(embedColsNames[1]);
      }

      embedXCol.init((i) => embeddings[i] ? embeddings[i][0] : undefined);
      embedYCol.init((i) => embeddings[i] ? embeddings[i][1] : undefined);
      const progress = (_nEpoch / epochsLength * 100);
      pg.update(progress, `Running sequence space ... ${progress.toFixed(0)}%`);
    }
    const embedColsNames = getEmbeddingColsNames(table);
    const withoutEmptyValues = DG.DataFrame.fromColumns([macroMolecule]).clone();
    const emptyValsIdxs = removeEmptyStringRows(withoutEmptyValues, macroMolecule);

    const chemSpaceParams: ISequenceSpaceParams = {
      seqCol: withoutEmptyValues.col(macroMolecule.name)!,
      methodName: methodName,
      similarityMetric: similarityMetric,
      embedAxesNames: embedColsNames,
      options: {...options, sparseMatrixThreshold: sparseMatrixThreshold ?? 0.8,
        usingSparseMatrix: table.rowCount > 20000},
    };

    const allowedRowCount = methodName === DimReductionMethods.UMAP ? 100000 : 15000;
    // number of rows which will be processed relatively fast
    const fastRowCount = methodName === DimReductionMethods.UMAP ? 5000 : 2000;
    if (table.rowCount > allowedRowCount) {
      grok.shell.warning(`Too many rows, maximum for sequence space is ${allowedRowCount}`);
      return;
    }

    async function getSeqSpace() {
      let resolveF: Function | null = null;

      const sub = grok.events.onViewerClosed.subscribe((args) => {
        const v = args.args.viewer as unknown as DG.Viewer<any>;
        if (v?.getOptions()?.look?.title && scatterPlot?.getOptions()?.look?.title &&
          v?.getOptions()?.look?.title === scatterPlot?.getOptions()?.look?.title) {
          grok.events.fireCustomEvent(DIMENSIONALITY_REDUCER_TERMINATE_EVENT, {});
          sub.unsubscribe();
          resolveF?.();
          pg.close();
        }
      });
      const sequenceSpaceResPromise = new Promise<ISequenceSpaceResult | undefined>(async (resolve) => {
        resolveF = resolve;
        const res = await getSequenceSpace(chemSpaceParams,
          options?.[BYPASS_LARGE_DATA_WARNING] ? undefined : progressFunc);
        resolve(res);
      });
      const sequenceSpaceRes = await sequenceSpaceResPromise;
      pg.close();
      sub.unsubscribe();
      return sequenceSpaceRes ? processResult(sequenceSpaceRes) : sequenceSpaceRes;
    }

    if (table.rowCount > fastRowCount && !options?.[BYPASS_LARGE_DATA_WARNING]) {
      ui.dialog().add(ui.divText(`Sequence space analysis might take several minutes.
    Do you want to continue?`))
        .onOK(async () => {
          await getSeqSpace();
        })
        .onCancel(() => { pg.close(); })
        .show();
    } else {
      return await getSeqSpace();
    }

    function processResult(sequenceSpaceRes: ISequenceSpaceResult): DG.ScatterPlotViewer | undefined {
      const embeddings = sequenceSpaceRes.coordinates;
      for (const col of embeddings) {
        const listValues = col.toList();
        emptyValsIdxs.forEach((ind: number) => listValues.splice(ind, 0, null));
        let embedCol = table.columns.byName(col.name);
        if (!embedCol) {
          embedCol = DG.Column.float(col.name, listValues.length);
          table.columns.add(embedCol);
        }
        embedCol.init((i) => listValues[i]);
      //table.columns.add(DG.Column.float(col.name, table.rowCount).init((i) => listValues[i]));
      }
      if (plotEmbeddings) {
        if (!scatterPlot) {
          scatterPlot = grok.shell
            .tableView(table.name)
            .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Sequence space'});
        }
        return scatterPlot;
      }
    }
  } catch (e) {
    console.error(e);
    pg.close();
  }


  /*   const encodedCol = encodeMonomers(macroMolecule);
  if (!encodedCol)
    return;
  const embedColsNames = getEmbeddingColsNames(table);
  const withoutEmptyValues = DG.DataFrame.fromColumns([encodedCol]).clone();
  const emptyValsIdxs = removeEmptyStringRows(withoutEmptyValues, encodedCol);

  const chemSpaceParams = {
    seqCol: withoutEmptyValues.col(encodedCol.name)!,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: embedColsNames
  };
  const sequenceSpaceRes = await sequenceSpace(chemSpaceParams);
  const embeddings = sequenceSpaceRes.coordinates;
  for (const col of embeddings) {
    const listValues = col.toList();
    emptyValsIdxs.forEach((ind: number) => listValues.splice(ind, 0, null));
    table.columns.add(DG.Column.fromList('double', col.name, listValues));
  }
  let sp;
  if (plotEmbeddings) {
    for (const v of grok.shell.views) {
      if (v.name === table.name)
        sp = (v as DG.TableView).scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Sequence space'});
    }
  } */
}

//top-menu: Bio | Convert | To Atomic Level...
//name: To Atomic Level
//description: Converts sequences to molblocks
//input: dataframe df [Input data table]
//input: column macroMolecule {semType: Macromolecule}
//input: bool nonlinear=false { description: Slower mode for cycling/branching HELM structures }
export async function toAtomicLevel(df: DG.DataFrame, macroMolecule: DG.Column, nonlinear: boolean): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Converting to atomic level ...');
  sequenceToMolfile(df, macroMolecule, nonlinear);
  pi.close();
}

//top-menu: Bio | Analyze | MSA...
//name: MSA
//description: Performs multiple sequence alignment
//tags: bio, panel
export function multipleSequenceAlignmentDialog(): void {
  multipleSequenceAlignmentUI();
}

//name: Multiple Sequence Alignment
//description: Multiple sequence alignment
//tags: bio
//input: column sequenceCol {semType: Macromolecule}
//input: column clustersCol
//output: column result
export async function alignSequences(sequenceCol: DG.Column<string> | null = null,
  clustersCol: DG.Column | null = null): Promise<DG.Column<string>> {
  return multipleSequenceAlignmentUI({col: sequenceCol, clustersCol});
}

//top-menu: Bio | Analyze | Composition
//name: Composition Analysis
//description: Visualizes sequence composition on a WebLogo plot
//meta.icon: files/icons/composition-analysis.svg
//output: viewer result
export async function compositionAnalysis(): Promise<void> {
  // Higher priority for columns with MSA data to show with WebLogo.
  const tv = grok.shell.tv;
  const df = tv.dataFrame;
  //@ts-ignore
  const colList: DG.Column[] = df.columns.toList().filter((col) => {
    if (col.semType != DG.SEMTYPE.MACROMOLECULE)
      return false;

    const _colUH = UnitsHandler.getOrCreate(col);
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
    const selectedCol = colList.find((c) => { return UnitsHandler.getOrCreate(c).isMsa(); });
    const colInput: DG.InputBase = ui.choiceInput(
      'Column', selectedCol ? selectedCol.name : colListNames[0], colListNames);
    ui.dialog({
      title: 'Composition Analysis',
      helpUrl: '/help/domains/bio/macromolecules.md#composition-analysis',
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
  } else {
    col = colList[0];
  }

  if (!col)
    return;

  await handler(col);
}

// 2023-05-17 Representations does not work at BioIT
// //name: Representations
// //tags: panel, widgets
// //input: cell macroMolecule {semType: Macromolecule}
// //output: widget result
// export async function peptideMolecule(macroMolecule: DG.Cell): Promise<DG.Widget> {
//   const monomersLibFile = await _package.files.readAsText(HELM_CORE_LIB_FILENAME);
//   const monomersLibObject: any[] = JSON.parse(monomersLibFile);
//
//   return representationsWidget(macroMolecule, monomersLibObject);
// }

//name: importFasta
//description: Opens FASTA file
//tags: file-handler
//meta.ext: fasta, fna, ffn, faa, frn, fa, fst
//input: string fileContent
//output: list tables
export function importFasta(fileContent: string): DG.DataFrame [] {
  const ffh = new FastaFileHandler(fileContent);
  return ffh.importFasta();
}

//name: importBam
//description: Opens Bam file
//tags: file-handler
//meta.ext: bam, bai
//input: string fileContent
//output: list tables
export function importBam(fileContent: string): DG.DataFrame [] {
  console.log(fileContent);
  return [];
}

//top-menu: Bio | Convert | Notation...
//name: convertDialog
export function convertDialog() {
  const col = getMacromoleculeColumns()[0];
  convert(col);
}

//name: monomerCellRenderer
//tags: cellRenderer
//meta.cellType: Monomer
//meta.columnTags: quality=Monomer
//output: grid_cell_renderer result
export function monomerCellRenderer(): MonomerCellRenderer {
  return new MonomerCellRenderer();
}

//name: testDetectMacromolecule
//input: string path {choices: ['Demo:Files/', 'System:AppData/']}
//output: dataframe result
export async function testDetectMacromolecule(path: string): Promise<DG.DataFrame> {
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
          //  `semType: ${semType}, units: ${col.getTag(DG.TAGS.UNITS)}`);
          // console.warn('file: "' + fileInfo.path + '", semType: "' + semType + '", ' +
          //   'units: "' + col.getTag(DG.TAGS.UNITS) + '"');

          res.push({
            file: fileInfo.path, result: 'detected', column: col.name,
            message: `units: ${col.getTag(DG.TAGS.UNITS)}`,
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

//top-menu: Bio | Convert | Split to Monomers...
//name: Split to Monomers
//input: dataframe table
//input: column sequence { semType: Macromolecule }
//output: dataframe result
//editor: Bio:SplitToMonomersEditor
export async function splitToMonomersTopMenu(table: DG.DataFrame, sequence: DG.Column): Promise<void> {
  await splitToMonomersUI(table, sequence);
}

//name: Bio: getHelmMonomers
//input: column sequence {semType: Macromolecule}
export function getHelmMonomers(sequence: DG.Column<string>): string[] {
  const uh = UnitsHandler.getOrCreate(sequence);
  const stats = uh.stats;
  return Object.keys(stats.freq);
}


//name: Sequence Similarity Search
//tags: viewer
//meta.icon: files/icons/sequence-similarity-viewer.svg
//output: viewer result
export function similaritySearchViewer(): SequenceSimilarityViewer {
  return new SequenceSimilarityViewer();
}

//top-menu: Bio | Search | Similarity Search
//name: similaritySearch
//description: Finds similar sequences
//output: viewer result
export function similaritySearchTopMenu(): void {
  const view = (grok.shell.v as DG.TableView);
  const viewer = view.addViewer('Sequence Similarity Search');
  view.dockManager.dock(viewer, 'down');
}

//name: Sequence Diversity Search
//tags: viewer
//meta.icon: files/icons/sequence-diversity-viewer.svg
//output: viewer result
export function diversitySearchViewer(): SequenceDiversityViewer {
  return new SequenceDiversityViewer();
}

//top-menu: Bio | Search | Diversity Search
//name: diversitySearch
//description: Finds the most diverse sequences
//output: viewer result
export function diversitySearchTopMenu() {
  const view = (grok.shell.v as DG.TableView);
  const viewer = view.addViewer('Sequence Diversity Search');
  view.dockManager.dock(viewer, 'down');
}

//name: SearchSubsequenceEditor
//tags: editor
//input: funccall call
export function searchSubsequenceEditor(call: DG.FuncCall) {
  const columns = getMacromoleculeColumns();
  if (columns.length === 1)
    call.func.prepare({macromolecules: columns[0]}).call(true);
  else
    new SubstructureSearchDialog(columns);
}

//top-menu: Bio | Search | Subsequence Search ...
//name: Subsequence Search
//input: column macromolecules
//editor: Bio:SearchSubsequenceEditor
export function SubsequenceSearchTopMenu(macromolecules: DG.Column): void {
  grok.shell.tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
    type: 'Bio:bioSubstructureFilter',
    column: macromolecules.name,
    columnName: macromolecules.name,
  });
  grok.shell.tv.grid.scrollToCell(macromolecules, 0);
}

//top-menu: Bio | Calculate | Identity...
//name: Identity Scoring
//description: Adds a column with fraction of matching monomers
//input: dataframe table [Table containing Macromolecule column]
//input: column macromolecules {semType: Macromolecule} [Sequences to score]
//input: string reference [Sequence, matching column format]
//output: column scores
export async function sequenceIdentityScoring(
  table: DG.DataFrame, macromolecule: DG.Column, reference: string
): Promise<DG.Column<number>> {
  const scores = calculateScores(table, macromolecule, reference, SCORE.IDENTITY);
  return scores;
}

//top-menu: Bio | Calculate | Similarity...
//name: Similarity Scoring
//description: Adds a column with similarity scores, calculated as sum of monomer fingerprint similarities
//input: dataframe table [Table containing Macromolecule column]
//input: column macromolecules {semType: Macromolecule} [Sequences to score]
//input: string reference [Sequence, matching column format]
//output: column scores
export async function sequenceSimilarityScoring(
  table: DG.DataFrame, macromolecule: DG.Column, reference: string
): Promise<DG.Column<number>> {
  const scores = calculateScores(table, macromolecule, reference, SCORE.SIMILARITY);
  return scores;
}

//name: saveAsFasta
//description: As FASTA...
//tags: fileExporter
export function saveAsFasta() {
  saveAsFastaUI();
}

//name: Bio Substructure Filter
//description: Substructure filter for macromolecules
//tags: filter
//output: filter result
//meta.semType: Macromolecule
export function bioSubstructureFilter(): BioSubstructureFilter {
  return new BioSubstructureFilter();
}

// -- Test apps --

//name: webLogoLargeApp
export async function webLogoLargeApp(): Promise<void> {
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

//name: getRegionApp
export async function getRegionApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('getRegion ...');
  try {
    const urlParams = new URLSearchParams(window.location.search);
    const app = new GetRegionApp(urlParams, 'getRegionApp');
    await app.init();
  } finally {
    pi.close();
  }
}

//name: getRegionHelmApp
export async function getRegionHelmApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('getRegion ...');
  try {
    const urlParams = new URLSearchParams(window.location.search);
    const df = await _package.files.readCsv('data/sample_HELM_empty_vals.csv');
    const app = new GetRegionApp(urlParams, 'getRegionHelmApp');
    await app.init({df: df, colName: 'HELM'});
  } finally {
    pi.close();
  }
}

// -- Handle context menu --

///name: addCopyMenu
//input: object cell
//input: object menu
export function addCopyMenu(cell: DG.Cell, menu: DG.Menu): void {
  addCopyMenuUI(cell, menu);
}

// -- Demo --

// demoBio01
//name: demoBioSimilarityDiversity
//meta.demoPath: Bioinformatics | Similarity, Diversity
//description: Sequence similarity tracking and evaluation dataset diversity
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Similarity,%20Diversity
//meta.isDemoScript: True
export async function demoBioSimilarityDiversity(): Promise<void> {
  await demoBio01UI();
}

// demoBio01a
//name:demoBioSequenceSpace
//meta.demoPath: Bioinformatics | Sequence Space
//description: Exploring sequence space of Macromolecules, comparison with hierarchical clustering results
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Sequence%20Space
//meta.isDemoScript: True
export async function demoBioSequenceSpace(): Promise<void> {
  await demoBio01aUI();
}

// demoBio01b
//name: demoBioActivityCliffs
//meta.demoPath: Bioinformatics | Activity Cliffs
//description: Activity Cliffs analysis on Macromolecules data
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Activity%20Cliffs
//meta.isDemoScript: True
export async function demoBioActivityCliffs(): Promise<void> {
  await demoBio01bUI();
}

// demoBio03
//name: demoBioAtomicLevel
//meta.demoPath: Bioinformatics | Atomic Level
//description: Atomic level structure of Macromolecules
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Atomic%20Level
//meta.isDemoScript: True
export async function demoBioAtomicLevel(): Promise<void> {
  await demoBio03UI();
}

// demoBio05
//name: demoBioHelmMsaSequenceSpace
//meta.demoPath: Bioinformatics | Helm, MSA, Sequence Space
//description: MSA and composition analysis on Helm data
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Helm,%20MSA,%20Sequence%20Space
//meta.isDemoScript: True
export async function demoBioHelmMsaSequenceSpace(): Promise<void> {
  await demoBio05UI();
}

//name: enumeratorColumnChoice
//input: dataframe df [Input data table]
//input: column macroMolecule
export async function enumeratorColumnChoice(df: DG.DataFrame, macroMolecule: DG.Column): Promise<void> {
  _setPeptideColumn(macroMolecule);
  await grok.data.detectSemanticTypes(df);
}

//name: PolyTool
//input: column molColumn {semType: Macromolecule}
//tags: panel, exclude-actions-panel
//output: widget result
export function getEnumeratorWidget(molColumn: DG.Column): DG.Widget {
  return _getEnumeratorWidget(molColumn);
}


//top-menu: Bio | Convert | SDF to JSON Library...
//name: SDF to JSON Library
//input: dataframe table
export async function sdfToJsonLib(table: DG.DataFrame) {
  const _jsonMonomerLibrary = createJsonMonomerLibFromSdf(table);
  const jsonMonomerLibrary = JSON.stringify(_jsonMonomerLibrary);
  DG.Utils.download(`${table.name}.json`, jsonMonomerLibrary);
}
