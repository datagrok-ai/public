/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

export const _package = new DG.Package();

import {MacromoleculeDifferenceCellRenderer, MonomerCellRenderer} from './utils/cell-renderer';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {runKalign, testMSAEnoughMemory} from './utils/multiple-sequence-alignment';
import {SequenceAlignment, Aligned} from './seq_align';
import {getEmbeddingColsNames, sequenceSpace, sequenceSpaceByFingerprints} from './analysis/sequence-space';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {createLinesGrid, createPropPanelElement, createTooltipElement, getChemSimilaritiesMarix, getSimilaritiesMarix} from './analysis/sequence-activity-cliffs';
import {createJsonMonomerLibFromSdf, encodeMonomers, getMolfilesFromSeq} from '@datagrok-libraries/bio/src/utils/monomer-utils';
import {HELM_CORE_LIB_FILENAME} from '@datagrok-libraries/bio/src/utils/const';
import {getMacroMol} from './utils/atomic-works';
import {MacromoleculeSequenceCellRenderer} from './utils/cell-renderer';
import {convert} from './utils/convert';
import {getMacroMolColumnPropertyPanel, representationsWidget} from './widgets/representations';
import {TAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule'
import {_toAtomicLevel} from '@datagrok-libraries/bio/src/utils/to-atomic-level';
import {FastaFileHandler} from '@datagrok-libraries/bio/src/utils/fasta-handler';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils';
import {
  generateManySequences,
  generateLongSequence,
  performanceTest
} from './tests/test-sequnces-generators';

import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import * as C from './utils/constants';
import {SequenceSimilarityViewer} from './analysis/sequence-similarity-viewer';
import {SequenceDiversityViewer} from './analysis/sequence-diversity-viewer';
import {invalidateMols, MONOMERIC_COL_TAGS, substructureSearchDialog} from './substructure-search/substructure-search';
import {saveAsFastaUI} from './utils/save-as-fasta';
import {BioSubstructureFilter} from './widgets/bio-substructure-filter';
import { getMonomericMols } from './calculations/monomerLevelMols';
import { delay } from '@datagrok-libraries/utils/src/test';

//tags: init
export async function initBio() {
}

//name: fastaSequenceCellRenderer
//tags: cellRenderer
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=fasta
//output: grid_cell_renderer result
export function fastaSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
  return new MacromoleculeSequenceCellRenderer();
}

//name: Sequence Renderer
//input: column molColumn {semType: Macromolecule}
//tags: panel
//output: widget result
export function macroMolColumnPropertyPanel(molColumn: DG.Column): DG.Widget {
  return getMacroMolColumnPropertyPanel(molColumn);
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


function checkInputColumnUi(
  col: DG.Column, name: string, allowedNotations: string[] = [], allowedAlphabets: string[] = []
): boolean {
  const [res, msg]: [boolean, string] = checkInputColumn(col, name, allowedNotations, allowedAlphabets);
  if (!res)
    grok.shell.warning(msg);
  return res;
}

export function checkInputColumn(
  col: DG.Column, name: string, allowedNotations: string[] = [], allowedAlphabets: string[] = []
): [boolean, string] {
  let res: boolean = true;
  let msg: string = '';

  const uh = new bio.UnitsHandler(col);
  if (col.semType !== DG.SEMTYPE.MACROMOLECULE) {
    grok.shell.warning(name + ' analysis is allowed for Macromolecules semantic type');
    res = false;
  } else {
    const notation: string = uh.notation;
    if (allowedNotations.length > 0 &&
      !allowedNotations.some((n) => notation.toUpperCase() == (n.toUpperCase()))
    ) {
      const notationAdd = allowedNotations.length == 0 ? 'any notation' :
        (`notation${allowedNotations.length > 1 ? 's' : ''} ${allowedNotations.map((n) => `"${n}"`).join(', ')} `);
      msg = `${name} + ' analysis is allowed for Macromolecules with notation ${notationAdd}.`;
      res = false;
    } else if (!uh.isHelm()) {
      // alphabet is not specified for 'helm' notation
      const alphabet: string = uh.alphabet;
      if (
        allowedAlphabets.length > 0 &&
        !allowedAlphabets.some((a) => alphabet.toUpperCase() == (a.toUpperCase()))
      ) {
        const alphabetAdd = allowedAlphabets.length == 0 ? 'any alphabet' :
          (`alphabet${allowedAlphabets.length > 1 ? 's' : ''} ${allowedAlphabets.map((a) => `"${a}"`).join(', ')}.`);
        msg = `${name} + ' analysis is allowed for Macromolecules with alphabet ${alphabetAdd}.`;
        res = false;
      }
    }
  }

  return [res, msg];
}

//name: sequenceAlignment
//input: string alignType {choices: ['Local alignment', 'Global alignment']}
// eslint-disable-next-line max-len
//input: string alignTable {choices: ['AUTO', 'NUCLEOTIDES', 'BLOSUM45', 'BLOSUM50','BLOSUM62','BLOSUM80','BLOSUM90','PAM30','PAM70','PAM250','SCHNEIDER','TRANS']}
//input: double gap
//input: string seq1
//input: string seq2
//output: object res
export function sequenceAlignment(alignType: string, alignTable: string, gap: number, seq1: string, seq2: string) {
  const toAlign = new SequenceAlignment(seq1, seq2, gap, alignTable);
  const res = alignType == 'Local alignment' ? toAlign.smithWaterman() : toAlign.needlemanWunch();
  return res;
}

//name: WebLogo
//description: WebLogo viewer
//tags: viewer, panel
//output: viewer result
export function webLogoViewer() {
  return new bio.WebLogoViewer();
}

//name: VdRegions
//description: V-Domain regions viewer
//tags: viewer, panel
//output: viewer result
export function vdRegionViewer() {
  return new VdRegionsViewer();
}

//top-menu: Bio | Sequence Activity Cliffs...
//name: Sequence Activity Cliffs
//description: detect activity cliffs
//input: dataframe table [Input data table]
//input: column macroMolecule {semType: Macromolecule}
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
export async function activityCliffs(df: DG.DataFrame, macroMolecule: DG.Column, activities: DG.Column,
  similarity: number, methodName: string): Promise<DG.Viewer | undefined> {
  if (!checkInputColumnUi(macroMolecule, 'Activity Cliffs'))
    return;
  const axesNames = getEmbeddingColsNames(df);
  const options = {
    'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005},
  };
  const tags = {
    'units': macroMolecule.getTag(DG.TAGS.UNITS),
    'aligned': macroMolecule.getTag(TAGS.aligned),
    'separator': macroMolecule.getTag(TAGS.separator),
    'alphabet': macroMolecule.getTag(TAGS.alphabet),
  };
  const sp = await getActivityCliffs(
    df,
    macroMolecule,
    null,
    axesNames,
    'Activity cliffs',
    activities,
    similarity,
    'Tanimoto',
    methodName,
    DG.SEMTYPE.MACROMOLECULE,
    tags,
    sequenceSpaceByFingerprints,
    getChemSimilaritiesMarix,
    createTooltipElement,
    createPropPanelElement,
    createLinesGrid,
    (options as any)[methodName]);
  return sp;
}

//top-menu: Bio | Sequence Space...
//name: Sequence Space
//input: dataframe table
//input: column macroMolecule { semType: Macromolecule }
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
//input: string similarityMetric { choices:["Tanimoto", "Asymmetric", "Cosine", "Sokal"] }
//input: bool plotEmbeddings = true
export async function sequenceSpaceTopMenu(table: DG.DataFrame, macroMolecule: DG.Column, methodName: string,
  similarityMetric: string = 'Tanimoto', plotEmbeddings: boolean): Promise<DG.Viewer | undefined> {
  //delay is required for initial function dialog to close before starting invalidating of molfiles. Otherwise dialog is freezing
  await delay(10);
  if (!checkInputColumnUi(macroMolecule, 'Sequence space'))
    return;

  const embedColsNames = getEmbeddingColsNames(table);
  const withoutEmptyValues = DG.DataFrame.fromColumns([macroMolecule]).clone();
  const emptyValsIdxs = removeEmptyStringRows(withoutEmptyValues, macroMolecule);

  const chemSpaceParams = {
    seqCol: withoutEmptyValues.col(macroMolecule.name)!,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: embedColsNames
  };
  const sequenceSpaceRes = await sequenceSpaceByFingerprints(chemSpaceParams);
  const embeddings = sequenceSpaceRes.coordinates;
  for (const col of embeddings) {
    const listValues = col.toList();
    emptyValsIdxs.forEach((ind: number) => listValues.splice(ind, 0, null));
    table.columns.add(DG.Column.fromList('double', col.name, listValues));
  }
  if (plotEmbeddings) {
    return grok.shell
      .tableView(table.name)
      .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Sequence space'});
  };

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
};

//top-menu: Bio | To Atomic Level...
//name: To Atomic Level
//description: returns molfiles for each monomer from HELM library
//input: dataframe df [Input data table]
//input: column macroMolecule {semType: Macromolecule}
export async function toAtomicLevel(df: DG.DataFrame, macroMolecule: DG.Column): Promise<void> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }
  if (!checkInputColumnUi(macroMolecule, 'To Atomic Level'))
    return;
  const monomersLibFile = await _package.files.readAsText(HELM_CORE_LIB_FILENAME);
  const monomersLibObject: any[] = JSON.parse(monomersLibFile);
  _toAtomicLevel(df, macroMolecule, monomersLibObject);
}

//top-menu: Bio | MSA...
//name: MSA
//input: dataframe table
//input: column sequence { semType: Macromolecule, units: ['fasta'], alphabet: ['DNA', 'RNA', 'PT'] }
//output: column result
export async function multipleSequenceAlignmentAny(table: DG.DataFrame, sequence: DG.Column): Promise<DG.Column | null> {
  const func: DG.Func = DG.Func.find({package: 'Bio', name: 'multipleSequenceAlignmentAny'})[0];

  if (!checkInputColumnUi(sequence, 'MSA', ['fasta'], ['DNA', 'RNA', 'PT']))
    return null;

  const unUsedName = table.columns.getUnusedName(`msa(${sequence.name})`);
  const msaCol = await runKalign(sequence, false, unUsedName);
  table.columns.add(msaCol);

  // This call is required to enable cell renderer activation
  await grok.data.detectSemanticTypes(table);

  // const tv: DG.TableView = grok.shell.tv;
  // tv.grid.invalidate();
  return msaCol;
}

//name: Bio | MSA
//tags: bio, panel
//input: column sequence { semType: Macromolecule }
//output: column result
export async function panelMSA(sequence: DG.Column): Promise<DG.Column | null> {
  return multipleSequenceAlignmentAny(sequence.dataFrame, sequence);
}

//name: Composition Analysis
//top-menu: Bio | Composition Analysis
//output: viewer result
export async function compositionAnalysis(): Promise<void> {
  // Higher priority for columns with MSA data to show with WebLogo.
  const tv = grok.shell.tv;
  const df = tv.dataFrame;
  //@ts-ignore
  const colList: DG.Column[] = df.columns.toList().filter((col) => {
    if (col.semType != DG.SEMTYPE.MACROMOLECULE)
      return false;

    const colUH = new bio.UnitsHandler(col);
    // TODO: prevent for cyclic, branched or multiple chains in Helm
    return true;
  });

  const handler = async (col: DG.Column) => {
    if (!checkInputColumnUi(col, 'Composition'))
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
    const selectedCol = colList.find((c) => { return (new bio.UnitsHandler(c)).isMsa(); });
    const colInput: DG.InputBase = ui.choiceInput(
      'Column', selectedCol ? selectedCol.name : colListNames[0], colListNames);
    ui.dialog({
      title: 'Composition Analysis',
      helpUrl: '/help/domains/bio/macromolecules.md#composition-analysis'
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

//top-menu: Bio | Sdf to Json lib...
//name: sdfToJsonLib
//input: dataframe table
export async function sdfToJsonLib(table: DG.DataFrame) {
  const jsonMonomerLibrary = createJsonMonomerLibFromSdf(table);
}

//name: Representations
//tags: panel, widgets
//input: cell macroMolecule {semType: Macromolecule}
//output: widget result
export async function peptideMolecule(macroMolecule: DG.Cell): Promise<DG.Widget> {
  const monomersLibFile = await _package.files.readAsText(HELM_CORE_LIB_FILENAME);
  const monomersLibObject: any[] = JSON.parse(monomersLibFile);

  return representationsWidget(macroMolecule, monomersLibObject);
}

//name: importFasta
//description: Opens FASTA file
//tags: file-handler
//meta.ext: fasta, fna, ffn, faa, frn, fa, fst
//input: string fileContent
//output: list tables
export function importFasta(fileContent: string): DG.DataFrame [] {
  const ffh = new bio.FastaFileHandler(fileContent);
  return ffh.importFasta();
}

//name: Bio | Convert ...
//friendly-name: Bio | Convert
//tags: panel, bio
//input: column col {semType: Macromolecule}
export function convertPanel(col: DG.Column): void {
  convert(col);
}

//name: monomerCellRenderer
//tags: cellRenderer
//meta.cellType: Monomer
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
            message: `units: ${col.getTag(DG.TAGS.UNITS)}`
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

//name: Bio | Split to monomers
//tags: panel, bio
//input: column col {semType: Macromolecule}
export function splitToMonomers(col: DG.Column<string>): void {
  if (!col.getTag(bio.TAGS.aligned).includes(C.MSA))
    return grok.shell.error('Splitting is applicable only for aligned sequences');

  const tempDf = splitAlignedSequences(col);
  const originalDf = col.dataFrame;
  for (const tempCol of tempDf.columns) {
    const newCol = originalDf.columns.add(tempCol);
    newCol.semType = C.SEM_TYPES.MONOMER;
    newCol.setTag(DG.TAGS.CELL_RENDERER, C.SEM_TYPES.MONOMER);
    newCol.setTag(C.TAGS.ALPHABET, col.getTag(C.TAGS.ALPHABET));
  }
  grok.shell.tv.grid.invalidate();
}

//name: Bio: getHelmMonomers
//input: column sequence {semType: Macromolecule}
export function getHelmMonomers(sequence: DG.Column<string>): string[] {
  const stats = bio.getStats(sequence, 1, bio.splitterAsHelm);
  return Object.keys(stats.freq);
}


//name: SequenceSimilaritySearchViewer
//tags: viewer
//output: viewer result
export function similaritySearchViewer(): SequenceSimilarityViewer {
  return new SequenceSimilarityViewer();
}

//top-menu: Bio | Similarity Search...
//name: similaritySearch
//description: finds the most similar sequence
//output: viewer result
export function similaritySearchTopMenu(): void {
  const view = (grok.shell.v as DG.TableView);
  const viewer = view.addViewer('SequenceSimilaritySearchViewer');
  view.dockManager.dock(viewer, 'down');
}

//name: SequenceDiversitySearchViewer
//tags: viewer
//output: viewer result
export function diversitySearchViewer(): SequenceDiversityViewer {
  return new SequenceDiversityViewer();
}

//top-menu: Bio | Diversity Search...
//name: diversitySearch
//description: finds the most diverse molecules
//output: viewer result
export function diversitySearchTopMenu() {
  const view = (grok.shell.v as DG.TableView);
  const viewer = view.addViewer('SequenceDiversitySearchViewer');
  view.dockManager.dock(viewer, 'down');
}

//name: Bio | Substructure search ...
//tags: panel, bio
//input: column col {semType: Macromolecule}
export function bioSubstructureSearch(col: DG.Column): void {
  substructureSearchDialog(col);
}

//name: saveAsFasta
//description: As FASTA...
//tags: fileExporter
export function saveAsFasta() {
  saveAsFastaUI();
}

//name: BioSubstructureFilter
//description: Substructure filter for macromolecules
//tags: filter
//output: filter result
//meta.semType: Macromolecule
export function bioSubstructureFilter(): BioSubstructureFilter {
  return new BioSubstructureFilter();
}
