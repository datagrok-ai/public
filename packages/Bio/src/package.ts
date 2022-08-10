/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

import {MacromoleculeDifferenceCellRenderer, MonomerCellRenderer} from './utils/cell-renderer';
import {WebLogo, SeqColStats} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {runKalign, testMSAEnoughMemory} from './utils/multiple-sequence-alignment';
import {SequenceAlignment, Aligned} from './seq_align';
import {Nucleotides} from '@datagrok-libraries/bio/src/nucleotides';
import {Aminoacids} from '@datagrok-libraries/bio/src/aminoacids';
import {getEmbeddingColsNames, sequenceSpace} from './utils/sequence-space';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {sequenceGetSimilarities, drawTooltip} from './utils/sequence-activity-cliffs';
import {createJsonMonomerLibFromSdf, encodeMonomers, getMolfilesFromSeq, HELM_CORE_LIB_FILENAME} from './utils/utils';
import {getMacroMol} from './utils/atomic-works';
import {MacromoleculeSequenceCellRenderer} from './utils/cell-renderer';
import {convert} from './utils/convert';
import {representationsWidget} from './widgets/representations';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {FastaFileHandler} from '@datagrok-libraries/bio/src/utils/fasta-handler';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils'


//name: fastaSequenceCellRenderer
//tags: cellRenderer
//meta.cellType: Sequence
//meta.columnTags: units=fasta
//output: grid_cell_renderer result
export function fastaSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
  return new MacromoleculeSequenceCellRenderer();
}

//name: separatorSequenceCellRenderer
//tags: cellRenderer
//meta.cellType: Sequence
//meta.columnTags: units=separator
//output: grid_cell_renderer result
export function separatorSequenceCellRenderer(): MacromoleculeSequenceCellRenderer {
  return new MacromoleculeSequenceCellRenderer();
}

function checkInputColumn(col: DG.Column, name: string,
  allowedNotations: string[] = [], allowedAlphabets: string[] = []): boolean {
  const notation: string = col.getTag(DG.TAGS.UNITS);
  const alphabet: string = col.getTag('alphabet')
  if (col.semType !== DG.SEMTYPE.MACROMOLECULE) {
    grok.shell.warning(name + ' analysis is allowed for Macromolecules semantic type');
    return false;
  } else if (
    (allowedAlphabets.length > 0 &&
      !allowedAlphabets.some((a) => alphabet.toUpperCase() == (a.toUpperCase()))) ||
    (allowedNotations.length > 0 &&
      !allowedNotations.some((n) => notation.toUpperCase() == (n.toUpperCase())))
  ) {
    const notationAdd = allowedNotations.length == 0 ? 'any notation' :
      (`notation${allowedNotations.length > 1 ? 's' : ''} ${allowedNotations.map((n) => `"${n}"`).join(', ')} `);
    const alphabetAdd = allowedNotations.length == 0 ? 'any alphabet' :
      (`alphabet${allowedAlphabets.length > 1 ? 's' : ''} ${allowedAlphabets.map((a) => `"${a}"`).join(', ')}.`);

    grok.shell.warning(name + ' analysis is allowed for Macromolecules with ' + notationAdd + ' and ' + alphabetAdd);
    return false;
  }

  return true;
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
  return new WebLogo();
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
  similarity: number, methodName: string): Promise<void> {
  if (!checkInputColumn(macroMolecule, 'Activity Cliffs'))
    return;
  const encodedCol = encodeMonomers(macroMolecule);
  if (!encodedCol)
    return;
  const axesNames = getEmbeddingColsNames(df);
  const options = {
    'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005},
  };
  const units = macroMolecule!.tags[DG.TAGS.UNITS];
  await getActivityCliffs(
    df,
    macroMolecule,
    encodedCol,
    axesNames,
    'Activity cliffs',
    activities,
    similarity,
    'Levenshtein',
    methodName,
    DG.SEMTYPE.MACROMOLECULE,
    units,
    sequenceSpace,
    sequenceGetSimilarities,
    drawTooltip,
    (options as any)[methodName]);
}

//top-menu: Bio | Sequence Space...
//name: Sequence Space
//input: dataframe table
//input: column macroMolecule { semType: Macromolecule }
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
//input: string similarityMetric { choices:["Levenshtein", "Tanimoto"] }
//input: bool plotEmbeddings = true
export async function sequenceSpaceTopMenu(table: DG.DataFrame, macroMolecule: DG.Column, methodName: string,
  similarityMetric: string = 'Levenshtein', plotEmbeddings: boolean): Promise<void> {
  if (!checkInputColumn(macroMolecule, 'Activity Cliffs'))
    return;
  const encodedCol = encodeMonomers(macroMolecule);
  if (!encodedCol)
    return;
  const embedColsNames = getEmbeddingColsNames(table);
  const withoutEmptyValues = DG.DataFrame.fromColumns([macroMolecule]).clone();
  const emptyValsIdxs = removeEmptyStringRows(withoutEmptyValues, encodedCol);

  const chemSpaceParams = {
    seqCol: withoutEmptyValues.col(macroMolecule.name)!,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: embedColsNames
  };
  const sequenceSpaceRes = await sequenceSpace(chemSpaceParams);
  const embeddings = sequenceSpaceRes.coordinates;
  for (const col of embeddings) {
      const listValues = col.toList();
      emptyValsIdxs.forEach((ind: number) => listValues.splice(ind, 0, null));
      table.columns.add(DG.Column.fromFloat32Array(col.name, listValues));
  }   
  if (plotEmbeddings) {
    for (const v of grok.shell.views) {
      if (v.name === table.name)
        (v as DG.TableView).scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Sequence space'});
    }
  }
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
  if (!checkInputColumn(macroMolecule, 'To Atomic Level'))
    return;

  let currentView: DG.TableView;
  for (const view of grok.shell.tableViews) {
    if (df.name === view.name)
      currentView = view;
  }

  // Some hack to activate Chem Molecule rendering
  const file2 = await _package.files.readAsText('tests/sar-small.csv');
  const df2 = DG.DataFrame.fromCsv(file2);
  const v2 = grok.shell.addTableView(df2);
  setTimeout(() => {
    grok.shell.closeTable(df2);
    v2.close();
    grok.shell.v = currentView;
  }, 100);

  const monomersLibFile = await _package.files.readAsText(HELM_CORE_LIB_FILENAME);
  const monomersLibObject: any[] = JSON.parse(monomersLibFile);
  const atomicCodes = getMolfilesFromSeq(macroMolecule, monomersLibObject);
  const result = await getMacroMol(atomicCodes!);

  const col = DG.Column.fromStrings('regenerated', result);
  col.semType = DG.SEMTYPE.MOLECULE;
  col.tags[DG.TAGS.UNITS] = 'molblock';
  df.columns.add(col, true);
}


//top-menu: Bio | MSA...
//name: MSA
//input: dataframe table
//input: column sequence { semType: Macromolecule }
//output: column result
export async function multipleSequenceAlignmentAny(table: DG.DataFrame, col: DG.Column): Promise<DG.Column | null> {
  if (!checkInputColumn(col, 'MSA', ['fasta'], ['DNA', 'RNA', 'PT']))
    return null;

  const unUsedName = table.columns.getUnusedName(`msa(${col.name})`);
  const msaCol = await runKalign(col, false, unUsedName);
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
export async function panelMSA(col: DG.Column): Promise<DG.Column | null> {
  return multipleSequenceAlignmentAny(col.dataFrame, col);
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

    const colUH = new UnitsHandler(col);
    // TODO: prevent for cyclic, branched or multiple chains in Helm
    return true;
  });

  const handler = async (col: DG.Column) => {
    if (!checkInputColumn(col, 'Composition'))
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
    const colInput: DG.InputBase = ui.choiceInput('Column', colListNames[0], colListNames);
    ui.dialog({
      title: 'R-Groups Analysis',
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
  const ffh = new FastaFileHandler(fileContent);
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

//name: MacromoleculeDifferenceCellRenderer
//tags: cellRenderer
//meta.cellType: MacromoleculeDifference
//output: grid_cell_renderer result
export function macromoleculeDifferenceCellRenderer(): MacromoleculeDifferenceCellRenderer {
  return new MacromoleculeDifferenceCellRenderer();
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
          //   'units: "' + col.getTag('units') + '"');

          res.push({
            file: fileInfo.path, result: 'detected', column: col.name,
            message: `units: ${col.getTag('units')}`
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
