/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as C from './utils/constants';

import {
  AlignedSequenceCellRenderer,
  AlignedSequenceDifferenceCellRenderer,
  AminoAcidsCellRenderer,
} from './utils/cell-renderer';
import {StackedBarChart} from './viewers/stacked-barchart-viewer';

import {analyzePeptidesWidget} from './widgets/analyze-peptides';
import {PeptideSimilaritySpaceWidget} from './utils/peptide-similarity-space';
import {manualAlignmentWidget} from './widgets/manual-alignment';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {peptideMoleculeWidget, getMolecule} from './widgets/peptide-molecule';
import {runKalign, testMSAEnoughMemory} from './utils/multiple-sequence-alignment';
import {msaWidget} from './widgets/multiple-sequence-alignment';
import {PeptideSpaceViewer} from './viewers/peptide-space-viewer';

export const _package = new DG.Package();
let currentTable: DG.DataFrame;
let alignedSequenceColumn: DG.Column;

async function main(chosenFile: string): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  const path = _package.webRoot + 'files/' + chosenFile;
  const peptides = await grok.data.loadTable(path);
  peptides.name = 'Peptides';
  const view = grok.shell.addTableView(peptides);
  view.name = 'PeptidesView';
  grok.shell.windows.showProperties = true;

  pi.close();
}

//name: Peptides
//tags: app
export async function Peptides(): Promise<void> {
  const wikiLink = ui.link('wiki', 'https://github.com/datagrok-ai/public/blob/master/help/domains/bio/peptides.md');
  const textLink = ui.inlineText(['For more details, see our ', wikiLink, '.']);

  const appDescription = ui.info(
    [
      ui.list([
        '- automatic recognition of peptide sequences',
        '- native integration with tons of Datagrok out-of-the box features (visualization, filtering, clustering, ' +
        'multivariate analysis, etc)',
        '- custom rendering in the spreadsheet',
        '- interactive logo plots',
        '- rendering residues',
        '- structure-activity relationship:',
        ' ',
        'a) highlighting statistically significant changes in activity in the [position, monomer] spreadsheet',
        'b) for the specific [position, monomer], visualizing changes of activity distribution (specific monomer in ' +
        'this position vs rest of the monomers in this position)',
        'c) interactivity',
      ]),
    ],
    'Use and analyse peptide sequence data to support your research:',
  );

  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  grok.shell.newView('Peptides', [
    appDescription,
    ui.info([textLink]),
    ui.divH([
      ui.button('Simple demo', () => main('aligned.csv'), ''),
      ui.button('Complex demo', () => main('aligned_2.csv'), ''),

    ]),
  ]);
}

//name: Peptides
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export async function peptidesPanel(col: DG.Column): Promise<DG.Widget> {
  if (!(col.temp['isAnalysisApplicable'] ?? true))
    return new DG.Widget(ui.divText('Analysis is not applicable'));

  [currentTable, alignedSequenceColumn] = getOrDefine(col.dataFrame, col);
  return analyzePeptidesWidget(currentTable, alignedSequenceColumn);
}

//name: peptide-sar-viewer
//description: Peptides SAR Viewer
//tags: viewer
//output: viewer result
export function sar(): SARViewer {
  return new SARViewer();
}

//name: peptide-sar-viewer-vertical
//description: Peptides Vertical SAR Viewer
//tags: viewer
//output: viewer result
export function sarVertical(): SARViewerVertical {
  return new SARViewerVertical();
}

//name: peptide-space-viewer
//description: Peptide Space Viewer
//tags: viewer
//output: viewer result
export function peptideSpace(): PeptideSpaceViewer {
  return new PeptideSpaceViewer();
}

//name: StackedBarchart Widget
//tags: panel, widgets
//input: column col {semType: aminoAcids}
//output: widget result
export async function stackedBarchartWidget(col: DG.Column): Promise<DG.Widget> {
  const viewer = await col.dataFrame.plot.fromType('StackedBarChartAA');
  const panel = ui.divH([viewer.root]);
  return new DG.Widget(panel);
}

//name: Peptide Molecule
//tags: panel, widgets
//input: string peptide {semType: alignedSequence}
//output: widget result
export async function peptideMolecule(peptide: string): Promise<DG.Widget> {
  [currentTable, alignedSequenceColumn] = getOrDefine();
  return peptideMoleculeWidget(peptide, currentTable);
}

//name: Peptide Molecule
//tags: panel, widgets
//input: string _aar {semType: aminoAcids}
//output: widget result
export async function peptideMolecule2(_aar: string): Promise<DG.Widget> {
  [currentTable, alignedSequenceColumn] = getOrDefine();
  const peptide = alignedSequenceColumn.get(currentTable.currentRowIdx);
  return peptideMolecule(peptide);
}

//name: StackedBarChartAA
//tags: viewer
//output: viewer result
export function stackedBarChart(): DG.JsViewer {
  return new StackedBarChart();
}

//name: alignedSequenceCellRenderer
//tags: cellRenderer
//meta.cellType: alignedSequence
//output: grid_cell_renderer result
export function alignedSequenceCellRenderer(): AlignedSequenceCellRenderer {
  return new AlignedSequenceCellRenderer();
}

//name: aminoAcidsCellRenderer
//tags: cellRenderer
//meta.cellType: aminoAcids
//output: grid_cell_renderer result
export function aminoAcidsCellRenderer(): AminoAcidsCellRenderer {
  return new AminoAcidsCellRenderer();
}

//name: Manual Alignment
//tags: panel, widgets
//input: string _monomer {semType: aminoAcids}
//output: widget result
export function manualAlignment(_monomer: string): DG.Widget {
  [currentTable, alignedSequenceColumn] = getOrDefine();
  //TODO: recalculate Molfile and Molecule panels on sequence update
  return manualAlignmentWidget(alignedSequenceColumn, currentTable);
}

//name: Peptide Space
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export async function peptideSpacePanel(col: DG.Column): Promise<DG.Widget> {
  [currentTable, alignedSequenceColumn] = getOrDefine(col.dataFrame, col);
  const widget = new PeptideSimilaritySpaceWidget(col, grok.shell.v as DG.TableView);
  return widget.draw();
}

//name: Molfile
//tags: panel, widgets
//input: string peptide { semType: alignedSequence }
//output: widget result
export async function peptideMolfile(peptide: string): Promise<DG.Widget> {
  [currentTable, alignedSequenceColumn] = getOrDefine();
  const smiles = getMolecule(peptide, alignedSequenceColumn.tags[C.TAGS.SEPARATOR]);
  return grok.functions.call('Chem:molfile', {'smiles': smiles}) as Promise<DG.Widget>;
}

//name: Molfile
//tags: panel, widgets
//input: string _aar { semType: aminoAcids }
//output: widget result
export async function peptideMolfile2(_aar: string): Promise<DG.Widget> {
  [currentTable, alignedSequenceColumn] = getOrDefine();
  const peptide = alignedSequenceColumn.get(currentTable.currentRowIdx);
  return peptideMolfile(peptide);
}

//name: Multiple sequence alignment
//tags: panel
//input: column col {semType: alignedSequence}
//output: dataframe result
export async function multipleSequenceAlignment(col: DG.Column): Promise<DG.DataFrame> {
  return msaWidget(col);
}

//name: Multiple sequence alignment for any column
//input: dataframe table
//input: column col
//output: dataframe result
export async function multipleSequenceAlignmentAny(table: DG.DataFrame, col: DG.Column): Promise<DG.DataFrame> {
  const msaCol = await runKalign(col, false);
  table.columns.add(msaCol);
  return table;
}

//name: Test multiple sequence alignment for any column
//input: dataframe _table
//input: column col
//output: column result
export async function runTestMSAEnoughMemory(_table: DG.DataFrame, col: DG.Column<string>): Promise<DG.Column<string>> {
  await testMSAEnoughMemory(col);
  return col;
}

//name: Get Peptides Structure
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export function getPeptidesStructure(col: DG.Column): DG.Widget {
  const getButtonTooltip = 'Retrieves peptides structure from customer database by special id column';
  const getButton = ui.button('Get structure', async () => {
    const progress = DG.TaskBarProgressIndicator.create('Getting structure...');
    try {
      const params = {peptidesTable: col.dataFrame};
      const result = await grok.functions.call('Customerextensions:getPeptidesStructure', params);
      const text = result ? 'Structure retreived' : 'Structure retreivial is not possible';
      grok.shell.info(text);
    } catch (e) {
      console.warn(e);
    } finally {
      progress.close();
    }
  }, getButtonTooltip);
  return new DG.Widget(getButton);
}

//name: alignedSequenceDifferenceCellRenderer
//tags: cellRenderer
//meta.cellType: alignedSequenceDifference
//output: grid_cell_renderer result
export function alignedSequenceDifferenceCellRenderer(): AlignedSequenceDifferenceCellRenderer {
  return new AlignedSequenceDifferenceCellRenderer();
}

function getOrDefine(dataframe?: DG.DataFrame, column?: DG.Column | null): [DG.DataFrame, DG.Column] {
  dataframe ??= grok.shell.t;
  column ??= dataframe.columns.bySemType(C.SEM_TYPES.ALIGNED_SEQUENCE);
  if (column === null)
    throw new Error('Table does not contain aligned sequence columns');

  return [dataframe, column];
}
