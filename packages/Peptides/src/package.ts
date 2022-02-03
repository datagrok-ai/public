/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  AlignedSequenceCellRenderer,
  AlignedSequenceDifferenceCellRenderer,
  AminoAcidsCellRenderer,
} from './utils/cell-renderer';
import {Logo} from './viewers/logo-viewer';
import {StackedBarChart} from './viewers/stacked-barchart-viewer';

import {analyzePeptidesWidget} from './widgets/analyze-peptides';
import {PeptideSimilaritySpaceWidget} from './utils/peptide-similarity-space';
import {manualAlignmentWidget} from './widgets/manual-alignment';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {peptideMoleculeWidget, getMolecule} from './widgets/peptide-molecule';
import {SubstViewer} from './viewers/subst-viewer';
import {runKalign, testMSAEnoughMemory} from './utils/multiple-sequence-alignment';

export const _package = new DG.Package();
let tableGrid: DG.Grid;
let currentDf: DG.DataFrame;
let alignedSequenceCol: DG.Column;
let view: DG.TableView;

async function main(chosenFile: string) {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  const path = _package.webRoot + 'files/' + chosenFile;
  const peptides = (await grok.data.loadTable(path));
  peptides.name = 'Peptides';
  peptides.setTag('dataType', 'peptides');
  const view = grok.shell.addTableView(peptides);
  tableGrid = view.grid;
  view.name = 'PeptidesView';
  grok.shell.windows.showProperties = true;

  pi.close();
}

//name: Peptides
//tags: app
export async function Peptides() {
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
  if (col.getTag('isAnalysisApplicable') === 'false')
    return new DG.Widget(ui.divText('Analysis is not applicable'));

  view = (grok.shell.v as DG.TableView);
  tableGrid = view.grid;
  currentDf = col.dataFrame;
  alignedSequenceCol = col;
  return await analyzePeptidesWidget(col, view, tableGrid, currentDf);
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

//name: substitution-analysis-viewer
//description: Substitution Analysis Viewer
//tags: viewer
//output: viewer result
export function subst(): SubstViewer {
  return new SubstViewer();
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
  return await peptideMoleculeWidget(peptide);
}

//name: Peptide Molecule
//tags: panel, widgets
//input: string aar {semType: aminoAcids}
//output: widget result
export async function peptideMolecule2(aar: string): Promise<DG.Widget> {
  const peptide = alignedSequenceCol.get(currentDf.currentRowIdx);
  return await peptideMolecule(peptide);
}

//name: StackedBarChartAA
//tags: viewer
//output: viewer result
export function stackedBarChart(): DG.JsViewer {
  return new StackedBarChart();
}

//name: alignedSequenceCellRenderer
//tags: cellRenderer, cellRenderer-alignedSequence
//meta-cell-renderer-sem-type: alignedSequence
//output: grid_cell_renderer result
export function alignedSequenceCellRenderer() {
  return new AlignedSequenceCellRenderer();
}

//name: aminoAcidsCellRenderer
//tags: cellRenderer, cellRenderer-aminoAcids
//meta-cell-renderer-sem-type: aminoAcids
//output: grid_cell_renderer result
export function aminoAcidsCellRenderer() {
  return new AminoAcidsCellRenderer();
}

//name: peptide-logo-viewer
//tags: viewer, panel
//output: viewer result
export function logov() {
  return new Logo();
}

//name: Manual Alignment
//tags: panel, widgets
//input: string monomer {semType: aminoAcids}
//output: widget result
export function manualAlignment(monomer: string) {
  //TODO: recalculate Molfile and Molecule panels on sequence update
  return manualAlignmentWidget(alignedSequenceCol, currentDf);
}

//name: Peptide Space
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export async function peptideSpacePanel(col: DG.Column): Promise<DG.Widget> {
  const widget = new PeptideSimilaritySpaceWidget(col, view ?? grok.shell.v);
  return await widget.draw();
}

//name: Molfile
//tags: panel, widgets
//input: string peptide { semType: alignedSequence }
//output: widget result
export async function peptideMolfile(peptide: string): Promise<DG.Widget> {
  const smiles = getMolecule(peptide);
  return await grok.functions.call('Chem:molfile', {'smiles': smiles});
}

//name: Molfile
//tags: panel, widgets
//input: string aar { semType: aminoAcids }
//output: widget result
export async function peptideMolfile2(aar: string): Promise<DG.Widget> {
  const peptide = alignedSequenceCol.get(currentDf.currentRowIdx);
  return await peptideMolfile(peptide);
}

//name: Multiple sequence alignment
//tags: panel
//input: column col {semType: alignedSequence}
//output: dataframe result
export async function multipleSequenceAlignment(col: DG.Column): Promise<DG.DataFrame> {
  const msaCol = await runKalign(col, true);
  const table = col.dataFrame;
  table.columns.add(msaCol);
  return table;
}

//name: Multiple sequence alignment for any column
//tags: panel
//input: column col
//output: dataframe result
export async function multipleSequenceAlignmentAny(col: DG.Column): Promise<DG.DataFrame> {
  const msaCol = await runKalign(col, false);
  const table = col.dataFrame;
  table.columns.add(msaCol);
  return table;
}

//name: Test multiple sequence alignment for any column
export async function runTestMSAEnoughMemory(col: DG.Column) {
  await testMSAEnoughMemory(col);
  return col;
}

//name: Substitution
//tags: panel, widgets
//input: dataframe table {semType: Substitution}
//output: widget result
export async function peptideSubstitution(table: DG.DataFrame): Promise<DG.Widget> {
  if (!table) {
    return new DG.Widget(ui.label('No substitution'));
  }
  const peptideLength = 17;
  const initialCol: DG.Column = table.columns.byName('Initial');
  const substitutedCol: DG.Column = table.columns.byName('Substituted');
  const substCounts = [];
  let cnt = 0;

  for (let i = 0; i < initialCol.length; ++i) {
    const initialPeptide: string = initialCol.get(i);
    const substPeptide: string = substitutedCol.get(i);
    const concat = initialPeptide + '#' + substPeptide;

    initialCol.set(i, concat);

    const initialAminos = initialPeptide.split('-');
    const substAminos = substPeptide.split('-');

    for (let j = 0; j < peptideLength; ++j) {
      if (initialAminos[j] != substAminos[j])
        substCounts[cnt++] = j;
    }
  }

  const countCol = DG.Column.fromInt32Array('substCounts', substCounts as unknown as Int32Array);
  const df = DG.DataFrame.fromColumns([countCol]);
  const barchart = df.plot.histogram({value: 'substCounts'});
  if (barchart) {
    barchart.root.style.width = '200px';
    barchart.root.style.marginLeft = '30px';
  }

  initialCol.semType = 'alignedSequenceDifference';
  initialCol.name = 'Substitution';
  table.columns.remove('Substituted');
  return new DG.Widget(ui.div([table.plot.grid().root]));
}

//name: alignedSequenceDifferenceCellRenderer
//tags: cellRenderer, cellRenderer-alignedSequenceDifference
//meta-cell-renderer-sem-type: alignedSequenceDifference
//output: grid_cell_renderer result
export function alignedSequenceDifferenceCellRenderer() {
  return new AlignedSequenceDifferenceCellRenderer();
}
