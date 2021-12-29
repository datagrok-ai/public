/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  AlignedSequenceCellRenderer,
  AminoAcidsCellRenderer,
} from './utils/cell-renderer';
import {Logo} from './viewers/logo-viewer';
import {StackedBarChart} from './viewers/stacked-barchart-viewer';

import {analyzePeptidesWidget} from './widgets/analyze-peptides';
import {PeptideSimilaritySpaceWidget} from './utils/peptide-similarity-space';
import {manualAlignmentWidget} from './widgets/manual-alignment';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {peptideMoleculeWidget} from './widgets/peptide-molecule';
import {SpiralPlot} from './viewers/spiral-plot';
import {SubstViewer} from './viewers/subst-viewer';

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

//name: Peptides App
//tags: app
export function Peptides() {
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
      ui.button('Open peptide sequences demonstration set', () => main('aligned.csv'), ''),
      ui.button('Open complex case demo', () => main('aligned_2.csv'), ''),
    ]),
  ]);
}

//name: Peptides
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export async function peptidesPanel(col: DG.Column): Promise<DG.Widget> {
  if (col.getTag('isAnalysisApplicable') === 'false') {
    return new DG.Widget(ui.divText('Analysis is not applicable'));
  }
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

//name: Spiral Plot
////input: dataframe table
////input: column activity
//tags: viewer, panel
//output: viewer result
export async function spiralPlot(): Promise<DG.Viewer> {//(table: DG.DataFrame, activity: DG.Column) {
// Read as dataframe
  const table = await grok.data.files.openTable('Demo:TestJobs:Files:DemoFiles/bio/peptides.csv');
  const activity = await table.columns.addNewCalculated('-log10(Activity)', '0-Log10(${Activity})');
  view = grok.shell.addTableView(table);
  return view.addViewer(SpiralPlot.fromTable(table, {valuesColumnName: activity.name}));
}
