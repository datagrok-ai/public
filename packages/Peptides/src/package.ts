/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

import {SARViewerBase} from './peptide-sar-viewer/sar-viewer';
import {
  AlignedSequenceCellRenderer,
  AminoAcidsCellRenderer,
} from './utils/cell-renderer';
import {Logo} from './peptide-logo-viewer/logo-viewer';
import {StackedBarChart, addViewerToHeader} from './stacked-barchart/stacked-barchart-viewer';
import {ChemPalette} from './utils/chem-palette';
// import { tTest, uTest } from './utils/misc';

import {DimensionalityReducer} from '@datagrok-libraries/utils/src/reduce_dimensionality';

import {getSequenceMolecularWeight} from './peptide-sar-viewer/molecular_measure';

export const _package = new DG.Package();
let tableGrid: DG.Grid;
let currentDf: DG.DataFrame;
let alignedSequenceCol: DG.Column;
let view: DG.TableView;

async function main(chosenFile: string) {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  //let peptides =
  //  await grok.data.loadTable('https://datagrok.jnj.com/p/ejaeger.il23peptideidp5562/il-23_peptide_idp-5562');
  const path = _package.webRoot + 'files/' + chosenFile;
  const peptides = (await grok.data.loadTable(path));
  peptides.name = 'Peptides';
  peptides.setTag('dataType', 'peptides');
  const view = grok.shell.addTableView(peptides);
  tableGrid = view.grid;
  // peptides.onSemanticTypeDetecting.subscribe((_: any) => {
  //   const regexp = new RegExp(/^([^-^\n]*-){2,49}(\w|\(|\))+$/);
  //   for (const col of peptides.columns) {
  //     col.semType = DG.Detector.sampleCategories(col, (s: any) => regexp.test(s.trim())) ? 'alignedSequence' : null;
  //     if (col.semType == 'alignedSequence') {
  //       expandColumn(col, tableGrid, (ent)=>{
  //         const subParts:string[] = ent.split('-');
  //         // eslint-disable-next-line no-unused-vars
  //         const [text, _] = processSequence(subParts);
  //         let textSize = 0;
  //         text.forEach((aar)=>{
  //           textSize += aar.length;
  //         });
  //         return textSize;
  //       });
  //     }
  //   }
  // });

  view.name = 'PeptidesView';

  grok.shell.windows.showProperties = true;

  pi.close();
}

//name: Peptides
//tags: app
export function Peptides() {
  const wikiLink = ui.link('wiki', 'https://github.com/datagrok-ai/public/blob/master/help/domains/bio/peptides.md');
  const textLink = ui.inlineText(['For more details, see our ', wikiLink, '.']);

  const appDescription = ui.info(
    [
      // ui.divText('\n To start the application :', {style: {'font-weight': 'bolder'}}),
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

  const annotationViewerDiv = ui.div();

  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  const mainDiv = ui.div();
  grok.shell.newView('Peptides', [
    appDescription,
    ui.info([textLink]),
    ui.div([
      ui.block25([
        ui.button('Open peptide sequences demonstration set', () => main('aligned.csv'), ''),
        ui.button('Open complex case demo', () => main('aligned_2.csv'), ''),
      ]),
      ui.block75([annotationViewerDiv]),
    ]),
    mainDiv,
  ]);
}

//name: peptideSimilaritySpace
//input: dataframe table
//input: column alignedSequencesColumn {semType: alignedSequence}
//input: string method {choices: ['TSNE', 'SPE', 'PSPE']} = 'TSNE'
//input: string measure {choices: ['Levenshtein', 'Jaro-Winkler']} = 'Levenshtein'
//input: int cyclesCount = 100
//output: graphics
export function peptideSimilaritySpace(
  table: DG.DataFrame,
  alignedSequencesColumn: DG.Column,
  method: string,
  measure: string,
  cyclesCount: number,
  activityColumnName: string,
) {
  const axesNames = ['~X', '~Y', 'MW'];

  const reducer = new DimensionalityReducer(
    alignedSequencesColumn.toList(),
    method,
    measure,
    {cycles: cyclesCount},
  );
  const embcols = reducer.transform(true);
  const columns = Array.from(embcols, (v, k) => (DG.Column.fromFloat32Array(axesNames[k], v)));

  const sequences = alignedSequencesColumn.toList();
  const mw: Float32Array = new Float32Array(sequences.length).fill(0);

  let currentSequence;
  for (let i = 0; i < sequences.length; ++i) {
    currentSequence = sequences[i];
    mw[i] = currentSequence == null ? 0 : getSequenceMolecularWeight(currentSequence);
  }

  columns.push(DG.Column.fromFloat32Array('MW', mw));

  const edf = DG.DataFrame.fromColumns(columns);

  // Add new axes.
  for (const axis of axesNames) {
    const col = table.col(axis);

    if (col == null) {
      table.columns.insert(edf.getCol(axis));
    } else {
      table.columns.replace(col, edf.getCol(axis));
    }
  }

  const view = (grok.shell.v as DG.TableView);

  const viewer = view.addViewer(DG.VIEWER.SCATTER_PLOT, {x: '~X', y: '~Y', color: activityColumnName, size: 'MW'});
  // (viewer as DG.ScatterPlotViewer).zoom(
  //   edf.getCol('~X').min,
  //   edf.getCol('~Y').min,
  //   edf.getCol('~X').max,
  //   edf.getCol('~Y').max,
  // );
  return viewer;
}

//name: Peptides
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export async function analyzePeptides(col: DG.Column): Promise<DG.Widget> {
  view = (grok.shell.v as DG.TableView);
  tableGrid = view.grid;
  currentDf = tableGrid.dataFrame!;
  let tempCol = null;
  for (const column of currentDf.columns.numerical) {
    tempCol = column.type === DG.TYPE.FLOAT ? column : null;
  }
  const defaultColumn: DG.Column = currentDf.col('activity') || currentDf.col('IC50') || tempCol;
  const histogramHost = ui.div([]);

  let hist: DG.Viewer;

  const activityScalingMethod = ui.choiceInput(
    'Activity scaling',
    'none',
    ['none', 'lg', '-lg'],
    async (currentMethod: string) => {
      const currentActivityCol = activityColumnChoice.value.name;
      const tempDf = currentDf.clone(currentDf.filter, [currentActivityCol]);
      switch (currentMethod) {
      case 'lg':
        await tempDf.columns.addNewCalculated('scaledActivity', 'Log10(${' + currentActivityCol + '})');
        break;
      case '-lg':
        await tempDf.columns.addNewCalculated('scaledActivity', '-1*Log10(${' + currentActivityCol + '})');
        break;
      default:
        await tempDf.columns.addNewCalculated('scaledActivity', '${' + currentActivityCol + '}');
        break;
      }
      hist = tempDf.plot.histogram({
        filteringEnabled: false,
        valueColumnName: 'scaledActivity',
        legendVisibility: 'Never',
        showXAxis: true,
        showColumnSelector: false,
        showRangeSlider: false,
      // bins: b,
      });
      histogramHost.lastChild?.remove();
      histogramHost.appendChild(hist.root);
    });
  activityScalingMethod.setTooltip('Function to apply for each value in activity column');

  const activityScalingMethodState = function(_: any) {
    activityScalingMethod.enabled =
      activityColumnChoice.value && DG.Stats.fromColumn(activityColumnChoice.value, currentDf.filter).min > 0;
    activityScalingMethod.fireChanged();
  };
  const activityColumnChoice = ui.columnInput(
    'Activity column',
    currentDf,
    defaultColumn,
    activityScalingMethodState,
  );
  activityColumnChoice.fireChanged();
  activityScalingMethod.fireChanged();

  const startBtn = ui.button('Launch SAR', async () => {
    if (activityColumnChoice.value.type === DG.TYPE.FLOAT) {
      const options = {
        'activityColumnColumnName': activityColumnChoice.value.name,
        'activityScalingMethod': activityScalingMethod.value,
      };
      for (let i = 0; i < tableGrid.columns.length; i++) {
        const col = tableGrid.columns.byIndex(i);
        if (col &&
            col.name &&
            col.name != 'IC50'&&
            col.column?.semType != 'aminoAcids'
        ) {
          //@ts-ignore
          tableGrid.columns.byIndex(i)?.visible = false;
        }
      }

      alignedSequenceCol = col;

      const sarViewer = view.addViewer('peptide-sar-viewer', options);
      const peptideSpaceViewer = peptideSimilaritySpace(
        currentDf,
        col,
        'TSNE',
        'Levenshtein',
        100,
        `${activityColumnChoice}Scaled`,
      );
      view.dockManager.dock(peptideSpaceViewer, 'down');
      view.dockManager.dock(sarViewer, 'right');

      const StackedBarchartProm = currentDf.plot.fromType('StackedBarChartAA');
      addViewerToHeader(tableGrid, StackedBarchartProm);
    } else {
      grok.shell.error('The activity column must be of floating point number type!');
    }
  });

  const viewer = await currentDf.plot.fromType('peptide-logo-viewer');

  return new DG.Widget(
    ui.divV([viewer.root, ui.inputs([activityColumnChoice, activityScalingMethod]), startBtn, histogramHost]),
  );
}

//name: peptide-sar-viewer
//description: Peptides SAR Viewer
//tags: viewer
//output: viewer result
export function sar(): SARViewerBase {
  return new SARViewerBase();
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
//input: string pep {semType: alignedSequence}
//output: widget result
export async function pepMolGraph(pep: string): Promise<DG.Widget> {
  const split = pep.split('-');
  const mols = [];
  for (let i = 1; i < split.length - 1; i++) {
    if (split[i] in ChemPalette.AASmiles) {
      const aar = ChemPalette.AASmiles[split[i]];
      mols[i] = aar.substr(0, aar.length - 1);
    } else if (!split[i] || split[i] == '-') {
      mols[i] = '';
    } else {
      return new DG.Widget(ui.divH([]));
    }
  }
  console.error(mols);
  console.error(mols.join('') + 'COOH');
  const sketch = grok.chem.svgMol(mols.join('') + 'O');
  const panel = ui.divH([sketch]);
  return new DG.Widget(panel);
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
  //TODO: update viewers right when the changes get applied
  const sequenceInput = ui.textInput('', alignedSequenceCol.get(currentDf.currentRowIdx));
  (sequenceInput.input as HTMLElement).style.height = '50px';
  (sequenceInput.input as HTMLElement).style.overflow = 'hidden';

  const applyChangesBtn = ui.button('Apply', () => {
    alignedSequenceCol.set(currentDf.currentRowIdx, sequenceInput.value);
  });

  const resetBtn = ui.button(
    ui.iconFA('redo'),
    () => sequenceInput.value = alignedSequenceCol.get(currentDf.currentRowIdx),
    'Reset',
  );
  $(resetBtn).addClass('dt-snippet-editor-icon dt-reset-icon');

  return new DG.Widget(ui.divV([resetBtn, sequenceInput.root, applyChangesBtn], 'dt-textarea-box'));
}

//name: testPeptideSimilaritySpace
export async function testPeptideSimilaritySpace() {
  const df = await grok.data.files.openTable('Demo:TestJobs:Files:DemoFiles/bio/peptides.csv');
  grok.shell.addTableView(df);
  peptideSimilaritySpace(df, df.getCol('AlignedSequence'), 'PSPE', 'Levenshtein', 100, 'Activity');
}
