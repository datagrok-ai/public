/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SARViewer} from './peptide-sar-viewer/sar-viewer';
import {AlignedSequenceCellRenderer, AminoAcidsCellRenderer} from './utils/cell-renderer';
import {Logo} from './peptide-logo-viewer/logo-viewer';
import {StackedBarChart, addViewerToHeader} from './stacked-barchart/stacked-barchart-viewer';
import {ChemPalette} from './utils/chem-palette';
// import { tTest, uTest } from './utils/misc';

export const _package = new DG.Package();
let tableGrid: DG.Grid;

async function main(chosenFile : string) {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  //let peptides =
  //  await grok.data.loadTable('https://datagrok.jnj.com/p/ejaeger.il23peptideidp5562/il-23_peptide_idp-5562');
  const path = _package.webRoot + 'files/' + chosenFile;
  const peptides = (await grok.data.loadTable(path));
  peptides.name = 'Peptides';
  peptides.setTag('dataType', 'peptides');
  const view = grok.shell.addTableView(peptides);
  tableGrid = view.grid;
  peptides.onSemanticTypeDetecting.subscribe((_) => {
    const regexp = new RegExp('^((.+)?-){5,49}(\\w|\\(|\\))+$');
    for (const col of peptides.columns) {
      col.semType = DG.Detector.sampleCategories(col, (s) => regexp.test(s)) ? 'alignedSequence' : null;
      if (col.semType == 'alignedSequence') {
        let maxWidth = 0;
        col.categories.forEach((seq:string) =>{
          maxWidth = maxWidth < seq.length ? seq.length : maxWidth;
        });
        tableGrid.columns.byName(col.name)!.width = maxWidth * 15;
      }
    }
  },
  );

  view.name = 'PeptidesView';

  grok.shell.windows.showProperties = true;

  pi.close();
}

//name: Peptides
//tags: app
export function Peptides() {
  const appDescription = ui.info(
    [
      ui.span(['For more details see LINK ']),
      ui.divText('\n To start the application :', {style: {'font-weight': 'bolder'}}),
      ui.divText('Select the corresponding .csv table with peptide sequences'),
    ], 'Transform peptide sequence data to research insights',
  );
  const annotationViewerDiv = ui.div();

  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  const mainDiv = ui.div();
  grok.shell.newView('Peptides', [
    appDescription,
    ui.h2('Choose .csv file'),
    ui.div([
      ui.block25([
        ui.button('Open simple case demo', () => main('aligned.csv'), ''),
        ui.button('Open complex case demo', () => main('aligned_2.csv'), ''),
      ]),
      ui.block75([annotationViewerDiv]),
    ]),
    mainDiv,
  ]);

  // tests test lol
  // const a1 = [13.3, 6.0, 20.0, 8.0, 14.0, 19.0, 18.0, 25.0, 16.0, 24.0, 15.0, 1.0, 15.0];
  // const a2 = [22.0, 16.0, 21.7, 21.0, 30.0, 26.0, 12.0, 23.2, 28.0, 23.0];
  // console.log(uTest(a1, a2));
  // console.log(tTest(a1, a2));
}

//name: Peptides
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export async function analyzePeptides(col: DG.Column): Promise<DG.Widget> {
  // let defaultColumn: DG.Column | null = col;
  let tempCol = null;
  for (const column of col.dataFrame.columns.numerical) {
    column.type === DG.TYPE.FLOAT ? tempCol = column : null;
  }
  const defaultColumn: DG.Column = col.dataFrame.col('activity') || col.dataFrame.col('IC50') || tempCol;

  const activityColumnChoice = ui.columnInput('Activity column', col.dataFrame, defaultColumn);
  const activityScalingMethod = ui.choiceInput('Activity scaling', 'none', ['none', 'lg', '-lg']);
  // it doesn't make sense to do apply log to values less than 0
  activityColumnChoice.onChanged((_: any) => {
    activityScalingMethod.enabled = activityColumnChoice.value && activityColumnChoice.value.min > 0;
    activityScalingMethod.setTooltip('Function to apply for each value in activity column');
  });
  activityColumnChoice.fireChanged();

  const startBtn = ui.button('Launch SAR', async () => {
    if (activityColumnChoice.value.type === DG.TYPE.FLOAT) {
      const options = {
        'activityColumnColumnName': activityColumnChoice.value.name,
        'activityScalingMethod': activityScalingMethod.value,
      };
      (grok.shell.v as DG.TableView).addViewer('peptide-sar-viewer', options);
      const widgProm = col.dataFrame.plot.fromType('StackedBarChartAA');
      addViewerToHeader(tableGrid, widgProm);
      // tableGrid.dataFrame!.columns.names().forEach((name:string)=>{
      //   col = tableGrid.dataFrame!.columns.byName(name);
      //   if (col.semType == 'aminoAcids') {
      //     let maxLen = 0;
      //     col.categories.forEach( (ent:string)=>{
      //       if ( ent.length> maxLen) {
      //         maxLen = ent.length;
      //       }
      //     });
      //     tableGrid.columns.byName(name)!.width = maxLen*10;
      //   }
      // });
    } else {
      grok.shell.error('The activity column must be of floating point number type!');
    }
  });

  const viewer = await col.dataFrame.plot.fromType('peptide-logo-viewer');

  return new DG.Widget(ui.divV([viewer.root, ui.inputs([activityColumnChoice, activityScalingMethod]), startBtn]));
}

//name: peptide-sar-viewer
//description: Peptides SAR Viewer
//tags: viewer
//output: viewer result
export function sar(): SARViewer {
  return new SARViewer();
}

//name: StackedBarchart Widget
//tags: panel, widgets
//input: column col {semType: aminoAcids}
//output: widget result

export async function stackedBarchartWidget(col:DG.Column):Promise<DG.Widget> {
  const viewer = await col.dataFrame.plot.fromType('StackedBarChartAA');
  const panel = ui.divH([viewer.root]);
  return new DG.Widget(panel);
}

//name: Peptide Molecule
//tags: panel, widgets
//input: string pep {semType: alignedSequence}
//output: widget result

export async function pepMolGraph(pep:string):Promise<DG.Widget> {
  const split = pep.split('-');
  const mols = [];
  for (let i = 1; i < split.length - 1; i++) {
    if (split[i] in ChemPalette.AASmiles ) {
      const aar = ChemPalette.AASmiles[split[i]]
      mols[i] = aar.substr(0, aar.length-1);
    } else if (!split[i]||split[i]=='-') {
      mols[i] = '';
    } else {
      return new DG.Widget(ui.divH([]));
    }
  }
  console.error(mols);
  console.error(mols.join('')+'COOH');
  const sketch = grok.chem.svgMol(mols.join(''));
  const panel = ui.divH([sketch]);
  return new DG.Widget(panel);
}

//name: StackedBarChartAA
//description: Creates an awesome viewer
//tags: viewer
//output: viewer result
export function stackedBarChart():DG.JsViewer {
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
