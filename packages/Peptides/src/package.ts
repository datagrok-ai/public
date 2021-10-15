/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {SARViewer} from './peptide-sar-viewer/sar-viewer';
import {AlignedSequenceCellRenderer, AminoAcidsCellRenderer} from './utils/cell-renderer';
import {Logo} from './peptide-logo-viewer/logo-viewer';
import {StackedBarChart, addViewerToHeader} from './stacked-barchart/stacked-barchart-viewer';
import {ChemPalette} from './utils/chem-palette';
import {describe} from './peptide-sar-viewer/describe';
import {SARViewerVertical} from './peptide-sar-viewer/sar-viewer-vertical';
// import { tTest, uTest } from './utils/misc';

export const _package = new DG.Package();
let tableGrid: DG.Grid;

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
  peptides.onSemanticTypeDetecting.subscribe((_) => {
    const regexp = new RegExp('^((.+)?-){5,49}(\\w|\\(|\\))+$');
    for (const col of peptides.columns) {
      col.semType = DG.Detector.sampleCategories(col, (s) => regexp.test(s)) ? 'alignedSequence' : null;
      if (col.semType == 'alignedSequence') {
        setTimeout(()=>{
          let maxWidth = 0;
          for (const s of col.categories) {
            const subParts:string[] = s.split('-');
            const simplified = !subParts.some((amino, index)=>{
              return amino.length>1&&index!=0&&index!=subParts.length-1;
            });
            const text:string[] = [];
            subParts.forEach((amino: string, index) => {
              if (index < subParts.length) {
                const gap = simplified?'':' ';
                amino += `${amino?'':'-'}${gap}`;
              }
              text.push(amino);
            });
            let textSize = 0;
            text.forEach((aar)=>{
              textSize += aar.length;
            });
            maxWidth = maxWidth < textSize ? textSize: maxWidth;
          }
          tableGrid.col(col.name)!.width = Math.min(maxWidth * 10, 1000);
        }, 500);
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
  const textLink = ui.div();
  textLink.innerHTML = `For more details, see our <a href="https://github.com/datagrok-ai/public/blob/master/help/domains/bio/peptides.md">[wiki]</a>.`;

  const appDescription = ui.info(
    [
      // ui.divText('\n To start the application :', {style: {'font-weight': 'bolder'}}),
      ui.divText(
        ' - automatic recognition of peptide sequences\n' +
      ' - native integration with tons of Datagrok out-of-the box features (visualization, filtering, clustering, multivariate analysis, etc)\n' +
      ' - custom rendering in the spreadsheet\n' +
      ' - interactive logo plots\n' +
      ' - rendering residues\n' +
      ' - structure-activity relationship:\n \n' +
      'a) highlighting statistically significant changes in activity in the [position, monomer] spreadsheet\n' +
      'b) for the specific [position, monomer], visualizing changes of activity distribution (specific monomer in this position vs rest of the monomers in this position)\n' +
      'c) interactivity\n \t',
      ),
    ],

    `Use and analyse peptide sequence data to support your research:\n`,
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
    //ui.h2('Choose .csv file'),
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

  const activityScalingMethod = ui.choiceInput('Activity scaling', 'none', ['none', 'lg', '-lg']);
  activityScalingMethod.setTooltip('Function to apply for each value in activity column');

  const activityScalingMethodState = function(_: any) {
    activityScalingMethod.enabled =
      activityColumnChoice.value && DG.Stats.fromColumn(activityColumnChoice.value, col.dataFrame.filter).min > 0;
  };
  const activityColumnChoice = ui.columnInput(
    'Activity column',
    col.dataFrame,
    defaultColumn,
    activityScalingMethodState,
  );
  activityColumnChoice.fireChanged();

  const startBtn = ui.button('Launch SAR', async () => {
    if (activityColumnChoice.value.type === DG.TYPE.FLOAT) {
      const options = {
        'activityColumnColumnName': activityColumnChoice.value.name,
        'activityScalingMethod': activityScalingMethod.value,
      };
      // @ts-ignore
      for (let i = 0; i < tableGrid.columns.length; i++) {
        const col = tableGrid.columns.byIndex(i);
        if (col &&
            col.name &&
            col.name != 'IC50'&&
            col.column?.semType != 'aminoAcids') {
          // @ts-ignore
          tableGrid.columns.byIndex(i)?.visible = false;
        }
      }

      //await describe(col.dataFrame, activityColumnChoice.value.name, activityScalingMethod.value, false, tableGrid);

      // @ts-ignore


      const viewer = DG.Viewer.fromType('peptide-sar-viewer', tableGrid.table, options);
      (grok.shell.v as DG.TableView).addViewer(viewer);
      const refNode = (grok.shell.v as DG.TableView).dockManager.dock(viewer, 'right');

      const hist = DG.Viewer.fromType('peptide-sar-viewer-vertical', tableGrid.table, options);
      (grok.shell.v as DG.TableView).addViewer(hist);
      (grok.shell.v as DG.TableView).dockManager.dock(hist, DG.DOCK_TYPE.DOWN, refNode);

      // (grok.shell.v as DG.TableView).addViewer('peptide-sar-viewer', options);
      // (grok.shell.v as DG.TableView).addViewer('peptide-sar-viewer-vertical', options);
      // @ts-ignore
      //view.dockManager.dock(ui.divText('bottom'), 'down');

      // @ts-ignore
      //console.error(sarViewer.view.dockNode);

      const StackedBarchartProm = col.dataFrame.plot.fromType('StackedBarChartAA');
      addViewerToHeader(tableGrid, StackedBarchartProm);

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

//name: peptide-sar-viewer-vertical
//description: Peptides SAR Viewer
//tags: viewer
//output: viewer result
export function sarVertical(): SARViewerVertical {
  return new SARViewerVertical();
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
//description: Creates an awesome viewer
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
