/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SARViewer} from './peptide-sar-viewer/sar-viewer';
import {AlignedSequenceCellRenderer} from './utils/cell-renderer'

export const _package = new DG.Package();

async function main() {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  //let peptides = await grok.data.loadTable('https://datagrok.jnj.com/p/ejaeger.il23peptideidp5562/il-23_peptide_idp-5562')
  const path = _package.webRoot + 'files/' + 'aligned.csv';
  const peptides = (await grok.data.loadTable(path));

  const view = grok.shell.addTableView(peptides);
  view.name = 'PeptidesView';

  pi.close();
}

//name: Peptides
//tags: app
export function Peptides() {
  const appDescription = ui.info(
    [
      ui.span(['For more deatails see LINK ']),
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
        ui.button('Open demo', () => main(), ''),
      ]),
      ui.block75([annotationViewerDiv]),
    ]),
    mainDiv,
  ]);
}

//name: SAR Viewer Help
//tags: panel, widget
//input: column _ {semType: alignedSequence}
//output: widget result
export function SARViewerHelp(_: DG.Column): DG.Widget {
  const helpStr = 
  "Circle size in the viewer is based on ratio and color is based on MAD.\n\
  \n\
  Statistics:\n\
  Count - number of peptides containing AAR at position\n\
  MAD - mean absolute deviation\n\
  Q1 - first quartile of activity\n\
  Median - median of activity\n\
  Q3 - third quartile of activity\n\
  IQR - interquartile range\n\
  CQV - coefficient of quartile variation (quartile coefficient of dispersion)\n\
  Ratio - share of peptides containing AAR at position\n";
  const div = ui.divV(helpStr.split("\n").map((line) => {
    return ui.divText(line);
  }))
  return new DG.Widget(div);
}

//name: Analyze Peptides
//tags: panel, widget
//input: column col {semType: alignedSequence}
//output: widget result
export function analyzePeptides(col: DG.Column): DG.Widget {
  const activityColumnChoice = ui.columnInput('Activity column', col.dataFrame, col);
  const activityScalingMethod = ui.choiceInput('Activity scaling', 'none', ['none', 'lg', '-lg']);
  const showHistogram = ui.boolInput('Show histogram', false);
  
  const startBtn = ui.button('Start', async () => {
    if (activityColumnChoice.value.type === DG.TYPE.FLOAT) {
      const peptidesView = grok.shell.v;
      if (peptidesView.type === DG.VIEW_TYPE.TABLE_VIEW) {
        const options = {
          'activityColumnColumnName': activityColumnChoice.value.name,
          'activityScalingMethod': activityScalingMethod.value,
          'showHistogram': showHistogram.value
        };
        (<DG.TableView>peptidesView).addViewer('peptide-sar-viewer', options);
      }
    } else {
      grok.shell.error("The activity column must be of double type!");
    }
  });
  return new DG.Widget(ui.divV([activityColumnChoice.root, activityScalingMethod.root, showHistogram.root, startBtn]));
}

//name: peptide-sar-viewer
//description: Peptides SAR Viewer
//tags: viewer
//output: viewer result
export function sar(): SARViewer {
  return new SARViewer();
}

//name: alignedSequenceCellRenderer
//tags: cellRenderer, cellRenderer-alignedSequence
//meta-cell-renderer-sem-type: alignedSequence
//output: grid_cell_renderer result
export function alignedSequenceCellRenderer() {
  return new AlignedSequenceCellRenderer();
}
