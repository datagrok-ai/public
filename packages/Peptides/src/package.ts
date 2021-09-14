/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {splitAlignedPeptides} from './splitAligned';
import {describe} from './SARViewer/describe';

export const _package = new DG.Package();

async function main() {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  //let peptides = await grok.data.loadTable('https://datagrok.jnj.com/p/ejaeger.il23peptideidp5562/il-23_peptide_idp-5562')
  const path = _package.webRoot + 'files/' + 'aligned.csv';
  const peptides = (await grok.data.loadTable(path));

  const view = grok.shell.addTableView(peptides);

  const grid = await describe(peptides, 'Activity');

  if (grid !== null) {
    view.addViewer(grid);
  }
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
    ], 'Transform peptide sequence data to research insights'
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
        ui.button('Open demo', () => main(), '')
      ]),
      ui.block75([annotationViewerDiv]),
    ]),
    mainDiv,
  ]);
}

//name: Split Sequence
//input: column peptideColumn {semType: alignedSequence}
//tags: panel
//output: widget result
export function splitAlignedSequence(peptideColumn: DG.Column) {
  grok.shell.addTableView(splitAlignedPeptides(peptideColumn));
}
