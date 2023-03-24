/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {analyzePeptidesUI} from './widgets/peptides';
import {PeptideSimilaritySpaceWidget} from './utils/peptide-similarity-space';
import {manualAlignmentWidget} from './widgets/manual-alignment';
import {MonomerPosition, MostPotentResiduesViewer} from './viewers/sar-viewer';

import {PeptideSpaceViewer} from './viewers/peptide-space-viewer';
import {LogoSummary} from './viewers/logo-summary';
import {MonomerWorks} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {PeptidesModel} from './model';

let monomerWorks: MonomerWorks | null = null;

export const _package = new DG.Package();

export function getMonomerWorks(): MonomerWorks | null {
  return monomerWorks;
}

//tags: init
export async function initPeptides(): Promise<void> {
  monomerWorks ??= new MonomerWorks(await grok.functions.call('Bio:getBioLib'));
}

async function openDemoData(chosenFile: string): Promise<void> {
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
export function Peptides(): void {
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
      ui.button('Simple demo', () => openDemoData('aligned.csv'), ''),
      ui.button('Complex demo', () => openDemoData('aligned_2.csv'), ''),
      ui.button('HELM demo', () => openDemoData('aligned_3.csv'), ''),
    ]),
  ]);
}

//top-menu: Bio | SAR | Peptides...
//name: Bio Peptides
export function peptidesDialog(): DG.Dialog {
  const analyzeObject = analyzePeptidesUI(grok.shell.t);
  const dialog = ui.dialog('Analyze Peptides').add(analyzeObject.host).onOK(async () => {
    const startSuccess = analyzeObject.callback();
    if (!startSuccess)
      dialog.show();
  });
  return dialog.show();
}

//name: Peptides
//tags: panel, widgets
//input: column col {semType: Macromolecule}
//output: widget result
export function peptidesPanel(col: DG.Column): DG.Widget {
  const analyzeObject = analyzePeptidesUI(col.dataFrame, col);
  return new DG.Widget(analyzeObject.host);
}

//name: peptide-sar-viewer
//description: Peptides SAR Viewer
//tags: viewer
//meta.icon: files/icons/peptide-sar-viewer.svg
//output: viewer result
export function sar(): MonomerPosition {
  return new MonomerPosition();
}

//name: peptide-sar-viewer-vertical
//description: Peptides Vertical SAR Viewer
//tags: viewer
//meta.icon: files/icons/peptide-sar-vertical-viewer.svg
//output: viewer result
export function sarVertical(): MostPotentResiduesViewer {
  return new MostPotentResiduesViewer();
}

//name: logo-summary-viewer
//tags: viewer
//meta.icon: files/icons/logo-summary-viewer.svg
//output: viewer result
export function logoSummary(): LogoSummary {
  return new LogoSummary();
}

//name: peptide-space-viewer
//description: Peptide Space Viewer
//tags: viewer
//meta.icon: files/icons/peptide-space-viewer.svg
//output: viewer result
export function peptideSpace(): PeptideSpaceViewer {
  return new PeptideSpaceViewer();
}

//name: Manual Alignment
//tags: panel, widgets
//input: string _monomer {semType: Monomer}
//output: widget result
export function manualAlignment(_monomer: string): DG.Widget {
  //TODO: recalculate Molfile and Molecule panels on sequence update
  const df = grok.shell.t;
  const model: PeptidesModel | null = df?.temp[PeptidesModel.modelName];
  if (!model)
    return new DG.Widget(ui.divText('Manual alignment works with peptides analysis'));

  const col = df.getCol(model.settings.sequenceColumnName!);
  return manualAlignmentWidget(col, df);
}

//name: Peptide Space
//tags: panel, widgets
//input: column col {semType: Macromolecule}
//output: widget result
export async function peptideSpacePanel(col: DG.Column): Promise<DG.Widget> {
  const widget = new PeptideSimilaritySpaceWidget(col, grok.shell.v as DG.TableView);
  return widget.draw();
}
