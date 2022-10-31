/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as C from './utils/constants';

import {analyzePeptidesWidget} from './widgets/peptides';
import {PeptideSimilaritySpaceWidget} from './utils/peptide-similarity-space';
import {manualAlignmentWidget} from './widgets/manual-alignment';
import {MutationCliffsViewer, MostPotentResiduesViewer} from './viewers/sar-viewer';

import {PeptideSpaceViewer} from './viewers/peptide-space-viewer';
import {LogoSummary} from './viewers/logo-summary';

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
      ui.button('HELM demo', () => main('aligned_3.csv'), ''),
    ]),
  ]);
}

//top-menu: Bio | Peptides...
//name: Bio Peptides
export async function peptidesDialog(): Promise<DG.Dialog> {
  const dialog = ui.dialog().add(await analyzePeptidesWidget(grok.shell.t));
  return dialog.show();
}

//name: Peptides
//tags: panel, widgets
//input: column col {semType: Macromolecule}
//output: widget result
export async function peptidesPanel(col: DG.Column): Promise<DG.Widget> {
  [currentTable, alignedSequenceColumn] = getOrDefine(col.dataFrame, col);
  return analyzePeptidesWidget(currentTable, alignedSequenceColumn);
}

//name: peptide-sar-viewer
//description: Peptides SAR Viewer
//tags: viewer
//output: viewer result
export function sar(): MutationCliffsViewer {
  return new MutationCliffsViewer();
}

//name: peptide-sar-viewer-vertical
//description: Peptides Vertical SAR Viewer
//tags: viewer
//output: viewer result
export function sarVertical(): MostPotentResiduesViewer {
  return new MostPotentResiduesViewer();
}

//name: logo-summary-viewer
//tags: viewer
//output: viewer result
export function logoSummary(): LogoSummary {
  return new LogoSummary();
}

//name: peptide-space-viewer
//description: Peptide Space Viewer
//tags: viewer
//output: viewer result
export function peptideSpace(): PeptideSpaceViewer {
  return new PeptideSpaceViewer();
}

//name: Manual Alignment
//tags: panel, widgets
//input: string _monomer {semType: Monomer}
//output: widget result
export function manualAlignment(_monomer: string): DG.Widget {
  [currentTable, alignedSequenceColumn] = getOrDefine();
  //TODO: recalculate Molfile and Molecule panels on sequence update
  return manualAlignmentWidget(alignedSequenceColumn, currentTable);
}

//name: Peptide Space
//tags: panel, widgets
//input: column col {semType: Macromolecule}
//output: widget result
export async function peptideSpacePanel(col: DG.Column): Promise<DG.Widget> {
  [currentTable, alignedSequenceColumn] = getOrDefine(col.dataFrame, col);
  const widget = new PeptideSimilaritySpaceWidget(col, grok.shell.v as DG.TableView);
  return widget.draw();
}

//name: Get Peptides Structure
//tags: panel, widgets
//input: column col {semType: Macromolecule}
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

function getOrDefine(dataframe?: DG.DataFrame, column?: DG.Column | null): [DG.DataFrame, DG.Column] {
  dataframe ??= grok.shell.t;
  // column ??= dataframe.columns.bySemType(C.SEM_TYPES.MACROMOLECULE);
  column ??= dataframe.getCol(C.COLUMNS_NAMES.MACROMOLECULE);
  if (column === null)
    throw new Error('Table does not contain aligned sequence columns');

  return [dataframe, column];
}
