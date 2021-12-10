/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

function nglPdbPanelWidget(pdbId: string) {
  const host = ui.div([], 'd4-ngl-viewer');
  //@ts-ignore
  const stage = new NGL.Stage(host);
  handleResize(host, stage);
  stage.loadFile(`rcsb://${pdbId}`, {defaultRepresentation: true});
  grok.shell.v.append(host);
  //return DG.Widget.fromRoot(host);
}

function handleResize(host: HTMLElement, stage: any) {
  const canvas = host.querySelector('canvas');

  function resize() {
    canvas!.width = Math.floor(canvas!.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(canvas!.clientHeight * window.devicePixelRatio);
    stage.handleResize();
  }

  ui.onSizeChanged(host).subscribe((_) => resize());
  resize();
}

async function fetchPDBAsString(pdbid: string, RCSB = 'https://files.rcsb.org/download/$.pdb') {
  const link = RCSB.replace('$', pdbid.toLowerCase());
  const response = await fetch(link);
  const pdbText = await response.text();
  return pdbText;
}

async function main() {
  //const body = await fetchPDBAsString('2v0a');
  //console.log(body);
  nglPdbPanelWidget('2v0a');
}

//name: PDB viewer app
//tags: app
export function PDBViewer() {
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
  grok.shell.newView('PDB viewer', [
    appDescription,
    ui.info([textLink]),
    ui.div([
      ui.block25([
        ui.button('Open peptide sequences demonstration set', () => main(), ''),
        ui.button('Open complex case demo', () => main(), ''),
      ]),
      ui.block75([annotationViewerDiv]),
    ]),
    mainDiv,
  ]);
}
