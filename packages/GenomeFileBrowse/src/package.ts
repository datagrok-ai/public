/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getGenomeViewer} from './genome-file-browse-viewer-wrapper'
export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

async function readGenomeConfig(fileContent: string, fileName: string): Promise<DG.View> {
  const root = ui.div([], { classes: 'd4-spectrum' });
  const spectrumComponent = getGenomeViewer(root, fileContent);
  const spectrumHtmlView = DG.View.fromRoot(spectrumComponent.host);
  spectrumHtmlView.name = fileName;
  spectrumComponent.loadData(fileContent);

  spectrumHtmlView.subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
    const cv = grok.shell.v;
    if ('id' in cv && cv.id && spectrumHtmlView.id && cv.id === spectrumHtmlView.id) {
      spectrumComponent.loadData(fileContent);
      spectrumHtmlView.root.appendChild(spectrumComponent.host);
    }
  }))
  return spectrumHtmlView;
}

//tags: fileViewer
//input: file file
//output: view v
export async function previewNMRFromDX(fileData: DG.FileInfo): Promise<DG.View> {
  let viewerHost  = ui.div();
  getGenomeViewer(viewerHost, await fileData.readAsString());
  let view = DG.View.fromRoot(viewerHost);
  grok.shell.addView(view);
  return view;
}


//name: spectraDataCheck
//input: string fileData
//output: bool res
export async function spectraDataCheck(fileData: string) {
  return !fileData.includes('NMR SPECTRUM');
}

//name: previewSpectraData
//tags: fileViewer
//meta.fileViewer: json
//meta.fileViewerCheck: SpectraViewer:spectraDataCheck
//input: file file
//output: view v
export async function previewSpectraData(fileData: DG.FileInfo): Promise<DG.View> {
  return await readGenomeConfig(await fileData.readAsString(), fileData.name);
}

//name: importSpectraData
//tags: file-handler
//meta.ext: json
//meta.fileViewerCheck: SpectraViewer:spectraDataCheck
//input: string fileString
//output: list v
export async function importSpectraData(fileString: string) {
  const res = await readGenomeConfig(fileString, 'Genome Browse');
  const newView = grok.shell.addView(res);
  newView != res && res.subs.forEach((s) => newView.subs.push(s));
  return [];
}