/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getGenomeBrowserComponent } from './genome-file-browse-viewer-wrapper'
export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

async function genomeBrowser(data: string, fileName: string = 'GenomeBrowser'): Promise<DG.View> {
  let viewerHost = ui.div();
  getGenomeBrowserComponent(viewerHost, data);
  viewerHost.classList.add('testClass');
  let view = DG.View.fromRoot(viewerHost);
  view.name = fileName;
  return view;
}

//name: previewGenomeBrowse
//tags: fileViewer
//meta.fileViewer: json
//meta.fileViewerCheck: GenomeBrowser:checkGenomeConfig
//meta.fileViewerCheck: GenomeBrowser:checkGenomeConfig
//input: file file
//output: view v
export async function previewGenomeFileBrowse(fileData: DG.FileInfo): Promise<DG.View> {

  return DG.View.fromViewAsync(async () => await genomeBrowser( await fileData.readAsString(), fileData.name));
}