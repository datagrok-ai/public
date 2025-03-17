/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getGenomeFileBrowserComponent} from './genome-file-browse-viewer-wrapper'
export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: genomeFileBrowse
//input: file file
//output: view v
export async function genomeFileBrowse(fileData: DG.FileInfo): Promise<DG.View> {
  let viewerHost  = ui.div();
  getGenomeFileBrowserComponent(viewerHost, await fileData.readAsString());
  let view = DG.View.fromRoot(viewerHost);
  view.name = fileData.name;
  grok.shell.addView(view);
  return view;
}


//tags: fileViewer
//meta.fileViewer: json
//meta.fileViewerCheck: GenomeFileBrowse:checkGenomeConfig
//input: file file
//output: view v
export async function previewGenomeFileBrowse(fileData: DG.FileInfo): Promise<DG.View>{
  let viewerHost  = ui.div();
  getGenomeFileBrowserComponent(viewerHost, await fileData.readAsString());
  let view = DG.View.fromRoot(viewerHost);
  view.name = fileData.name;
  grok.shell.addView(view);
  return view;
}
  

//name: genomeFileBrowseHadnler
//input: string bytes
//output: list tables
//tags: file-handler
//meta.fileViewerCheck: GenomeFileBrowse:checkGenomeConfig
//meta.ext: json
export async function genomeFileBrowseHadnler(bytes: string) {
  let viewerHost  = ui.div();
  getGenomeFileBrowserComponent(viewerHost, bytes);
  let view = DG.View.fromRoot(viewerHost);
  view.name = 'GenomeFileBrowse';
  grok.shell.addView(view);
  return [];
}


//name: checkGenomeConfig
//input: string content
//output: bool result
export function checkGenomeConfig(content: string) {
  let jsonObj  = JSON.parse(content);
  return content.length < 1000000 && jsonObj.tracks && (jsonObj.assemblies || jsonObj.assembly);//megabyte
}
