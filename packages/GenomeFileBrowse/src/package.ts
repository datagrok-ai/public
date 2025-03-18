/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getGenomeFileBrowserComponent } from './genome-file-browse-viewer-wrapper'
export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

async function genomeFileBrowse(data: string, fileName : string = 'GenomeFileBrowse'): Promise<DG.View> {
  let viewerHost = ui.div();
  getGenomeFileBrowserComponent(viewerHost, data);
  let view = DG.View.fromRoot(viewerHost);
  view.name = fileName;
  return view;
}


//name: previewGenomeFileBrowse
//tags: fileViewer
//meta.fileViewer: json
//meta.fileViewerCheck: GenomeFileBrowse:checkGenomeConfig
//input: file file
//output: view v
export async function previewGenomeFileBrowse(fileData: DG.FileInfo): Promise<DG.View> {
  let view = await genomeFileBrowse( await fileData.readAsString(), fileData.name);
  grok.shell.addView(view);
  return view;
}


//name: genomeFileBrowseHadnler
//input: string bytes
//output: list tables
//tags: file-handler
//meta.fileViewerCheck: GenomeFileBrowse:checkGenomeConfig
//meta.ext: json
export async function genomeFileBrowseHadnler(data: string) {
  let view = await genomeFileBrowse(data);
  grok.shell.addView(view);
  return [];
}


//name: checkGenomeConfig
//input: string content
//output: bool result
export function checkGenomeConfig(content: string) { 
  if(content.length < 1000000)
    return false;
  let jsonObj = JSON.parse(content);
  const uriRegex = new RegExp('"uri"');
  const uriHttpRegex = new RegExp('"uri"\\s*:\\s*"http');
  if(!!(jsonObj.tracks && (jsonObj.assemblies || jsonObj.assembly))){//megabyte and content checks 
    if (content.match(uriRegex)?.length === content.match(uriHttpRegex)?.length)
      return true;
    grok.shell.warning('Can not load the file as Genome File Browse viewer due there is uri existing without HTTP');
  }
  return false;
}
