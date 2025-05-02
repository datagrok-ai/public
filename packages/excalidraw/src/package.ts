/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { renderExcalidraw } from './react-thing';
import './excal.css';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Excalidraw
//tags: app
//meta.browsePath: Misc
//output: view v
export function excalidrawApp() {
  const div = renderExcalidraw();
  const v= DG.View.fromRoot(div.actualRoot);
  v.name = 'Excalidraw';
  return v;
}

//name: excalfileViewer
//description: Excalidraw viewer
//tags: fileViewer
//meta.fileViewer: excalidraw
//meta.ext: excalidraw
//input: file fileContent
//output: view v
export async function excalfileViewer(fileContent: DG.FileInfo) {
  if (!fileContent) {
    grok.shell.info('No file content provided');
    return;
  }
  const fileC = await fileContent.readAsString();
  if (!fileC) {
    grok.shell.info('No file content provided');
    return;
  }
  const json = JSON.parse(fileC);
  if (!json) {
    grok.shell.info('No JSON content provided');
    return;
  }
  const div = renderExcalidraw();
  const v = DG.View.fromRoot(div.actualRoot);
  v.name = fileContent.name ?? 'Excalidraw';
  setTimeout(() => {
    const api = div.getApi();
    api.updateScene(json);

  }, 500);
  return v;
}