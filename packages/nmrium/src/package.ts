/* Do not change these import lines to match external modules in webpack configuration */
import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './bp.css';
import '@blueprintjs/icons/lib/css/blueprint-icons.css';
import { getNMRiumView, loadNMRiumData } from './utils';

export const _package = new DG.Package();
async function addNmriumView(extension: string, bytes: string) {
  const {v, nmriumId} = getNMRiumView();
  grok.shell.addView(v);
  await loadNMRiumData(`file.${extension}`, bytes, nmriumId);
  return [];
}

//name: jdxFileHandler
//input: string bytes
//output: list tables
//tags: file-handler
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: jdx
export async function jdxFileHandler(bytes: string) {
  return await addNmriumView('jdx', bytes);
}

//name: jdxFileHandler
//input: string bytes
//output: list tables
//tags: file-handler
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: dx
export async function dxFileHandler(bytes: string) {
  return await addNmriumView('dx', bytes);
}

//name: nmriumFileHandler
//input: string bytes
//output: list tables
//tags: file-handler
//meta.ext: nmrium
export async function nmriumFileHandler(bytes: string) {
  return await addNmriumView('nmrium', bytes);
}

//tags: fileViewer
//meta.fileViewer: nmrium
//input: file file
//output: view v
export function previewNMRData(fileData: DG.FileInfo): DG.View {
  const {v, nmriumId} = getNMRiumView();
  (async () => {
    const dataString = await fileData.readAsString();
    await loadNMRiumData(fileData.name, dataString, nmriumId);
  })();

  return v;
}

//tags: fileViewer
//meta.fileViewer: dx, jdx
//input: file file
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//output: view v
export function previewNMRFromDX(fileData: DG.FileInfo): DG.View {
  return previewNMRData(fileData);
}

//name: checkNmriumJdx
//input: string content
//output: bool result
export function checkNmriumJdx(content: string) {
  return content.includes('NMR SPECTRUM');
}


