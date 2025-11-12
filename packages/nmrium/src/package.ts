/* Do not change these import lines to match external modules in webpack configuration */
import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './bp.css';
import '@blueprintjs/icons/lib/css/blueprint-icons.css';
import { getNMRiumView, loadNMRiumData } from './utils';
export * from './package.g';
export const _package = new DG.Package();
async function addNmriumView(extension: string, bytes: string) {
  const {v, nmriumId} = getNMRiumView();
  grok.shell.addView(v);
  await loadNMRiumData(`file.${extension}`, bytes, nmriumId);
  return [];
}

export class PackageFunctions{
  @grok.decorators.fileHandler({
    fileViewerCheck: 'Nmrium:checkNmriumJdx',
    ext: 'jdx',
    outputs: [{name: 'tables', type: 'list'}]
  })
  static async jdxFileHandler(bytes: string) {
  
    return await addNmriumView('jdx', bytes);
  }


  @grok.decorators.fileHandler({
    fileViewerCheck: 'Nmrium:checkNmriumJdx',
    ext: 'dx',
    name: 'jdxFileHandler',
    outputs: [{name: 'tables', type: 'list'}]
  })
  static async dxFileHandler(bytes: string) {
  
    return await addNmriumView('dx', bytes);
  }


  @grok.decorators.fileHandler({ext: 'nmrium', outputs: [{name: 'tables', type: 'list'}]})
  static async nmriumFileHandler(bytes: string) {
    return await addNmriumView('nmrium', bytes);
  }


  @grok.decorators.fileViewer({fileViewer: 'nmrium'})
  static previewNMRData(file: DG.FileInfo): DG.View {
    const {v, nmriumId} = getNMRiumView();
    (async () => {
      const dataString = await file.readAsString();
      await loadNMRiumData(file.name, dataString, nmriumId);
    })();

    return v;
  }


  @grok.decorators.fileViewer({
    fileViewer: 'dx, jdx',
    fileViewerCheck: 'Nmrium:checkNmriumJdx',
  })
  static previewNMRFromDX(file: DG.FileInfo): DG.View {
    return PackageFunctions.previewNMRData(file);
  }


  @grok.decorators.func()
  static checkNmriumJdx (content: string) : boolean{
    return content.includes('NMR SPECTRUM');
  }
}

