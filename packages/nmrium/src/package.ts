/* Do not change these import lines to match external modules in webpack configuration */
import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getNMRiumComponent } from "../nmrium-wrapper/nmrium-react-wrapper/src/NMRiumWrapper";
import { NMRiumEvents } from "../nmrium-wrapper/nmrium-react-wrapper/src/NMRiumWrapper";
import '@blueprintjs/core/lib/css/blueprint.css';
import '@blueprintjs/icons/lib/css/blueprint-icons.css';
import { delay } from '@datagrok-libraries/utils/src/test';

export const _package = new DG.Package();

// //name: jdxFileHandler
// //input: list bytes
// //output: list tables
// //tags: file-handler
// //meta.ext: jdf
// export async function jdfFileHandler(bytes: string) {
//   let root = await addNmriumView();
//   const blob = new Blob([bytes], { type: 'text/plain' });
//   const file = new File([blob], "file.jdf", { type: "text/plain" });
//   await delay(100);
//   NMRiumEvents.trigger('load', {
//     data: [file],
//     type: 'file', wrapper: root.children[0]?.id
//   });
//   return [];
// }

//name: jdxFileHandler
//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: jdx
export async function jdxFileHandler(bytes: string) {
  let root = await addNmriumView();
  const blob = new Blob([bytes], { type: 'text/plain' });
  const file = new File([blob], "file.jdx", { type: "text/plain" });
  await delay(100);
  NMRiumEvents.trigger('load', {
    data: [file],
    type: 'file', wrapper: root.children[0]?.id
  });
  return [];
}

//name: jdxFileHandler
//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: dx
export async function dxFileHandler(bytes: string) {
  let root = await addNmriumView();
  const blob = new Blob([bytes], { type: 'text/plain' });
  const file = new File([blob], "file.dx", { type: "text/plain" });
  await delay(100);
  NMRiumEvents.trigger('load', {
    data: [file],
    type: 'file', wrapper: root.children[0]?.id
  });
  return [];
}

//name: nmriumFileHandler
//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: nmrium
export async function nmriumFileHandler(bytes: string) {
  let root = await addNmriumView();
  const blob = new Blob([bytes], { type: 'text/plain' });
  const file = new File([blob], "file.nmrium", { type: "text/plain" });
  await delay(100);
  NMRiumEvents.trigger('load', {
    data: [file],
    type: 'file', wrapper: root.children[0]?.id
  });
  return [];
}

//tags: fileViewer
//meta.fileViewer: dx, jdx, nmrium
//input: file file
//output: view v
export function  previewNMRData (fileData: DG.FileInfo): DG.View {
  const root = ui.div();
  root.style.width = '100%';
  root.style.height = '100%';
  const v = DG.View.fromRoot(root);
  v.name = 'NMRium';
  getNMRiumComponent(root as any)

  setTimeout(async () => {
    const blob = new Blob([await fileData.readAsString()], { type: 'text/plain' });
    const file = new File([blob], fileData.name, { type: "text/plain" });
    NMRiumEvents.trigger('load', {
      data: [file],
      type: 'file', wrapper: root.children[0]?.id
    });
  }, 100)
  return v;
}

//name: addNmriumView
export async function addNmriumView() {
  const root = ui.div();
  root.style.width = '100%';
  root.style.height = '100%';
  const v = DG.View.fromRoot(root);
  v.name = 'NMRium';
  grok.shell.addView(v);

  getNMRiumComponent(root as any)
  return root;
}