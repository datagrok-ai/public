/* Do not change these import lines to match external modules in webpack configuration */
import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getNMRiumComponent} from "../nmrium-wrapper/nmrium-react-wrapper/src/NMRiumWrapper";
import {NMRiumEvents} from "../nmrium-wrapper/nmrium-react-wrapper/src/NMRiumWrapper";
import '@blueprintjs/core/lib/css/blueprint.css';
import '@blueprintjs/icons/lib/css/blueprint-icons.css';

export const _package = new DG.Package();

//name: info
//tags: autostart
export function info() {
  grok.shell.info(_package.webRoot);
  // @ts-ignore
  window.NMRiumEvents = NMRiumEvents;
}


//name: addNmriumView
export function addNmriumView() {
  const root = ui.div();
  root.style.width = '100%';
  root.style.height = '100%';
  const v = DG.View.fromRoot(root);
  v.name = 'NMRium';
  grok.shell.addView(v);


  // function getNMRiumComponent(root: HTMLElement) {
  //   const props = {};
  //   const component = React.createElement(NMRiumWrapper, props, null);
  //   const root1 = ReactDOM.createRoot(root);
  //   root1.render(component);
  //   return root;
  // }
  // so that at the time of creation nmriumHost is already available and has some size
  setTimeout(() => getNMRiumComponent(root as any), 1000);
}