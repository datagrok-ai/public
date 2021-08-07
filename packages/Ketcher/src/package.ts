/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Editor } from "ketcher-react";
import { RemoteStructServiceProvider } from 'ketcher-core'
import React from 'react';
import ReactDOM from 'react-dom';
import buildKetcherAsync from "ketcher-react/dist/script";

export let _package = new DG.Package();

const structServiceProvider = new RemoteStructServiceProvider(
  'foo_api_path'
  //process.env.REACT_APP_API_PATH!
)

//name: ketch
export async function ketch() {

  let props = {
    //staticResourcesUrl: process.env.PUBLIC_URL,
    staticResourcesUrl: _package.webRoot!,
    structServiceProvider: structServiceProvider
  };

  // Andrew: this works but I don't know how to set smiles, or receive events
  let host = ui.div([], { style: {height: '500px', width: '500px'}});
  let component = React.createElement(Editor, props, null);
  ReactDOM.render(component, host);

  // Andrew: this seems to be what we need but I can't get it to build
  // let ketcher = await buildKetcherAsync({
  //   element: host,
  //   structServiceProvider: structServiceProvider,
  //   staticResourcesUrl: _package.webRoot!
  // });
  //ketcher.setMolecule('Fc1cc(Cl)ccc1Br');

  ui.dialog()
    .add(host)
    .show();
}