import {chem} from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as DG from 'datagrok-api/dg';
import {RemoteStructServiceProvider} from "ketcher-core";
import React from "react";
import {Editor} from "ketcher-react";
import ReactDOM from "react-dom";
import {_package} from "./package";

const structServiceProvider = new RemoteStructServiceProvider(
  'foo_api_path'
  //process.env.REACT_APP_API_PATH!
)

export class KetcherSketcher extends chem.SketcherBase {
  constructor() {
    super();

    let props = {
      //staticResourcesUrl: process.env.PUBLIC_URL,
      staticResourcesUrl: _package.webRoot!,
      structServiceProvider: structServiceProvider
    };

    // Andrew: this works but I don't know how to set smiles, or receive events
    let host = ui.div([], { style: {height: '500px', width: '500px'}});
    let component = React.createElement(Editor, props, null);
    ReactDOM.render(component, host);

    this.root.appendChild(host);

    // Andrew: this seems to be what we need but I can't get it to build
    // let ketcher = await buildKetcherAsync({
    //   element: host,
    //   structServiceProvider: structServiceProvider,
    //   staticResourcesUrl: _package.webRoot!
    // });
    //ketcher.setMolecule('Fc1cc(Cl)ccc1Br');
  }
}