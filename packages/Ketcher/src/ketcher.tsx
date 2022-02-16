import * as React from "react";
import * as ReactDOM from "react-dom";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import * as DG from "datagrok-api/dg";
import {_package} from "./package";
import {Editor} from "ketcher-react";
import {StandaloneStructServiceProvider} from "ketcher-standalone";
import {Ketcher} from "ketcher-core";
import "ketcher-react/dist/index.css";
import "./editor.css";

let sketcherId = 0;

export class KetcherSketcher extends grok.chem.SketcherBase {
  declare _ketcher: Ketcher;

  constructor() {
    super();
    let structServiceProvider = new StandaloneStructServiceProvider();
 
    let props = {
      staticResourcesUrl: !_package.webRoot
        ? ""
        : _package.webRoot.substring(0, _package.webRoot.length - 1),
      structServiceProvider: structServiceProvider,
      errorHandler: (message: string) => 
      {
        console.log('Skecther error', message);
      },
      onInit: (ketcher: Ketcher) => {
        this._ketcher = ketcher;
      },
    };

    let host = ui.div([], { style: { width: "fit-content" } });
 
    let component = React.createElement(Editor, props, null);
    ReactDOM.render(component, host);
    let sketcherConentDiv =  document.querySelectorAll("div.ui-div > div.grok-sketcher.ui-box");
    if (sketcherConentDiv[0])
    {
      sketcherConentDiv[0].setAttribute("style", "width: fit-content; height: fit-content;");
    }
    this.root.appendChild(host);
  }

  async init() {
    let id = `ketcher-${sketcherId++}`;
    this.root.id = id;
    this.onChanged.next(null);
  }
 
  async smilesToMol(smiles: string) {
    let molBlock: string = await grok.functions.call("Chem:convertMolecule", { 'molecule': smiles, 'from': 'smiles', 'to': 'molblock'});
    return molBlock;
  }

  detach() {
    super.detach();
  }

  async getSmiles(): Promise<string> {
    return await this._ketcher?.getSmiles();
  }

  setSmiles(smiles: string) {
    this.smilesToMol(smiles).then((molBlock) => {
      try {
        this._ketcher?.setMolecule(molBlock).then(() => {});
      } catch (e) {
        console.log(e);
        return; // fallback value
      }
    });
  }

  async getMolFile(): Promise<string> {
     return await this._ketcher?.getMolfile();
  }

  setMolFile(molfile: string) {
      try {
        this._ketcher?.setMolecule(molfile).then(()=>{})
      } catch (e) {
        console.log(e);
        return; // fallback value
    }
  }
}
