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
import {_jsThen} from "datagrok-api/src/utils";

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
      errorHandler: (message: string) => {
        console.log("Skecther error", message);
      },
      onInit: (ketcher: Ketcher) => {
        this._ketcher = ketcher;
        (this._ketcher.editor as any).subscribe("change", async (e: any) => {
          await this.setSmilesSmartsMolfile();
          this.onChanged.next(null);
        });
        this._ketcher.editor.zoom(0.5);
      },
    };

    let host = ui.div([], { style: { width: "500px", height: "400px" } });
    host.style.setProperty('overflow', 'hidden', 'important');

    let component = React.createElement(Editor, props, null);
    ReactDOM.render(component, host);

    this.root.appendChild(host);
  }

  async init() {
    let id = `ketcher-${sketcherId++}`;
    this.root.id = id;
    this.onChanged.next(null);
  }

  get supportedExportFormats() {
    return ["smiles", "mol", "smarts"];
  }

  get smiles() {
    return this._smiles;
  }

  set smiles(smiles) {
    this.setSmiles(smiles);
  }

  async getSmiles(): Promise<string> {
    return await this._ketcher?.getSmiles();
  }

  async setSmiles(smiles: string) {
    this.setKetcherMolecule(smiles);
  }

  async getMolFile(): Promise<string> {
    return await this._ketcher?.getMolfile();
  }

  get molFile() {
    return this._molFile;
  }

  set molFile(molfile: string) {
    this.setMolFile(molfile);
  }

  async setMolFile(molfile: string) {
    this.setKetcherMolecule(molfile);
  }

  get smarts() {
    return this._smarts;
  }

  set smarts(smarts: string) {
    this.setSmarts(smarts);
  }

  async getSmarts(): Promise<string> {
    return this._smarts;
  }

  async setSmarts(smarts: string) {
    this.setKetcherMolecule(smarts);
  }

  setKetcherMolecule(molecule: string) {
    try {
      this._ketcher?.setMolecule(molecule);
    } catch (e) {
      console.log(e);
      return;
    }
  }

  async setSmilesSmartsMolfile(){
    this._smiles = await this.getSmiles();
    this._molFile = await this.getMolFile();
    this._smarts = await this.getSmarts();
  }

  detach() {
    super.detach();
  }

}
