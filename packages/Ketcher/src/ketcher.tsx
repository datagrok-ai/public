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
import { chem } from "datagrok-api/grok";

let sketcherId = 0;

export class KetcherSketcher extends grok.chem.SketcherBase {

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
        this._sketcher = ketcher;
        (this._sketcher.editor as any).subscribe("change", async (_: any) => {
          this.host!._smiles = await this._sketcher!.getSmiles();
          this.host!._molfile = this.host!.molFileUnits === grok.chem.MOLV2000 ?
            await this._sketcher!.getMolfile('v2000'): await this._sketcher!.getMolfile('v3000');
          this.onChanged.next(null);
        });
        this._sketcher.editor.zoom(0.5);
      },
    };

    let host = ui.div([], { style: { width: "500px", height: "400px" } });
    host.style.setProperty('overflow', 'hidden', 'important');

    let component = React.createElement(Editor, props, null);
    ReactDOM.render(component, host);

    this.root.appendChild(host);
  }

  async init(host: chem.Sketcher) {
    this.host = host;
    let id = `ketcher-${sketcherId++}`;
    this.root.id = id;
    this.onChanged.next(null);
  }

  get supportedExportFormats() {
    return ["smiles", "mol", "smarts"];
  }

  get smiles() {
    return this._sketcher ? this.host!._smiles : this.host!.getSmiles();
  }

  set smiles(smiles) {
    this.setKetcherMolecule(smiles);
  }

  get molFile() {
    return this._sketcher ? this.host!._molfile : this.host!.getMolFile();
  }

  set molFile(molfile: string) {
    this.setKetcherMolecule(molfile);
  }

  get molV3000() {
    return this._sketcher ? this.host!._molfile : this.host!.getMolFile();
  }

  set molV3000(molfile: string) {
    this.setKetcherMolecule(molfile);
  }

  async getSmarts(): Promise<string> {
    return this._sketcher ? await this._sketcher.getSmarts() : await this.host!.getSmarts();
  }

  set smarts(smarts: string) {
    this.setKetcherMolecule(smarts);
  }

  setKetcherMolecule(molecule: string) {
    try {
      this._sketcher?.setMolecule(molecule);
    } catch (e) {
      console.log(e);
      return;
    }
  }

  detach() {
    super.detach();
  }

}
