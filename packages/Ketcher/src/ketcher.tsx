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
  declare _molFile: string;

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
        (this._ketcher.editor as any).subscribe("change", (e: any) => {

          try {
            this._ketcher.getMolfile().then((molfile) => {
              this._molFile = molfile;
            });
          } catch (ex) {
            console.log(ex);
          }
          this.onChanged.next(null);
        });
      },
    };

    let host = ui.div([], { style: { width: "700px", height: "500px" } });

    let component = React.createElement(Editor, props, null);
    ReactDOM.render(component, host);
    let sketcherConentDiv = document.querySelectorAll(
      "div.ui-div > div.grok-sketcher.ui-box"
    );
    if (sketcherConentDiv[0]) {
      sketcherConentDiv[0].setAttribute(
        "style",
        "width: fit-content; height: fit-content;"
      );
    }
    this.root.appendChild(host);
  }

  async init() {
    this._molFile = "";
    let id = `ketcher-${sketcherId++}`;
    this.root.id = id;
    this.onChanged.next(null);
  }

  async smilesToMol(smiles: string): Promise<string> {
    return await grok.functions.call("Chem:convertMolecule", {
      molecule: smiles,
      from: "smiles",
      to: "molblock",
    });
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
        return;
      }
    });
  }

  async getMolFile(): Promise<string> {
    return await this._ketcher?.getMolfile();
  }

  setMolFile(molfile: string) {
    this._molFile = molfile;
    try {
      this._ketcher?.setMolecule(molfile).then(() => {});
    } catch (e) {
      console.log(e);
      return;
    }
  }

  get supportedExportFormats() {
    return ["smiles", "mol"];
  }

  get smiles() {
    return this._molFile;
  }

  set smiles(smiles) {
    this.smilesToMol(smiles).then((molFile) => this.setMolFile(molFile));
  }

  get molFile() {
    return this._molFile;
  }

  set molFile(molfile: string) {
    this.setMolFile(molfile);
  }
}
