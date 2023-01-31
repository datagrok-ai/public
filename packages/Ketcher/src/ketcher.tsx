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
import { KETCHER_MOLV2000, KETCHER_MOLV3000 } from "./constants";

let sketcherId = 0;

export class KetcherSketcher extends grok.chem.SketcherBase {

  _smiles: string | null = null;
  _molV2000: string | null = null;
  _molV3000: string | null = null;
  _smarts: string | null = null;
  _sketcher: Ketcher | null = null;

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
          this._smiles = await this._sketcher!.getSmiles();
          this._molV2000 = await this._sketcher!.getMolfile(KETCHER_MOLV2000);
          this._molV3000 = await this._sketcher!.getMolfile(KETCHER_MOLV3000);
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
  }

  get supportedExportFormats() {
    return ["smiles", "mol", "smarts"];
  }

  get smiles() {
    return this._smiles === null ?
    this._molV2000 !== null ? DG.chem.convert(this._molV2000, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles) :
    this._molV3000 !== null ? DG.chem.convert(this._molV3000, DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smiles) :
    this._smarts !== null ? DG.chem.smilesFromSmartsWarning() : '' : this._smiles;
  }

  set smiles(smiles: string) {
    this._smiles = smiles;
    this._molV2000 = null;
    this._molV3000 = null;
    this._smarts = null;
    this.setKetcherMolecule(smiles);
  }

  get molFile() {
    return this._molV2000 === null ?
      this._molV3000 !== null ? DG.chem.convert(this._molV3000, DG.chem.Notation.V3KMolBlock, DG.chem.Notation.MolBlock) :
      this._smiles !== null ? DG.chem.convert(this._smiles, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock) :
      this._smarts !== null ? DG.chem.convert(this._smarts, DG.chem.Notation.Smarts, DG.chem.Notation.MolBlock) : '' : this._molV2000;
  }

  set molFile(molfile: string) {
    this._molV2000 = molfile;
    this._molV3000 = null;
    this._smarts = null;
    this._smiles = null
    this.setKetcherMolecule(molfile);
  }

  get molV3000() {
    return this._molV3000 === null ?
      this._molV2000 !== null ? DG.chem.convert(this._molV2000, DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock) :
      this._smiles !== null ? DG.chem.convert(this._smiles, DG.chem.Notation.Smiles, DG.chem.Notation.V3KMolBlock) :
      this._smarts !== null ? DG.chem.convert(this._smarts, DG.chem.Notation.Smarts, DG.chem.Notation.V3KMolBlock) : '' : this._molV3000;
  }

  set molV3000(molfile: string) {
    this._molV3000 = molfile;
    this._molV2000 = null;
    this._smarts = null;
    this._smiles = null;
    this.setKetcherMolecule(molfile);
  }

  async getSmarts(): Promise<string> {
    return this._sketcher ? await this._sketcher.getSmarts() : this._smarts ?? '';
  }

  set smarts(smarts: string) {
    this._smarts = smarts;
    this._smiles = null;
    this._molV2000 = null;
    this._molV3000 = null;
    this.setKetcherMolecule(smarts);
  }

  get isInitialized() {
    return this._sketcher !== null;
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
