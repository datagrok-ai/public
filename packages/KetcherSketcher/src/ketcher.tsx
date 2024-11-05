import * as React from "react";
import * as ReactDOM from 'react-dom/client';
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import * as DG from "datagrok-api/dg";
import {_package} from "./package";
import {Editor} from "ketcher-react";
import * as PolymerKetcher from "ketcher-macromolecules";
import {StandaloneStructServiceProvider} from "ketcher-standalone";
import {Ketcher} from "ketcher-core";
import "ketcher-react/dist/index.css";
import "../css/editor.css";
import { chem } from "datagrok-api/grok";
import { KETCHER_MOLV2000, KETCHER_MOLV3000, KETCHER_WINDOW_OBJECT } from "./constants";
import { ketcherMacromoleculeMode, setShowPolymerEditor, ToggleControl } from "./modeControl";
import { Subscription } from "rxjs";

let sketcherId = 0;

export class KetcherSketcher extends grok.chem.SketcherBase {

  _smiles: string | null = null;
  _molV2000: string | null = null;
  _molV3000: string | null = null;
  _smarts: string | null = null;
  _sketcher: Ketcher | null = null;
  ketcherHost: HTMLDivElement;
  componentKetcher;
  componentPolymerKetcher;
  togglerComponent;
  ketcherModeSub: Subscription;
  switchedMode = false;
  compRoot;

  constructor() {
    super();
    let structServiceProvider = new StandaloneStructServiceProvider();

    this.togglerComponent = React.createElement(ToggleControl);
    this.ketcherModeSub = ketcherMacromoleculeMode.subscribe((value: boolean) => {
      this.switchedMode = true;
      value ? this.compRoot.render(this.componentPolymerKetcher) : this.compRoot.render(this.componentKetcher);
    })

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
        // grok.dapi.userDataStorage.getValue(KETCHER_OPTIONS, KETCHER_USER_STORAGE, true).then((opts: string) => {
        //   if (opts) {
        //     this._sketcher?.editor.setOptions(opts);
        //   }
        // });
        //@ts-ignore
        window[KETCHER_WINDOW_OBJECT] = ketcher;
        this.setMoleculeFromHost();
        (this._sketcher.editor as any).subscribe("change", async (_: any) => {
          try {
            this._smiles = await this._sketcher!.getSmiles();
          } catch { //in case we are working with smarts - getSmiles() will fail with exception
            this._smiles = null;
          }
          this._molV2000 = await this._sketcher!.getMolfile(KETCHER_MOLV2000);
          this._molV3000 = await this._sketcher!.getMolfile(KETCHER_MOLV3000);
          this.onChanged.next(null);
        });
      },
      togglerComponent: this.togglerComponent
    };

    this.ketcherHost = ui.div([], 'ketcher-host');

    this.componentKetcher = React.createElement(Editor, props, null);
    this.compRoot = ReactDOM.createRoot(this.ketcherHost);   
    this.componentPolymerKetcher = React.createElement(PolymerKetcher.Editor, {
      togglerComponent: this.togglerComponent
    });
    this.compRoot.render(this.componentKetcher);
    this.root.appendChild(this.ketcherHost);
  }

  resize() {
    if (this.switchedMode) {
      this.switchedMode = false;
      return;
    }
    if (!this.ketcherHost.classList.contains('ketcher-resizing'))
      this.ketcherHost.classList.add('ketcher-resizing');
  }

  async init(host: chem.Sketcher) {
    this.host = host;
    let id = `ketcher-${sketcherId++}`;
    this.root.id = id;
    if (this.host.isResizing)
      this.ketcherHost.classList.add('ketcher-resizing');
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

  setMoleculeFromHost(): void {
    if (this.host) {
      if (this.host!._molfile !== null) {
        if (this.host!.molFileUnits === chem.Notation.MolBlock)
          this.molFile = this.host!._molfile;
        if (this.host!.molFileUnits === chem.Notation.V3KMolBlock)
          this.molV3000 = this.host!._molfile;
        return;
      }
      if (this.host!._smiles !== null) {
        this.smiles = this.host!._smiles;
        return;
      }
      if (this.host!._smarts !== null) {
        this.smarts = this.host!._smarts;
        return;
      }
    }
  }

  detach() {
   // grok.dapi.userDataStorage.postValue(KETCHER_OPTIONS, KETCHER_USER_STORAGE, JSON.stringify(this._sketcher?.editor.options()), true);
    super.detach();
    this.compRoot.unmount();
    this.ketcherModeSub.unsubscribe();
    setShowPolymerEditor(false);
    
  }

}
