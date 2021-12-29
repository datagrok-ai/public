import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import { MiscMethods } from "./misc.js"
import { PdbEntry } from '../pdb-entry.js';


// export let _package = new DG.Package();

export class NglAspect {

  stage: any;
  entry: PdbEntry;
  colorScheme: any;
  nglHost: HTMLElement;
  repChoice: DG.InputBase;
  selection: any;

  public async init(view: DG.TableView, entry: PdbEntry, colorScheme: {}, nglHost: HTMLElement, repChoice: DG.InputBase, twinSelections: {}) {

    this.colorScheme = colorScheme;
    let col_background = this.colorScheme["col_background"];
    let col_helix = this.colorScheme[""]

    nglHost.style.backgroundColor = col_background;
    view.box = true;
    this.entry = entry;
    this.repChoice = repChoice;
    this.selection = twinSelections;
    //@ts-ignore
    this.stage = new NGL.Stage(nglHost);
    this.render(true);
    this.nglResize(nglHost);
  }

  public async render(reload: boolean) {
    let switchObj = this.selection;

    let col_chain = this.colorScheme["col_chain"];
    let col_helix = this.colorScheme["col_helix"];
    let col_highlight = this.colorScheme["col_highlight"];

    //highlights in NGL
    let scheme_buffer: string[][] = [];

    Object.keys(switchObj).forEach((keyChain) => {
      Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

        let index = keyFtStart;
        if (switchObj[keyChain][keyFtStart]['state'] === true) {
          scheme_buffer.push([col_highlight, `${index} and :${keyChain}`]);
        }
      });
    });

    let schemeId;

    this.entry.entities[0].chains.forEach(chain => {
      let str_buffer = '';
      Object.keys(chain.tracks).forEach(track => {
        if(track === "HELIX_P"){
          let index1 = chain.tracks[track][0];
          let index2 = chain.tracks[track].slice(-1);
          str_buffer = str_buffer + ` or ${index1}-${index2} and :${chain.id}`;
        }
      });

      scheme_buffer.push([col_helix, str_buffer]);
      scheme_buffer.push([col_chain, `* and :${chain.id}`]);
    

    });

    //@ts-ignore
    schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);



    if (reload) {
      this.stage.removeAllComponents();
      await this.loadPdb(this.entry.body, this.repChoice, { color: schemeId });

      //recovering ball and stick residual at cartoon view
      if (this.repChoice.value === "cartoon") {
        Object.keys(switchObj).forEach((keyChain) => {
          Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

            let index = keyFtStart;

            let selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
            let addedRepresentation = this.stage.compList[0].addRepresentation("ball+stick", { sele: selection });
            switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;

            if (switchObj[keyChain][keyFtStart]['state'] !== true) {
              addedRepresentation.setVisibility(false);
            }
          });
        });
      }
    } else {
      this.stage.compList[0].addRepresentation(this.repChoice.value, { color: schemeId });

      if (this.repChoice.value === "cartoon") {
        Object.keys(switchObj).forEach((keyChain) => {
          Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

            let index = keyFtStart;

            if (typeof switchObj[keyChain][keyFtStart]['representation'] === 'undefined') {
              let selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
              let addedRepresentation = this.stage.compList[0].addRepresentation("ball+stick", { sele: selection });
              switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;
            }

            if (switchObj[keyChain][keyFtStart]['state'] !== true) {
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(false);
            } else {
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(true);
            }
          });
        });
      }
    }
  }

  // load the 3D model
  private async loadPdb(pdbStr: string, repChoice:DG.InputBase, schemeObj: any) {
    var stringBlob = new Blob([pdbStr], { type: 'text/plain' });
    await this.stage.loadFile(stringBlob, { ext: "pdb" }).then(function (o: any) {
      o.addRepresentation(repChoice.value, schemeObj);
      o.autoView();
    });
  }

  // viewer resize
  private _resize(host: HTMLElement) {
    let canvas = host.querySelector('canvas');
    canvas!.width = Math.floor(host.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(host.clientHeight * window.devicePixelRatio);
    this.stage.handleResize();
  }

  private nglResize(host: HTMLElement) {
    ui.onSizeChanged(host).subscribe((_) => this._resize(host));
    this._resize(host);
  }
}