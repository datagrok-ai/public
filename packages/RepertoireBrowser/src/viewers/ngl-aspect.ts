import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import { MiscMethods } from "./misc"

export class NglAspect {

  stage: any;
  schemeObj: any;
  pdbStr: string;
  json: any;
  colorScheme: any;
  nglHost: HTMLElement;
  repChoice: DG.InputBase;
  cdrScheme: DG.InputBase;
  paratopes: DG.InputBase;
  selection: any;

  public async init(view: DG.TableView, pdbStr: string, json: any, colorScheme: any, nglHost: HTMLElement,
    repChoice: DG.InputBase, cdrScheme: DG.InputBase, paratopes: DG.InputBase, twinSelections: any) {

    this.colorScheme = colorScheme;
    let col_background = this.colorScheme["col_background"];

    nglHost.style.backgroundColor = col_background;
    view.box = true;
    this.pdbStr = pdbStr;
    this.json = json;
    this.repChoice = repChoice;
    this.cdrScheme = cdrScheme;
    this.paratopes = paratopes;
    this.selection = twinSelections;
    //@ts-ignore
    this.stage = new NGL.Stage(nglHost);
    this.render(true);
    this.nglResize(nglHost);
  }

  public async render(reload: boolean) {
    let switchObj = this.selection;

    let col_heavy_chain = this.colorScheme["col_heavy_chain"];
    let col_light_chain = this.colorScheme["col_light_chain"];
    let col_cdr = this.colorScheme["col_cdr"];
    let col_para = this.colorScheme["col_para"];
    let col_partopes_low = this.colorScheme["col_partopes_low"]; //col_para in rgb
    let col_partopes_high = this.colorScheme["col_partopes_high"];
    let col_highlight = (this.cdrScheme.value === 'default' || this.paratopes.value === true) ?
      this.colorScheme["col_highlight"] : this.colorScheme["col_highlight_cdr"];

    //highlights in NGL
    let scheme_buffer: string[][] = [];

    Object.keys(switchObj).forEach((keyChain) => {
      Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

        let index = keyChain === "H" ? this.json.map_H[keyFtStart] : this.json.map_L[keyFtStart];
        if (switchObj[keyChain][keyFtStart]['state'] === true) {
          scheme_buffer.push([col_highlight, `${index} and :${keyChain}`]);
        }
      });
    });

    let schemeId;
    if (this.paratopes.value === true) {
      let palette = MiscMethods.interpolateColors(col_partopes_low, col_partopes_high, 100);
      Object.keys(this.json.parapred_predictions).forEach((chain) => {
        Object.keys(this.json.parapred_predictions[chain]).forEach((index) => {

          let nindex = chain === "H" ? this.json.map_H[index] : this.json.map_L[index];

          scheme_buffer.push([
            palette[Math.round(this.json.parapred_predictions[chain][index] * 100)],
            `${nindex} and :${chain}`
          ]);
        })

      })
      scheme_buffer.push([col_para, "* and :H"]);
      scheme_buffer.push([col_para, "* and :L"]);
      //@ts-ignore
      schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
    }
    else {
      if (this.cdrScheme.value === 'default') {
        scheme_buffer.push([col_heavy_chain, "* and :H"]);
        scheme_buffer.push([col_light_chain, "* and :L"]);
        //@ts-ignore
        schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
      } else {
        Object.keys(this.json.cdr_ranges).forEach((str) => {
          if (str.includes(this.cdrScheme.value + '_CDRH')) {
            let str_buffer = '';
            for (let i = 0; i < Object.keys(this.json.cdr_ranges[str]).length; i++) {
              let nindex1 = this.json.map_H[this.json.cdr_ranges[str][i][0]];
              let nindex2 = this.json.map_H[this.json.cdr_ranges[str][i][1]];

              str_buffer = str_buffer + ` or ${nindex1}-${nindex2} and :H`;
            }
            str_buffer = str_buffer.slice(4);
            scheme_buffer.push([col_cdr, str_buffer]);
            scheme_buffer.push([col_heavy_chain, "* and :H"]);

          } else if (str.includes(this.cdrScheme.value + '_CDRL')) {
            let str_buffer = ''
            for (let i = 0; i < Object.keys(this.json.cdr_ranges[str]).length; i++) {
              let nindex1 = this.json.map_L[this.json.cdr_ranges[str][i][0]];
              let nindex2 = this.json.map_L[this.json.cdr_ranges[str][i][1]];

              str_buffer = str_buffer + ` or ${nindex1}-${nindex2} and :L`;
            }
            str_buffer = str_buffer.slice(4);
            scheme_buffer.push([col_cdr, str_buffer]);
            scheme_buffer.push([col_light_chain, "* and :L"]);
          }
        });
        //@ts-ignore
        schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
      }
    }

    if (reload) {
      this.stage.removeAllComponents();
      await this.loadPdb(this.pdbStr, this.repChoice, { color: schemeId });

      //recovering ball and stick residual at cartoon view
      if (this.repChoice.value === "cartoon") {
        Object.keys(switchObj).forEach((keyChain) => {
          Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

            let index = keyChain === "H" ? this.json.map_H[keyFtStart] : this.json.map_L[keyFtStart];

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

            let index = keyChain === "H" ? this.json.map_H[keyFtStart] : this.json.map_L[keyFtStart];

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
  private async loadPdb(pdbStr: string, repChoice: DG.InputBase, schemeObj: any) {
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