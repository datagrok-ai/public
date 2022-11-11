import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MiscMethods} from './misc';
import {PdbType} from '../utils/data-loader';

export class NglAspect {
  stage: any;
  schemeObj: any;
  pdbStr: PdbType;
  json: any;
  colorScheme: any;
  nglHost: HTMLElement;
  repChoice: DG.InputBase;
  cdrScheme: DG.InputBase;
  paratopes: DG.InputBase;
  selection: any;

  public async init(
    view: DG.TableView, pdbStr: PdbType, json: any, colorScheme: any, nglHost: HTMLElement,
    repChoice: DG.InputBase, cdrScheme: DG.InputBase, paratopes: DG.InputBase, twinSelections: any
  ) {
    this.colorScheme = colorScheme;
    nglHost.style.backgroundColor = this.colorScheme['colBackground'];
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
    const switchObj = this.selection;

    const colHeavyChain = this.colorScheme['colHeavyChain'];
    const colLightChain = this.colorScheme['colLightChain'];
    const colCdr = this.colorScheme['colCdr'];
    const colPara = this.colorScheme['colPara'];
    const colParatopesLow = this.colorScheme['colParatopesLow']; //col_para in rgb
    const colParatopesHigh = this.colorScheme['colParatopesHigh'];
    const colHighlight = (this.cdrScheme.value === MiscMethods.NoSchemeItem || this.paratopes.value === true) ?
      this.colorScheme['colHighlight'] : this.colorScheme['colHighlightCdr'];

    //highlights in NGL
    const schemeBuffer: string[][] = [];

    Object.keys(switchObj).forEach((keyChain) => {
      Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {
        const index = keyChain === 'H' ? this.json.map_H[keyFtStart] : this.json.map_L[keyFtStart];
        if (switchObj[keyChain][keyFtStart]['state'] === true)
          schemeBuffer.push([colHighlight, `${index} and :${keyChain}`]);
      });
    });

    let schemeId;
    if (this.paratopes.value === true) {
      const palette = MiscMethods.interpolateColors(colParatopesLow, colParatopesHigh, 100);
      Object.keys(this.json.parapred_predictions).forEach((chain) => {
        Object.keys(this.json.parapred_predictions[chain]).forEach((index) => {
          const nindex = chain === 'H' ? this.json.map_H[index] : this.json.map_L[index];
          schemeBuffer.push([
            palette[Math.round(this.json.parapred_predictions[chain][index] * 100)],
            `${nindex} and :${chain}`]);
        });
      });
      schemeBuffer.push([colPara, '* and :H']);
      schemeBuffer.push([colPara, '* and :L']);
      //@ts-ignore
      schemeId = NGL.ColormakerRegistry.addSelectionScheme(schemeBuffer);
    } else {
      if (this.cdrScheme.value === MiscMethods.NoSchemeItem) {
        schemeBuffer.push([colHeavyChain, '* and :H']);
        schemeBuffer.push([colLightChain, '* and :L']);
        //@ts-ignore
        schemeId = NGL.ColormakerRegistry.addSelectionScheme(schemeBuffer);
      } else {
        Object.keys(this.json.cdr_ranges).forEach((str) => {
          if (str.includes(this.cdrScheme.value + '_CDRH')) {
            let strBuffer = '';
            for (let i = 0; i < Object.keys(this.json.cdr_ranges[str]).length; i++) {
              const nindex1 = this.json.map_H[this.json.cdr_ranges[str][i][0]];
              const nindex2 = this.json.map_H[this.json.cdr_ranges[str][i][1]];

              strBuffer = strBuffer + ` or ${nindex1}-${nindex2} and :H`;
            }
            strBuffer = strBuffer.slice(4);
            schemeBuffer.push([colCdr, strBuffer]);
            schemeBuffer.push([colHeavyChain, '* and :H']);
          } else if (str.includes(this.cdrScheme.value + '_CDRL')) {
            let strBuffer = '';
            for (let i = 0; i < Object.keys(this.json.cdr_ranges[str]).length; i++) {
              const nindex1 = this.json.map_L[this.json.cdr_ranges[str][i][0]];
              const nindex2 = this.json.map_L[this.json.cdr_ranges[str][i][1]];

              strBuffer = strBuffer + ` or ${nindex1}-${nindex2} and :L`;
            }
            strBuffer = strBuffer.slice(4);
            schemeBuffer.push([colCdr, strBuffer]);
            schemeBuffer.push([colLightChain, '* and :L']);
          }
        });
        //@ts-ignore
        schemeId = NGL.ColormakerRegistry.addSelectionScheme(schemeBuffer);
      }
    }

    if (reload || this.stage.compList.length === 0) {
      this.stage.removeAllComponents();
      await this.loadPdb(this.pdbStr, this.repChoice, {color: schemeId});

      //recovering ball and stick residual at cartoon view
      if (this.repChoice.value === 'cartoon') {
        Object.keys(switchObj).forEach((keyChain) => {
          Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {
            const index = keyChain === 'H' ? this.json.map_H[keyFtStart] : this.json.map_L[keyFtStart];

            const selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
            const addedRepresentation = this.stage.compList[0].addRepresentation('ball+stick', {sele: selection});
            switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;

            if (switchObj[keyChain][keyFtStart]['state'] !== true)
              addedRepresentation.setVisibility(false);
          });
        });
      }
    } else {
      this.stage.compList[0].addRepresentation(this.repChoice.value, {color: schemeId});

      if (this.repChoice.value === 'cartoon') {
        Object.keys(switchObj).forEach((keyChain) => {
          Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {
            const index = keyChain === 'H' ? this.json.map_H[keyFtStart] : this.json.map_L[keyFtStart];

            if (typeof switchObj[keyChain][keyFtStart]['representation'] === 'undefined') {
              const selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
              const addedRepresentation = this.stage.compList[0].addRepresentation('ball+stick', {sele: selection});
              switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;
            }

            if (switchObj[keyChain][keyFtStart]['state'] !== true)
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(false);
            else
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(true);
          });
        });
      }
    }
  }

  // load the 3D model
  private async loadPdb(pdbStr: PdbType, repChoice: DG.InputBase, schemeObj: any) {
    const stringBlob = new Blob([<string>pdbStr], {type: 'text/plain'});
    await this.stage.loadFile(stringBlob, {ext: 'pdb'}).then(function(o: any) {
      o.addRepresentation(repChoice.value, schemeObj);
      o.autoView();
    });
  }

  // viewer resize
  private _resize(host: HTMLElement) {
    const canvas = host.querySelector('canvas');
    canvas!.width = Math.floor(host.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(host.clientHeight * window.devicePixelRatio);
    this.stage.handleResize();
  }

  private nglResize(host: HTMLElement) {
    ui.onSizeChanged(host).subscribe((_) => this._resize(host));
    this._resize(host);
  }
}
