import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
//@ts-ignore
import mutcodes from "../externalData/mutcodes.json";
import { MiscMethods } from "./misc.js"
import { PdbEntry } from '../pdb-entry';

export class PvizAspect {
  ngl: any;
  pviz: any;
  pVizParams: any;

  entry: PdbEntry;
  jsonObs: any;
  colorScheme: any;
  hosts: { [key: string]: HTMLDivElement; };
  selection: { [key: string]: {} };

  async init(entry: PdbEntry, colorScheme: any, pVizHosts: { [key: string]: HTMLDivElement; }, twinSelections: any) {

    //@ts-ignore
    this.pviz = window.pviz;
    this.pVizParams = {};
    this.entry = entry

    this.pVizParams.seq = {};
    this.entry.entities[0].chains.forEach(chain => {
      this.pVizParams.seq[chain.id] = entry.entities[0].sequence;
    });

    this.colorScheme = colorScheme;
    this.hosts = pVizHosts;

    this.selection = twinSelections;

    this.helixMapping();

    this.entry.entities[0].chains.forEach(async (chain) => {
      await this.render(chain.id);
      await this.resize(chain.id);
    });
  }

  public async render(chain: string) {

    let host = this.hosts[chain];

    //@ts-ignore
    if ($(host).width() !== 0) {
      let seq = this.pVizParams.seq[chain];

      let seqEntry = new this.pviz.SeqEntry({
        sequence: seq
      });
      new this.pviz.SeqEntryAnnotInteractiveView({
        model: seqEntry,
        collapsible: true,
        el: host
      }).render();

      let switchObj = this.selection;
      let pVizParams = this.pVizParams;
      let pv = this;

      //adding all features
      seqEntry.addFeatures(pVizParams.helixMap[chain].helixFeatureMap);

      //gradient coloring
      //if (this.paratopes.value === true) { this.applyGradient(this.pVizParams.parMap[chain].par_color_obj); }

      //await this.color(chain);

      //mouse over handlers

      //mouse click handlers
    }
  }

  helixMapping() {

    let helixMap: {[key: string]: {}} = {};

    this.entry.entities[0].chains.forEach(async (chain) => {
      let helixFeatureMap: any[] = [];
      Object.keys(chain.tracks).forEach(track => {
        if (track === "HELIX_P") {
          helixFeatureMap.push({
            category: 'Alpha helix',
            type: 'alpha',
            start: chain.tracks[track][0],
            end: chain.tracks[track].slice(-1),
            text: '',
            improbable: true
          })
        }
      });

      helixMap[chain.id] = { helixFeatureMap: helixFeatureMap };
    });

    this.pVizParams.helixMap = (helixMap);
  }

  // async color(chosenTracksChain) {

  //   let switchObj = this.selection;
  //   let pVizParams = this.pVizParams;

  //   Object.keys(switchObj).forEach((keyChain) => {

  //     let selectorStr = 'g.feature.data.D rect.feature';
  //     let denElements = document.querySelectorAll(selectorStr);
  //     let denNumbers = pVizParams.denMap[keyChain].den_el_obj;
  //     Object.keys(switchObj[keyChain]).forEach((keyPosition) => {

  //       let position = parseInt(keyPosition);

  //       if (keyChain === chosenTracksChain) {
  //         //densities
  //         if (denNumbers.indexOf(position) !== -1) {
  //           if (switchObj[keyChain][keyPosition]['state'] === false) {
  //             //@ts-ignore
  //             denElements[denNumbers.indexOf(position)].style.fill = pVizParams.denMap[keyChain].den_color_obj['D'][denNumbers.indexOf(position)];
  //           } else {
  //             //@ts-ignore
  //             denElements[denNumbers.indexOf(position)].style.fill = 'black';
  //           }
  //         }

  //         //predicted PTMS
  //         let listsPredictedPtms = [];

  //         this.ptmChoices.forEach(ptm => {
  //           let selectorStrPTM = 'g.feature.' + ptm.replace(" ", "_") + ' rect.feature';
  //           let elPTM = document.querySelectorAll(selectorStrPTM);
  //           let el_lstPTM = pVizParams.ptmMap[chosenTracksChain].ptm_el_obj[mutcodes[ptm.replace(" ", "_")]];
  //           listsPredictedPtms[ptm] = [elPTM, el_lstPTM];
  //         });

  //         Object.keys(listsPredictedPtms).forEach(ptm => {
  //           let elPTM = listsPredictedPtms[ptm][0];
  //           let el_lstPTM = listsPredictedPtms[ptm][1];
  //           if (typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1) {
  //             if (switchObj[keyChain][position]['state'] === false) {
  //               elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.ptmMap[keyChain].ptm_color_obj[mutcodes[ptm.replace(" ", "_")]][el_lstPTM.indexOf(position)];
  //             } else {
  //               elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
  //             }
  //           }
  //         });

  //         //motif PTMS
  //         let listsMotifsPtms = [];

  //         this.motChoices.forEach(ptm => {
  //           let selectorStrPTM = 'g.feature.' + mutcodes[ptm.replaceAll(" ", "_")] + "."
  //             + ptm.replaceAll(" ", "_").replaceAll(")", "\\)").replaceAll("(", "\\(") + ' rect.feature';
  //           let elPTM = document.querySelectorAll(selectorStrPTM);
  //           let el_lstPTM = pVizParams.motMap[chosenTracksChain].mot_el_obj[mutcodes[ptm.replace(" ", "_")]];
  //           listsMotifsPtms[ptm] = [elPTM, el_lstPTM];
  //         });

  //         Object.keys(listsMotifsPtms).forEach(ptm => {
  //           let elPTM = listsMotifsPtms[ptm][0];
  //           let el_lstPTM = listsMotifsPtms[ptm][1];
  //           if (typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1) {
  //             if (switchObj[keyChain][position]['state'] === false) {
  //               elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.motMap[keyChain].mot_color_obj[mutcodes[ptm.replace(" ", "_")]][el_lstPTM.indexOf(position)];
  //             } else {
  //               elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
  //             }
  //           }
  //         });

  //         //observed PTMS
  //         let listsObservedPtms = [];

  //         this.obsChoices.forEach(ptm => {
  //           let selectorStrPTM = 'g.feature.' + ptm.replaceAll(" ", "_") + "."
  //             + ptm.replaceAll(" ", "_") + ' rect.feature';
  //           let elPTM = document.querySelectorAll(selectorStrPTM);
  //           let el_lstPTM = pVizParams.obsMap[chosenTracksChain].obs_el_obj[ptm.replace(" ", "_")];
  //           listsObservedPtms[ptm] = [elPTM, el_lstPTM];
  //         });

  //         Object.keys(listsObservedPtms).forEach(ptm => {
  //           let elPTM = listsObservedPtms[ptm][0];
  //           let el_lstPTM = listsObservedPtms[ptm][1];
  //           if (typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1) {
  //             if (switchObj[keyChain][position]['state'] === false) {
  //               elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.obsMap[keyChain].obs_color_obj[ptm.replace(" ", "_")][el_lstPTM.indexOf(position)];
  //             } else {
  //               elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
  //             }
  //           }
  //         });
  //       }
  //     });
  //   });
  // }

  // applyGradient(gradient_obj) {
  //   Object.keys(gradient_obj).forEach((ptm_track) => {
  //     let selectorStr = 'g.feature.' + ptm_track + ' rect.feature';
  //     let el = document.querySelectorAll(selectorStr);
  //     for (let i = 0; i < el.length; i++) {
  //       //@ts-ignore
  //       el[i].style.fill = gradient_obj[ptm_track][i];
  //     }
  //   })
  // }

  // resize handle
  private async resize(chain: string) {
    let host = this.hosts[chain];

    ui.onSizeChanged(host).subscribe(async (_) => {
      await this.render(chain)
    });
  }
}

