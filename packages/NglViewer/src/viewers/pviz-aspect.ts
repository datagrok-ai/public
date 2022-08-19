import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
//@ts-ignore
import mutcodes from "../externalData/mutcodes.json";
import { MiscMethods } from "./misc.js"

export class PvizAspect {
  ngl: any;
  pviz: any;
  pVizParams: any;

  entry: string;
  jsonObs: any;
  colorScheme: any;
  hosts: { [key: string]: HTMLDivElement; };
  selection: { [chain: string]: {[key: string]: any} };
  chains: any;

  async init(entry: string, colorScheme: any, pVizHosts: { [key: string]: HTMLDivElement; }, twinSelections: any, chains: any) {

    //@ts-ignore
    this.pviz = window.pviz;
    this.pVizParams = {};
    this.entry = entry

    this.pVizParams.seq = {};
    this.chains = chains;
    this.chains.forEach(chain => {
      //this.pVizParams.seq[chain] = entry.entities[0].sequence;
    });

    this.colorScheme = colorScheme;
    this.hosts = pVizHosts;

    this.selection = twinSelections;

    //this.helixMapping();
    this.siteMapping();
    this.chains.forEach(async (chain: string) => {
      await this.render(chain);
      await this.resize(chain);
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
      //seqEntry.addFeatures(pVizParams.helixMap[chain].helixFeatureMap);
      seqEntry.addFeatures(pVizParams.siteMap[chain].featureMap);


      //gradient coloring
      //if (this.paratopes.value === true) { this.applyGradient(this.pVizParams.parMap[chain].par_color_obj); }

      //await this.color(chain);

      //mouse over handlers

      //mouse click handlers
      this.pviz.FeatureDisplayer.addClickCallback(['active'], async function (ft: any) {
        if (switchObj[chain][ft.start] === undefined) {
          switchObj[chain][ft.start] = {};
          switchObj[chain][ft.start]['state'] = true;
        } else {
          switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
        }
        grok.events.fireCustomEvent("selectionChanged", null);
        await pv.color(chain);
      });
    }
  }

  helixMapping() {

    let helixMap: { [key: string]: {} } = {};

    // this.chains.forEach(async (chain: string) => {
    //   let helixFeatureMap: any[] = [];
    //   Object.keys(chain.tracks).forEach(track => {
    //     if (track === "HELIX_P") {
    //       helixFeatureMap.push({
    //         category: 'Alpha helix',
    //         type: 'alpha',
    //         start: chain.tracks[track][0],
    //         end: chain.tracks[track].slice(-1),
    //         text: '',
    //         improbable: true
    //       })
    //     }
    //   });

    //   helixMap[chain.id] = { helixFeatureMap: helixFeatureMap };
    // });

    this.pVizParams.helixMap = (helixMap);
  }

  siteMapping() {
    let siteMap: { [key: string]: any } = {};
    let chains = ["A", "B"];
    chains.forEach((chain) => {
      let featureMap: any[] = [];
      let siteArray: number[] = [24, 25, 26, 27, 28, 29, 48, 49, 50, 51, 52];


      siteArray.forEach(point => {
        featureMap.push({
          groupSet: 'Active site',
          category: 'Active',
          type: 'active',
          start: point,
          end: point,
          text: '',
          improbable: true
        })
      })
      siteMap[chain] = { featureMap: featureMap, elObj: siteArray };
    })

    this.pVizParams.siteMap = (siteMap);
  }

  async color(chosenTracksChain: string) {

    let switchObj = this.selection;
    let pVizParams = this.pVizParams;

    Object.keys(switchObj).forEach((keyChain) => {

      let selectorStr = 'g.feature.data.active rect.feature';
      let siteElements = document.querySelectorAll(selectorStr);
      let siteNumbers = pVizParams.siteMap[keyChain].elObj;
      Object.keys(switchObj[keyChain]).forEach((keyPosition) => {

        let position = parseInt(keyPosition);

        if (keyChain === chosenTracksChain) {
          //densities
          if (siteNumbers.indexOf(position) !== -1) {
            if (switchObj[keyChain][keyPosition]['state'] === false) {
              //@ts-ignore
              siteElements[siteNumbers.indexOf(position)].style.fill = 'red';
            } else {
              //@ts-ignore
              siteElements[siteNumbers.indexOf(position)].style.fill = 'black';
            }
          }
        }
      });
    });
  }

  applyGradient(gradient_obj: any) {
    Object.keys(gradient_obj).forEach((ptm_track) => {
      let selectorStr = 'g.feature.' + ptm_track + ' rect.feature';
      let el = document.querySelectorAll(selectorStr);
      for (let i = 0; i < el.length; i++) {
        //@ts-ignore
        el[i].style.fill = gradient_obj[ptm_track][i];
      }
    })
  }

  // resize handle
  private async resize(chain: string) {
    let host = this.hosts[chain];

    ui.onSizeChanged(host).subscribe(async (_) => {
      await this.render(chain)
    });
  }
}

