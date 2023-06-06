import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//@ts-ignore
import {PdbEntry} from '../pdb-entry';

export class PvizAspect {
  ngl: any;
  pviz: any;
  pVizParams: any;

  entry: PdbEntry;
  jsonObs: any;
  colorScheme: any;
  hosts: { [key: string]: HTMLDivElement; };
  selection: { [chain: string]: {[key: string]: any} };

  async init(entry: PdbEntry, colorScheme: any, pVizHosts: { [key: string]: HTMLDivElement; }, twinSelections: any) {
    //@ts-ignore
    this.pviz = window.pviz;
    this.pVizParams = {};
    this.entry = entry;

    this.pVizParams.seq = {};
    this.entry.entities[0].chains.forEach((chain) => {
      this.pVizParams.seq[chain.id] = entry.entities[0].sequence;
    });

    this.colorScheme = colorScheme;
    this.hosts = pVizHosts;

    this.selection = twinSelections;

    //this.helixMapping();
    this.siteMapping();
    this.entry.entities[0].chains.forEach(async (chain) => {
      await this.render(chain.id);
      await this.resize(chain.id);
    });
  }

  public async render(chain: string) {
    const host = this.hosts[chain];

    //@ts-ignore
    if ($(host).width() !== 0) {
      const seq = this.pVizParams.seq[chain];

      const seqEntry = new this.pviz.SeqEntry({
        sequence: seq,
      });
      new this.pviz.SeqEntryAnnotInteractiveView({
        model: seqEntry,
        collapsible: true,
        el: host,
      }).render();

      const switchObj = this.selection;
      const pVizParams = this.pVizParams;
      const pv = this;

      //adding all features
      //seqEntry.addFeatures(pVizParams.helixMap[chain].helixFeatureMap);
      seqEntry.addFeatures(pVizParams.siteMap[chain].featureMap);


      //gradient coloring
      //if (this.paratopes.value === true) { this.applyGradient(this.pVizParams.parMap[chain].par_color_obj); }

      //await this.color(chain);

      //mouse over handlers

      //mouse click handlers
      this.pviz.FeatureDisplayer.addClickCallback(['active'], async function(ft: any) {
        if (switchObj[chain][ft.start] === undefined) {
          switchObj[chain][ft.start] = {};
          switchObj[chain][ft.start]['state'] = true;
        } else {
          switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state'];
        }
        grok.events.fireCustomEvent('selectionChanged', null);
        await pv.color(chain);
      });
    }
  }

  helixMapping() {
    const helixMap: { [key: string]: {} } = {};

    this.entry.entities[0].chains.forEach(async (chain) => {
      const helixFeatureMap: any[] = [];
      Object.keys(chain.tracks).forEach((track) => {
        if (track === 'HELIX_P') {
          helixFeatureMap.push({
            category: 'Alpha helix',
            type: 'alpha',
            start: chain.tracks[track][0],
            end: chain.tracks[track].slice(-1),
            text: '',
            improbable: true,
          });
        }
      });

      helixMap[chain.id] = {helixFeatureMap: helixFeatureMap};
    });

    this.pVizParams.helixMap = (helixMap);
  }

  siteMapping() {
    const siteMap: { [key: string]: any } = {};
    const chains = ['A', 'B'];
    chains.forEach((chain) => {
      const featureMap: any[] = [];
      const siteArray: number[] = [24, 25, 26, 27, 28, 29, 48, 49, 50, 51, 52];


      siteArray.forEach((point) => {
        featureMap.push({
          groupSet: 'Active site',
          category: 'Active',
          type: 'active',
          start: point,
          end: point,
          text: '',
          improbable: true,
        });
      });
      siteMap[chain] = {featureMap: featureMap, elObj: siteArray};
    });

    this.pVizParams.siteMap = (siteMap);
  }

  async color(chosenTracksChain: string) {
    const switchObj = this.selection;
    const pVizParams = this.pVizParams;

    Object.keys(switchObj).forEach((keyChain) => {
      const selectorStr = 'g.feature.data.active rect.feature';
      const siteElements = document.querySelectorAll(selectorStr);
      const siteNumbers = pVizParams.siteMap[keyChain].elObj;
      Object.keys(switchObj[keyChain]).forEach((keyPosition) => {
        const position = parseInt(keyPosition);

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

  applyGradient(gradientObj: any) {
    Object.keys(gradientObj).forEach((ptmTrack) => {
      const selectorStr = 'g.feature.' + ptmTrack + ' rect.feature';
      const el = document.querySelectorAll(selectorStr);
      for (let i = 0; i < el.length; i++) {
        //@ts-ignore
        el[i].style.fill = gradientObj[ptmTrack][i];
      }
    });
  }

  // resize handle
  private async resize(chain: string) {
    const host = this.hosts[chain];

    ui.onSizeChanged(host).subscribe(async (_) => {
      await this.render(chain);
    });
  }
}

