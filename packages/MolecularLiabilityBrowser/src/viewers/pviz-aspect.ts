import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {MiscMethods} from './misc';
import {
  ColorSchemeType,
  DataLoader,
  ObsFeatureType,
  IFeature,
  JsonType,
  MutcodesDataType,
  PtmFeatureType,
  PvizMotType,
  PvizType,
  SeqEntry,
  PvizParamType,
  CdrFeatureType,
  PvizCdrType,
  PvizParType,
  ParFeatureType,
  DenFeatureType,
  PointType,
  PvizPtmType, PvizDenType
} from '../utils/data-loader';
import {InputBase} from 'datagrok-api/dg';

export class PvizAspect {
  private dataLoader: DataLoader;

  get mutcodes(): MutcodesDataType { return this.dataLoader.mutcodes; };

  ngl: any;
  pviz: PvizType;
  pVizParams: PvizParamType;

  json: JsonType;
  jsonObs: any;
  colorScheme: ColorSchemeType;
  hostH: any;
  hostL: any;
  ptmChoices: string[];
  motChoices: string[];
  obsChoices: string[];
  cdrScheme: InputBase;
  paratopes: InputBase;
  ptmProb: number;
  selection: any;

  constructor(dataLoader: DataLoader) {
    this.dataLoader = dataLoader;
  }

  async init(
    json: JsonType, jsonObs: Object, colorScheme: ColorSchemeType,
    pVizHostH: HTMLDivElement, pVizHostL: HTMLDivElement,
    ptmChoices: string[], ptmMotifChoices: string[], ptmObsChoices: string[],
    cdrScheme: InputBase, paratopes: InputBase, ptmProb: number, twinSelections: { [chain: string]: Object }
  ) {
    //@ts-ignore
    this.pviz = window.pviz;
    this.pVizParams = {
      seq: {'H': json.heavy_seq, 'L': json.light_seq},
      cdrMap: {}, parMap: {}, denMap: {}, ptmMap: {}, motMap: {}, obsMap: {}
    };

    this.json = json;
    this.jsonObs = jsonObs;
    this.colorScheme = colorScheme;
    this.hostH = pVizHostH;
    this.hostL = pVizHostL;

    this.selection = twinSelections;

    this.cdrMapping(cdrScheme);
    this.parMapping(paratopes);
    this.denMapping();
    this.ptmMapping(ptmChoices, ptmProb);
    this.motMapping(ptmMotifChoices, ptmProb);
    if (this.jsonObs !== null) this.obsMapping(ptmObsChoices);

    await this.render('H');

    await this.resize('H');
    await this.resize('L');
  }

  public async render(chain: string) {
    const host = chain == 'H' ? this.hostH : this.hostL;

    //@ts-ignore
    if ($(host).width() !== 0) {
      const seq = this.pVizParams.seq[chain];

      const seqEntry: SeqEntry = new this.pviz.SeqEntry({sequence: seq});
      new this.pviz.SeqEntryAnnotInteractiveView({
        model: seqEntry,
        collapsible: true,
        el: host,
      }).render();

      const ptmCodes = Object.keys(this.pVizParams.ptmMap[chain].ptmColorObj);
      const motCodes = Object.keys(this.pVizParams.motMap[chain].motColorObj);

      let obsCodes = null;
      if (this.jsonObs !== null)
        obsCodes = Object.keys(this.pVizParams.obsMap[chain].obsColorObj);

      const switchObj = this.selection;
      const pVizParams = this.pVizParams;
      const pv = this;

      //adding all features
      seqEntry.addFeatures(pVizParams.cdrMap[chain].cdrFeatureMap);
      if (this.paratopes.value === true)
        seqEntry.addFeatures(this.pVizParams.parMap[chain].parFeatureMap);

      seqEntry.addFeatures(this.pVizParams.denMap[chain].denFeatureMap);
      seqEntry.addFeatures(this.pVizParams.ptmMap[chain].ptmFeatureMap);
      seqEntry.addFeatures(this.pVizParams.motMap[chain].motFeatureMap);

      //gradient coloring
      if (this.paratopes.value === true) this.applyGradient(this.pVizParams.parMap[chain].parColorObj);
      if (this.jsonObs !== null) seqEntry.addFeatures(this.pVizParams.obsMap[chain].obsFeatureMap);

      this.applyGradient(this.pVizParams.denMap[chain].denColorObj);
      this.applyGradient(this.pVizParams.ptmMap[chain].ptmColorObj);
      this.applyGradient(this.pVizParams.motMap[chain].motColorObj);
      if (this.jsonObs !== null) this.applyGradient(this.pVizParams.obsMap[chain].obsColorObj);

      await this.color(chain);

      //mouse over handlers
      this.pviz.FeatureDisplayer.addMouseoverCallback(['P'], async function(ft: IFeature) {
        const selectorStr = 'g.feature.data.P.Paratope_predictions rect.feature';
        const elList = document.querySelectorAll<SVGRectElement>(selectorStr);
        const elPosList = pVizParams.parMap[chain].parElObj;
        const probLst = pVizParams.parMap[chain].parProbObj;
        const el: SVGRectElement = elList[elPosList.indexOf((ft.start).toString())];
        const prob = probLst[elPosList.indexOf((ft.start).toString())];
        const elClientRect = el.getBoundingClientRect();
        ui.tooltip.show(
          ui.span([`Probability: ${prob.toFixed(2)}`]),
          elClientRect.left + 10,
          elClientRect.top + 10);
      }).addMouseoutCallback(['P'], function(ft: IFeature) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(['D'], async function(ft: IFeature) {
        const selectorStr = 'g.feature.data.D rect.feature';
        const elList = document.querySelectorAll<SVGRectElement>(selectorStr);
        const elPosList: number[] = pVizParams.denMap[chain].denElObj;
        const probList: number[] = pVizParams.denMap[chain].denProbObj;
        const ptmList = pVizParams.denMap[chain].denPtmArr;
        const el: SVGRectElement = elList[elPosList.indexOf(ft.start)];
        const prob: number = probList[elPosList.indexOf(ft.start)];
        const ptmArPointList = ptmList[elPosList.indexOf(ft.start)];
        const ptmsStrList: string[] = Array<string>(2);

        for (let i = 0; i < ptmArPointList.length; i++) {
          const ptmArPoint = ptmArPointList[i];
          const ptmArProb: string = (ptmArPoint[1] > 1 ? ptmArPoint[1] / 100 : ptmArPoint[1]).toFixed(2);
          ptmsStrList[i] = ptmArPoint[0].replace('_', ' ') + ' probability  ~' + ptmArProb;
        }
        const ptmsStr = ptmsStrList.join('\n');
        const elClientRect = el.getBoundingClientRect();

        ui.tooltip.show(
          ui.divText(`Probability: ${prob.toFixed(2)}\n${ptmsStr}`),
          elClientRect.left + 10,
          elClientRect.top + 10);
      }).addMouseoutCallback(['D'], function(ft: IFeature) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(ptmCodes, async (ft: IFeature) => {
        const selectorStr = 'g.feature.' + this.mutcodes[ft.category.replaceAll(' ', '_')] + '.' +
          ft.category.replaceAll(' ', '_').replaceAll(')', '\\)').replaceAll('(', '\\(') + ' rect.feature';
        const elList = document.querySelectorAll<SVGRectElement>(selectorStr);
        const elPosList: number[] = pVizParams.ptmMap[chain].ptmElObj[this.mutcodes[ft.category.replaceAll(' ', '_')]];
        const el: SVGRectElement = elList[elPosList.indexOf(ft.start)];
        const elClientRect = el.getBoundingClientRect();

        ui.tooltip.show(
          ui.span([`${ft.category}`]),
          elClientRect.left + 10,
          elClientRect.top + 10);
      }).addMouseoutCallback(ptmCodes, function(ft: IFeature) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(motCodes, async (ft: IFeature) => {
        const selectorStr = 'g.feature.' + this.mutcodes[ft.category.replaceAll(' ', '_')] + '.' +
          ft.category.replaceAll(' ', '_').replaceAll(')', '\\)').replaceAll('(', '\\(') + ' rect.feature';
        const elList = document.querySelectorAll<SVGRectElement>(selectorStr);
        const elPosList = pVizParams.motMap[chain].motElObj[this.mutcodes[ft.category.replaceAll(' ', '_')]];
        const el: SVGRectElement = elList[elPosList.indexOf(ft.start)];
        const elClientRect = el.getBoundingClientRect();

        ui.tooltip.show(
          ui.span([`${ft.category}`]),
          elClientRect.left + 10,
          elClientRect.top + 10);
      }).addMouseoutCallback(motCodes, function(ft: IFeature) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(obsCodes, async function(ft: IFeature) {
        const selectorStr = 'g.feature.' + ft.category.replaceAll(' ', '_') + '.' +
          ft.category.replaceAll(' ', '_') + ' rect.feature';

        const elList = document.querySelectorAll<SVGRectElement>(selectorStr);
        const elPosList: number[] = pVizParams.obsMap[chain].obsElObj[ft.category];
        const types: string[] = pVizParams.obsMap[chain].obsProbObj[ft.start][0];
        const probabilities: number[] = pVizParams.obsMap[chain].obsProbObj[ft.start][1];
        const el: SVGRectElement = elList[elPosList.indexOf(ft.start)];

        let ptmsStr = '';

        for (let i = 0; i < types.length; i++) {
          const prob = probabilities[i] == 0 ? 'NA' : (probabilities[i] == 0.01 ? 0 : probabilities[i].toFixed(2));
          ptmsStr += '\n' + types[i] + ' ~' + prob;
        }
        const elClientRect = el.getBoundingClientRect();

        ui.tooltip.show(
          ui.divText(`Probability: ${ptmsStr}`),
          elClientRect.left + 10,
          elClientRect.top + 10);
      }).addMouseoutCallback(obsCodes, function(ft: IFeature) {
        ui.tooltip.hide();
      });

      const featureClickCallback = async (ft: IFeature) => {
        if (switchObj[chain][ft.start] === undefined) {
          switchObj[chain][ft.start] = {};
          switchObj[chain][ft.start]['state'] = true;
        } else {
          switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state'];
        }

        grok.events.fireCustomEvent('selectionChanged', null);
        await pv.color(chain);
      };
      //mouse click handlers
      this.pviz.FeatureDisplayer.addClickCallback(['D'], featureClickCallback);
      this.pviz.FeatureDisplayer.addClickCallback(ptmCodes, featureClickCallback);
      this.pviz.FeatureDisplayer.addClickCallback(motCodes, featureClickCallback);
      this.pviz.FeatureDisplayer.addClickCallback(obsCodes, featureClickCallback);
    }
  }

  async color(chosenTracksChain: string) {
    const switchObj = this.selection;
    const pVizParams = this.pVizParams;

    Object.keys(switchObj).forEach((keyChain: string) => {
      const selectorStr = 'g.feature.data.D rect.feature';
      const denElements = document.querySelectorAll(selectorStr) as NodeListOf<HTMLElement>;
      const denNumbers: number[] = pVizParams.denMap[keyChain].denElObj;
      Object.keys(switchObj[keyChain]).forEach((keyPosition) => {
        const position = parseInt(keyPosition);

        if (keyChain === chosenTracksChain) {
          //densities
          if (denNumbers.indexOf(position) !== -1) {
            if (switchObj[keyChain][keyPosition]['state'] === false) {
              denElements[denNumbers.indexOf(position)].style.fill =
                pVizParams.denMap[keyChain].denColorObj['D'][denNumbers.indexOf(position)];
            } else {
              denElements[denNumbers.indexOf(position)].style.fill = 'black';
            }
          }

          //predicted PTMS
          const listsPredictedPtms: { [ptm: string]: any } = {};

          this.ptmChoices.forEach((ptm: string) => {
            const selectorStrPTM = 'g.feature.' + ptm.replace(' ', '_') + ' rect.feature';
            const elPTM = document.querySelectorAll(selectorStrPTM);
            const elLstPTM = pVizParams.ptmMap[chosenTracksChain].ptmElObj[this.mutcodes[ptm.replace(' ', '_')]];
            listsPredictedPtms[ptm] = [elPTM, elLstPTM];
          });

          Object.keys(listsPredictedPtms).forEach((ptm) => {
            const elPTM = listsPredictedPtms[ptm][0];
            const elLstPTM = listsPredictedPtms[ptm][1];
            if (typeof elLstPTM !== 'undefined' && elLstPTM.indexOf(position) !== -1) {
              if (switchObj[keyChain][position]['state'] === false) {
                elPTM[elLstPTM.indexOf(position)].style.fill =
                  pVizParams.ptmMap[keyChain]
                    .ptmColorObj[this.mutcodes[ptm.replace(' ', '_')]][elLstPTM.indexOf(position)];
              } else {
                elPTM[elLstPTM.indexOf(position)].style.fill = 'black';
              }
            }
          });

          //motif PTMS
          const listsMotifsPtms: { [ptm: string]: any } = {};

          this.motChoices.forEach((ptm: string) => {
            const selectorStrPTM = 'g.feature.' + this.mutcodes[ptm.replaceAll(' ', '_')] + '.' +
              ptm.replaceAll(' ', '_').replaceAll(')', '\\)').replaceAll('(', '\\(') + ' rect.feature';
            const elPTM = document.querySelectorAll(selectorStrPTM);
            const elLstPTM = pVizParams.motMap[chosenTracksChain].motElObj[this.mutcodes[ptm.replace(' ', '_')]];
            listsMotifsPtms[ptm] = [elPTM, elLstPTM];
          });

          Object.keys(listsMotifsPtms).forEach((ptm) => {
            const elPTM = listsMotifsPtms[ptm][0];
            const elLstPTM = listsMotifsPtms[ptm][1];
            if (typeof elLstPTM !== 'undefined' && elLstPTM.indexOf(position) !== -1) {
              if (switchObj[keyChain][position]['state'] === false) {
                elPTM[elLstPTM.indexOf(position)].style.fill =
                  pVizParams.motMap[keyChain]
                    .motColorObj[this.mutcodes[ptm.replace(' ', '_')]][elLstPTM.indexOf(position)];
              } else {
                elPTM[elLstPTM.indexOf(position)].style.fill = 'black';
              }
            }
          });

          //observed PTMS
          const listsObservedPtms: { [ptm: string]: any } = [];

          if (this.jsonObs !== null) {
            this.obsChoices.forEach((ptm) => {
              const selectorStrPTM = 'g.feature.' + ptm.replaceAll(' ', '_') + '.' +
                ptm.replaceAll(' ', '_') + ' rect.feature';
              const elPTM = document.querySelectorAll(selectorStrPTM);
              const elLstPTM = pVizParams.obsMap[chosenTracksChain].obsElObj[ptm.replace(' ', '_')];
              listsObservedPtms[ptm] = [elPTM, elLstPTM];
            });

            Object.keys(listsObservedPtms).forEach((ptm) => {
              const elPTM = listsObservedPtms[ptm][0];
              const elLstPTM = listsObservedPtms[ptm][1];
              if (typeof elLstPTM !== 'undefined' && elLstPTM.indexOf(position) !== -1) {
                if (switchObj[keyChain][position]['state'] === false) {
                  elPTM[elLstPTM.indexOf(position)].style.fill =
                    pVizParams.obsMap[keyChain].obsColorObj[ptm.replace(' ', '_')][elLstPTM.indexOf(position)];
                } else {
                  elPTM[elLstPTM.indexOf(position)].style.fill = 'black';
                }
              }
            });
          }
        }
      });
    });
  }

  // mapping objects for sequence rendering
  cdrMapping(cdrScheme: InputBase) {
    this.cdrScheme = cdrScheme;
    const cdrMap: { [chain: string]: PvizCdrType } = {};
    const chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain) => {
      const cdrFeatureMap: CdrFeatureType[] = [];
      const cdrRanges: PointType[] = this.json.cdr_ranges[this.cdrScheme.value + '_CDR' + chain + '_ranges'];
      if (this.cdrScheme.value !== MiscMethods.NoSchemeItem) {
        Object.values(cdrRanges).forEach((range) => {
          cdrFeatureMap.push({
            category: 'CDR region',
            type: 'CDR',
            start: range[0],
            end: range[1],
            text: '',
            improbable: true,
          });
        });
      }
      cdrMap[chain] = {cdrFeatureMap: cdrFeatureMap};
    });

    this.pVizParams.cdrMap = (cdrMap);
  }

  parMapping(paratopes: InputBase) {
    this.paratopes = paratopes;
    const parMap: { [chain: string]: PvizParType } = {};
    const chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain: string) => {
      const parFeatureMap: ParFeatureType[] = [];
      const parColorArr: string[] = [];
      const parElObj: string[] = [];
      const parProbObj: number[] = [];
      const palette: string[] = MiscMethods.interpolateColors('(255, 255, 255)', '(255, 0, 255)', 100);

      Object.keys(this.json.parapred_predictions[chain]).forEach((index) => {
        parFeatureMap.push({
          category: 'Paratope predictions',
          type: 'P',
          start: parseInt(index),
          end: parseInt(index),
          text: '',
          improbable: true,
        });
        parColorArr.push(palette[Math.round(this.json.parapred_predictions[chain][index] * 100)]);
        parElObj.push(index);
        parProbObj.push(this.json.parapred_predictions[chain][index]);
      });

      parMap[chain] = {
        parFeatureMap: parFeatureMap,
        parColorObj: {'P': parColorArr},
        parElObj: parElObj,
        parProbObj: parProbObj
      };
    });

    this.pVizParams.parMap = (parMap);
  };

  denMapping() {
    const denMap: { [chain: string]: PvizDenType } = {};
    const chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain: string) => {
      const denFeatureList: number[] = [];
      let denColorArr = new Array(this.pVizParams.seq[chain].length).fill(-1);
      let denPtmArr = new Array(this.pVizParams.seq[chain].length).fill([]);
      const palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

      Object.keys(this.json.ptm_predictions[chain]).forEach((ptm: string) => {
        if (this.json.ptm_predictions[chain][ptm][0][1] <= 1) {
          this.json.ptm_predictions[chain][ptm].forEach((point: PointType) => {
            if (!(denFeatureList.includes(point[0])))
              denFeatureList.push(point[0]);

            denColorArr[point[0]] = denColorArr[point[0]] == -1 ?
              (point[1] > 1 ? point[1] / 100 : point[1]) :
              1 - (1 - denColorArr[point[0]]) * (1 - (point[1] > 1 ? point[1] / 100 : point[1]));
            denPtmArr[point[0]] = denPtmArr[point[0]].concat([[ptm, point[1]]]);
          });
        }
      });

      denFeatureList.sort((a, b) => a - b);
      const denElObj = denFeatureList.slice();
      const denFeatureMap: DenFeatureType[] = denFeatureList.map(function(ft) {
        return {
          groupSet: 'Predicted PTM density',
          category: '',
          type: 'D',
          start: ft,
          end: ft,
          text: '',
          improbable: true,
        };
      });

      denColorArr = denColorArr.filter((x) => {
        return x !== -1;
      });

      denPtmArr = denPtmArr.filter((x) => {
        return x.length > 0;
      });

      const denProbObj = denColorArr.slice();

      for (let i = 0; i < denColorArr.length; i++)
        denColorArr[i] = palette[Math.round(denColorArr[i] * 4)];

      denMap[chain] = {
        denFeatureMap: denFeatureMap,
        denColorObj: {'D': denColorArr},
        denElObj: denElObj,
        denProbObj: denProbObj,
        denPtmArr: denPtmArr
      };
    });

    this.pVizParams.denMap = (denMap);
  }

  ptmMapping(ptmChoices: string[], prob: number) {
    this.ptmChoices = ptmChoices;
    this.ptmProb = prob;

    const ptmMap: { [chain: string]: PvizPtmType } = {};
    const chains: string[] = Object.keys(this.pVizParams.seq);
    chains.forEach((chain: string) => {
      const ptmFeatureMap: PtmFeatureType[] = [];
      const ptmColorObj: { [ptm: string]: string[] } = {};
      const ptmElObj: { [ptm: string]: number[] } = {};
      const ptmProbObj: { [ptm: string]: number[] } = {};
      const palette: string[] = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

      this.ptmChoices.forEach((ptm: string) => {
        const ptmArray: PointType[] = this.json.ptm_predictions[chain][ptm.replace(' ', '_')];
        if (ptmArray !== undefined) {
          const ptmColorArr: string[] = [];
          const ptmElArr: number[] = [];
          const ptmProbArr: number[] = [];
          ptmArray.forEach((point: PointType) => {
            if (point[1] > prob) {
              ptmFeatureMap.push({
                groupSet: 'Predicted PTMs',
                category: ptm,
                type: this.mutcodes[ptm.replace(' ', '_')],
                start: point[0],
                end: point[0],
                text: '',
                improbable: true,
              });
              ptmColorArr.push(palette[point[1] > 1 ? Math.round(point[1] * 4) / 100 : Math.round(point[1] * 4)]);
              ptmElArr.push(point[0]);
              ptmProbArr.push(point[1]);
            }
          });
          if (ptmColorArr.length > 0) {
            ptmColorObj[this.mutcodes[ptm.replace(' ', '_')]] = ptmColorArr;
            ptmElObj[this.mutcodes[ptm.replace(' ', '_')]] = ptmElArr;
            ptmProbObj[this.mutcodes[ptm.replace(' ', '_')]] = ptmProbArr;
          }
        }
      });
      ptmMap[chain] = {
        ptmFeatureMap: ptmFeatureMap, ptmColorObj: ptmColorObj,
        ptmElObj: ptmElObj, ptmProbObj: ptmProbObj,
      };
    });

    this.pVizParams.ptmMap = (ptmMap);
  }

  motMapping(ptmChoices: string[], prob: number) {
    this.motChoices = ptmChoices;
    const motMap: { [chain: string]: PvizMotType } = {};
    const chains: string[] = Object.keys(this.pVizParams.seq);
    chains.forEach((chain: string) => {
      const ptmFeatureMap: PtmFeatureType[] = [];
      const ptmColorObj: { [ptm: string]: string[] } = {};
      const ptmElObj: { [ptm: string]: number[] } = {};
      const ptmProbObj: { [ptm: string]: number[] } = {};
      const palette: string[] = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

      ptmChoices.forEach((ptm: string) => {
        const ptmArray = this.json.ptm_predictions[chain][ptm.replaceAll(' ', '_')];
        if (ptmArray !== undefined) {
          const ptmColorArr: string[] = [];
          const ptmElArr: number[] = [];
          const ptmProbArr: number[] = [];
          ptmArray.forEach((point: number[]) => {
            if (point[1] > prob) {
              ptmFeatureMap.push({
                groupSet: 'Predicted PTMs',
                category: ptm,
                type: this.mutcodes[ptm.replaceAll(' ', '_')],
                start: point[0],
                end: point[0],
                text: '',
                improbable: true,
              });
              ptmColorArr.push(palette[point[1] > 1 ? Math.round(point[1] * 4) / 100 : Math.round(point[1] * 4)]);
              ptmElArr.push(point[0]);
              ptmProbArr.push(0.99);
            }
          });
          if (ptmColorArr.length > 0) {
            ptmColorObj[this.mutcodes[ptm.replaceAll(' ', '_')]] = ptmColorArr;
            ptmElObj[this.mutcodes[ptm.replaceAll(' ', '_')]] = ptmElArr;
            ptmProbObj[this.mutcodes[ptm.replaceAll(' ', '_')]] = ptmProbArr;
          }
        }
      });
      motMap[chain] = {
        motFeatureMap: ptmFeatureMap, motColorObj: ptmColorObj,
        motElObj: ptmElObj, motProbObj: ptmProbObj,
      };
    });

    this.pVizParams.motMap = (motMap);
  }

  obsMapping(ptmChoices: string[]) {
    this.obsChoices = ptmChoices;
    const obsMap: { [chain: string]: any } = {};
    const chains = Object.keys(this.pVizParams.seq);
    const palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

    chains.forEach((chain) => {
      const obsFeatureMap: ObsFeatureType[] = [];
      const obsColorObj: { [ptm: string]: any[] } = {};
      const obsElObj: { [ptm: string]: any[] } = {};
      const obsProbObj: { [pos: string]: [string[], number[]] } = {};

      ptmChoices.forEach((ptm: string) => {
        const ptmTree = this.jsonObs[chain][ptm.replace(' ', '_')];
        if (ptmTree !== undefined) {
          const obsColorArr: any[] = [];
          const obsElArr: any[] = [];
          const obsTypesArr = [];
          const obsTypesProbsArr = [];

          Object.keys(ptmTree).forEach((type) => {
            const point = this.jsonObs[chain][ptm.replace(' ', '_')][type];
            if (!obsElArr.includes(point[0]))
              obsElArr.push(point[0]);
          });

          obsElArr.forEach((position) => {
            const types: string[] = [];
            const typesProbs: number[] = [];
            let prob = -1;
            Object.keys(ptmTree).forEach((type) => {
              const point: number[] = this.jsonObs[chain][ptm][type];
              if (point[0] === position) {
                types.push(type);
                typesProbs.push(point[1]);
                const addProb = point[1] == 0.01 ? 0 : point[1];
                prob = prob == -1 ? addProb / 100 : 1 - (1 - prob) * (1 - addProb / 100);
              }
            });

            obsColorArr.push(palette[Math.round(prob * 4)]);
            obsTypesArr.push(types);
            obsTypesProbsArr.push(typesProbs);

            obsFeatureMap.push({
              groupSet: 'Observed PTMs',
              category: ptm,
              type: ptm,
              start: position,
              end: position,
              text: ptm,
              improbable: true,
            });

            obsColorObj[ptm.replace(' ', '_')] = obsColorArr;
            obsElObj[ptm.replace(' ', '_')] = obsElArr;
            obsProbObj[position] = [types, typesProbs];
          });
        }
      });

      obsMap[chain] = {
        obsFeatureMap: obsFeatureMap, obsColorObj: obsColorObj,
        obsElObj: obsElObj, obsProbObj: obsProbObj,
      };
    });

    this.pVizParams.obsMap = (obsMap);
  }

  applyGradient(gradientObj: { [ptmTrack: string]: string[] }) {
    Object.keys(gradientObj).forEach((ptmTrack) => {
      const selectorStr = 'g.feature.' + ptmTrack + ' rect.feature';
      const elList = document.querySelectorAll<SVGRectElement>(selectorStr);
      for (let i = 0; i < elList.length; i++)
        elList[i].style.fill = gradientObj[ptmTrack][i];
    });
  }

  // resize handle
  private async resize(chain: string) {
    const host = chain == 'H' ? this.hostH : this.hostL;

    ui.onSizeChanged(host).subscribe(async (_) => {
      await this.render(chain);
    });
  }
}
