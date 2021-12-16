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

  json: any;
  jsonObs: any;
  colorScheme: any;
  hostH: any;
  hostL: any;
  ptmChoices: any;
  motChoices: any;
  obsChoices: any;
  cdrScheme: any;
  paratopes: any;
  ptmProb: any;
  selection: any;

  async init(json, jsonObs, colorScheme, pVizHostH, pVizHostL,
    ptmChoices, ptmMotifChoices, ptmObsChoices, cdrScheme, paratopes, ptmProb, twinSelections) {

    //@ts-ignore
    this.pviz = window.pviz;
    this.pVizParams = {};
    this.pVizParams.seq = { 'H': json.heavy_seq, 'L': json.light_seq }

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
    if (this.jsonObs !== null) { this.obsMapping(ptmObsChoices); };

    await this.render('H');

    await this.resize('H');
    await this.resize('L');
    //@ts-ignore
    //return MiscMethods.setDockSize(view, inputs.nglNode, inputs.sequenceTabs, paratopes);
  }

  public async render(chain: string) {

    let host = chain == "H" ? this.hostH : this.hostL;

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

      let ptmCodes = Object.keys(this.pVizParams.ptmMap[chain].ptm_color_obj);
      let motCodes = Object.keys(this.pVizParams.motMap[chain].mot_color_obj);

      let obsCodes = null;
      if (this.jsonObs !== null) {
        obsCodes = Object.keys(this.pVizParams.obsMap[chain].obs_color_obj);
      }
      let switchObj = this.selection;
      let pVizParams = this.pVizParams;
      let pv = this;

      //adding all features
      seqEntry.addFeatures(pVizParams.cdrMap[chain].cdr_feature_map);
      if (this.paratopes.value === true) {
        seqEntry.addFeatures(this.pVizParams.parMap[chain].par_feature_map);
      }
      seqEntry.addFeatures(this.pVizParams.denMap[chain].den_feature_map);
      seqEntry.addFeatures(this.pVizParams.ptmMap[chain].ptm_feature_map);
      seqEntry.addFeatures(this.pVizParams.motMap[chain].mot_feature_map);

      //gradient coloring
      if (this.paratopes.value === true) { this.applyGradient(this.pVizParams.parMap[chain].par_color_obj); }
      if (this.jsonObs !== null) { seqEntry.addFeatures(this.pVizParams.obsMap[chain].obs_feature_map); }

      this.applyGradient(this.pVizParams.denMap[chain].den_color_obj);
      this.applyGradient(this.pVizParams.ptmMap[chain].ptm_color_obj);
      this.applyGradient(this.pVizParams.motMap[chain].mot_color_obj);
      if (this.jsonObs !== null) { this.applyGradient(this.pVizParams.obsMap[chain].obs_color_obj) };

      await this.color(chain);

      //mouse over handlers
      this.pviz.FeatureDisplayer.addMouseoverCallback(['P'], async function (ft) {

        let selectorStr = 'g.feature.data.P.Paratope_predictions rect.feature';
        let el = document.querySelectorAll(selectorStr);
        let el_lst = pVizParams.parMap[chain].par_el_obj;
        let prob_lst = pVizParams.parMap[chain].par_prob_obj;
        //@ts-ignore
        el = el[el_lst.indexOf((ft.start).toString())];
        let prob = prob_lst[el_lst.indexOf((ft.start).toString())];

        ui.tooltip.show(
          ui.span([`Probability: ${prob.toFixed(2)}`]),
          //@ts-ignore
          el.getBoundingClientRect().left + 10,
          //@ts-ignore
          el.getBoundingClientRect().top + 10
        );

      }).addMouseoutCallback(['P'], function (ft) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(['D'], async function (ft) {

        let selectorStr = 'g.feature.data.D rect.feature';
        let el = document.querySelectorAll(selectorStr);
        let el_lst = pVizParams.denMap[chain].den_el_obj;
        let prob_lst = pVizParams.denMap[chain].den_prob_obj;
        let ptm_list = pVizParams.denMap[chain].den_ptm_arr;
        //@ts-ignore
        el = el[el_lst.indexOf(ft.start)];
        let prob = prob_lst[el_lst.indexOf(ft.start)];
        let ptmsArPoint = ptm_list[el_lst.indexOf(ft.start)];
        let ptmsStr = "";

        for (let i = 0; i < ptmsArPoint.length; i++) {
          ptmsStr += "\n" + ptmsArPoint[i][0].replace("_", " ") + " probability  ~" + (ptmsArPoint[i][1] > 1 ? ptmsArPoint[i][1] / 100 : ptmsArPoint[i][1]).toFixed(2);
        }

        ui.tooltip.show(
          ui.divText(`Probability: ${prob.toFixed(2)}${ptmsStr}`),
          //@ts-ignore
          el.getBoundingClientRect().left + 10,
          //@ts-ignore
          el.getBoundingClientRect().top + 10
        );

      }).addMouseoutCallback(['D'], function (ft) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(ptmCodes, async function (ft) {

        let selectorStr = 'g.feature.' + mutcodes[ft.category.replaceAll(" ", "_")] + "."
          + ft.category.replaceAll(" ", "_").replaceAll(")", "\\)").replaceAll("(", "\\(") + ' rect.feature';
        let el = document.querySelectorAll(selectorStr);
        let el_lst = pVizParams.ptmMap[chain].ptm_el_obj[mutcodes[ft.category.replaceAll(" ", "_")]];
        //@ts-ignore
        el = el[el_lst.indexOf(ft.start)];

        ui.tooltip.show(
          ui.span([`${ft.category}`]),
          //@ts-ignore
          el.getBoundingClientRect().left + 10,
          //@ts-ignore
          el.getBoundingClientRect().top + 10
        );

      }).addMouseoutCallback(ptmCodes, function (ft) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(motCodes, async function (ft) {

        let selectorStr = 'g.feature.' + mutcodes[ft.category.replaceAll(" ", "_")] + "."
          + ft.category.replaceAll(" ", "_").replaceAll(")", "\\)").replaceAll("(", "\\(") + ' rect.feature';
        let el = document.querySelectorAll(selectorStr);
        let el_lst = pVizParams.motMap[chain].mot_el_obj[mutcodes[ft.category.replaceAll(" ", "_")]];
        //@ts-ignore
        el = el[el_lst.indexOf(ft.start)];

        ui.tooltip.show(
          ui.span([`${ft.category}`]),
          //@ts-ignore
          el.getBoundingClientRect().left + 10,
          //@ts-ignore
          el.getBoundingClientRect().top + 10
        );

      }).addMouseoutCallback(motCodes, function (ft) {
        ui.tooltip.hide();
      });

      this.pviz.FeatureDisplayer.addMouseoverCallback(obsCodes, async function (ft) {
        let selectorStr = 'g.feature.' + ft.category.replaceAll(" ", "_") + "."
          + ft.category.replaceAll(" ", "_") + ' rect.feature';

        let el = document.querySelectorAll(selectorStr);
        let el_lst = pVizParams.obsMap[chain].obs_el_obj[ft.category];
        let types = pVizParams.obsMap[chain].obs_prob_obj[ft.start][0];
        let probabilities = pVizParams.obsMap[chain].obs_prob_obj[ft.start][1];

        //@ts-ignore
        el = el[el_lst.indexOf(ft.start)];

        let ptmsStr = "";

        for (let i = 0; i < types.length; i++) {
          let prob = probabilities[i] == 0 ? "NA" : (probabilities[i] == 0.01 ? 0 : probabilities[i].toFixed(2));
          ptmsStr += "\n" + types[i] + " ~" + prob;
        }

        ui.tooltip.show(
          ui.divText(`Probability: ${ptmsStr}`),
          //@ts-ignore
          el.getBoundingClientRect().left + 10,
          //@ts-ignore
          el.getBoundingClientRect().top + 10
        );

      }).addMouseoutCallback(obsCodes, function (ft) {
        ui.tooltip.hide();
      });

      //mouse click handlers
      this.pviz.FeatureDisplayer.addClickCallback(['D'], async function (ft) {
        if (switchObj[chain][ft.start] === undefined) {
          switchObj[chain][ft.start] = {};
          switchObj[chain][ft.start]['state'] = true;
        } else {
          switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
        }
        grok.events.fireCustomEvent("selectionChanged", null);
        await pv.color(chain);
      });

      this.pviz.FeatureDisplayer.addClickCallback(ptmCodes, async function (ft) {
        if (switchObj[chain][ft.start] === undefined) {
          switchObj[chain][ft.start] = {};
          switchObj[chain][ft.start]['state'] = true;
        } else {
          switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
        }
        grok.events.fireCustomEvent("selectionChanged", null);
        await pv.color(chain);
      });

      this.pviz.FeatureDisplayer.addClickCallback(motCodes, async function (ft) {
        if (switchObj[chain][ft.start] === undefined) {
          switchObj[chain][ft.start] = {};
          switchObj[chain][ft.start]['state'] = true;
        } else {
          switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
        }
        grok.events.fireCustomEvent("selectionChanged", null);
        await pv.color(chain);
      });

      this.pviz.FeatureDisplayer.addClickCallback(obsCodes, async function (ft) {
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

  async color(chosenTracksChain) {

    let switchObj = this.selection;
    let pVizParams = this.pVizParams;

    Object.keys(switchObj).forEach((keyChain) => {

      let selectorStr = 'g.feature.data.D rect.feature';
      let denElements = document.querySelectorAll(selectorStr);
      let denNumbers = pVizParams.denMap[keyChain].den_el_obj;
      Object.keys(switchObj[keyChain]).forEach((keyPosition) => {

        let position = parseInt(keyPosition);

        if (keyChain === chosenTracksChain) {
          //densities
          if (denNumbers.indexOf(position) !== -1) {
            if (switchObj[keyChain][keyPosition]['state'] === false) {
              //@ts-ignore
              denElements[denNumbers.indexOf(position)].style.fill = pVizParams.denMap[keyChain].den_color_obj['D'][denNumbers.indexOf(position)];
            } else {
              //@ts-ignore
              denElements[denNumbers.indexOf(position)].style.fill = 'black';
            }
          }

          //predicted PTMS
          let listsPredictedPtms = [];

          this.ptmChoices.forEach(ptm => {
            let selectorStrPTM = 'g.feature.' + ptm.replace(" ", "_") + ' rect.feature';
            let elPTM = document.querySelectorAll(selectorStrPTM);
            let el_lstPTM = pVizParams.ptmMap[chosenTracksChain].ptm_el_obj[mutcodes[ptm.replace(" ", "_")]];
            listsPredictedPtms[ptm] = [elPTM, el_lstPTM];
          });

          Object.keys(listsPredictedPtms).forEach(ptm => {
            let elPTM = listsPredictedPtms[ptm][0];
            let el_lstPTM = listsPredictedPtms[ptm][1];
            if (typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1) {
              if (switchObj[keyChain][position]['state'] === false) {
                elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.ptmMap[keyChain].ptm_color_obj[mutcodes[ptm.replace(" ", "_")]][el_lstPTM.indexOf(position)];
              } else {
                elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
              }
            }
          });

          //motif PTMS
          let listsMotifsPtms = [];

          this.motChoices.forEach(ptm => {
            let selectorStrPTM = 'g.feature.' + mutcodes[ptm.replaceAll(" ", "_")] + "."
              + ptm.replaceAll(" ", "_").replaceAll(")", "\\)").replaceAll("(", "\\(") + ' rect.feature';
            let elPTM = document.querySelectorAll(selectorStrPTM);
            let el_lstPTM = pVizParams.motMap[chosenTracksChain].mot_el_obj[mutcodes[ptm.replace(" ", "_")]];
            listsMotifsPtms[ptm] = [elPTM, el_lstPTM];
          });

          Object.keys(listsMotifsPtms).forEach(ptm => {
            let elPTM = listsMotifsPtms[ptm][0];
            let el_lstPTM = listsMotifsPtms[ptm][1];
            if (typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1) {
              if (switchObj[keyChain][position]['state'] === false) {
                elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.motMap[keyChain].mot_color_obj[mutcodes[ptm.replace(" ", "_")]][el_lstPTM.indexOf(position)];
              } else {
                elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
              }
            }
          });

          //observed PTMS
          let listsObservedPtms = [];

          this.obsChoices.forEach(ptm => {
            let selectorStrPTM = 'g.feature.' + ptm.replaceAll(" ", "_") + "."
              + ptm.replaceAll(" ", "_") + ' rect.feature';
            let elPTM = document.querySelectorAll(selectorStrPTM);
            let el_lstPTM = pVizParams.obsMap[chosenTracksChain].obs_el_obj[ptm.replace(" ", "_")];
            listsObservedPtms[ptm] = [elPTM, el_lstPTM];
          });

          Object.keys(listsObservedPtms).forEach(ptm => {
            let elPTM = listsObservedPtms[ptm][0];
            let el_lstPTM = listsObservedPtms[ptm][1];
            if (typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1) {
              if (switchObj[keyChain][position]['state'] === false) {
                elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.obsMap[keyChain].obs_color_obj[ptm.replace(" ", "_")][el_lstPTM.indexOf(position)];
              } else {
                elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
              }
            }
          });
        }
      });
    });
  }

  // mapping objects for sequence rendering
  cdrMapping(cdrScheme) {
    this.cdrScheme = cdrScheme;
    let cdrMap = {}
    let chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain) => {
      let cdr_feature_map = [];
      let cdr_ranges = this.json.cdr_ranges[this.cdrScheme.value + '_CDR' + chain + '_ranges']
      if (this.cdrScheme.value !== 'default') {
        Object.values(cdr_ranges).forEach((range) => {
          cdr_feature_map.push({
            category: 'CDR region',
            type: 'CDR',
            start: range[0],
            end: range[1],
            text: '',
            improbable: true
          })
        })
      }
      cdrMap[chain] = { cdr_feature_map: cdr_feature_map }
    })

    this.pVizParams.cdrMap = (cdrMap);
  }

  parMapping(paratopes) {
    this.paratopes = paratopes;
    let parMap = {}
    let chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain) => {
      let par_feature_map = [];
      let par_color_arr = [];
      let par_el_obj = [];
      let par_prob_obj = [];
      let palette = MiscMethods.interpolateColors('(255, 255, 255)', '(255, 0, 255)', 100);

      Object.keys(this.json.parapred_predictions[chain]).forEach((index) => {
        par_feature_map.push({
          category: 'Paratope predictions',
          type: 'P',
          start: index,
          end: index,
          text: '',
          improbable: true
        })
        par_color_arr.push(palette[Math.round(this.json.parapred_predictions[chain][index] * 100)]);
        par_el_obj.push(index);
        par_prob_obj.push(this.json.parapred_predictions[chain][index]);
      })

      parMap[chain] = { par_feature_map: par_feature_map, par_color_obj: { 'P': par_color_arr }, par_el_obj: par_el_obj, par_prob_obj: par_prob_obj };
    })

    this.pVizParams.parMap = (parMap)
  }

  denMapping() {
    let denMap = {}
    let chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain) => {
      let den_feature_map = [];
      let den_color_arr = new Array(this.pVizParams.seq[chain].length).fill(-1);
      let den_ptm_arr = new Array(this.pVizParams.seq[chain].length).fill([]);
      let palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

      Object.keys(this.json.ptm_predictions[chain]).forEach((ptm) => {
        if (this.json.ptm_predictions[chain][ptm][0][1] <= 1) {
          this.json.ptm_predictions[chain][ptm].forEach(point => {
            if (!(den_feature_map.includes(point[0]))) {
              den_feature_map.push(point[0]);
            }
            den_color_arr[point[0]] = den_color_arr[point[0]] == -1 ? (point[1] > 1 ? point[1] / 100 : point[1]) : 1 - (1 - den_color_arr[point[0]]) * (1 - (point[1] > 1 ? point[1] / 100 : point[1]));
            den_ptm_arr[point[0]] = den_ptm_arr[point[0]].concat([[ptm, point[1]]]);
          })
        }
      })

      den_feature_map.sort((a, b) => a - b);
      let den_el_obj = den_feature_map.slice();
      den_feature_map = den_feature_map.map(function (ft) {
        return {
          groupSet: 'Predicted PTM density',
          category: '',
          type: 'D',
          start: ft,
          end: ft,
          text: '',
          improbable: true
        }
      });

      den_color_arr = den_color_arr.filter((x) => {
        return x !== -1
      });

      den_ptm_arr = den_ptm_arr.filter((x) => {
        return x.length > 0
      });

      let den_prob_obj = den_color_arr.slice();

      for (let i = 0; i < den_color_arr.length; i++) {
        den_color_arr[i] = palette[Math.round(den_color_arr[i] * 4)];
      }

      denMap[chain] = { den_feature_map: den_feature_map, den_color_obj: { 'D': den_color_arr }, den_el_obj: den_el_obj, den_prob_obj: den_prob_obj, den_ptm_arr: den_ptm_arr };
    })

    this.pVizParams.denMap = (denMap)
  }

  ptmMapping(ptmChoices, prob) {
    this.ptmChoices = ptmChoices;
    this.ptmProb = prob;

    let ptmMap = {}
    let chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain) => {
      let ptm_feature_map = [];
      let ptm_color_obj = {};
      let ptm_el_obj = {};
      let ptm_prob_obj = {};
      let palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

      this.ptmChoices.forEach(ptm => {
        let ptm_array = this.json.ptm_predictions[chain][ptm.replace(" ", "_")];
        if (ptm_array !== undefined) {

          let ptm_color_arr = [];
          let ptm_el_arr = [];
          let ptm_prob_arr = [];
          ptm_array.forEach(point => {
            if (point[1] > prob) {

              ptm_feature_map.push({
                groupSet: 'Predicted PTMs',
                category: ptm,
                type: mutcodes[ptm.replace(" ", "_")],
                start: point[0],
                end: point[0],
                text: '',
                improbable: true
              })
              ptm_color_arr.push(palette[point[1] > 1 ? Math.round(point[1] * 4) / 100 : Math.round(point[1] * 4)]);
              ptm_el_arr.push(point[0]);
              ptm_prob_arr.push(point[1]);
            }
          })
          if (ptm_color_arr.length > 0) {
            ptm_color_obj[mutcodes[ptm.replace(" ", "_")]] = ptm_color_arr;
            ptm_el_obj[mutcodes[ptm.replace(" ", "_")]] = ptm_el_arr;
            ptm_prob_obj[mutcodes[ptm.replace(" ", "_")]] = ptm_prob_arr;
          }
        }
      })
      ptmMap[chain] = {
        ptm_feature_map: ptm_feature_map, ptm_color_obj: ptm_color_obj,
        ptm_el_obj: ptm_el_obj, ptm_prob_obj: ptm_prob_obj
      };
    })

    this.pVizParams.ptmMap = (ptmMap);
  }

  motMapping(ptmChoices, prob) {
    this.motChoices = ptmChoices;
    let motMap = {}
    let chains = Object.keys(this.pVizParams.seq);
    chains.forEach((chain) => {
      let ptm_feature_map = [];
      let ptm_color_obj = {};
      let ptm_el_obj = {};
      let ptm_prob_obj = {};
      let palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

      ptmChoices.forEach(ptm => {
        let ptm_array = this.json.ptm_predictions[chain][ptm.replaceAll(" ", "_")];
        if (ptm_array !== undefined) {

          let ptm_color_arr = [];
          let ptm_el_arr = [];
          let ptm_prob_arr = [];
          ptm_array.forEach(point => {
            if (point[1] > prob) {

              ptm_feature_map.push({
                groupSet: 'Predicted PTMs',//'Motif PTMs',
                category: ptm,
                type: mutcodes[ptm.replaceAll(" ", "_")],
                start: point[0],
                end: point[0],
                text: '',
                improbable: true
              })
              ptm_color_arr.push(palette[point[1] > 1 ? Math.round(point[1] * 4) / 100 : Math.round(point[1] * 4)]);
              ptm_el_arr.push(point[0]);
              ptm_prob_arr.push(0.99);
            }
          })
          if (ptm_color_arr.length > 0) {
            ptm_color_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_color_arr;
            ptm_el_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_el_arr;
            ptm_prob_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_prob_arr;
          }
        }
      })
      motMap[chain] = {
        mot_feature_map: ptm_feature_map, mot_color_obj: ptm_color_obj,
        mot_el_obj: ptm_el_obj, mot_prob_obj: ptm_prob_obj
      };
    })

    this.pVizParams.motMap = (motMap);
  }

  obsMapping(ptmChoices) {
    this.obsChoices = ptmChoices;
    let obsMap = {}
    let chains = Object.keys(this.pVizParams.seq);
    let palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

    chains.forEach((chain) => {
      let obs_feature_map = [];
      let obs_color_obj = {};
      let obs_el_obj = {};
      let obs_prob_obj = {};

      ptmChoices.forEach(ptm => {
        let ptm_tree = this.jsonObs.ptm_observed[chain][ptm.replace(" ", "_")];
        if (ptm_tree !== undefined) {

          let obs_color_arr = [];
          let obs_el_arr = [];
          let obs_types_arr = [];
          let obs_types_probs_arr = [];

          Object.keys(ptm_tree).forEach((type) => {
            let point = this.jsonObs.ptm_observed[chain][ptm.replace(" ", "_")][type];
            if (!obs_el_arr.includes(point[0])) {
              obs_el_arr.push(point[0]);
            }
          });

          obs_el_arr.forEach((position) => {
            let types = [];
            let types_probs = [];
            let prob = -1;
            Object.keys(ptm_tree).forEach((type) => {
              let point = this.jsonObs.ptm_observed[chain][ptm][type];
              if (point[0] === position) {
                types.push(type);
                types_probs.push(point[1]);
                let addProb = point[1] == 0.01 ? 0 : point[1];
                prob = prob == -1 ? addProb / 100 : 1 - (1 - prob) * (1 - addProb / 100);
              }
            });

            obs_color_arr.push(palette[Math.round(prob * 4)]);
            obs_types_arr.push(types);
            obs_types_probs_arr.push(types_probs);

            obs_feature_map.push({
              groupSet: 'Observed PTMs',
              category: ptm,
              type: ptm,
              start: position,
              end: position,
              text: ptm,
              improbable: true
            })

            obs_color_obj[ptm.replace(" ", "_")] = obs_color_arr;
            obs_el_obj[ptm.replace(" ", "_")] = obs_el_arr;
            obs_prob_obj[position] = [types, types_probs];
          });
        }
      });

      obsMap[chain] = {
        obs_feature_map: obs_feature_map, obs_color_obj: obs_color_obj,
        obs_el_obj: obs_el_obj, obs_prob_obj: obs_prob_obj
      };
    });

    this.pVizParams.obsMap = (obsMap);
  }

  applyGradient(gradient_obj) {
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
    let host = chain == "H" ? this.hostH : this.hostL;

    ui.onSizeChanged(host).subscribe(async (_) => {
      await this.render(chain)
    });
  }
}

