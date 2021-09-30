import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import json from "./TPP000153303.json";
import mutcodes from "./mutcodes.json";
import {MiscMethods} from "./misc.js"


export class PvizMethods {

    async init(view, inputs, ngl) {

        this.ngl = ngl
        this.pviz = window.pviz;
        this.pVizParams = {};
        this.pVizParams.seq = {'H':json.heavy_seq, 'L': json.light_seq}
        this.pVizParams.ptmMap = this.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value)
        this.pVizParams.denMap = this.ptmDenMapping()
        this.pVizParams.parMap = this.paratopeMapping()
        this.pVizParams.cdrMap = this.cdrMapping(inputs.cdr_scheme.value)

        await this.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value);

        await this.pvizResize(inputs.pViz_host_H, 'H', inputs.paratopes);
        await this.pvizResize(inputs.pViz_host_L, 'L', inputs.paratopes);
        MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);

    }

    // mapping objects for sequence rendering
    ptmMapping(ptm_choices, prob) {

        let ptmMap = {}
        let chains = Object.keys(this.pVizParams.seq);
        chains.forEach((chain) => {
            let ptm_feature_map = [];
            let ptm_color_obj = {};
            let ptm_el_obj = {};
            let ptm_prob_obj = {};
            let palette = MiscMethods.interpolateColors('(255, 255, 0)','(255, 0, 0)', 5);

            ptm_choices.forEach(ptm => {
                let ptm_array = json.ptm_predictions[chain][ptm];
                if (ptm_array !== undefined) {

                    let ptm_color_arr = [];
                    let ptm_el_arr = [];
                    let ptm_prob_arr = [];
                    ptm_array.forEach(point => {
                        if(point[1] > prob) {

                            ptm_feature_map.push({
                                groupSet: 'PTMs',
                                category : mutcodes[ptm],
                                type : mutcodes[ptm],
                                start : point[0],
                                end : point[0],
                                text : mutcodes[ptm],
                                improbable : true
                            })
                            ptm_color_arr.push(palette[Math.round(point[1]*4)])
                            ptm_el_arr.push(point[0]);
                            ptm_prob_arr.push(point[1]);
                        }
                    })
                    if (ptm_color_arr.length > 0) {
                        ptm_color_obj[mutcodes[ptm]] = ptm_color_arr;
                        ptm_el_obj[mutcodes[ptm]] = ptm_el_arr;
                        ptm_prob_obj[mutcodes[ptm]] = ptm_prob_arr;
                    }
                }
            })
            ptmMap[chain] = {ptm_feature_map: ptm_feature_map, ptm_color_obj: ptm_color_obj,
                ptm_el_obj: ptm_el_obj, ptm_prob_obj: ptm_prob_obj};
        })

        return(ptmMap);
    }

    ptmDenMapping() {

        let denMap = {}
        let chains = Object.keys(this.pVizParams.seq);
        chains.forEach((chain) => {
            let den_feature_map = [];
            let den_el_obj  = [];
            let den_prob_obj = [];
            let den_color_arr = new Array(this.pVizParams.seq[chain].length).fill(1);
            let palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

            Object.values(json.ptm_predictions[chain]).forEach((ptm_array) => {
                ptm_array.forEach(point => {
                    if (!(den_feature_map.includes(point[0]))) {
                        den_feature_map.push(point[0]);
                    }
                    den_color_arr[point[0]] = den_color_arr[point[0]] * (1 - point[1])
                    den_el_obj.push(point[0]);
                    den_prob_obj.push(point[1]);
                })
            })

            den_el_obj.sort(function(a, b) {return a - b; });;

            den_feature_map.sort((a, b) => a - b);
            den_feature_map = den_feature_map.map(function (ft) {
                return {
                    groupSet: 'PTM density',
                    category: '',
                    type: 'D',
                    start: ft,
                    end: ft,
                    text: '',
                    improbable: true
                }
            });

            den_color_arr = den_color_arr.filter((x) => {
                return x !== 1
            });
            for (let i = 0; i < den_color_arr.length; i++) {
                den_color_arr[i] = palette[Math.round((1 - den_color_arr[i]) * 4)];
            }

            denMap[chain] = {den_feature_map: den_feature_map, den_color_obj: {'D': den_color_arr}, den_el_obj: den_el_obj, den_prob_obj:den_prob_obj}
        })


        return (denMap)
    }

    paratopeMapping() {

        let parMap = {}
        let chains = Object.keys(this.pVizParams.seq);
        chains.forEach((chain) => {
            let par_feature_map = [];
            let par_color_arr = [];
            let par_el_obj  = [];
            let par_prob_obj = [];
            let palette = MiscMethods.interpolateColors('(255, 255, 255)', '(255, 0, 255)', 100);

            Object.keys(json.parapred_predictions[chain]).forEach((index) => {
                par_feature_map.push({
                    category: 'Paratope predictions',
                    type: 'P',
                    start: index,
                    end: index,
                    text: '',
                    improbable: true
                })
                par_color_arr.push(palette[Math.round(json.parapred_predictions[chain][index] * 100)]);
                par_el_obj.push(index);
                par_prob_obj.push(json.parapred_predictions[chain][index]);
            })

            parMap[chain] = {par_feature_map: par_feature_map, par_color_obj: {'P': par_color_arr}, par_el_obj : par_el_obj, par_prob_obj:par_prob_obj}
        })

        return(parMap)
    }

    cdrMapping(cdr_scheme) {

        let cdrMap = {}
        let chains = Object.keys(this.pVizParams.seq);
        chains.forEach((chain) => {
            let cdr_feature_map = [];
            let cdr_ranges = json.cdr_ranges[cdr_scheme + '_CDR' + chain + '_ranges']
            if (cdr_scheme !== 'default') {
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
            cdrMap[chain] = {cdr_feature_map: cdr_feature_map}
        })

        return (cdrMap)
    }

    applyGradient(gradient_obj) {
        Object.keys(gradient_obj).forEach((ptm_track) => {
            let selectorStr = 'g.feature.' + ptm_track + ' rect.feature';
            let el = document.querySelectorAll(selectorStr);
            for (let i = 0; i < el.length; i++) {
                el[i].style.fill = gradient_obj[ptm_track][i];
            }
        })
    }

    // main sequence rendering func
    async loadSequence(host, chain, paratopes) {

        if( $(host).width() !== 0) {
            let seq = this.pVizParams.seq[chain]

            let seqEntry = new this.pviz.SeqEntry({
                sequence: seq
            });
            new this.pviz.SeqEntryAnnotInteractiveView({
                model: seqEntry,
                collapsible: true,
                el: host
            }).render();

            let mod_codes = Object.keys(this.pVizParams.ptmMap[chain].ptm_color_obj);

            mod_codes.forEach((mod) => {
                this.pviz.FeatureDisplayer.trackHeightPerCategoryType[mod] = 1.5;
                this.pviz.FeatureDisplayer.setStrikeoutCategory(mod);
            });

            this.pviz.FeatureDisplayer.trackHeightPerCategoryType['P'] = 1.5;
            this.pviz.FeatureDisplayer.setStrikeoutCategory('P');
            

            let switchObj = {'H':{}, 'L':{}}
            let pVizParams = this.pVizParams;
            let stage = this.ngl.stage;
            this.pviz.FeatureDisplayer.addClickCallback (mod_codes, async function(ft) {

                let selectorStr = 'g.feature.' + ft.category + ' rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.ptmMap[chain].ptm_el_obj[ft.category];

                let sidechains = `${ft.start + 1} and :${chain} and (not backbone or .CA or (PRO and .N))`
                let r;
                if (switchObj[chain][ft.start] === undefined) {
                    r = stage.compList[0].addRepresentation("ball+stick", {sele: sidechains});
                    switchObj[chain][ft.start] = {};
                    switchObj[chain][ft.start]['state'] = false;
                    switchObj[chain][ft.start]['rep'] = r
                    el[el_lst.indexOf(ft.start)].style.fill = 'black';
                } else {
                    r = switchObj[chain][ft.start]['rep'];
                    r.setVisibility(switchObj[chain][ft.start]['state']);
                    if (switchObj[chain][ft.start]['state'] === false) {
                        el[el_lst.indexOf(ft.start)].style.fill = pVizParams.ptmMap[chain].ptm_color_obj[ft.category][el_lst.indexOf(ft.start)];
                    } else {
                        el[el_lst.indexOf(ft.start)].style.fill = 'black';
                    }
                    switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
                }
            })

            this.pviz.FeatureDisplayer.addMouseoverCallback(mod_codes, async function(ft) {

                let selectorStr = 'g.feature.' + ft.category + ' rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.ptmMap[chain].ptm_el_obj[ft.category];
                let prob_lst = pVizParams.ptmMap[chain].ptm_prob_obj[ft.category];
                el = el[el_lst.indexOf(ft.start)];
                let prob =  prob_lst[el_lst.indexOf(ft.start)];

                ui.tooltip.show(
                    ui.span([`${ft.category}: Pr ~${prob.toFixed(2)}`]),
                    el.getBoundingClientRect().left + 10,
                    el.getBoundingClientRect().top + 10
                );

            }).addMouseoutCallback(mod_codes, function(ft) {
                ui.tooltip.hide();
            });

            this.pviz.FeatureDisplayer.addMouseoverCallback(['P'], async function(ft) {

                let selectorStr = 'g.feature.data.P.Paratope_predictions rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.parMap[chain].par_el_obj;
                let prob_lst = pVizParams.parMap[chain].par_prob_obj;
                el = el[el_lst.indexOf((ft.start).toString())];
                let prob =  prob_lst[el_lst.indexOf((ft.start).toString())];
                
                ui.tooltip.show(
                    ui.span([`Probability: ${prob.toFixed(2)}`]),
                    el.getBoundingClientRect().left + 10,
                    el.getBoundingClientRect().top + 10
                );

            }).addMouseoutCallback(mod_codes, function(ft) {
                ui.tooltip.hide();
            });

            this.pviz.FeatureDisplayer.addMouseoverCallback(['D'], async function(ft) {

                let selectorStr = 'g.feature.data.D rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.denMap[chain].den_el_obj;
                let prob_lst = pVizParams.denMap[chain].den_prob_obj;
                el = el[el_lst.indexOf(ft.start)];
                let prob =  prob_lst[el_lst.indexOf(ft.start)];
                
                ui.tooltip.show(
                    ui.span([`Probability: ${prob.toFixed(2)}`]),
                    el.getBoundingClientRect().left + 10,
                    el.getBoundingClientRect().top + 10
                );

            }).addMouseoutCallback(mod_codes, function(ft) {
                ui.tooltip.hide();
            });

            if (paratopes === true) {
                seqEntry.addFeatures(this.pVizParams.parMap[chain].par_feature_map);
            }
            seqEntry.addFeatures(this.pVizParams.ptmMap[chain].ptm_feature_map);
            seqEntry.addFeatures(this.pVizParams.denMap[chain].den_feature_map);
            seqEntry.addFeatures(this.pVizParams.cdrMap[chain].cdr_feature_map);
            this.applyGradient(this.pVizParams.ptmMap[chain].ptm_color_obj);
            this.applyGradient(this.pVizParams.denMap[chain].den_color_obj);
            this.applyGradient(this.pVizParams.parMap[chain].par_color_obj);

        }
    }

    // resize handle
    async pvizResize(host, chain, paratopes) {
        ui.onSizeChanged(host).subscribe(async (_) => {
            await this.loadSequence(host, chain, paratopes.value)
        });
    }

}

