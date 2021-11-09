import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import mutcodes from "./mutcodes.json";
import {MiscMethods} from "./misc.js"


export class PvizMethods {

    async init(view, inputs, ngl, json) {

        this.ngl = ngl
        this.pviz = window.pviz;
        this.pVizParams = {};
        this.pVizParams.seq = {'H':json.heavy_seq, 'L': json.light_seq};
        this.pVizParams.ptmMap = this.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value, json);
        this.pVizParams.ptmMotifsMap = this.ptmMotifsMapping(inputs.ptm_motif_choices.value, inputs.ptm_prob.value, json)
        this.pVizParams.denMap = this.ptmDenMapping(json);
        this.pVizParams.parMap = this.paratopeMapping(json);
        this.pVizParams.cdrMap = this.cdrMapping(inputs.cdr_scheme.value, json);

        await this.loadSequence(inputs, 'H', json);

        await this.pvizResize(inputs, 'H', json);
        await this.pvizResize(inputs, 'L', json);
        MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
    }

    // mapping objects for sequence rendering
    ptmMapping(ptm_choices, prob, json) {
        let ptmMap = {}
        let chains = Object.keys(this.pVizParams.seq);
        chains.forEach((chain) => {
            let ptm_feature_map = [];
            let ptm_color_obj = {};
            let ptm_el_obj = {};
            let ptm_prob_obj = {};
            let palette = MiscMethods.interpolateColors('(255, 255, 0)','(255, 0, 0)', 5);

            ptm_choices.forEach(ptm => {
                let ptm_array = json.ptm_predictions[chain][ptm.replaceAll(" ", "_")];
                if (ptm_array !== undefined) {

                    let ptm_color_arr = [];
                    let ptm_el_arr = [];
                    let ptm_prob_arr = [];
                    ptm_array.forEach(point => {
                        if(point[1] > prob) {

                            ptm_feature_map.push({
                                groupSet: 'Predicted PTMs',
                                category : ptm,
                                type : mutcodes[ptm.replaceAll(" ", "_")],
                                start : point[0],
                                end : point[0],
                                text : ptm,
                                improbable : true
                            })
                            ptm_color_arr.push(palette[point[1] > 1 ? Math.round(point[1]*4)/100: Math.round(point[1]*4)]);
                            ptm_el_arr.push(point[0]);
                            ptm_prob_arr.push(point[1]);
                        }
                    })
                    if (ptm_color_arr.length > 0) {
                        ptm_color_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_color_arr;
                        ptm_el_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_el_arr;
                        ptm_prob_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_prob_arr;
                    }
                }
            })
            ptmMap[chain] = {ptm_feature_map: ptm_feature_map, ptm_color_obj: ptm_color_obj,
                ptm_el_obj: ptm_el_obj, ptm_prob_obj: ptm_prob_obj};
        })

        return(ptmMap);
    }

    ptmMotifsMapping(ptm_choices, prob, json) {
        let ptmMotifsMap = {}
        let chains = Object.keys(this.pVizParams.seq);
        chains.forEach((chain) => {
            let ptm_feature_map = [];
            let ptm_color_obj = {};
            let ptm_el_obj = {};
            let ptm_prob_obj = {};
            let palette = MiscMethods.interpolateColors('(255, 255, 0)','(255, 0, 0)', 5);

            ptm_choices.forEach(ptm => {
                let ptm_array = json.ptm_predictions[chain][ptm.replaceAll(" ", "_")];
                if (ptm_array !== undefined) {

                    let ptm_color_arr = [];
                    let ptm_el_arr = [];
                    let ptm_prob_arr = [];
                    ptm_array.forEach(point => {
                        if(point[1] > prob) {

                            ptm_feature_map.push({
                                groupSet: 'Predicted PTMs',//'Motif PTMs',
                                category : ptm,
                                type : mutcodes[ptm.replaceAll(" ", "_")],
                                start : point[0],
                                end : point[0],
                                text : ptm,
                                improbable : true
                            })
                            ptm_color_arr.push(palette[point[1] > 1 ? Math.round(point[1]*4)/100: Math.round(point[1]*4)]);
                            ptm_el_arr.push(point[0]);
                            ptm_prob_arr.push(point[1]);
                        }
                    })
                    if (ptm_color_arr.length > 0) {
                        ptm_color_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_color_arr;
                        ptm_el_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_el_arr;
                        ptm_prob_obj[mutcodes[ptm.replaceAll(" ", "_")]] = ptm_prob_arr;
                    }
                }
            })
            ptmMotifsMap[chain] = {ptm_feature_map: ptm_feature_map, ptm_color_obj: ptm_color_obj,
                ptm_el_obj: ptm_el_obj, ptm_prob_obj: ptm_prob_obj};
        })

        return(ptmMotifsMap);
    }

    ptmDenMapping(json) {
        let denMap = {}
        let chains = Object.keys(this.pVizParams.seq);
        chains.forEach((chain) => {
            //let length = chain === "H"? json.H_length : json.L_length;

            let den_feature_map = [];
            let den_color_arr = new Array(this.pVizParams.seq[chain].length).fill(-1);
            let den_ptm_arr = new Array(this.pVizParams.seq[chain].length).fill([]);
            let palette = MiscMethods.interpolateColors('(255, 255, 0)', '(255, 0, 0)', 5);

            Object.keys(json.ptm_predictions[chain]).forEach((ptm) => {

                if(json.ptm_predictions[chain][ptm][0][1] <=1 ){
                    json.ptm_predictions[chain][ptm].forEach(point => {

                        let index = point[0];
                        if (!(den_feature_map.includes(index))) {
                            den_feature_map.push(index);
                        }
                        den_color_arr[index] = den_color_arr[index] == -1 ? (point[1] > 1 ? point[1]/100: point[1]) : 1 - (1 - den_color_arr[index]) * (1 -  (point[1] > 1 ? point[1]/100: point[1]));
                        den_ptm_arr[index] = den_ptm_arr[index].concat([[ptm, point[1]]]);
                    })
                }
            });

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

            den_ptm_arr = den_ptm_arr.filter((x) =>{
                return x.length > 0
            });

            let den_prob_obj = den_color_arr.slice();
 
            for (let i = 0; i < den_color_arr.length; i++) {
                den_color_arr[i] = palette[Math.round(den_color_arr[i] * 4)];
            }

            denMap[chain] = {den_feature_map: den_feature_map, den_color_obj: {'D': den_color_arr}, den_el_obj:den_el_obj, den_prob_obj:den_prob_obj, den_ptm_arr:den_ptm_arr};
        })

        return (denMap)
    }

    paratopeMapping(json) {
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

            parMap[chain] = {par_feature_map: par_feature_map, par_color_obj: {'P': par_color_arr}, par_el_obj : par_el_obj, par_prob_obj:par_prob_obj};
        })

        return(parMap)
    }

    cdrMapping(cdr_scheme, json) {

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
    async loadSequence(inputs, chain, json, reLoad = false) {

        let host = chain == "H" ? inputs.pViz_host_H : inputs.pViz_host_L;

        if( $(host).width() !== 0) {
            let seq = this.pVizParams.seq[chain];

            let seqEntry = new this.pviz.SeqEntry({
                sequence: seq
            });
            new this.pviz.SeqEntryAnnotInteractiveView({
                model: seqEntry,
                collapsible: true,
                el: host
            }).render();

            let mod_codes = Object.keys(this.pVizParams.ptmMap[chain].ptm_color_obj);
            let mod_motifs_codes = Object.keys(this.pVizParams.ptmMotifsMap[chain].ptm_color_obj);

            mod_codes.forEach((mod) => {
                this.pviz.FeatureDisplayer.trackHeightPerCategoryType[mod] = 1.5;
                this.pviz.FeatureDisplayer.setStrikeoutCategory(mod);
            });

            mod_motifs_codes.forEach((mod) => {
                this.pviz.FeatureDisplayer.trackHeightPerCategoryType[mod] = 1.5;
                this.pviz.FeatureDisplayer.setStrikeoutCategory(mod);
            });

            this.pviz.FeatureDisplayer.trackHeightPerCategoryType['P'] = 1.5;
            this.pviz.FeatureDisplayer.setStrikeoutCategory('P');
            

            let switchObj = inputs.pVizNglRelation;
            let pVizParams = this.pVizParams;
            let pv = this;
            
            this.pviz.FeatureDisplayer.addClickCallback (mod_codes, async function(ft) {               
                if (switchObj[chain][ft.start] === undefined) {     
                    switchObj[chain][ft.start] = {};
                    switchObj[chain][ft.start]['state'] = true;
                } else {
                    switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
                }

                await pv.consistentlyColorpVizNGL(inputs, chain, json, reLoad);
            })

            this.pviz.FeatureDisplayer.addClickCallback (mod_motifs_codes, async function(ft) {               
                if (switchObj[chain][ft.start] === undefined) {
                    switchObj[chain][ft.start] = {};
                    switchObj[chain][ft.start]['state'] = true;
                } else {
                    switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
                }

                await pv.consistentlyColorpVizNGL(inputs, chain, json, reLoad);
            })

            this.pviz.FeatureDisplayer.addClickCallback (['D'], async function(ft) {
                if (switchObj[chain][ft.start] === undefined) {
                    switchObj[chain][ft.start] = {};
                    switchObj[chain][ft.start]['state'] = true;
                } else {
                    switchObj[chain][ft.start]['state'] = !switchObj[chain][ft.start]['state']
                }

                await pv.consistentlyColorpVizNGL(inputs, chain, json, reLoad);
            })

            this.pviz.FeatureDisplayer.addMouseoverCallback(mod_codes, async function(ft) {

                let selectorStr = 'g.feature.' + ft.category.replaceAll(" ", "_") + ' rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.ptmMap[chain].ptm_el_obj[mutcodes[ft.category.replaceAll(" ", "_")]];
                let prob_lst = pVizParams.ptmMap[chain].ptm_prob_obj[mutcodes[ft.category.replaceAll(" ", "_")]];
                el = el[el_lst.indexOf(ft.start)];
                let prob =  prob_lst[el_lst.indexOf(ft.start)];

                ui.tooltip.show(
                    ui.span([`${ft.category} probability ${prob.toFixed(2)}`]),
                    el.getBoundingClientRect().left + 10,
                    el.getBoundingClientRect().top + 10
                );

            }).addMouseoutCallback(mod_codes, function(ft) {
                ui.tooltip.hide();
            });

            this.pviz.FeatureDisplayer.addMouseoverCallback(mod_motifs_codes, async function(ft) {

                let selectorStr = 'g.feature.' + mutcodes[ft.category.replaceAll(" ", "_")] + "." 
                                        + ft.category.replaceAll(" ", "_").replaceAll(")", "\\)").replaceAll("(", "\\(");
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.ptmMotifsMap[chain].ptm_el_obj[mutcodes[ft.category.replaceAll(" ", "_")]];
                el = el[el_lst.indexOf(ft.start)];

                ui.tooltip.show(
                    ui.span([`${ft.category}`]),
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

            }).addMouseoutCallback(['P'], function(ft) {
                ui.tooltip.hide();
            });

            this.pviz.FeatureDisplayer.addMouseoverCallback(['D'], async function(ft) {

                let selectorStr = 'g.feature.data.D rect.feature';
                let el = document.querySelectorAll(selectorStr);
                let el_lst = pVizParams.denMap[chain].den_el_obj;
                let prob_lst = pVizParams.denMap[chain].den_prob_obj;
                let ptm_list =  pVizParams.denMap[chain].den_ptm_arr;
                el = el[el_lst.indexOf(ft.start)];
                let prob =  prob_lst[el_lst.indexOf(ft.start)];
                let ptmsArPoint = ptm_list[el_lst.indexOf(ft.start)];
                let ptmsStr = "";

                for(let i = 0; i < ptmsArPoint.length; i++){
                    ptmsStr += "\n" + ptmsArPoint[i][0].replaceAll("_", " ") + " probability  " + (ptmsArPoint[i][1] > 1 ? ptmsArPoint[i][1]/100: ptmsArPoint[i][1]).toFixed(2);
                }
                
                ui.tooltip.show(
                    ui.divText(`Probability: ${prob.toFixed(2)}${ptmsStr}`),
                    el.getBoundingClientRect().left + 10,
                    el.getBoundingClientRect().top + 10
                );

            }).addMouseoutCallback(['D'], function(ft) {
                ui.tooltip.hide();
            });

            if (inputs.paratopes.value === true) {
               seqEntry.addFeatures(this.pVizParams.parMap[chain].par_feature_map);
            }
            seqEntry.addFeatures(this.pVizParams.ptmMap[chain].ptm_feature_map);
            seqEntry.addFeatures(this.pVizParams.ptmMotifsMap[chain].ptm_feature_map);
            seqEntry.addFeatures(this.pVizParams.denMap[chain].den_feature_map);
            seqEntry.addFeatures(this.pVizParams.cdrMap[chain].cdr_feature_map);
            this.applyGradient(this.pVizParams.ptmMap[chain].ptm_color_obj);
            this.applyGradient(this.pVizParams.ptmMotifsMap[chain].ptm_color_obj);
            this.applyGradient(this.pVizParams.denMap[chain].den_color_obj);
            this.applyGradient(this.pVizParams.parMap[chain].par_color_obj);
            this.consistentlyColorpVizNGL(inputs, chain, json, reLoad);
        }
    }

    // resize handle
    async pvizResize(inputs, chain, json) {
        let host = chain == "H" ? inputs.pViz_host_H : inputs.pViz_host_L;

        ui.onSizeChanged(host).subscribe(async (_) => {
            await this.loadSequence(inputs, chain, json)
        });
    }

    async consistentlyColorpVizNGL(inputs, chosenTracksChain, json, reLoad){

        let switchObj = inputs.pVizNglRelation;
        let pVizParams = this.pVizParams;
        let colorScheme = inputs.colorScheme;
        let ngl = this.ngl

        let col_heavy_chain = colorScheme["col_heavy_chain"];
        let col_light_chain = colorScheme["col_light_chain"];
        let col_cdr = colorScheme["col_cdr"];
        let col_para = colorScheme["col_para"];
        let col_partopes_low = colorScheme["col_partopes_low"]; //col_para in rgb
        let col_partopes_high = colorScheme["col_partopes_high"];
        let col_highlight =  (inputs.cdr_scheme.value === 'default' || inputs.paratopes.value === true) ? 
                                                colorScheme["col_highlight"]: colorScheme["col_highlight_cdr"];

        //highlights in NGL
        let scheme_buffer = [];
               
        Object.keys(switchObj).forEach((keyChain) =>{
            Object.keys(switchObj[keyChain]).forEach((keyFtStart) =>{

                let index = keyChain === "H"? json.map_H[keyFtStart] : json.map_L[keyFtStart];
                if(switchObj[keyChain][keyFtStart]['state'] === true){
                    scheme_buffer.push([col_highlight, `${index} and :${keyChain}`]);
                }
            });
        });
        
        let schemeId;
        if(inputs.paratopes.value === true){
            let palette = MiscMethods.interpolateColors(col_partopes_low, col_partopes_high, 100);
            Object.keys(json.parapred_predictions).forEach((chain) => {
                Object.keys(json.parapred_predictions[chain]).forEach((index) => {

                    let nindex = chain === "H"? json.map_H[index] : json.map_L[index];

                    scheme_buffer.push([
                        palette[Math.round(json.parapred_predictions[chain][index] * 100)],
                        `${nindex} and :${chain}`
                    ]);
                })

            })
            scheme_buffer.push([col_para, "* and :H"]);
            scheme_buffer.push([col_para, "* and :L"]);
            schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
        }
        else{
            if (inputs.cdr_scheme.value === 'default') {
                scheme_buffer.push([col_heavy_chain, "* and :H"]);
                scheme_buffer.push([col_light_chain, "* and :L"]);
                schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
            } else{
                Object.keys(json.cdr_ranges).forEach((str) => {
                    if (str.includes(inputs.cdr_scheme.value + '_CDRH')) {
                        let str_buffer = '';
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            let nindex1 =  json.map_H[json.cdr_ranges[str][i][0]];
                            let nindex2 =  json.map_H[json.cdr_ranges[str][i][1]];

                            str_buffer = str_buffer + ` or ${nindex1}-${nindex2} and :H`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push([col_cdr, str_buffer]);
                        scheme_buffer.push([col_heavy_chain, "* and :H"]);

                    } else if (str.includes(inputs.cdr_scheme.value + '_CDRL')) {
                        let str_buffer = ''
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            let nindex1 =  json.map_L[json.cdr_ranges[str][i][0]];
                            let nindex2 =  json.map_L[json.cdr_ranges[str][i][1]];

                            str_buffer = str_buffer + ` or ${nindex1}-${nindex2} and :L`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push([col_cdr, str_buffer]);
                        scheme_buffer.push([col_light_chain, "* and :L"]);
                    }
                });
                schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
            }
        }

        if(reLoad){
            ngl.stage.removeAllComponents();
            await ngl.loadPdb(ngl.path, inputs.repChoice, {color: schemeId});

                //recovering ball and stick residual at cartoon view
            if(inputs.repChoice.value === "cartoon"){
                Object.keys(switchObj).forEach((keyChain) =>{
                    Object.keys(switchObj[keyChain]).forEach((keyFtStart) =>{

                        let index = keyChain === "H"? json.map_H[keyFtStart] : json.map_L[keyFtStart];

                        let selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
                        let addedRepresentation = ngl.stage.compList[0].addRepresentation("ball+stick", {sele: selection}); 
                        switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;
                                    
                        if(switchObj[keyChain][keyFtStart]['state'] !== true){
                            addedRepresentation.setVisibility(false);
                        }
                    });
                });
            }
        } else{
            ngl.stage.compList[0].addRepresentation(inputs.repChoice.value, {color: schemeId});

            if(inputs.repChoice.value === "cartoon"){
                Object.keys(switchObj).forEach((keyChain) =>{
                    Object.keys(switchObj[keyChain]).forEach((keyFtStart) =>{

                        let index = keyChain === "H"? json.map_H[keyFtStart] : json.map_L[keyFtStart];

                        if (typeof switchObj[keyChain][keyFtStart]['representation'] === 'undefined'){
                            let selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
                            let addedRepresentation = ngl.stage.compList[0].addRepresentation("ball+stick", {sele: selection}); 
                            switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;
                        }
       
                        if(switchObj[keyChain][keyFtStart]['state'] !== true){
                            switchObj[keyChain][keyFtStart]['representation'].setVisibility(false);
                        } else{
                            switchObj[keyChain][keyFtStart]['representation'].setVisibility(true);
                        }
                    });
                });
            }
        }
       
        //colors of selected pViz
        let selectorStr = 'g.feature.data.D rect.feature';
        let el = document.querySelectorAll(selectorStr);

        let lists_ptm = [];
        let lists_ptm_motifs = [];

        inputs.ptm_choices.value.forEach(ptm => {
            let selectorStrPTM = 'g.feature.' + ptm.replaceAll(" ", "_") + ' rect.feature';
            let elPTM = document.querySelectorAll(selectorStrPTM);
            let el_lstPTM = pVizParams.ptmMap[chosenTracksChain].ptm_el_obj[mutcodes[ptm.replaceAll(" ", "_")]];
            lists_ptm[ptm] =[elPTM, el_lstPTM];
        });

        inputs.ptm_motif_choices.value.forEach(ptm => {
            let selectorStrPTM = 'g.feature.' + mutcodes[ptm.replaceAll(" ", "_")] + "." 
                     + ptm.replaceAll(" ", "_").replaceAll(")", "\\)").replaceAll("(", "\\(");
            let elPTM = document.querySelectorAll(selectorStrPTM);
            let el_lstPTM = pVizParams.ptmMotifsMap[chosenTracksChain].ptm_el_obj[mutcodes[ptm.replaceAll(" ", "_")]];
            lists_ptm_motifs[ptm] =[elPTM, el_lstPTM];
        });
        
        Object.keys(switchObj).forEach((keyChain) =>{
            let el_lst = pVizParams.denMap[keyChain].den_el_obj;

            Object.keys(switchObj[keyChain]).forEach((keyFtStart) =>{

                let position = parseInt(keyFtStart);

                if (keyChain === chosenTracksChain)
                {
                    if (switchObj[keyChain][keyFtStart]['state'] === false) {
                        el[el_lst.indexOf(position)].style.fill = pVizParams.denMap[keyChain].den_color_obj['D'][el_lst.indexOf(position)];
                    } else {
                        el[el_lst.indexOf(position)].style.fill = 'black';
                    }

                    Object.keys(lists_ptm).forEach(ptm =>{
                        let elPTM = lists_ptm[ptm][0];
                        let el_lstPTM = lists_ptm[ptm][1];
                        if(typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1){
                            if (switchObj[keyChain][position]['state'] === false) {
                                elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.ptmMap[keyChain].ptm_color_obj[mutcodes[ptm.replaceAll(" ", "_")]][el_lstPTM.indexOf(position)];
                            } else {
                                elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
                            }
                        }
                    });

                    Object.keys(lists_ptm_motifs).forEach(ptm =>{
                        let elPTM = lists_ptm_motifs[ptm][0];
                        let el_lstPTM = lists_ptm_motifs[ptm][1];
                        if(typeof el_lstPTM !== 'undefined' && el_lstPTM.indexOf(position) !== -1){
                            if (switchObj[keyChain][position]['state'] === false) {
                                elPTM[el_lstPTM.indexOf(position)].style.fill = pVizParams.ptmMotifsMap[keyChain].ptm_color_obj[mutcodes[ptm.replaceAll(" ", "_")]][el_lstPTM.indexOf(position)];
                            } else {
                                elPTM[el_lstPTM.indexOf(position)].style.fill = 'black';
                            }
                        }
                    });
                }
            });
        });
    }
}

