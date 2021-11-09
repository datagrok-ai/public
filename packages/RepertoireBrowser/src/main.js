import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import {MiscMethods} from "./misc.js"
import {PvizMethods} from "./pViz.js"
import {MsaMethods} from "./msa.js"
import {NglMethods} from "./ngl.js"
import {RepertoireBrowserPanel} from "./repertoire-browser-panel.js"
import json from "./example.json";
import jsonNumbering from "./exampleNums.json";

import {_package} from "./package";

export class LaunchBrowser {

    async init(view) {

        let jsonInit = json;
        let jsonN = jsonNumbering;

        let hNumberingStr = jsonN.heavy_numbering;
        let lNumberingStr = jsonN.light_numbering;

        let hNumbering =[];
        let lNumbering =[];

        for(let i =0; i < hNumberingStr.length;i++)
            hNumbering.push(parseInt(hNumberingStr[i].replaceAll(" ", "")));
        
        for(let i =0; i < lNumberingStr.length;i++)
            lNumbering.push(parseInt(lNumberingStr[i].replaceAll(" ", "")));


        let jsonStr = jsonInit;
        jsonStr["map_H"] = hNumbering;
        jsonStr["map_L"] = lNumbering;
        //let jsonStr = {};

        // jsonStr["heavy_seq"] = jsonInit.heavy_seq;
        // jsonStr["light_seq"] = jsonInit.light_seq;
        // jsonStr["heavy_chain"] = jsonInit.heavy_chain;
        // jsonStr["light_chain"] = jsonInit.light_chain;

        // let ptmp ={};
        // let npp ={};
        // let cr = {};
        // let lengthH = 0;
        // let lengthL = 0;
        // Object.keys(jsonInit.ptm_predictions).forEach((chain) => {
        //     let add = {};
        //     Object.keys(jsonInit.ptm_predictions[chain]).forEach((ptm) => {
        //         let addPTM = [];

        //         for(let i =0; i < jsonInit.ptm_predictions[chain][ptm].length; i++){
        //             let oindex = jsonInit.ptm_predictions[chain][ptm][i][0];
        //             let nindex = chain ==="H"? hNumbering[oindex] : lNumbering[oindex];
        //             addPTM.push([nindex,jsonInit.ptm_predictions[chain][ptm][i][1]]);

        //             if(chain === "H")
        //                 lengthH = nindex>lengthH? nindex: lengthH;
        //             else
        //                 lengthL = nindex>lengthL? nindex: lengthL;
        //         }

        //         add[ptm] = addPTM;
        //     })
        //     ptmp[chain] = add;
        // });
        // Object.keys(jsonInit.parapred_predictions).forEach((chain) => {
        //     let add = {};
        //     Object.keys(jsonInit.parapred_predictions[chain]).forEach((index) => {
        //         let nindex = chain ==="H"? hNumbering[parseInt(index)] : lNumbering[parseInt(index)];
        //         add[nindex.toString()] = jsonInit.parapred_predictions[chain][index];
        //         if(chain === "H")
        //             lengthH = nindex>lengthH? nindex: lengthH;
        //         else
        //             lengthL = nindex>lengthL? nindex: lengthL;
        //     })
        //     npp[chain] = add;
        // });
        // Object.keys(jsonInit.cdr_ranges).forEach((range) => {
        //     let addCDR = [];

        //     for(let i =0; i < jsonInit.cdr_ranges[range].length; i++){
        //         let oindex1 = jsonInit.cdr_ranges[range][i][0];
        //         let nindex1 = range.includes("CDRH")? hNumbering[oindex1] : lNumbering[oindex1];

        //         let oindex2 = jsonInit.cdr_ranges[range][i][1];
        //         let nindex2 = range.includes("CDRH")? hNumbering[oindex2] : lNumbering[oindex2];
                
        //         addCDR.push([nindex1, nindex2]);

        //         if(range.includes("CDRH"))
        //             lengthH = nindex2>lengthH? nindex2: lengthH;
        //         else
        //             lengthL = nindex2>lengthL? nindex2: lengthL;
        //     }

        //     cr[range] = addCDR;
        // });

        // jsonStr["ptm_predictions"] = jsonInit.ptm_predictions;
        // jsonStr["parapred_predictions"] = jsonInit.parapred_predictions;
        // jsonStr["cdr_ranges"] = jsonInit.cdr_ranges;
        // jsonStr["ptm_predictions_pdb"] = ptmp;
        // jsonStr["parapred_predictions_pdb"] = npp;
        // jsonStr["cdr_ranges_pdb"] = cr;
        // jsonStr["H_length"] = lengthH;
        // jsonStr["L_length"] = lengthL;

        ///// MAIN BODY ////
        let inputs = new RepertoireBrowserPanel();
        await inputs.init(view, jsonStr);

        let ngl = new NglMethods();
        await ngl.init(view, inputs, jsonStr);

        let pViz = new PvizMethods();
        await pViz.init(view, inputs, ngl, jsonStr);

        // let msa = new MsaMethods();
        // msa.init(view, inputs);

        grok.events.onTooltipShown.subscribe((args) => {
            if (args.args.context instanceof DG.Column) {
                switch(args.args.context.name){
                    case 'cdr length': 
                    args.args.element.innerHTML = "Complementarity-determining region length<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding, length is expressed in number of residuals<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;
                  
                    case 'surface cdr hydrophobicity': 
                    args.args.element.innerHTML = "Complementarity-determining region surface hydrophobicity<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding, hydrophobicity is expressed as polar surface in A<sup>2</sup><br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;

                    case 'positive cdr charge': 
                    args.args.element.innerHTML = "Complementarity-determining region positive charge<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;

                    case 'negative cdr charge': 
                    args.args.element.innerHTML = "Complementarity-determining region negative charge<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;

                    case 'sfvcsp': 
                    args.args.element.innerHTML = "Structural Fv Charge Symmetry Parameter<br>" +
                    "is a measure of aggregate-inducing electrostatic attraction of the entire variable region<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;
                  
                    default:
                    break;
                }
            }
        });

        inputs.repChoice.onChanged(async () => {

            //ngl.stage.removeAllComponents();
            //let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes, inputs.colorScheme);
            //await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);

            await pViz.loadSequence(inputs, 'H', jsonStr, true)
            await pViz.loadSequence(inputs, 'L', jsonStr, true)
        });

        inputs.cdr_scheme.onChanged(async () => {

            // ngl.stage.removeAllComponents();
            // let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes, inputs.colorScheme);
            // await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);

            pViz.pVizParams.cdrMap = pViz.cdrMapping(inputs.cdr_scheme.value, jsonStr)
            await pViz.loadSequence(inputs, 'H', jsonStr)
            await pViz.loadSequence(inputs, 'L', jsonStr)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.paratopes.onChanged(async () => {

            // ngl.stage.removeAllComponents();
            // let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes, inputs.colorScheme);
            // await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);

            await pViz.loadSequence(inputs, 'H', jsonStr)
            await pViz.loadSequence(inputs, 'L', jsonStr)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.ptm_choices.onChanged(async () => {

            pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value, jsonStr )
            await pViz.loadSequence(inputs, 'H', jsonStr)
            await pViz.loadSequence(inputs, 'L', jsonStr)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.ptm_motif_choices.onChanged(async () => {

            pViz.pVizParams.ptmMotifsMap = pViz.ptmMotifsMapping(inputs.ptm_motif_choices.value, inputs.ptm_prob.value, jsonStr)
            await pViz.loadSequence(inputs, 'H', jsonStr)
            await pViz.loadSequence(inputs, 'L', jsonStr)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.ptm_prob.onChanged(async () => {

            pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value, jsonStr)
            await pViz.loadSequence(inputs, 'H', jsonStr)
            await pViz.loadSequence(inputs, 'L', jsonStr)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        //inputs.msaContentChoice.onChanged(() => { msa.drawAlignments(); });

        //DG.debounce(view.table.onCurrentRowChanged, 200).subscribe(() => { msa.drawAlignments(); });

    }

}