import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import {MiscMethods} from "./misc.js"
import {PvizMethods} from "./pViz.js"
import {MsaMethods} from "./msa.js"
import {NglMethods} from "./ngl.js"
import {RepertoireBrowserPanel} from "./repertoire-browser-panel.js"


import {_package} from "./package";

export class LaunchBrowser {

    async init(view) {


        ///// MAIN BODY ////
        let inputs = new RepertoireBrowserPanel();
        await inputs.init(view);

        let ngl = new NglMethods();
        await ngl.init(view, inputs);

        let pViz = new PvizMethods();
        await pViz.init(view, inputs, ngl);

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

            ngl.stage.removeAllComponents();
            let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes, inputs.colorScheme);
            await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);
        });

        inputs.cdr_scheme.onChanged(async () => {

            ngl.stage.removeAllComponents();
            let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes, inputs.colorScheme);
            await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);

            pViz.pVizParams.cdrMap = pViz.cdrMapping(inputs.cdr_scheme.value)
            await pViz.loadSequence(inputs, 'H')
            await pViz.loadSequence(inputs, 'L')

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.paratopes.onChanged(async () => {

            ngl.stage.removeAllComponents();
            let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes, inputs.colorScheme);
            await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);

            await pViz.loadSequence(inputs, 'H')
            await pViz.loadSequence(inputs, 'L')

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.ptm_choices.onChanged(async () => {

            pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value)
            await pViz.loadSequence(inputs, 'H')
            await pViz.loadSequence(inputs, 'L')

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.ptm_prob.onChanged(async () => {

            pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value)
            await pViz.loadSequence(inputs, 'H')
            await pViz.loadSequence(inputs, 'L')

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        //inputs.msaContentChoice.onChanged(() => { msa.drawAlignments(); });

        DG.debounce(view.table.onCurrentRowChanged, 200).subscribe(() => { msa.drawAlignments(); });

    }

}