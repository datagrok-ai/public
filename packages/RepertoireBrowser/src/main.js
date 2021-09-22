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


        inputs.repChoice.onChanged(async () => {

            ngl.stage.removeAllComponents();
            let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes);
            await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);
        });

        inputs.cdr_scheme.onChanged(async () => {

            ngl.stage.removeAllComponents();
            let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes);
            await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);

            pViz.pVizParams.cdrMap = pViz.cdrMapping(inputs.cdr_scheme.value)
            await pViz.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value)
            await pViz.loadSequence(inputs.pViz_host_L, 'L', inputs.paratopes.value)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.paratopes.onChanged(async () => {

            ngl.stage.removeAllComponents();
            let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes);
            await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);

            await pViz.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value)
            await pViz.loadSequence(inputs.pViz_host_L, 'L', inputs.paratopes.value)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.ptm_choices.onChanged(async () => {

            pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value)
            await pViz.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value)
            await pViz.loadSequence(inputs.pViz_host_L, 'L', inputs.paratopes.value)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.ptm_prob.onChanged(async () => {

            pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value)
            await pViz.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value)
            await pViz.loadSequence(inputs.pViz_host_L, 'L', inputs.paratopes.value)

            MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        });

        inputs.msaContentChoice.onChanged(() => { msa.drawAlignments(); });

        DG.debounce(view.table.onCurrentRowChanged, 200).subscribe(() => { msa.drawAlignments(); });

    }

}