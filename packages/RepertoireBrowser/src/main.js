import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import json from "./TPP000153303.json";
import MiscMethods from "./misc.js"
import PvizMethods from "./pViz.js"
import MsaMethods from "./msa.js"
import NglMethods from "./ngl.js"
import RepertoireBrowserPanel from "./repertoire-browser-panel.js"


import {_package} from "./package";

export class LaunchBrowser {

    async init(view) {


        ///// MAIN BODY ////
        let inputs = new RepertoireBrowserPanel();
        await inputs.init(view);

        let ngl = new NglMethods();
        await ngl.init(view, inputs);

        let pViz = new PvizMethods();
        await pViz.init(inputs, ngl);

        let msa = new MsaMethods();
        msa.init(view, inputs);

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

        // inputs.paratopes.onChanged(async () => {
        //
        //     ngl.stage.removeAllComponents();
        //     let schemeObj = ngl.CDR3(inputs.cdr_scheme, inputs.paratopes);
        //     await ngl.loadPdb(ngl.path, inputs.repChoice, schemeObj);
        //
        //     await pViz.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value)
        //     await pViz.loadSequence(inputs.pViz_host_L, 'L', inputs.paratopes.value)
        //
        //     MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        // });
        //
        // inputs.ptm_choices.onChanged(async () => {
        //
        //     pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value)
        //     await pViz.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value)
        //     await pViz.loadSequence(inputs.pViz_host_L, 'L', inputs.paratopes.value)
        //
        //     MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        // });
        //
        // inputs.ptm_prob.onChanged(async () => {
        //
        //     pViz.pVizParams.ptmMap = pViz.ptmMapping(inputs.ptm_choices.value, inputs.ptm_prob.value)
        //     await pViz.loadSequence(inputs.pViz_host_H, 'H', inputs.paratopes.value)
        //     await pViz.loadSequence(inputs.pViz_host_L, 'L', inputs.paratopes.value)
        //
        //     MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
        // });


        inputs.msaContentChoice.onChanged(() => { msa.drawAlignments(); });
        DG.debounce(view.table.onCurrentRowChanged, 200).subscribe(() => { msa.drawAlignments(); });

        // let seqHeavyCol = table.col('sequence_alignment_heavy');
        // let seqLightCol = table.col('sequence_alignment_light');
        // let germHeavyCol = table.col('germline_alignment_heavy');
        // let germLightCol = table.col('germline_alignment_light');
        //
        // let seqAlignHeavyAACol = table.col('sequence_alignment_aa_heavy');
        // let seqAlignLightAACol = table.col('sequence_alignment_aa_light');
        // let germAlignHeavyAACol = table.col('germline_alignment_aa_heavy');
        // let germAlignLightAACol = table.col('germline_alignment_aa_light');
        //
        // let vStartHeavy = table.col('v_alignment_start_heavy');
        // let dStartHeavy = table.col('d_alignment_start_heavy');
        // let jStartHeavy = table.col('j_alignment_start_heavy');
        // let vEndHeavy = table.col('v_alignment_end_heavy');
        // let dEndHeavy = table.col('d_alignment_end_heavy');
        // let jEndHeavy = table.col('j_alignment_end_heavy');
        //
        // let vStartLight = table.col('v_alignment_start_light');
        // let dStartLight = table.col('d_alignment_start_light');
        // let jStartLight = table.col('j_alignment_start_light');
        // let vEndLight = table.col('v_alignment_end_light');
        // let dEndLight = table.col('d_alignment_end_light');
        // let jEndLight = table.col('j_alignment_end_light');


        // const msaOpts = {
        //     el: msa_host_L,
        //     vis: {
        //         conserv: false,
        //         overviewbox: false,
        //         consensus: true,
        //         seqlogo: true,
        //         scaleslider: false,
        //     },
        //     conf: {
        //         dropImport: true
        //     },
        //     bootstrapMenu: false,
        //     zoomer: {
        //         boxRectHeight: 1,
        //         boxRectWidth: 1,
        //         labelNameLength: 110,
        //         labelFontsize: 10,
        //         labelIdLength: 20
        //     }
        // };
        // const msaL = new msa.msa(msaOpts);
        // msaOpts.el = msa_host_H;
        // const msaH = new msa.msa(msaOpts);
        // const gffParser = msa.io.gff;
        // DG.debounce(table.onCurrentRowChanged, 200).subscribe(drawAlignments);
        // msaMethods.drawAlignments();


        ngl.nglResize(inputs.ngl_host, ngl.stage);
        await pViz.pvizResize(inputs.pViz_host_H, 'H', inputs.paratopes);
        await pViz.pvizResize(inputs.pViz_host_L, 'L', inputs.paratopes);
        //
        // MiscMethods.setDockSize(view, inputs.ngl_node, inputs.sequence_tabs, inputs.paratopes);
    }

}