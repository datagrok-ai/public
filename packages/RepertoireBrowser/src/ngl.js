import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import json from "./example.json";
import {MiscMethods} from "./misc.js"
import {_package} from "./package";

// export let _package = new DG.Package();

export class NglMethods {

    async init(view, inputs, json) {

        let colorScheme = inputs.colorScheme;
        let col_background = colorScheme["col_background"];

        inputs.ngl_host.style.backgroundColor = col_background;
        view.box = true;
        this.stage = new NGL.Stage(inputs.ngl_host);
        this.path = _package.webRoot + 'pdbfiles/' + 'example.pdb';

        this.schemeObj = this.CDR3(inputs.cdr_scheme, inputs.paratopes, json, colorScheme);

        await this.loadPdb(this.path, inputs.repChoice, this.schemeObj);
        this.nglResize(inputs.ngl_host);
    }

    // ---- NGL ----
    // create a color scheme for CDR3 regions
    CDR3(cdr_scheme, paratopes, json, colorScheme) {

        let col_heavy_chain = colorScheme["col_heavy_chain"];
        let col_light_chain = colorScheme["col_light_chain"];
        let col_cdr = colorScheme["col_cdr"];
        let col_para = colorScheme["col_para"];
        let col_partopes_low = colorScheme["col_partopes_low"]; //col_para in rgb
        let col_partopes_high = colorScheme["col_partopes_high"];


        let schemeId;

        if (paratopes.value === true) {
            let palette = MiscMethods.interpolateColors(col_partopes_low, col_partopes_high, 100);
            let selectionScheme = [];
            Object.keys(json.parapred_predictions).forEach((chain) => {
                Object.keys(json.parapred_predictions[chain]).forEach((index) => {

                    selectionScheme.push([
                        palette[Math.round(json.parapred_predictions[chain][index] * 100)],
                        `${index} and :${chain}`
                    ]);
                })

            })
            selectionScheme.push([col_para, "* and :H"]);
            selectionScheme.push([col_para, "* and :L"]);
            schemeId = NGL.ColormakerRegistry.addSelectionScheme(selectionScheme);
        } else {
            if (cdr_scheme.value === 'default') {
                schemeId = NGL.ColormakerRegistry.addSelectionScheme([
                    [col_heavy_chain, "* and :H"],
                    [col_light_chain, "* and :L"]
                ]);
            } else {
                let scheme_buffer = [];
                Object.keys(json.cdr_ranges).forEach((str) => {
                    if (str.includes(cdr_scheme.value + '_CDRH')) {
                        let str_buffer = '';
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${json.cdr_ranges[str][i][0]}-${json.cdr_ranges[str][i][1]} and :H`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push([col_cdr, str_buffer]);
                        scheme_buffer.push([col_heavy_chain, "* and :H"]);

                    } else if (str.includes(cdr_scheme.value + '_CDRL')) {
                        let str_buffer = ''
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${json.cdr_ranges[str][i][0]}-${json.cdr_ranges[str][i][1]} and :L`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push([col_cdr, str_buffer]);
                        scheme_buffer.push([col_light_chain, "* and :L"]);
                    }
                });
                schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
            }
        }
        return {color: schemeId};
    }

    // load the 3D model
    async loadPdb(bytes, repChoice, schemeObj) {
        await this.stage.loadFile(bytes).then(function (o) {
            o.addRepresentation(repChoice.value, schemeObj);
            o.autoView();
        });
    }

    // viewer resize
    _resize(host) {
        let canvas = host.querySelector('canvas');
        canvas.width = Math.floor(host.clientWidth * window.devicePixelRatio);
        canvas.height = Math.floor(host.clientHeight * window.devicePixelRatio);
        this.stage.handleResize();
    }

    nglResize(host) {
        ui.onSizeChanged(host).subscribe((_) => this._resize(host));
        this._resize(host);
    }

}