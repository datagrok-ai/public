import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import json from "./TPP000153303.json";
import {MiscMethods} from "./misc.js"
import {_package} from "./package";

// export let _package = new DG.Package();

export class NglMethods {

    async init(view, inputs) {

        let col_background = 'white';

        inputs.ngl_host.style.backgroundColor = col_background;
        view.box = true;
        this.stage = new NGL.Stage(inputs.ngl_host);
        this.path = _package.webRoot + 'pdbfiles/' + 'TPP000153303.pdb';
        this.schemeObj = this.CDR3(inputs.cdr_scheme, inputs.paratopes);

        await this.loadPdb(this.path, inputs.repChoice, this.schemeObj);
        this.nglResize(inputs.ngl_host);

    }

    // ---- NGL ----
    // create a color scheme for CDR3 regions
    CDR3(cdr_scheme, paratopes) {
  
        let col_heavy_chain = '#0069a7';
        let col_light_chain = '#f1532b';
        let col_cdr = '#45d145';
        let col_para = '#b0c4de';
        let col_partopes_low = '(176,196,222)'; //col_para in rgb
        let col_partopes_high = '(255, 0, 255)';


        let schemeId;
        let baseH = col_heavy_chain;
        let baseL = col_light_chain;

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
                    [baseH, "* and :H"],
                    [baseL, "* and :L"]
                ]);
            } else {
                let scheme_buffer = [];
                Object.keys(json.cdr_ranges).forEach((str) => {
                    if (str.includes(cdr_scheme.value + '_CDRH')) {
                        let str_buffer = ''
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${json.cdr_ranges[str][i][0]}-${json.cdr_ranges[str][i][1]} and :H`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push([col_cdr, str_buffer]);
                        scheme_buffer.push([baseH, "* and :H"]);

                    } else if (str.includes(cdr_scheme.value + '_CDRL')) {
                        let str_buffer = ''
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${json.cdr_ranges[str][i][0]}-${json.cdr_ranges[str][i][1]} and :L`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push([col_cdr, str_buffer]);
                        scheme_buffer.push([baseL, "* and :L"]);
                    }
                });
                schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
            }
        }
        return {color: schemeId};
    }

    // load the 3D model
    async loadPdb(bytes, repChoice, schemeObj) {
        this.stage.loadFile(bytes).then(function (o) {
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