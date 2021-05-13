import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import json from "./TPP000153303.json";
import MiscMethods from "./misc.js"
import {_package} from "./package";

// export let _package = new DG.Package();

export default class NglMethods {

    async init(view, inputs) {

        // inputs.ngl_host = ui.div([],'d4-ngl-viewer');
        inputs.ngl_host.style.backgroundColor ='black';
        view.box = true;
        this.stage = new NGL.Stage(inputs.ngl_host);
        this.path = _package.webRoot + 'pdbfiles/' + 'TPP000153303.pdb';
        this.schemeObj = this.CDR3(inputs.cdr_scheme, inputs.paratopes);

        await this.loadPdb(this.path, inputs.repChoice, this.schemeObj);
    }

    // ---- NGL ----
    // create a color scheme for CDR3 regions
    CDR3(cdr_scheme, paratopes) {
        let schemeId;
        let baseH = 'darkblue';
        let baseL = 'darkred';

        if (paratopes.value === true) {
            let palette = MiscMethods.interpolateColors('(255, 255, 255)', '(255, 0, 255)', 100);
            let selectionScheme = [];
            Object.keys(json.parapred_predictions).forEach((chain) => {
                Object.keys(json.parapred_predictions[chain]).forEach((index) => {
                    selectionScheme.push([
                        palette[Math.round(json.parapred_predictions[chain][index] * 100)],
                        `${index} and :${chain}`
                    ]);
                })
            })
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
                        scheme_buffer.push(["limegreen", str_buffer]);
                        scheme_buffer.push([baseH, "* and :H"]);

                    } else if (str.includes(cdr_scheme.value + '_CDRL')) {
                        let str_buffer = ''
                        for (let i = 0; i < Object.keys(json.cdr_ranges[str]).length; i++) {
                            str_buffer = str_buffer + ` or ${json.cdr_ranges[str][i][0]}-${json.cdr_ranges[str][i][1]} and :L`;
                        }
                        str_buffer = str_buffer.slice(4);
                        scheme_buffer.push(["limegreen", str_buffer]);
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

    nglResize(host, stage) {
        let canvas = host.querySelector('canvas');

        function resize() {
            canvas.width = Math.floor(host.clientWidth * window.devicePixelRatio);
            canvas.height = Math.floor(host.clientHeight * window.devicePixelRatio);
            stage.handleResize();
        }

        ui.onSizeChanged(host).subscribe((_) => resize());
        resize();
    }

}