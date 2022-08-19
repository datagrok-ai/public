import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export class MiscMethods {

    // processes JSON to derive scheme names
    static extract_schemes(json: any) {
        let raw_scheme_names = Object.keys(json.cdr_ranges);
        let schemes_lst = ['default'];
        raw_scheme_names.forEach((str) => {
            let strArr = str.split('_')
            if (schemes_lst.includes(strArr[0]) === false) {
                schemes_lst.push(strArr[0]);
            }
        })
        return schemes_lst;
    }

    // color interpolation
    static interpolateColors(color1: string, color2: string, steps: number) {

        function interpolateColor(color1: number[], color2: number[], factor: number) {
            if (arguments.length < 3) {
                factor = 0.5;
            }
            var result = color1.slice();
            for (var i = 0; i < 3; i++) {
                result[i] = Math.round(result[i] + factor * (color2[i] - color1[i]));
            }
            let hex_col = "#" + ((1 << 24) + (result[0] << 16) + (result[1] << 8) + result[2]).toString(16).slice(1);
            return hex_col;
        };

        var stepFactor = 1 / (steps - 1),
            interpolatedColorArray = [];

        let colors1 = color1.match(/\d+/g)!.map(Number);
        let colors2 = color2.match(/\d+/g)!.map(Number);

        for (var i = 0; i < steps; i++) {
            interpolatedColorArray.push(interpolateColor(colors1, colors2, stepFactor * i));
        }

        return interpolatedColorArray;
    }

    // ---- Resizing ----
    static setDockSize(view: DG.TableView, node:any, nodeContent: any) {

        let rootNodeHeight = view.dockManager.rootNode.container.containerElement.clientHeight;
        let newHeight = 0;
        //@ts-ignore
        newHeight = $("#feature-viewer").outerHeight(true) + 55;
        newHeight = 1 / (rootNodeHeight / newHeight);

        //@ts-ignore
        return view.dockManager.dock(nodeContent, 'down', node, 'Sequence', newHeight.toFixed(2));
    }
}
