/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import map from "./maps.json";

export let _package = new DG.Package();

//name: toScript
//top-menu: convert | toScript
export async function toScript() {

    // recursive multi-level object merging
    function assignOnlyIntersection(target, source) {
        Object.keys(target).forEach(key => {
            (target[key] == null) && delete target[key]
            if (target[key] in source) {
                target[key] = source[target[key]];
            } else if (typeof target[key] === 'object') {
                assignOnlyIntersection(target[key], source)
            }
        })
    }

    // dynamic default string trimming
    function dynamicReplace(stRing, optionsObj, map) {
        Object.keys(optionsObj).forEach(key => {
            (optionsObj[key] == null) && delete optionsObj[key]
            if (typeof optionsObj[key] === 'object' &&
                Object.keys(optionsObj[key]).every(elem => elem != '0')) {
                stRing = dynamicReplace(stRing, optionsObj[key], map);
            } else {
                stRing = stRing.replace("!(" + key + ")",map[key]);
                stRing = stRing.split("!(" + key + ")").join(optionsObj[key]);
            }
        })
        return stRing;
    }

    // main slicing function
    async function strReplace(optionsObj) {

        // extract default string and viewer type
        let stRing;
        let paramsMap;
        if (optionsObj.type === 'Trellis plot') {
            stRing = map.plotScripts[optionsObj.look.innerViewerLook.type] +
                " + facet_grid(!(xColumnNames)!(yColumnNames))"; // need to append 'Look'
            paramsMap = map.additionalOps[optionsObj.look.innerViewerLook.type];
        } else {
            stRing = map.plotScripts[optionsObj.type];
            paramsMap = map.additionalOps[optionsObj.type];
        }

        // Replace misc grok codes with R analogues
        assignOnlyIntersection(optionsObj, map.miscCodes);

        // trim the default code string
        stRing = dynamicReplace(stRing, optionsObj, paramsMap);
        stRing = stRing.replace(/!\([^)]*\) */g, "");

        // add a print statement
        stRing = stRing + "\nprint(plt)"
        return stRing;
    }

    // test viewer + options
    let view = grok.shell.addTableView(grok.data.demo.demog());

    // Scatter plot
    // let plot = view.scatterPlot({
    //     x: 'height',
    //     y: 'weight',
    //     color: 'race'
    // });
    // plot.setOptions({
    //     showRegressionLine: true,
    //     markerType: 'square'
    // });

    // // Histogram
    // let plot = view.histogram({
    //     value: 'weight'
    // });
    // plot.setOptions({
    //     bins : 40
    // });

    // // Bar chart
    // let plot = view.barChart({
    //     split: 'study',
    //     value: 'age'
    // });
    // plot.setOptions({
    //     valueAggrType : "skew"
    // });

    // // Box plot
    // let plot = view.boxPlot({
    //     value : 'height',
    //     category: 'site',
    //     markerColorColumnName: 'sex'
    //     // binColorColumnName: 'height'
    // })

    // Correlation plot
    let plot = view.corrPlot({
        xs: ['age', 'weight', 'height'],
        ys: ['age', 'weight', 'height'],
    });

    // let plot = view.lineChart();


    // collect viewer properties
    let options = JSON.parse(plot.getOptions());

    //choose and slice the string
    let rCode = await strReplace(options);

    // output code in a dialogue window
    // let input = document.createElement('TEXTAREA');
    // input.value = rCode
    // ui.dialog('OUTPUT SCRIPT')
    //     .add(input)
    //     .onOK(() => { grok.shell.info('OK!'); })
    //     .showModal(true);

    ui.dialog('OUTPUT SCRIPT')
        .add(view.root)
        .onOK(() => { grok.shell.info('OK!'); })
        .showModal(true);

    // run the generated script in R
    view.addViewer('Scripting Viewer', {
        script: map.header + rCode
    });
}
