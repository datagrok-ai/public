/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import map from "./maps.json";

export let _package = new DG.Package();

//name: toScript
//top-menu: convert | toScript
export async function toScript() {
    
    //main slicing function
    async function strReplace(optionsObj) {

        grok.shell.info(Object.entries(optionsObj.look));

        let stRing = map.plotScripts[optionsObj.type];
        let paramsMap = map.additionalOps[optionsObj.type];

        //decompose getOptions() output
        let toRemove = Object.keys(optionsObj.look);
        let toInsert = Object.values(optionsObj.look);

        //decompose custom mappings
        let mapKeys = Object.keys(paramsMap);
        let mapVals = Object.values(paramsMap);

        //trim the generalized code version
        let i;
        for (i = 0; i < mapKeys.length; i++) {

            if (toRemove.includes(mapKeys[i])) {
                stRing = stRing.replace("." + mapKeys[i],mapVals[i]);
            } else {
                stRing = stRing.replace("." + mapKeys[i],"");
            }

        }

        //fill in the actual parameters
        for (i = 0; i < toRemove.length; i++) {

            //replace grok specific toInsert values with R codes
            if (Object.keys(map.miscCodes).includes(toInsert[i])) {
                toInsert[i] = map.miscCodes[toInsert[i]];
            }

            //replace missing values with NULL
            if (toInsert[i] == null) {
                toInsert[i] = 'NULL';
            }

            //replace all string parameter markers with corresponding values
            // stRing = stRing.replace(/toRemove[i]/g, toInsert[i]);
            stRing = stRing.split("." + toRemove[i]).join(toInsert[i]);
        }

        stRing = stRing + "\nprint(plt)"
        return stRing;
    }

    //test viewer + options
    let view = grok.shell.addTableView(grok.data.demo.demog());

    // //Scatter plot
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

    // Bar chart
    // let plot = view.barChart({
    //     split: 'study',
    //     value: 'age'
    // });
    // plot.setOptions({
    //     valueAggrType : "skew"
    // });

    // //Box plot
    // let plot = view.boxPlot({
    //     value : 'height',
    //     category: 'site',
    //     markerColorColumnName: 'sex'
    //     // binColorColumnName: 'height'
    // })

    //Correlation plot
    let plot = view.corrPlot({
        xs: ['age', 'weight', 'height'],
        ys: ['age', 'weight', 'height'],
    });


    //collect viewer properties
    let options = JSON.parse(plot.getOptions());

    //choose and slice the string
    let rCode = await strReplace(options);

    //output code in a dialogue window
    let input = document.createElement('TEXTAREA');
    input.value = rCode
    ui.dialog('Output script')
        .add(input)
        .onOK(() => { grok.shell.info('OK!'); })
        .show();

    //run the generated script in R
    view.addViewer('Scripting Viewer', {
        script: map.header + rCode
    });
}
