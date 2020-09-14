/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import map from "./maps.json";

export let _package = new DG.Package();

grok.events.onContextMenu.subscribe((args) => {
    if (args.args.context instanceof DG.Viewer) {
        args.args.menu.item('to R script', async () => {

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
            function pltGenerate(options) {

                let plt;
                if (options.type === 'Scatter plot') {
                    plt = DG.Viewer.scatterPlot(args.args.context.table, options.look);
                } else if (options.type === 'Histogram') {
                    plt = DG.Viewer.histogram(args.args.context.table, options.look);
                } else if (options.type === 'Bar chart') {
                    plt = DG.Viewer.barChart(args.args.context.table, options.look);
                } else if (options.type === 'Box plot') {
                    plt = DG.Viewer.boxPlot(args.args.context.table, options.look);
                } else if (options.type === 'Correlation plot') {
                    plt = DG.Viewer.correlationPlot(args.args.context.table, options.look);
                } else if (options.type === 'Line chart') {
                    plt = DG.Viewer.lineChart(args.args.context.table, options.look);
                }

                return plt;
            }
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

            let options = JSON.parse(args.args.context.getOptions());
            let plt = pltGenerate(options);
            let rCode = await strReplace(options);

            let block =
                $(ui.splitV([
                    ui.textArea(rCode),
                    ui.splitH([
                        plt,
                        DG.Viewer.fromType('Scripting Viewer',args.args.context.table,{script: map.header + rCode})])
                ])).css('flex-grow', '1');

            ui.dialog('OUTPUT SCRIPT').add(block[0]).showModal(true);
        });
    }
});

//name: exportFunc
//tags: autostart
export function export_123() {}
