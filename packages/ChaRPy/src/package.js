/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import mapR from "./mapR.json";
import mapPy from "./mapPy.json";

export let _package = new DG.Package();

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
function dynamicReplace(colsList, stRing, optionsObj, map) {
    Object.keys(optionsObj).forEach(key => {
        (optionsObj[key] == null) && delete optionsObj[key];
        if (key.includes('ColumnName')) {
            colsList.push("'" + optionsObj[key] + "'");
        }
        if (typeof optionsObj[key] === 'object' &&
            Object.keys(optionsObj[key]).every(elem => elem != '0')) {
            stRing = dynamicReplace(colsList, stRing, optionsObj[key], map)[0];
        } else {
            stRing = stRing.replace("!(" + key + ")",map[key]);
            stRing = stRing.split("!(" + key + ")").join(optionsObj[key]);
        }
    })
    return [stRing, colsList];
}

grok.events.onContextMenu.subscribe((args) => {
    if (args.args.context instanceof DG.Viewer) {
        args.args.menu.item('to R script', async () => {

            async function strReplace(optionsObj,map) {

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
                let colsList = [];
                stRing = dynamicReplace(colsList, stRing, optionsObj, paramsMap)[0];
                stRing = stRing.replace(/!\([^)]*\) */g, "");

                // add a print statement
                stRing = stRing + "\nprint(plt)"
                return stRing;
            }

            let options = JSON.parse(args.args.context.getOptions());
            let view = DG.Viewer.fromType(options.type, args.args.context.table, options.look);
            let rCode = await strReplace(options,mapR);

            let block =
                $(ui.splitV([
                    ui.textArea(rCode),
                    ui.splitH([
                        view,
                        DG.Viewer.fromType('Scripting Viewer',args.args.context.table,{script: mapR.header + rCode})])
                ])).css('flex-grow', '1');

            ui.dialog('OUTPUT SCRIPT').add(block[0]).showModal(true);
        });
    }
});

grok.events.onContextMenu.subscribe((args) => {
    if (args.args.context instanceof DG.Viewer) {
        args.args.menu.item('to Python script', async () => {

            async function strReplace(optionsObj,map) {

                // extract default string and viewer type
                let pyString = map.plotScripts[optionsObj.type];
                let paramsMap = map.additionalOps[optionsObj.type];

                // Replace misc grok codes with R analogues
                assignOnlyIntersection(optionsObj, map.miscCodes);

                // trim the default code string
                let colsList = [];
                let dynamicOut = dynamicReplace(colsList,pyString, optionsObj, paramsMap);
                pyString = dynamicOut[0];
                colsList = dynamicOut[1];
                pyString = pyString.replace("!(colsList)", colsList);
                pyString = pyString.replace(/!\([^)]*\) */g, "");

                // add a print statement
                return pyString;
            }

            let options = JSON.parse(args.args.context.getOptions());
            let view = DG.Viewer.fromType(options.type, args.args.context.table, options.look);
            let pyCode = await strReplace(options,mapPy);

            let block =
                $(ui.splitV([
                    ui.textArea(pyCode),
                    ui.splitH([
                        view,
                        DG.Viewer.fromType('Scripting Viewer',args.args.context.table,{script: mapPy.header + pyCode})])
                ])).css('flex-grow', '1');

            ui.dialog('OUTPUT SCRIPT').add(block[0]).showModal(true);
        });
    }
});

//name: exportFunc
//tags: autostart
export function export_123() {}
