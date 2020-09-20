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
            // stRing = stRing.replace("!(" + key + ")",map[key]);
            stRing = stRing.split("!(" + key + ")").join(map[key]);
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

                    stRing = map.plotScripts[optionsObj.look.innerViewerLook['#type'].toLowerCase().
                        replace('look','')] + " + facet_grid(!(xColumnNames)!(yColumnNames))";
                    paramsMap = map.additionalOps[optionsObj.look.innerViewerLook['#type'].toLowerCase().
                        replace('look','')];

                } else {
                    stRing = map.plotScripts[optionsObj.type.toLowerCase().replace(' ','')];
                    paramsMap = map.additionalOps[optionsObj.type.toLowerCase().replace(' ','')];
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
            let viewLeft = DG.Viewer.fromType(options.type,
              args.args.context.table, options.look);
            let rCode = await strReplace(options, mapR);
            let viewRight = DG.Viewer.fromType('Scripting Viewer',
              args.args.context.table, {script: mapR.header + rCode});

            let block =
                $(ui.splitV([
                    ui.textArea(rCode),
                    ui.splitH([
                        viewLeft,
                        viewRight])
                ])).css('flex-grow', '1');

            ui.dialog('OUTPUT SCRIPT').onClose(() => {
              viewLeft.dataFrame = new DG.DataFrame();
              viewRight.dataFrame = new DG.DataFrame();
            }).add(block[0]).showModal(true);
        });
    }
});

grok.events.onContextMenu.subscribe((args) => {
    if (args.args.context instanceof DG.Viewer) {
        args.args.menu.item('to Python script', async () => {

            async function strReplace(optionsObj,map) {

                let pyString;
                let paramsMap;
                if (optionsObj.type === 'Trellis plot') {

                    pyString = map.plotScripts[optionsObj.look.innerViewerLook['#type'].toLowerCase().
                    replace('look','')] ;
                    paramsMap = map.additionalOps[optionsObj.look.innerViewerLook['#type'].toLowerCase().
                    replace('look','')];

                } else {
                    pyString = map.plotScripts[optionsObj.type.toLowerCase().replace(' ','')];
                    paramsMap = map.additionalOps[optionsObj.type.toLowerCase().replace(' ','')];
                }

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
            let viewLeft = DG.Viewer.fromType(options.type,
              args.args.context.table, options.look);
            let pyCode = await strReplace(options, mapPy);
            let viewRight = DG.Viewer.fromType('Scripting Viewer',
              args.args.context.table, {script: mapPy.header + pyCode});

            let block =
                $(ui.splitV([
                    ui.textArea(pyCode),
                    ui.splitH([
                        viewLeft,
                        viewRight])
                ])).css('flex-grow', '1');

            ui.dialog('OUTPUT SCRIPT').onClose(() => {
              viewLeft.dataFrame = new DG.DataFrame();
              viewRight.dataFrame = new DG.DataFrame();
            }).add(block[0]).showModal(true);
        });
    }
});

//name: exportFunc
//tags: autostart
export function toScriptInit() {}
