import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
// import * as DG from "datagrok-api/dg";
import json from "./TPP000153303.json";

export default class MiscMethods {

    // processes JSON to derive scheme names
    static extract_schemes() {
        let raw_scheme_names = Object.keys(json.cdr_ranges);
        let schemes_lst = ['default'];
        raw_scheme_names.forEach((str) => {
            str = str.split('_')
            if (schemes_lst.includes(str[0]) === false) {
                schemes_lst.push(str[0]);
            }
        })
        return schemes_lst;
    }

    // color interpolation
    static interpolateColors(color1, color2, steps) {

        function interpolateColor(color1, color2, factor) {
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

        color1 = color1.match(/\d+/g).map(Number);
        color2 = color2.match(/\d+/g).map(Number);

        for(var i = 0; i < steps; i++) {
            interpolatedColorArray.push(interpolateColor(color1, color2, stepFactor * i));
        }

        return interpolatedColorArray;
    }

    // ---- Resizing ----
    static setDockSize(view, node, nodeContent) {
        let nodeContentHeight = 0;
        let rootNodeHeight = view.dockManager.rootNode.container.containerElement.clientHeight;
        let newHeight = 0;
        newHeight = $("#feature-viewer").outerHeight(true) + 55;
        newHeight = 1/(rootNodeHeight/newHeight);
        // newHeight = Math.ceil(newHeight*100)/100;

        return view.dockManager.dock(nodeContent, 'down', node, 'Sequence', newHeight.toFixed(2));
    }

    // ---- Save/Load ----
    // selection saving
    // async static save_load(table, acc) {
    //
    //     async function saveSelectedRows(table, uniqueId, connection, fileToSave) {
    //
    //         let indexes = table
    //             .groupBy([`${uniqueId.value}`])
    //             .whereRowMask(table.selection)
    //             .aggregate();
    //
    //         let data = `${table.toString()}\n${uniqueId.stringValue}\n${indexes.col(0).toList()}`;
    //         await grok.dapi.files.writeAsText(`${connection}${fileToSave.value}.txt`, data);
    //     }
    //
    //     async function loadSelectedRows(table, connection, savedFilesList) {
    //
    //         let res = await grok.dapi.files.readAsText(`${connection}${savedFilesList.value}`);
    //         res = res.split("\n");
    //         let uniqueColumnName = res[1];
    //         let values = res[2].split(',');
    //         console.log(values);
    //         table.rows.select((row) => values.includes(row[`${uniqueColumnName}`]));
    //
    //     }
    //
    //
    //     let fileToSave = ui.stringInput('FileName', 'filename');
    //     let connection = 'Demo:TestJobs:Files:DemoFiles/';
    //
    //     let files = await grok.dapi.files.list(connection, false, '');
    //     files = files.map((e) => e.path);
    //     let savedFilesList = await ui.choiceInput('Saved Rows', ' ', files)
    //     // let uniqueId = ui.columnInput('Unique id column', table, table.col('tenx_barcode'));
    //     let uniqueId = ui.stringInput('Unique id column', 'tenx_barcode');
    //
    //
    //     let saveDialog = () => {
    //         ui.dialog('Save rows to file')
    //             .add(uniqueId).add(fileToSave)
    //             .onOK(() => saveSelectedRows(table, uniqueId, connection, fileToSave)).show();
    //     };
    //
    //     let loadDialog = async () => {
    //         let files = await grok.dapi.files.list(connection, false, '');
    //         files = files.map((e) => e.path);
    //         let savedFilesList = await ui.choiceInput('Saved Rows', ' ', files)
    //
    //         ui.dialog('Load rows from file')
    //             .add(savedFilesList)
    //             .onOK(() => loadSelectedRows(table, connection, savedFilesList)).show();
    //     };
    //
    //
    //     let saveRowsButton = ui.button('SAVE');
    //     saveRowsButton.addEventListener("click", saveDialog);
    //
    //     let loadRowsButton = ui.button('LOAD')
    //     loadRowsButton.addEventListener("click", loadDialog);
    //
    //     acc.addPane('Save/Load', () => ui.divH([saveRowsButton, loadRowsButton]));
    // }

}