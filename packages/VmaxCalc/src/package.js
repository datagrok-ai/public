/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: vmax
export async function Vmax() {

    async function calculate(dataTable, farea, pressure, vbatch, tbatch){
        let f = await grok.functions.eval("Vmaxcalc:VmaxR");

        let call = f.prepare({
            'test_data': dataTable,
            'test_area': farea.value,
            'test_pressure': pressure.value,
            'vbatch': vbatch.value,
            'tbatch': tbatch.value
        });

        await call.call();
        let render_data = call.getParamValue('proc_data');
        let reg_params = call.getParamValue('reg_params');
        return [render_data, reg_params];
    }

    function regParamsToTable(t) {
        let string = `${t.get('vars',0)}: ${Math.round(t.get('vals',0)*Math.pow(10,10))/Math.pow(10,10)}`
        for (let i = 1; i <= t.rowCount - 1; i++) {
            string = string + `\n${t.get('vars',i)}: ${Math.round(t.get('vals',i)*Math.pow(10,10))/Math.pow(10,10)}`;
        }
        return string;
    }

    let v = ui.dialog('Model information');

    let tableName = ui.choiceInput('Table', null, grok.shell.tableNames);
    let dataTable = grok.shell.tableByName(tableName.value);

    let farea = ui.floatInput('Filter area (cm^2)',3.5);
    let pressure = ui.floatInput('Pressure (psi)',10);
    let vbatch = ui.floatInput('Vbatch (L)',25);
    vbatch.setTooltip('Desired batch size to process');
    let tbatch = ui.floatInput('tbatch (hr)',0.5);
    tbatch.setTooltip('Desitred process time');

    let inputs = ui.div([ui.inputs([tableName, farea, pressure, vbatch, tbatch])]);
    v.add(inputs).onOK(async ()=> {
        let results = await calculate(dataTable, farea, pressure, vbatch, tbatch);
        let results_view = grok.shell.addTableView(results[0]);
        results_view.scatterPlot({
            x: 'time (hr)',
            y: 't/V (hr/(L/m2))',
            markerDefaultSize: 8,
            showRegressionLine: true
        });

        let acc = results_view.toolboxPage.accordion;
        acc.addPane('Regression parameters',
            () => ui.div([DG.Viewer.grid(results[1]).root]),
            // () => ui.divText([regParamsToTable(results[1])]),
            true, acc.panes[0])
    }).show();

}
