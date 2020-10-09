/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

async function importPy(data,samplingFreq,signalType){

    let f = await grok.functions.eval("Biosensors:importPyphysio");

    let call = f.prepare({
        'ecg_data':data,
        'fsamp':samplingFreq.value,
        'signalType':signalType.value
    });

    await call.call();
    return call.getParamValue('plt');
}

async function filterPy(data,samplingFreq,signalType,filter1,filterparams){

    let f;
    let call;
    if (filter1.value === 'IIR') {

        f = await grok.functions.eval("Biosensors:filterIirPyphysio");

        call = f.prepare({
            'ecg_data':data,
            'fsamp':samplingFreq.value,
            'signalType':signalType.value,
            'fp': filterparams[0].value,
            'fs': filterparams[1].value,
            'ftype': filterparams[2].value,
        });

    }
    if (filter1.value === 'normalize') {

        f = await grok.functions.eval("Biosensors:filterNormPyphysio");

        call = f.prepare({
            'ecg_data':data,
            'fsamp':samplingFreq.value,
            'signalType':signalType.value,
            'norm_method': filterparams[0].value,
        });

    }
    if (filter1.value === 'resample') {

        f = await grok.functions.eval("Biosensors:filterResPyphysio");

        call = f.prepare({
            'ecg_data':data,
            'fsamp':samplingFreq.value,
            'signalType':signalType.value,
            'fout': filterparams[0].value,
            'kind': filterparams[1].value
        });

    }
    await call.call();
    return call.getParamValue('plt');
}


//name: pipelineDemo
export async function pipelineDemo() {

    function paramSelector(x) {

        let methodparams;
        if (x === 'IIR') {
            let fp = ui.intInput('fp', 45);
            let fs = ui.intInput('fs', 50);
            let ftype = ui.choiceInput('ftype','',['ellip']);
            methodparams = [fp,fs,ftype];

        }
        if (x === 'normalize') {

            let normMethod = ui.choiceInput('norm_method','',['standard']);
            methodparams = [normMethod];
        }
        if (x === 'resample') {

            let fout = ui.intInput('fout', 4096);
            let kind = ui.choiceInput('kind','',['cubic']);
            methodparams = [fout,kind];
        }
        return methodparams;
    }

    let v = ui.dialog('DEMO PIPELINE');

    /////////////////////////////////////////

    // Import
    let tableName = ui.choiceInput('Table', null, grok.shell.tableNames);
    let dataTable = grok.shell.tableByName(tableName.value);

    let signalClass = ui.choiceInput('Signal class','',['Evenly signal','Unevenly signal']);
    signalClass.setTooltip('Dependent on the nature of oscilations Periodic/Aperiodic');

    let samplingFreq = ui.intInput('Sampling frequency', 1024);
    samplingFreq.setTooltip('Number of samples per second');

    let signalType = ui.choiceInput('Signal type','',['ecg','eda']);
    signalType.setTooltip('Nature of the physiological signal');

    // Filter
    let filter1 = ui.choiceInput('1st filter', '', ['IIR','normalize','resample']);
    // let filter2 = ui.choiceInput('2nd filter', '', ['IIR','normalize','resample']);
    // let filter3 = ui.choiceInput('3rd filter', '', ['IIR','normalize','resample']);

    // Information extraction
    let infoType = ui.choiceInput('To extract', '', ['Beat from ECG']);

    // Indicators
    let indicator = ui.choiceInput('Indicator preset', '', ['HRV']);

    /////////////////////////////////////////

    // Import containers
    let containerIMPORT = ui.div();
    let containerOGplot = ui.div();

    // Filter containers
    let containerFILTER1 = ui.div();
    let paramsContainer1 = ui.div();
    let accFILTER1 = ui.accordion();

    // let containerFILTER2 = ui.div();
    // let paramsContainer2 = ui.div();
    // let accFILTER2 = ui.accordion();
    //
    // let containerFILTER3 = ui.div();
    // let paramsContainer3 = ui.div();
    // let accFILTER3 = ui.accordion();
    let containerFLplot = ui.div();

    // Information extraction containers
    let containerINFO = ui.div();
    let containerINFplot = ui.div();

    // Indicator containers
    let containerIndicator = ui.div();

    /////////////////////////////////////////

    // Import dialogue
    let importInputs = ui.inputs([tableName, signalClass, samplingFreq, signalType]);
    containerIMPORT.appendChild(importInputs);
    containerOGplot.appendChild(ui.bigButton('PLOT ORIGINAL',async () => {

        let plotOG = await importPy(dataTable,samplingFreq,signalType);
        let tableView = grok.shell.getTableView(tableName.value);
        let node1 = tableView.dockManager.dock(plotOG, 'right', null, 'Original plot');

    }));

    // Filter dialogue
    let filterparams;
    let filterInputs1 = ui.inputs([filter1]);
    containerFILTER1.appendChild(filterInputs1);
    filter1.onChanged(function () {
        $(paramsContainer1).empty()
        filterparams = paramSelector(filter1.value)
        paramsContainer1.appendChild(ui.inputs(filterparams));
    });
    accFILTER1.addPane('parameters', () => paramsContainer1);

    // let filterInputs2 = ui.inputs([filter2]);
    // containerFILTER2.appendChild(filterInputs2);
    // filter2.onChanged(function () {
    //     $(paramsContainer2).empty()
    //     let filterparams2 = paramSelector(filter2.value)
    //     paramsContainer2.appendChild(ui.inputs(filterparams2));
    // });
    // accFILTER2.addPane('parameters', () => paramsContainer2);
    //
    // let filterInputs3 = ui.inputs([filter3]);
    // containerFILTER3.appendChild(filterInputs3);
    // filter3.onChanged(function () {
    //     $(paramsContainer3).empty()
    //     let filterparams3 = paramSelector(filter3.value)
    //     paramsContainer3.appendChild(ui.inputs(filterparams3));
    // });
    // accFILTER3.addPane('parameters', () => paramsContainer3);

    containerFLplot.appendChild(ui.bigButton('PLOT FILTERED',async () => {

        let plotFL = await filterPy(dataTable,samplingFreq,signalType,filter1,filterparams);
        let tableView = grok.shell.getTableView(tableName.value);
        let node2 = tableView.dockManager.dock(plotFL, 'right', null, 'Filtered plot');

    }));

    // Information extraction dialogue
    let infoInputs = ui.inputs([infoType]);
    containerINFO.appendChild(infoInputs);
    containerINFplot.appendChild(ui.bigButton('EXTRACT INFO',async () => {}));

    // Indicators dialogue
    let indicatorInputs = ui.inputs([indicator]);
    containerIndicator.appendChild(indicatorInputs);


    /////////////////////////////////////////

    v.add(containerIMPORT);
    v.add(containerOGplot);
    v.add(containerFILTER1);
    v.add(accFILTER1);
    // v.add(containerFILTER2);
    // v.add(accFILTER2);
    // v.add(containerFILTER3);
    // v.add(containerFLplot);
    // v.add(accFILTER3);
    v.add(containerFLplot);
    v.add(containerINFO);
    v.add(containerINFplot);
    v.add(containerIndicator).onOK(()=> {}).show();

}
