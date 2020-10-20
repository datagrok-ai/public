/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

async function importPy(data,samplingFreq,signalType){

    let f = await grok.functions.eval("Pyphysio:importPyphysio");

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

        f = await grok.functions.eval("Pyphysio:filterIirPyphysio");

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

        f = await grok.functions.eval("Pyphysio:filterNormPyphysio");

        call = f.prepare({
            'ecg_data':data,
            'fsamp':samplingFreq.value,
            'signalType':signalType.value,
            'norm_method': filterparams[0].value,
        });

    }
    if (filter1.value === 'resample') {

        f = await grok.functions.eval("Pyphysio:filterResPyphysio");

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

async function  applyFilter(data,samplingFreq,signalType, paramsT){

    let f = await grok.functions.eval("Pyphysio:filterPyphysio");

    let call = f.prepare({
        'ecg_data':data,
        'fsamp':samplingFreq.value,
        'signalType':signalType.value,
        'paramsT':paramsT
    });

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

    // Information extraction
    let infoType = ui.choiceInput('To extract', '', ['Beat from ECG']);

    // Indicators
    let indicator = ui.choiceInput('Indicator preset', '', ['HRV']);

    /////////////////////////////////////////

    // Import containers
    let containerIMPORT = ui.div();
    let containerOGplot = ui.div();

    // Filter containers
    let filterButton = ui.div();
    let containerFILTER = ui.div();
    let paramsContainer = ui.div();
    let accFILTER = ui.accordion();
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
    let filtersLST = [];
    let allParams = [];
    let filterInputs = ui.inputs(filtersLST);
    containerFILTER.appendChild(filterInputs);

    let i = 0;
    filterButton.appendChild(ui.button('ADD FILTER',async () => {

        var str ="filter" + i + " = ui.choiceInput('filter №' + (i + 1), '', ['IIR','normalize','resample'])";
        eval(str);

        filtersLST[i] = eval('filter' + i);
        let filterInputs1 = ui.inputs(filtersLST);

        containerFILTER.replaceChild(filterInputs1,filterInputs);
        filterInputs = filterInputs1;

        eval('filter' + i).onChanged(function () {
            $(paramsContainer).empty();
            let val = eval('filter' + (i-1)).value;
            allParams[i-1] = paramSelector(val);
            paramsContainer.appendChild(ui.inputs(allParams[i-1]));
        })
        i++;
    }));
    accFILTER.addPane('parameters', () => paramsContainer)

    containerFLplot.appendChild(ui.bigButton('PLOT FILTERED',async () => {

        let paramsT = DG.DataFrame.create(filtersLST.length);
        paramsT.columns.addNew('filter', 'string');

        for(let j=0; j<filtersLST.length; j++) {
            paramsT.columns.byName('filter').set(j,filtersLST[j].value);

            Object.keys(allParams[j]).forEach(key => {
                if(!paramsT.columns.names().includes(key)) {
                    // definitely needs reworking
                    if(typeof(allParams[j][key].value) === 'number'){
                        paramsT.columns.addNew(key,'int');
                    } else {
                        paramsT.columns.addNew(key,typeof(allParams[j][key].value));
                    }
                }
                paramsT.columns.byName(key).set(j,allParams[j][key].value);
            })
        }

        let plotFL = await applyFilter(dataTable,samplingFreq,signalType,paramsT);
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
    v.add(containerFILTER);
    v.add(accFILTER);
    v.add(filterButton);
    v.add(containerFLplot);
    v.add(containerINFO);
    v.add(containerINFplot);
    v.add(containerIndicator).onOK(()=> {}).show();

}
