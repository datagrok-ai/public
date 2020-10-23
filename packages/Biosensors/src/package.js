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
    return call.getParamValue('fig');
}

async function  applyFilter(data,samplingFreq,signalType, paramsT){

    let f = await grok.functions.eval("Biosensors:filtersPyphysio");

    let call = f.prepare({
        'ecg_data':data,
        'fsamp':samplingFreq.value,
        'signalType':signalType.value,
        'paramsT':paramsT
    });

    await call.call();
    return call.getParamValue('fig');
}

async function extractInfo(data,samplingFreq,signalType, paramsT,infoType){

    let f = await grok.functions.eval("Biosensors:infoPyphysio");

    let call = f.prepare({
        'ecg_data':data,
        'fsamp':samplingFreq.value,
        'signalType':signalType.value,
        'paramsT':paramsT,
        'info': infoType.value
    });

    await call.call();
    return call.getParamValue('fig');
}

async function toIndicators(data,samplingFreq,signalType,paramsT,infoType,indicator){

    let f = await grok.functions.eval("Biosensors:indicatorsPyphysio");

    let call = f.prepare({
        'ecg_data': data,
        'fsamp': samplingFreq.value,
        'signalType': signalType.value,
        'paramsT': paramsT,
        'info': infoType.value,
        'preset': indicator.value
    });

    await call.call();
    return call.getParamValue('FD_HRV_df');
}

function paramsToTable(filtersLST,allParams){

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
    return paramsT;
}

// Table preprocessing function
// Created a new truncated dataframe
// Input: colsFilter (type: list), list of columns to keep
//        table (type: dataframe), original dataframe
// Output: t (type: dataframe), new truncated dataframe
function tableTrim(cols, table){
    let l = [];
    for (let j = 0; j < cols.length; j++) {
        if (cols[j] !== '') {
            l.push(table.columns.byName(cols[j]));
        }
    }
    let t = DG.DataFrame.fromColumns(l);
    return t;
}


//name: Biosensors
//tags: panel, widgets
//input: dataframe table
//condition: analysisCondition(table)
export function Biosensors(table){

    let v = ui.dialog('DEMO PIPELINE');

    let column = ui.columnsInput('Biosensor', table);
    column.setTooltip('choose one column of biosensor data');

    let samplingFreq = ui.intInput('Sampling frequency', 2048);
    samplingFreq.setTooltip('Number of samples per second');

    let containerImport = ui.div();
    containerImport.appendChild(ui.inputs([column,samplingFreq]));

    v.add(containerImport).onOK(() => {
        grok.shell.addTableView(tableTrim(column.value,table));
    }).show();

}




//name: pipelineDemo
export async function pipelineDemo() {

    function paramSelector(x) {

        let methodparams;
        if (x === 'IIR') {
            let fp = ui.intInput('fp', 45);
            let fs = ui.intInput('fs', 50);
            let ftype = ui.choiceInput('ftype','',['ellip']);
            methodparams = {'fp':fp,'fs':fs,'ftype':ftype};

        }
        if (x === 'normalize') {

            let normMethod = ui.choiceInput('norm_method','',['standard']);
            methodparams = {'normMethod' : normMethod};
        }
        if (x === 'resample') {

            let fout = ui.intInput('fout', 4096);
            let kind = ui.choiceInput('kind','',['cubic']);
            methodparams = {'fout':fout,'kind':kind};
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

    let samplingFreq = ui.intInput('Sampling frequency', 2048);
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
    let tableView = grok.shell.getTableView(tableName.value);

    /////////////////////////////////////////

    // Import dialogue
    let importInputs = ui.inputs([tableName, signalClass, samplingFreq, signalType]);
    containerIMPORT.appendChild(importInputs);
    let node1;
    containerOGplot.appendChild(ui.bigButton('PLOT ORIGINAL',async () => {

        let plotOG = await importPy(dataTable,samplingFreq,signalType);
        node1 = tableView.dockManager.dock(plotOG, 'right', null, 'Original plot');

    }));


    // Filter dialogue
    let filtersLST = [];
    let allParams = [];
    let paramsT;
    let filterInputs = ui.inputs(filtersLST);
    containerFILTER.appendChild(filterInputs);
    let i = 0;
    filterButton.appendChild(ui.button('ADD FILTER',async () => {

        filtersLST[i] = ui.choiceInput('filter №' + (i+1), '', ['IIR','normalize','resample']);
        let filterInputs1 = ui.inputs(filtersLST);

        containerFILTER.replaceChild(filterInputs1,filterInputs);
        filterInputs = filterInputs1;

        filtersLST[i].onChanged(function () {
            $(paramsContainer).empty();
            let val = filtersLST[i-1].value;
            allParams[i-1] = paramSelector(val);
            paramsContainer.appendChild(ui.inputs(Object.values(allParams[i-1])));
        })
        i++;
    }));
    accFILTER.addPane('parameters', () => paramsContainer)
    let node2;
    containerFLplot.appendChild(ui.bigButton('PLOT FILTERED',async () => {

        paramsT = paramsToTable(filtersLST,allParams);
        let plotFL = await applyFilter(dataTable,samplingFreq,signalType,paramsT);
        node2 = tableView.dockManager.dock(plotFL, 'fill', node1, 'Filtered plot');

    }));


    // Information extraction dialogue
    let infoInputs = ui.inputs([infoType]);
    containerINFO.appendChild(infoInputs);
    let node3;
    containerINFplot.appendChild(ui.bigButton('EXTRACT INFO',async () => {

        paramsT = paramsToTable(filtersLST,allParams);
        let plotInfo = await extractInfo(dataTable,samplingFreq,signalType,paramsT,infoType);
        node3 = tableView.dockManager.dock(plotInfo, 'fill', node2, infoType.value);

    }));


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
    v.add(containerIndicator).onOK(async ()=> {

        paramsT = paramsToTable(filtersLST,allParams);
        let indicatorDf = await toIndicators(dataTable,samplingFreq,signalType,paramsT,infoType,indicator);
        let newView = grok.shell.addTableView(indicatorDf);

    }).show();

}

