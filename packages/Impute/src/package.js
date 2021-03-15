/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//Preprocessing and metadata collection function
//Passes the dataframe to a python script, which generates plots and collects additional metadata
//Inputs:  data (type: dataframe) - table with missing values
//Outputs: plt1 (type: HTML element) - NA correlation plot
//         plt2 (type: HTML element) - matrix plot
//         plt3 (type: HTML element) - cluster plot
async function cleanMeta(data) {
    grok.shell.info('Preprocessing data . . .');
    let f = await grok.functions.eval("Impute:cleanMetaImpl");

    let call = f.prepare({
        'data': data
    });

    await call.call();
    let plt1 = call.getParamValue('naPlot');
    let plt2 = call.getParamValue('matrixPlot');
    let plt3 = call.getParamValue('clusterPlot')
    return [plt1, plt2, plt3];
}

//Main imputation function
//Calls the correct imputation method with previously specified hyperparameters
//Inputs:  data (type: dataframe) - processed table containing missing values
//         method (type: string) - imputation method selected by the user
//         methodparams (type: list) - list of chosen hyperparameters
//         columns (type: list) - list of columns to be imputed
//Outputs: completeData (type: dataframe) - table of completed columns
async function imputeWithMethod(data,method,methodparams,columns){
    let completeData;

    //mixed imputation with mice
    if (method === 'mice::mice') {
        completeData = await grok.functions.call('Impute:miceImpl',
            {
                'data':data,
                'columns':columns,
                'maxIter':methodparams[0].value
            });
    }

    //Hmisc aregImpute
    if (method === 'Hmisc::aregImpute:pmm' || method === 'Hmisc::aregImpute:regression'
        || method === 'Hmisc::aregImpute:normpmm') {
        completeData = await grok.functions.call('Impute:aregImputeImpl',
            {
                'data':data,
                'columns':columns,
                'nk': methodparams[0].value,
                'tlinear':methodparams[1].value,
                'type': method,
                'pmmtype': methodparams[2].value,
                'match': methodparams[3].value,
                'btMethod': methodparams[4].value,
                'burnin': methodparams[5].value
            });
    }

    //VIM kNN
    if (method === 'VIM::kNN') {
        completeData = await grok.functions.call('Impute:knnImpl',
            {
                'data':data,
                'columns':columns,
                'k':methodparams[0].value
            });
    }

    //Factorial Analysis for Mixed Data"
    if (method === 'missMDA::FAMD') {
        completeData = await grok.functions.call('Impute:imputeFAMDImpl',
            {
                'data':data,
                'columns':columns,
                'ncp':methodparams[0].value,
                'method':methodparams[1].value,
                'regCoeff':methodparams[2].value
            });
    }

    //Principal Component Analysis
    if (method === 'missMDA::PCA') {
        completeData = await grok.functions.call('Impute:imputePCAImpl',
            {
                'data':data,
                'columns':columns,
                'method':methodparams[0].value,
                'regCoeff':methodparams[1].value
            });
    }

    //Multiple Correspondence Analysis
    if (method === 'missMDA::MCA') {
        completeData = await grok.functions.call('Impute:imputeMCAImpl',
            {
                'data':data,
                'columns':columns,
                'ncp':methodparams[0].value,
                'method':methodparams[1].value,
                'regCoeff':methodparams[2].value
            });
    }

    //pcaMethods
    if (method === 'pcaMethods::nipals' || method === 'pcaMethods::ppca'
        || method === 'pcaMethods::bpca' || method === 'pcaMethods::nlpca') {
        completeData = await grok.functions.call('Impute:pcaMethodsImpl',
            {
                'data':data,
                'columns':columns,
                'scaling':methodparams[0].value,
                'method':method
            });
    }

    //random forest imputation with missForest
    if (method === 'missForest::missForest'){
        completeData = await grok.functions.call('Impute:missForestImpl',
            {
                'data':data,
                'columns':columns,
                'maxiter':methodparams[0].value,
                'ntree':methodparams[1].value,
                'decreasing':methodparams[2].value,
                'replace':methodparams[3].value,
                'parallelize':methodparams[4].value
            });
    }

    return(completeData);
}

//Column filtering function
//Makes a new table from specified columns
//Inputs:  data (type: dataframe) - original table
//         colsList (type: list) - names of columns to be kept
//Outputs: t (type: dataframe) - trimmed table
function tableFilter(data,colsList){
    let l = [];
    colsList.forEach(col => l.push(data.columns.byName(col)));
    let t = DG.DataFrame.fromColumns(l);
    return(t);
}

//Warnings generator function
//Alerts user if too many values per column are missing, displays the columns and their stats in a summary table,
//allows user to choose which columns are to be removed, parses column inputs.
//Inputs:  data (type: dataframe) - original processed table with missing values
//         columnsImpute (type: Object) - columns to be filled in
//         columnsData (type: Object) - inferences are made from these columns
//         varRemoveThreshold (type: Object) - minimum allowed fraction of missing values per column
//         summaryTable (type: Object) - empty Object filled in by the warnings functions
//         warningsContainer (type: ui.div()) - container for the info in summaryTable
//Outputs: selectedCols (type: list) - list of all columns selected by the user
//         columnsImpute (type: list) - list of columns to be imputed
//         columnsData (type: list) - list of columns to be used for making inferences
function warnings(data,columnsImpute,columnsData,varRemoveThreshold,summaryTable,warningsContainer){

    columnsImpute = columnsImpute.value.map((col) => col.name).join().split(',');
    columnsData = columnsData.value.map((col) => col.name).join().split(',');
    let selectedCols = [...new Set(columnsData.concat(columnsImpute))].filter(el => el !== '');
    selectedCols.forEach(col => {

        if (data.columns.byName(col).stats.missingValueCount /
            data.columns.byName(col).stats.totalCount > varRemoveThreshold.value) {

            let tempStats = [];
            tempStats[0] = data.columns.byName(col).stats.missingValueCount;
            tempStats[1] = data.columns.byName(col).stats.missingValueCount /
                data.columns.byName(col).stats.totalCount;
            tempStats[2] = ui.boolInput('',true);
            summaryTable[col] = tempStats;
        }
    });

    let tb = ui.table(Object.entries(summaryTable), (item, idx) =>
            [`${item[0]}:`, `${item[1][0]} (${item[1][1]})`, item[1][2].root],
        ['column', 'missing', 'keep']);

    $(warningsContainer).empty()
    warningsContainer.appendChild(ui.h2('Columns summary'));
    warningsContainer.appendChild(tb);

    return([selectedCols,columnsImpute,columnsData]);

}

//Column list filtering function
//Makes sure columns deselected by the user are removed from the list
//Inputs:  selectedCols (type: list) - list of all columns selected by the user
//         summaryTable (type: Object) - Object containing info on which columns to drop
//Outputs: keepCols (type: list) - filtered column list
function colsSieve(selectedCols,summaryTable){
    let dropCols = [];
    Object.keys(summaryTable).forEach(key => {
        if(summaryTable[key][2].value === false){
            dropCols.push(key);
        }
    })
    let keepCols = selectedCols.filter(x => !dropCols.includes(x));
    return(keepCols);
}

//Parameter selection function
//Chooses hyperparameter inputs based on the selected algorithm
//Inputs:  method (type: string) - imputation method selected by the user
//Outputs: methodparams (type: list) - list of ui.inputs related to the selected method
function paramSelector(method) {

    let methodparams;
    if (method === 'mice::mice') {
        let p1 = ui.intInput('Imputation iterations', 20);
        p1.setTooltip('number of imputation iterations');
        methodparams  = [p1];
    }

    if (method === 'Hmisc::aregImpute:pmm' || method === 'Hmisc::aregImpute:regression' || method === 'Hmisc::aregImpute:normpmm') {

        let p1 = ui.intInput('Number of knots',0);
        p1.setTooltip('number of knots to use for continuous variables')

        let p2 = ui.boolInput('Tlinear transform',true);
        p1.setTooltip('set to FALSE to allow a target variable to have a nonlinear left-hand-side transformation')

        let p3 = ui.choiceInput('PMM type',2,[1,2,3]);
        p1.setTooltip('type of matching to be used for predictive mean matching')

        let p4 = ui.choiceInput('Match','weighted',['weighted', 'closest', 'kclosest']);
        p1.setTooltip('donor observation sampling method')

        let p5 = ui.choiceInput('Boot method','simple',['simple','approximate bayesian']);
        p1.setTooltip('')

        let p6 = ui.intInput('Burnin',10);
        p1.setTooltip('more burn-in iterations may be required when multiple variables are missing on the same observations')

        methodparams = [p1, p2, p3, p4, p5, p6];

    }

    if (method === 'VIM::kNN') {
        let p1 = ui.intInput('Nearest neighbours', 10);
        p1.setTooltip('number of nearest neighbours to be used by kNN');
        methodparams  = [p1];
    }

    if (method === 'missMDA::FAMD') {

        let p1 = ui.intInput('ncp', 9);
        p1.setTooltip('number of components used to predict the missing entries ~(ncol() - 2)');

        let p2 = ui.choiceInput('Reconstruction method','',['Regularized','EM']);

        let p3 = ui.floatInput('Regularization coefficient',1);
        p3.setTooltip('Used only if Regularized method is selected');
        methodparams  = [p1, p2, p3];
    }

    if (method === 'missMDA::PCA') {
        let p1 = ui.choiceInput('Reconstruction method','',['Regularized','EM']);

        let p2 = ui.floatInput('Regularization coefficient',1);
        p2.setTooltip('Used only if Regularized method is selected');
        methodparams  = [p1, p2];
    }

    if (method === 'missMDA::MCA') {

        let p1 = ui.intInput('ncp', 9);
        p1.setTooltip('number of components used to predict the missing entries ~(ncol() - 2)');

        let p2 = ui.choiceInput('Reconstruction method','',['Regularized','EM']);

        let p3 = ui.floatInput('Regularization coefficient',1);
        p3.setTooltip('Used only if Regularized method is selected');
        methodparams  = [p1, p2, p3];
    }

    if (method === 'pcaMethods::nipals' || x === 'pcaMethods::ppca' || x === 'pcaMethods::bpca' || x === 'pcaMethods::nlpca') {

        let p1 = ui.choiceInput('Scaling method','uv',['uv','vector','pareto']);
        p1.setTooltip('matrix is centered and normalized using the selected method');
        methodparams = [p1];
    }

    if (method === 'missForest::missForest'){
        let p1 = ui.intInput('maxiter',10);

        let p2 = ui.intInput('ntree',100);

        let p3 = ui.boolInput('decreasing',true);

        let p4 = ui.boolInput('replace',true);

        let p5 = ui.choiceInput('parallelize','no',['no','variables','forests']);

        methodparams = [p1, p2, p3, p4, p5];
    }
    return methodparams;
}

//Method selection function
//Offers a subset of algorithms suitable for a given table
//Inputs:  dfDtype (type: string) - table data type
//Outputs: methodLst (type: ui.choiceInput) - selection of algorithms to choose from
function methodSelector(dfDtype) {

    let methodLst = ui.choiceInput('Algorithm','none',['none']);
    if (dfDtype === 'continuous') {
        methodLst = ui.choiceInput('Algorithm','',['Hmisc::aregImpute:pmm',
            'Hmisc::aregImpute:regression', 'Hmisc::aregImpute:normpmm', 'VIM::kNN', 'mice::mice',
            'missMDA::PCA', 'missMDA::FAMD', 'pcaMethods::nipals',  'pcaMethods::ppca',
            'pcaMethods::bpca', 'pcaMethods::nlpca','missForest::missForest']);
    }

    if (dfDtype === 'categorical') {
        methodLst = ui.choiceInput('Algorithm','',['Hmisc::aregImpute:pmm',
            'VIM::kNN', 'mice::mice', 'missMDA::FAMD', 'missMDA::MCA','missForest::missForest']);
    }

    if (dfDtype === 'mixed') {
        methodLst = ui.choiceInput('Algorithm','',['Hmisc::aregImpute:pmm',
            'VIM::kNN', 'mice::mice', 'missMDA::FAMD','missForest::missForest']);
    }

    return methodLst;
}

//top-menu: ML | Impute | byMethod
export async function byMethod() {

    //Create dialogue window
    let v = ui.dialog('Missing Value Imputation');

    //region 1st Wizard window
    //region INPUTS
    //container0 inputs
    let tableName = ui.choiceInput('Table', null, grok.shell.tableNames);
    let dataTable = grok.shell.tableByName(tableName.value);
    let columnsImpute = ui.columnsInput('Impute', dataTable);
    let columnsData = ui.columnsInput('Impute from', dataTable);
    let varRemoveThreshold = ui.floatInput('Column NA threshold',0.05)
    let pltInputs0 = ui.inputs([tableName,columnsImpute,columnsData,varRemoveThreshold]);

    //create container0, add inputs
    let container0 = ui.div();
    container0.appendChild(pltInputs0);

    //refresh column inputs when changing table
    tableName.onChanged(function() {
        dataTable = grok.shell.tableByName(tableName.value);
        columnsImpute = ui.columnsInput('Impute', dataTable);
        columnsData = ui.columnsInput('Data', dataTable);
        let pltInputs1 = ui.inputs([tableName,columnsImpute,columnsData,varRemoveThreshold]);
        container0.replaceChild(pltInputs1,pltInputs0);
        pltInputs0 = pltInputs1;
    });
    //endregion

    //region WARNINGS
    //create warnings, parse column inputs
    let colsArr;
    let summaryTable;
    let warningsContainer = ui.div();
    columnsImpute.onChanged(function() {
        summaryTable = {};
        colsArr = warnings(dataTable,columnsImpute,columnsData,varRemoveThreshold,summaryTable,warningsContainer);
    });
    columnsData.onChanged(function() {
        summaryTable = {};
        colsArr = warnings(dataTable,columnsImpute,columnsData,varRemoveThreshold,summaryTable,warningsContainer);
    });
    varRemoveThreshold.onChanged(function() {
        summaryTable = {};
        colsArr = warnings(dataTable,columnsImpute,columnsData,varRemoveThreshold,summaryTable,warningsContainer);
    });
    //endregion

    //region PLOTTING
    //display analytic plots
    let plottingContainer = ui.div();
    plottingContainer.appendChild(ui.button('GENERATE PLOTS',async () => {

        //preprocess and extract metadata
        let filteredTable = tableFilter(dataTable,colsSieve(colsArr[0],summaryTable));
        let metaOut = await cleanMeta(filteredTable);

        grok.shell.info('Generating: NA correlation, dendrogram and matrix plots');

        let tableView = grok.shell.getTableView(tableName.value);
        tableView.resetLayout()
        let node1 = tableView.dockManager.dock(metaOut[0], 'right', null, 'naPlot');
        let node2 = tableView.dockManager.dock(metaOut[1], 'fill', node1, 'matrixPlot');
        let node3 = tableView.dockManager.dock(metaOut[2], 'fill', node1, 'clusterPlot');

    }));
    //endregion
    //endregion

    //region 2nd Wizard window
    //region INPUTS
    //container1 inputs
    let dfDtype = ui.choiceInput('Data type','',['continuous','categorical','mixed']);
    let method = methodSelector(dfDtype.value);
    let impInputs0 = ui.inputs( [dfDtype, method]);
    let methodparams;

    //create container1, add inputs
    let container1 = ui.div();
    let container2 = ui.div();
    let acc0 = ui.accordion();
    container1.appendChild(impInputs0);

    //refresh algorithms and parameters input when changing dfDtype
    dfDtype.onChanged(function () {
        $(container2).empty()
        method = methodSelector(dfDtype.value);
        let impInputs1 = ui.inputs([dfDtype, method]);
        container1.replaceChild(impInputs1, impInputs0);
        impInputs0 = impInputs1;
        method.onChanged(function () {
            $(container2).empty()
            methodparams = paramSelector(method.value);
            container2.appendChild(ui.inputs(methodparams));
        });
    });
    acc0.addPane('parameters', () => container2)
    //endregion

    //region ASSEMBLE & IMPUTE
    //add containers to dialogue window
    v.add(container0);
    v.add(warningsContainer);
    v.add(plottingContainer);
    v.add(container1);
    v.add(acc0).onOK(async ()=> {

        //progress indicator
        let pi = DG.TaskBarProgressIndicator.create('Imputing...');

        //display basic info
        grok.shell.info('Imputing with: ' + method.value);
        grok.shell.info(methodparams.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));

        //impute with selected method
        let filteredTable = tableFilter(dataTable,colsSieve(colsArr[0],summaryTable));
        let completeData = await imputeWithMethod(filteredTable,method.value,methodparams,
            colsSieve(colsArr[1],summaryTable));
        grok.shell.addTableView(completeData);
        pi.close();

    }).show();
    //endregion
    //endregion
}

grid.onCellPrepare()