/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//preprocessing and metadata collection function
async function cleanMeta(data,columns,varRemovalThreshold,indRemovalThreshold,naCoding,meta) {
    grok.shell.info('Preprocessing data . . .');
    let f = await grok.functions.eval("Impute:cleanMetaImpl");

    let call = f.prepare({
        'data': data,
        'columns': columns.value,
        'varRemovalThreshold': varRemovalThreshold.value,
        'indRemovalThreshold': indRemovalThreshold.value,
        'naCoding': naCoding.value,
        'meta': meta
    });

    if (meta === false) {

        await call.call();
        return call.getParamValue('cleanDf');


    } else {

        await call.call();
        let plt1 = call.getParamValue('naPlot');
        let plt2 = call.getParamValue('matrixPlot');
        let plt3 = call.getParamValue('clusterPlot')
        return [plt1, plt2, plt3];

    }

}

async function imputeWithMethod(data,method,methodparams){
    let completeData;

    //mixed imputation with mice
    if (method === 'mice::mice') {
        completeData = await grok.functions.call('Impute:miceImpl',
            {
                'data':data,
                'maxIter':methodparams[0].value
            });
    }

    //Hmisc aregImpute
    if (method === 'Hmisc::aregImpute:pmm' || method === 'Hmisc::aregImpute:regression'
        || method === 'Hmisc::aregImpute:normpmm') {
        completeData = await grok.functions.call('Impute:aregImputeImpl',
            {
                'data':data,
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
                'k':methodparams[0].value
            });
    }

    //Factorial Analysis for Mixed Data"
    if (method === 'missMDA::FAMD') {
        completeData = await grok.functions.call('Impute:imputeFAMDImpl',
            {
                'data':data,
                'catCols':methodparams[0].value,
                'ncp':methodparams[1].value,
                'method':methodparams[2].value,
                'regCoeff':methodparams[3].value
            });
    }

    //Principal Component Analysis
    if (method === 'missMDA::PCA') {
        completeData = await grok.functions.call('Impute:imputePCAImpl',
            {
                'data':data,
                'catCols':methodparams[0].value,
                'method':methodparams[1].value,
                'regCoeff':methodparams[2].value
            });
    }

    //Multiple Correspondence Analysis
    if (method === 'missMDA::MCA') {
        completeData = await grok.functions.call('Impute:imputeMCAImpl',
            {
                'data':data,
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
                'scaling':methodparams[0].value,
                'method':method
            });
    }

    //random forest imputation with missForest
    if (method === 'missForest::missForest'){
        completeData = await grok.functions.call('Impute:missForestImpl',
            {
                'data':data,
                'maxiter':methodparams[0].value,
                'ntree':methodparams[1].value,
                'decreasing':methodparams[2].value,
                'replace':methodparams[3].value,
                'parallelize':methodparams[4].value
            });
    }

    return(completeData);
}


//top-menu: ML | Impute | byMethod
export async function byMethod() {

    //parameter selection function
    function paramSelector(x) {

        if (x === 'mice::mice') {
            let p1 = ui.intInput('Imputation iterations', 20);
            p1.setTooltip('number of imputation iterations');
            methodparams  = [p1];
        }

        if (x === 'Hmisc::aregImpute:pmm' || x === 'Hmisc::aregImpute:regression' || x === 'Hmisc::aregImpute:normpmm') {

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

        if (x === 'VIM::kNN') {
            let p1 = ui.intInput('Nearest neighbours', 10);
            p1.setTooltip('number of nearest neighbours to be used by kNN');
            methodparams  = [p1];
        }

        if (x === 'missMDA::FAMD') {
            let p1 = ui.columnsInput('Categorical columns', dataTable);
            p1.setTooltip('columns to be treated as categorical variables');

            let p2 = ui.intInput('ncp', 9);
            p2.setTooltip('number of components used to predict the missing entries ~(ncol() - 2)');

            let p3 = ui.choiceInput('Reconstruction method','',['Regularized','EM']);

            let p4 = ui.floatInput('Regularization coefficient',1);
            p4.setTooltip('Used only if Regularized method is selected');
            methodparams  = [p1, p2, p3, p4];
        }

        if (x === 'missMDA::PCA') {
            let p1 = ui.columnsInput('Categorical columns', dataTable);
            p1.setTooltip('columns to be treated as categorical variables');

            let p2 = ui.choiceInput('Reconstruction method','',['Regularized','EM']);

            let p3 = ui.floatInput('Regularization coefficient',1);
            p3.setTooltip('Used only if Regularized method is selected');
            methodparams  = [p1, p2, p3];
        }

        if (x === 'missMDA::MCA') {

            let p1 = ui.intInput('ncp', 9);
            p1.setTooltip('number of components used to predict the missing entries ~(ncol() - 2)');

            let p2 = ui.choiceInput('Reconstruction method','',['Regularized','EM']);

            let p3 = ui.floatInput('Regularization coefficient',1);
            p3.setTooltip('Used only if Regularized method is selected');
            methodparams  = [p1, p2, p3];
        }

        if (x === 'pcaMethods::nipals' || x === 'pcaMethods::ppca' || x === 'pcaMethods::bpca' || x === 'pcaMethods::nlpca') {

            let p1 = ui.choiceInput('Scaling method','none',['none','uv','vector','pareto']);
            p1.setTooltip('matrix is centered and normalized using the selected method');
            methodparams = [p1];
        }

        if (x === 'missForest::missForest'){
            let p1 = ui.intInput('maxiter',10);

            let p2 = ui.intInput('ntree',100);

            let p3 = ui.boolInput('decreasing',true);

            let p4 = ui.boolInput('replace',true);

            let p5 = ui.choiceInput('parallelize','no',['no','variables','forests']);

            methodparams = [p1, p2, p3, p4, p5];
        }
        return methodparams;
    }

    function methodSelector(x) {

        let placeHolder = ui.choiceInput('Algorithm','none',['none']);
        if (x === 'continuous') {
            placeHolder = ui.choiceInput('Algorithm','',['Hmisc::aregImpute:pmm',
                'Hmisc::aregImpute:regression', 'Hmisc::aregImpute:normpmm', 'VIM::kNN', 'mice::mice',
                'missMDA::PCA', 'missMDA::FAMD', 'pcaMethods::nipals',  'pcaMethods::ppca',
                'pcaMethods::bpca', 'pcaMethods::nlpca','missForest::missForest']);
        }

        if (x === 'categorical') {
            placeHolder = ui.choiceInput('Algorithm','',['Hmisc::aregImpute:pmm',
                'VIM::kNN', 'mice::mice', 'missMDA::FAMD', 'missMDA::MCA','missForest::missForest']);
        }

        if (x === 'mixed') {
            placeHolder = ui.choiceInput('Algorithm','',['Hmisc::aregImpute:pmm',
                'VIM::kNN', 'mice::mice', 'missMDA::FAMD','missForest::missForest']);
        }

        return placeHolder;
    }

    let v = ui.dialog('Missing Value Imputation');

    //container0 inputs
    let tableName = ui.choiceInput('Table', null, grok.shell.tableNames);
    let dataTable = grok.shell.tableByName(tableName.value);
    let columns = ui.columnsInput('Columns', dataTable);

    let varRemovalThreshold = ui.floatInput('Column NA threshold (%)', 0.5);
    varRemovalThreshold.setTooltip('all columns with NA % above the threshold will be removed');

    let indRemovalThreshold = ui.floatInput('Row NA threshold (%)', 0.5);
    indRemovalThreshold.setTooltip('all rows with NA % above the threshold will be removed');

    let naCoding = ui.floatInput('NA coding', null);
    naCoding.setTooltip('manually convert anomalous values to NA')

    //container1 inputs
    let dfDtype = ui.choiceInput('Data type','',['continuous','categorical','mixed']);
    let method = methodSelector(dfDtype.value);

    //metadata collection switch
    let meta;
    let methodparams = [];

    //add containers
    let container0 = ui.div();
    let container1 = ui.div();
    let container2 = ui.div();

    //add accordeons
    let acc0 = ui.accordion();


    //data preprocessing
    let pltInputs0 = ui.inputs([tableName,columns,varRemovalThreshold,indRemovalThreshold,naCoding]);
    container0.appendChild(pltInputs0);
    tableName.onChanged(function() {
        dataTable = grok.shell.tableByName(tableName.value);
        columns = ui.columnsInput('Columns', dataTable);
        let pltInputs1 = ui.inputs([tableName,columns,varRemovalThreshold,indRemovalThreshold,naCoding]);
        container0.replaceChild(pltInputs1,pltInputs0);
        pltInputs0 = pltInputs1;
    });

    //plotting
    container0.appendChild(ui.button('GENERATE PLOTS',async () => {

        //preprocess and extract metadata
        meta = true;
        let metaOut = await cleanMeta(dataTable,columns,varRemovalThreshold,indRemovalThreshold,naCoding,meta);

        grok.shell.info('Generating: NA correlation, dendrogram and matrix plots');


        let tableView = grok.shell.getTableView(tableName.value);
        let node1 = tableView.dockManager.dock(metaOut[0], 'right', null, 'naPlot');
        let node2 = tableView.dockManager.dock(metaOut[1], 'fill', node1, 'matrixPlot');
        let node3 = tableView.dockManager.dock(metaOut[2], 'fill', node1, 'clusterPlot');

    }));


    //impute method selection
    let impInputs0 = ui.inputs( [dfDtype, method]);
    container1.appendChild(impInputs0);

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

    //add containers to dialogue
    v.add(container0);
    v.add(container1);
    v.add(acc0).onOK(async ()=> {

        //progress indicator
        let pi = DG.TaskBarProgressIndicator.create('Imputing...');

        //preprocess without extracting metadata and impute
        meta = false;
        let cleanOut = await cleanMeta(dataTable,columns,varRemovalThreshold,indRemovalThreshold,naCoding,meta);
        grok.shell.info('Imputing with: ' + method.value);
        grok.shell.info(methodparams.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));

        //imputation function
        let completeData = await imputeWithMethod(cleanOut,method.value,methodparams);
        grok.shell.addTableView(completeData);
        pi.close();

    }).show();
}

