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

    await call.call();
    let df1 = call.getParamValue('outputDF1');
    let df2 = call.getParamValue('outputDF2');

    return [df1, df2];
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
    if (method === 'Hmisc::aregImpute') {
        completeData = await grok.functions.call('Impute:aregImputeImpl',
            {
                'data':data,
                'burnin':methodparams[0].value
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

    return(completeData);
}


//top-menu: ML | Impute | Visualize
export async function byMethod() {

    //parameter selection function
    function paramSelector(x) {

        if (x === 'mice::mice') {
            let p1 = ui.intInput('Imputation iterations', 20);
            p1.setTooltip('number of imputation iterations');
            methodparams  = [p1];
        }

        if (x === 'Hmisc::aregImpute') {
            let p1 = ui.intInput('Burnin', 5);
            p1.setTooltip('amount of initial iterations to warm up the imputer');
            methodparams  = [p1];
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
        return methodparams;
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
    let method = ui.choiceInput('Algorithm','',['mice::mice','Hmisc::aregImpute',
        'VIM::kNN','missMDA::FAMD','missMDA::PCA','missMDA::MCA']);

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
    let inputs0 = ui.inputs([tableName,columns,varRemovalThreshold,indRemovalThreshold,naCoding]);
    container0.appendChild(inputs0);
    tableName.onChanged(function() {
        dataTable = grok.shell.tableByName(tableName.value);
        columns = ui.columnsInput('Columns', dataTable);
        let inputs1 = ui.inputs([tableName,columns,varRemovalThreshold,indRemovalThreshold,naCoding]);
        container0.replaceChild(inputs1,inputs0);
        inputs0 = inputs1;
    });

    //plotting
    container0.appendChild(ui.button('GENERATE PLOTS',async () => {

        //preprocess and extract metadata
        meta = true;
        let metaOut = await cleanMeta(dataTable,columns,varRemovalThreshold,indRemovalThreshold,naCoding,meta);

        grok.shell.info('Generating: NA correlation, dendrogram and matrix plots');
        grok.shell.addTableView(metaOut[0]);
        grok.shell.addTableView(metaOut[1]);

    }));


    //method parameter selection
    container1.appendChild(ui.inputs([method]));
    method.onChanged(function () {
        $(container2).empty()
        methodparams = paramSelector(method.value);
        container2.appendChild(ui.inputs(methodparams));
    });
    acc0.addPane('method parameters', () => container2)

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
        let completeData = await imputeWithMethod(cleanOut[0],method.value,methodparams);
        grok.shell.addTableView(completeData);
        pi.close();

    }).show();
}

