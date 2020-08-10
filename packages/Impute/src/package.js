/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//preprocessing and metadata collection function
async function cleanMeta(data,columns,varRemovalThreshold,indRemovalThreshold,meta) {
    grok.shell.info('Preprocessing data . . .');
    let cleanOut = await grok.functions.call('Impute:cleanMetaImpl',
        {
            'data': data,
            'columns': columns,
            'varRemovalThreshold': varRemovalThreshold.value,
            'indRemovalThreshold': indRemovalThreshold.value,
            'meta': meta
        });
    return(cleanOut);
}

async function imputeWithMethod(data,method,methodparams){
    let completeData;

    //mixed imputation with mice
    if (method === 'mice') {
        completeData = await grok.functions.call('Impute:miceImpl',
            {
                'data':data,
                'impNum':methodparams[0].value,
                'maxIter':methodparams[1].value,
                'whichImp':methodparams[2].value
            });
    }

    //Hmisc aregImpute
    if (method === 'aregImpute') {
        completeData = await grok.functions.call('Impute:aregImputeImpl',
            {
                'data':data,
                'burnin':methodparams[0].value
            });
    }

    if (method === 'KNN') {
        completeData = await grok.functions.call('Impute:knnImpl',
            {
                'data':data,
                'k':methodparams[0].value
            });
    }



    return(completeData);
}


//top-menu: ML | Impute | Visualize
//input: dataframe data [Input data table]
//input: column_list columns [list of all columns of interest]
export async function byMethod(data,columns) {

    //parameter selection function
    function paramSelector(x) {
        if (x === 'mice') {
            let p1 = ui.intInput('impNum', 5);
            let p2 = ui.intInput('maxIter', 20);
            let p3 = ui.intInput('whichImp', 1);
            methodparams  = [p1,p2,p3];
        }

        if (x === 'aregImpute') {
            let p1 = ui.intInput('burnin', 5);
            methodparams  = [p1];
        }

        if (x === 'KNN') {
            let p1 = ui.intInput('nearest neighbours', 10);
            methodparams  = [p1];
        }

        return methodparams;
    }

    let v = ui.dialog('Preprocessing parameters');

    //container0 inputs
    let varRemovalThreshold = ui.floatInput('Column NA max threshold (%)', 0.5);
    let indRemovalThreshold = ui.floatInput('Row NA max threshold(%)', 0.5);
    //acc0 input
    let plotOptions = ui.multiChoiceInput('Plotting options', ['NA correlation'], ['NA correlation', 'matrix plot', 'dendrogram']);

    //container1 inputs
    let method = ui.choiceInput('Algorithm','',['mice','aregImpute','KNN']);

    //metadata collection switch
    let meta;
    let methodparams = [];

    //add containers
    let container0 = ui.div();
    let container1 = ui.div();
    let container2 = ui.div();

    //add accordeons
    let acc0 = ui.accordion();
    let acc1 = ui.accordion();


    //data preprocessing
    container0.appendChild(ui.inputs([varRemovalThreshold, indRemovalThreshold]));
    acc0.addPane('additional plotting', () => ui.div([ui.inputs([plotOptions]), ui.button('PLOT',async () => {

        //preprocess and extract metadata
        meta = true;
        let metaOut = await cleanMeta(data,columns,varRemovalThreshold,indRemovalThreshold,meta);
        grok.shell.info('Generating: ' + plotOptions.value + ' plot(s)');

        // plotting in JS
        let view  =  grok.shell.addTableView(metaOut);

        view.corrPlot({
            xs: metaOut.columns.names(),
            ys: metaOut.columns.names(),
        });

    })]));


    //method parameter selection
    container1.appendChild(ui.inputs([method]));
    method.onChanged(function () {
        $(container2).empty()
        methodparams = paramSelector(method.value);
        container2.appendChild(ui.inputs(methodparams));
    });
    acc1.addPane('method parameters', () => container2)

    //add containers to dialogue
    v.add(container0);
    v.add(acc0);
    v.add(container1);
    v.add(acc1).onOK(async ()=> {

        //preprocess without extracting metadata and impute
        meta = false;
        let cleanOut = await cleanMeta(data,columns,varRemovalThreshold,indRemovalThreshold,meta);
        grok.shell.info('Imputing with: ' + method.value);
        grok.shell.info(methodparams.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));

        //imputation function
        let completeData = await imputeWithMethod(cleanOut,method.value,methodparams);
        grok.shell.addTableView(completeData);

    }).show();
}
