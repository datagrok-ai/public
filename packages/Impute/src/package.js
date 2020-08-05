/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//preprocessing and metadata collection function
async function _clean(data,columns,varRemovalThreshold,indRemovalThreshold,meta) {
    grok.shell.info('Preprocessing data . . .');
    let cleaned_corr = await grok.functions.call('Impute:CleanMetaImpl',
        {
            'data': data,
            'columns': columns,
            'varRemovalThreshold': varRemovalThreshold.value,
            'indRemovalThreshold': indRemovalThreshold.value,
            'meta': meta.value
        });
    return(cleaned_corr)
}


//top-menu: ML | Impute | Visualize
//input: dataframe data [Input data table]
//input: column_list columns [list of all columns of interest]
export async function byMethod(data,columns) {

    //parameter selection function
    function paramSelector(x) {
        if (x === 'PMM') {
            let p1 = ui.intInput('p1', 1);
            let p2 = ui.intInput('p2', 2);
            methodparams  = [p1,p2];
        }
        if (x === 'missForest') {
            let p1 = ui.intInput('f1', 3);
            let p2 = ui.intInput('f2', 4);
            methodparams  = [p1,p2];
        }
        return methodparams;
    }

    let v = ui.dialog('Preprocessing parameters');

    //container0 inputs
    let varRemovalThreshold = ui.floatInput('column NA max threshold (%)', 0.5);
    let indRemovalThreshold = ui.floatInput('row NA max threshold(%)', 0.5);
    //acc0 input
    let plotOptions = ui.multiChoiceInput('plotting options', ['NA correlation'], ['NA correlation', 'matrix plot', 'dendrogram']);

    //container1 inputs
    let method = ui.choiceInput('choose an algorithm','',['PMM','missForest']);

    //metadata collection switch
    let meta = false;
    let methodparams = [];

    //add containers
    let container0 = ui.div();
    let container1 = ui.div();
    let container2 = ui.div();

    //add accordeons
    let acc0 = ui.accordion();
    let acc1 = ui.accordion();


    //
    container0.appendChild(ui.inputs([varRemovalThreshold, indRemovalThreshold]));
    acc0.addPane('additional plotting', () => ui.div([ui.inputs([plotOptions]), ui.button('PLOT',() => {
        meta = true;
        _clean(data,columns,varRemovalThreshold,indRemovalThreshold,meta);
        grok.shell.info('Generating: ' + plotOptions.value + ' plot(s)');
    })]));


    container1.appendChild(ui.inputs([method]));
    method.onChanged(function () {
        $(container2).empty()
        container2.appendChild(ui.inputs(paramSelector(method.value)));
    });

    acc1.addPane('method parameters', () => container2)

    v.add(container0);
    v.add(acc0);
    v.add(container1);
    v.add(acc1).onOK(()=> {
        _clean(data,columns,varRemovalThreshold,indRemovalThreshold,meta);
        grok.shell.info('Imputing with: ' + method.value);
    }).show();

    grok.shell.addTableView(DF)
}