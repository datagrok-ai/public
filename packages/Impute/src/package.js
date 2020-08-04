/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

async function _clean(data,columns,varRemovalThreshold,indRemovalThreshold,meta) {
    grok.shell.info('Preprocessing data . . .');
    var cleaned = await grok.functions.call('Impute:CleanMetaImpl',
        {
            'data': data,
            'columns': columns,
            'varRemovalThreshold': varRemovalThreshold.value,
            'indRemovalThreshold': indRemovalThreshold.value,
            'meta': meta.value
        });
}


let container;

//top-menu: ML | Impute | Visualize
//input: dataframe data [Input data table]
//input: column_list columns [list of all columns of interest]
export async function ByMethod(data,columns) {

    let v = ui.dialog('Preprocessing parameters');

    let varRemovalThreshold = ui.floatInput('column NA max threshold (%)', 0.5);
    let indRemovalThreshold = ui.floatInput('row NA max threshold(%)', 0.5);
    let plotOptions = ui.multiChoiceInput('plotting options', ['NA correlation'], ['NA correlation', 'matrix plot', 'dendrogram']);
    let meta = false;

    container = ui.div();

    container.appendChild(ui.inputs([varRemovalThreshold, indRemovalThreshold]));
    container.append(ui.button('Preprocess',() => {
        meta = true;
        container.appendChild(ui.inputs([plotOptions]));
    }));


    v.add(container).onOK(()=> _clean(data,columns,varRemovalThreshold,indRemovalThreshold,meta)).show();
          ui.choiceInput()
}