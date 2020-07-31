/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//input: dataframe data [Input data table]
//input: column_list columns [list of all columns of interest]
//input: double VarRemovalThreshold = 0.5 [% missing per column]
//input: double IndRemovalThreshold = 0.5 [% missing per row]

export async function CleanData(data,columns,VarRemovalThreshold,IndRemovalThreshold) {

    let cleaned = await grok.functions.call('Impute:CleanDataImpl',
        {
            'data': data,
            'columns': columns,
            'VarRemovalThreshold': VarRemovalThreshold,
            'IndRemovalThreshold': IndRemovalThreshold
        });
    grok.shell.addTableView(cleaned);
    return cleaned;
}


