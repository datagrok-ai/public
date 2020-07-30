/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: TemplateJS
//input: dataframe data [Input data table]
//input: column_list columns [list of all columns of interest]
//input: double VarRemovalThreshold = 0.5 [% missingness per column]
//input: double IndRemovalThreshold = 0.5 [% missingness per row]
export async function FFF(data,columns,VarRemovalThreshold,IndRemovalThreshold) {

    console.log(data);
    let cleaned = await grok.functions.call('Impute:CleanData',
        {'data':data,'columns':columns,'VarRemovalThreshold':VarRemovalThreshold,'IndRemovalThreshold':IndRemovalThreshold});

    grok.shell.addTableView(cleaned);
    return cleaned;
}


