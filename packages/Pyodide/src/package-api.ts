import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {
  //file and blob
  export async function pyodideBlobInputOutput(blobInput: any): Promise<any> {
    return await grok.functions.call('Pyodide:PyodideBlobInputOutput', { blobInput });
  }

  export async function pyodideBool(bool_input: boolean): Promise<boolean> {
    return await grok.functions.call('Pyodide:PyodideBool', { bool_input });
  }

  //calc column
  export async function pyodideCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('Pyodide:PyodideCalcColumn', { x });
  }

  //column list input
  export async function pyodideColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('Pyodide:PyodideColumnList', { df, cols });
  }

  export async function dataFrameInputCallJS(input_df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('Pyodide:DataFrameInputCallJS', { input_df });
  }

  export async function dataFrameInputCallPy(input_df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('Pyodide:DataFrameInputCallPy', { input_df });
  }

  //datetime input/output
  export async function pyodideDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('Pyodide:PyodideDate', { input_datetime });
  }

  export async function datetimeInputCallJS(input_datetime: any): Promise<any> {
    return await grok.functions.call('Pyodide:DatetimeInputCallJS', { input_datetime });
  }

  export async function datetimeInputCallPy(input_datetime: any): Promise<any> {
    return await grok.functions.call('Pyodide:DatetimeInputCallPy', { input_datetime });
  }

  //datetime input/output
  export async function datetimeTestPy(input_datetime: any): Promise<any> {
    return await grok.functions.call('Pyodide:DatetimeTestPy', { input_datetime });
  }

  export async function pyodideDepsTest(bool_input: boolean): Promise<boolean> {
    return await grok.functions.call('Pyodide:PyodideDepsTest', { bool_input });
  }

  export async function pyodideDouble(double_input: number): Promise<number> {
    return await grok.functions.call('Pyodide:PyodideDouble', { double_input });
  }

  //dataframe input/output
  export async function editDFPy(input_df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('Pyodide:EditDFPy', { input_df });
  }

  export async function pyodideEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('Pyodide:PyodideEmptyDataframe', {});
  }

  //graphics output column input
  export async function pyodideGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('Pyodide:PyodideGraphics', { df, xName, yName });
  }

  export async function pyodideInt(integer_input: number): Promise<number> {
    return await grok.functions.call('Pyodide:PyodideInt', { integer_input });
  }

  //map input/output
  export async function pyodideMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('Pyodide:PyodideMap', { input_map, unique_key });
  }

  export async function simpleInputsCallJS(in1: boolean, in2: number, in3: number, in4: string): Promise<boolean> {
    return await grok.functions.call('Pyodide:SimpleInputsCallJS', { in1, in2, in3, in4 });
  }

  export async function simpleInputsCallPy(in1: boolean, in2: number, in3: number, in4: string): Promise<boolean> {
    return await grok.functions.call('Pyodide:SimpleInputsCallPy', { in1, in2, in3, in4 });
  }

  export async function simpleInputsPy(in1: boolean, in2: number, in3: number, in4: string): Promise<boolean> {
    return await grok.functions.call('Pyodide:SimpleInputsPy', { in1, in2, in3, in4 });
  }

  //df performance
  export async function pyodideSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('Pyodide:PyodideSingleDf', { df });
  }

  export async function pyodideString(string_input: string): Promise<string> {
    return await grok.functions.call('Pyodide:PyodideString', { string_input });
  }
}

export namespace Funcs {
  export async function initPyodide(): Promise<any> {
    return await grok.functions.call('Pyodide:InitPyodide', {});
  }

  export async function makeVectorCode(script: any): Promise<any> {
    return await grok.functions.call('Pyodide:MakeVectorCode', { script });
  }

  export async function pyodideLanguageHandler(scriptCall: any): Promise<any> {
    return await grok.functions.call('Pyodide:PyodideLanguageHandler', { scriptCall });
  }
}
