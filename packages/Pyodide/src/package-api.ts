import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  //file and blob
  export async function pyodideBlobInputOutput(blobInput: any): Promise<any> {
    return await grok.functions.call('@datagrok/pyodide:PyodideBlobInputOutput', { blobInput });
  }

  export async function pyodideBool(bool_input: boolean): Promise<boolean> {
    return await grok.functions.call('@datagrok/pyodide:PyodideBool', { bool_input });
  }

  //calc column
  export async function pyodideCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('@datagrok/pyodide:PyodideCalcColumn', { x });
  }

  //column list input
  export async function pyodideColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/pyodide:PyodideColumnList', { df, cols });
  }

  //datetime input/output
  export async function pyodideDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('@datagrok/pyodide:PyodideDate', { input_datetime });
  }

  export async function pyodideDouble(double_input: number): Promise<number> {
    return await grok.functions.call('@datagrok/pyodide:PyodideDouble', { double_input });
  }

  export async function pyodideEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/pyodide:PyodideEmptyDataframe', {});
  }

  //graphics output column input
  export async function pyodideGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/pyodide:PyodideGraphics', { df, xName, yName });
  }

  export async function pyodideInt(integer_input: number): Promise<number> {
    return await grok.functions.call('@datagrok/pyodide:PyodideInt', { integer_input });
  }

  //map input/output
  export async function pyodideMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('@datagrok/pyodide:PyodideMap', { input_map, unique_key });
  }

  //df performance
  export async function pyodideSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/pyodide:PyodideSingleDf', { df });
  }

  export async function pyodideString(string_input: string): Promise<string> {
    return await grok.functions.call('@datagrok/pyodide:PyodideString', { string_input });
  }
}

export namespace funcs {
  export async function initPyodide(): Promise<any> {
    return await grok.functions.call('@datagrok/pyodide:InitPyodide', {});
  }

  export async function makeVectorCode(script: any): Promise<any> {
    return await grok.functions.call('@datagrok/pyodide:MakeVectorCode', { script });
  }

  export async function pyodideLanguageHandler(scriptCall: any): Promise<any> {
    return await grok.functions.call('@datagrok/pyodide:PyodideLanguageHandler', { scriptCall });
  }
}
