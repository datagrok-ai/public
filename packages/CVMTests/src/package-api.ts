import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  export async function dummyPython(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:DummyPython', {});
  }

  export async function juliaParamsTest(i: number): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaParamsTest', { i });
  }

  //file and blob
  export async function juliaFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function juliaCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaCalcColumn', { x });
  }

  //column list input
  export async function juliaColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaColumnList', { df, cols });
  }

  export async function juliaIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaIntColumn', {});
  }

  //df input/output
  export async function juliaDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function juliaDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaDate', { input_datetime });
  }

  //Echo string script
  export async function juliaEcho(string_input: string): Promise<string> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaEcho', { string_input });
  }

  //df input/output
  export async function juliaEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaEmptyDataframe', {});
  }

  //file lines count
  export async function juliaLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function juliaGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function juliaSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function juliaMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaMap', { input_map, unique_key });
  }

  //df performance
  export async function juliaSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:JuliaSingleDf', { df });
  }

  export async function octaveParamsTest(i: number, d: number, b: boolean, s: string, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveParamsTest', { i, d, b, s, df });
  }

  //file and blob
  export async function octaveFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function octaveCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveCalcColumn', { x });
  }

  //column list input
  export async function octaveColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveColumnList', { df, cols });
  }

  //df input/output
  export async function octaveDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function octaveDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveDate', { input_datetime });
  }

  //Echo string script
  export async function octaveEcho(string_input: string): Promise<string> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveEcho', { string_input });
  }

  //Exponential Operator Test
  export async function octaveExponentialOperatorTest(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveExponentialOperatorTest', {});
  }

  //file lines count
  export async function octaveLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function octaveGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function octaveSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function octaveMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveMap', { input_map, unique_key });
  }

  //df performance
  export async function octaveSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveSingleDf', { df });
  }

  //df performance
  export async function octaveStdout(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:OctaveStdout', { df });
  }

  export async function pythonParamsTest(i: number, d: number, b: boolean, s: string, m: any, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonParamsTest', { i, d, b, s, m, df });
  }

  //file and blob
  export async function pythonFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function pythonCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonCalcColumn', { x });
  }

  //column list input
  export async function pythonColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonColumnList', { df, cols });
  }

  //df input/output
  export async function pythonDataframeGraphicsCached(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonDataframeGraphicsCached', { df });
  }

  export async function pythonIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonIntColumn', {});
  }

  //df input/output
  export async function pythonDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function pythonDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonDate', { input_datetime });
  }

  //Echo string script
  export async function pythonEcho(string_input: string): Promise<string> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonEcho', { string_input });
  }

  //df input/output
  export async function pythonEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonEmptyDataframe', {});
  }

  //anchors count in html
  export async function pythonAnchorsCount(html: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonAnchorsCount', { html });
  }

  //df input/output
  export async function pythonException(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonException', { df });
  }

  //image pixel count
  export async function imagePixelCount(fileInput: DG.FileInfo): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:ImagePixelCount', { fileInput });
  }

  //graphics output column input
  export async function pythonGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function pythonSimpleCached(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonSimpleCached', { integer_input, double_input, bool_input, string_input });
  }

  //primitive inputs and outputs test
  export async function pythonSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function pythonMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonMap', { input_map, unique_key });
  }

  //file lines count
  export async function pythonLinesCount(file: DG.FileInfo): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonLinesCount', { file });
  }

  //df performance
  export async function pythonSingleDfCached(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonSingleDfCached', { df });
  }

  //df performance
  export async function pythonSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonSingleDf', { df });
  }

  //choices test
  export async function pythonStringChoices(choices: string): Promise<string> {
    return await grok.functions.call('@datagrok/cvm-tests:PythonStringChoices', { choices });
  }

  export async function rParamsTest(i: number, d: number, b: boolean, s: string, dt: any, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:RParamsTest', { i, d, b, s, dt, df });
  }

  //file and blob
  export async function rFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('@datagrok/cvm-tests:RFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function rCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:RCalcColumn', { x });
  }

  //column list input
  export async function rColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:RColumnList', { df, cols });
  }

  export async function rIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:RIntColumn', {});
  }

  //df input/output
  export async function rDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:RDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function rDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:RDate', { input_datetime });
  }

  //file lines count
  export async function rLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:RLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function rGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:RGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function rSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('@datagrok/cvm-tests:RSimple', { integer_input, double_input, bool_input, string_input });
  }

  //df performance
  export async function rSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:RSingleDf', { df });
  }

  export async function rowsTest(): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/cvm-tests:RowsTest', {});
  }
}

export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/cvm-tests:Info', {});
  }
}
