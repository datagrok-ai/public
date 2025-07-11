import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  export async function dummyPython(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:DummyPython', {});
  }

  export async function juliaParamsTest(i: number): Promise<number> {
    return await grok.functions.call('CVMTests:JuliaParamsTest', { i });
  }

  //file and blob
  export async function juliaFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CVMTests:JuliaFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function juliaCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CVMTests:JuliaCalcColumn', { x });
  }

  //column list input
  export async function juliaColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:JuliaColumnList', { df, cols });
  }

  export async function juliaIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:JuliaIntColumn', {});
  }

  //df input/output
  export async function juliaDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:JuliaDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function juliaDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CVMTests:JuliaDate', { input_datetime });
  }

  //Echo string script
  export async function juliaEcho(string_input: string): Promise<string> {
    return await grok.functions.call('CVMTests:JuliaEcho', { string_input });
  }

  //df input/output
  export async function juliaEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:JuliaEmptyDataframe', {});
  }

  //file lines count
  export async function juliaLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('CVMTests:JuliaLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function juliaGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CVMTests:JuliaGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function juliaSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CVMTests:JuliaSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function juliaMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('CVMTests:JuliaMap', { input_map, unique_key });
  }

  //df performance
  export async function juliaSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:JuliaSingleDf', { df });
  }

  export async function octaveParamsTest(i: number, d: number, b: boolean, s: string, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('CVMTests:OctaveParamsTest', { i, d, b, s, df });
  }

  //file and blob
  export async function octaveFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CVMTests:OctaveFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function octaveCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CVMTests:OctaveCalcColumn', { x });
  }

  //column list input
  export async function octaveColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:OctaveColumnList', { df, cols });
  }

  //df input/output
  export async function octaveDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:OctaveDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function octaveDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CVMTests:OctaveDate', { input_datetime });
  }

  //Echo string script
  export async function octaveEcho(string_input: string): Promise<string> {
    return await grok.functions.call('CVMTests:OctaveEcho', { string_input });
  }

  //Exponential Operator Test
  export async function octaveExponentialOperatorTest(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:OctaveExponentialOperatorTest', {});
  }

  //file lines count
  export async function octaveLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('CVMTests:OctaveLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function octaveGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CVMTests:OctaveGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function octaveSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CVMTests:OctaveSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function octaveMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('CVMTests:OctaveMap', { input_map, unique_key });
  }

  //df performance
  export async function octaveSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:OctaveSingleDf', { df });
  }

  //df performance
  export async function octaveStdout(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:OctaveStdout', { df });
  }

  export async function pythonParamsTest(i: number, d: number, b: boolean, s: string, m: any, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('CVMTests:PythonParamsTest', { i, d, b, s, m, df });
  }

  //file and blob
  export async function pythonFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CVMTests:PythonFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function pythonCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CVMTests:PythonCalcColumn', { x });
  }

  //column list input
  export async function pythonColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonColumnList', { df, cols });
  }

  //df input/output
  export async function pythonDataframeGraphicsCached(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonDataframeGraphicsCached', { df });
  }

  export async function pythonIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonIntColumn', {});
  }

  //df input/output
  export async function pythonDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function pythonDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CVMTests:PythonDate', { input_datetime });
  }

  //Echo string script
  export async function pythonEcho(string_input: string): Promise<string> {
    return await grok.functions.call('CVMTests:PythonEcho', { string_input });
  }

  //df input/output
  export async function pythonEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonEmptyDataframe', {});
  }

  //anchors count in html
  export async function pythonAnchorsCount(html: string): Promise<number> {
    return await grok.functions.call('CVMTests:PythonAnchorsCount', { html });
  }

  //df input/output
  export async function pythonException(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonException', { df });
  }

  //image pixel count
  export async function imagePixelCount(fileInput: DG.FileInfo): Promise<number> {
    return await grok.functions.call('CVMTests:ImagePixelCount', { fileInput });
  }

  //graphics output column input
  export async function pythonGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CVMTests:PythonGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function pythonSimpleCached(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CVMTests:PythonSimpleCached', { integer_input, double_input, bool_input, string_input });
  }

  //primitive inputs and outputs test
  export async function pythonSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CVMTests:PythonSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function pythonMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('CVMTests:PythonMap', { input_map, unique_key });
  }

  //file lines count
  export async function pythonLinesCount(file: DG.FileInfo): Promise<number> {
    return await grok.functions.call('CVMTests:PythonLinesCount', { file });
  }

  //df performance
  export async function pythonSingleDfCached(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonSingleDfCached', { df });
  }

  //df performance
  export async function pythonSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:PythonSingleDf', { df });
  }

  //choices test
  export async function pythonStringChoices(choices: string): Promise<string> {
    return await grok.functions.call('CVMTests:PythonStringChoices', { choices });
  }

  export async function rParamsTest(i: number, d: number, b: boolean, s: string, dt: any, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('CVMTests:RParamsTest', { i, d, b, s, dt, df });
  }

  //file and blob
  export async function rFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CVMTests:RFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function rCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CVMTests:RCalcColumn', { x });
  }

  //column list input
  export async function rColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:RColumnList', { df, cols });
  }

  export async function rIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:RIntColumn', {});
  }

  //df input/output
  export async function rDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:RDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function rDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CVMTests:RDate', { input_datetime });
  }

  //file lines count
  export async function rLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('CVMTests:RLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function rGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CVMTests:RGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function rSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CVMTests:RSimple', { integer_input, double_input, bool_input, string_input });
  }

  //df performance
  export async function rSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:RSingleDf', { df });
  }

  export async function rowsTest(): Promise<DG.DataFrame> {
    return await grok.functions.call('CVMTests:RowsTest', {});
  }
}

export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('CVMTests:Info', {});
  }
}
