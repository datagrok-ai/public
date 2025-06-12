import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {
  export async function dummyPython(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:DummyPython', {});
  }

  export async function juliaParamsTest(i: number): Promise<number> {
    return await grok.functions.call('CvmTests:JuliaParamsTest', { i });
  }

  //file and blob
  export async function juliaFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CvmTests:JuliaFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function juliaCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CvmTests:JuliaCalcColumn', { x });
  }

  //column list input
  export async function juliaColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:JuliaColumnList', { df, cols });
  }

  export async function juliaIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:JuliaIntColumn', {});
  }

  //df input/output
  export async function juliaDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:JuliaDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function juliaDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CvmTests:JuliaDate', { input_datetime });
  }

  //Echo string script
  export async function juliaEcho(string_input: string): Promise<string> {
    return await grok.functions.call('CvmTests:JuliaEcho', { string_input });
  }

  //df input/output
  export async function juliaEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:JuliaEmptyDataframe', {});
  }

  //file lines count
  export async function juliaLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('CvmTests:JuliaLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function juliaGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CvmTests:JuliaGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function juliaSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CvmTests:JuliaSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function juliaMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('CvmTests:JuliaMap', { input_map, unique_key });
  }

  //df performance
  export async function juliaSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:JuliaSingleDf', { df });
  }

  export async function octaveParamsTest(i: number, d: number, b: boolean, s: string, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('CvmTests:OctaveParamsTest', { i, d, b, s, df });
  }

  //file and blob
  export async function octaveFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CvmTests:OctaveFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function octaveCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CvmTests:OctaveCalcColumn', { x });
  }

  //column list input
  export async function octaveColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:OctaveColumnList', { df, cols });
  }

  //df input/output
  export async function octaveDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:OctaveDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function octaveDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CvmTests:OctaveDate', { input_datetime });
  }

  //Echo string script
  export async function octaveEcho(string_input: string): Promise<string> {
    return await grok.functions.call('CvmTests:OctaveEcho', { string_input });
  }

  //Exponential Operator Test
  export async function octaveExponentialOperatorTest(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:OctaveExponentialOperatorTest', {});
  }

  //file lines count
  export async function octaveLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('CvmTests:OctaveLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function octaveGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CvmTests:OctaveGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function octaveSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CvmTests:OctaveSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function octaveMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('CvmTests:OctaveMap', { input_map, unique_key });
  }

  //df performance
  export async function octaveSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:OctaveSingleDf', { df });
  }

  //df performance
  export async function octaveStdout(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:OctaveStdout', { df });
  }

  export async function pythonParamsTest(i: number, d: number, b: boolean, s: string, m: any, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('CvmTests:PythonParamsTest', { i, d, b, s, m, df });
  }

  //file and blob
  export async function pythonFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CvmTests:PythonFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function pythonCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CvmTests:PythonCalcColumn', { x });
  }

  //column list input
  export async function pythonColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonColumnList', { df, cols });
  }

  //df input/output
  export async function pythonDataframeGraphicsCached(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonDataframeGraphicsCached', { df });
  }

  export async function pythonIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonIntColumn', {});
  }

  //df input/output
  export async function pythonDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function pythonDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CvmTests:PythonDate', { input_datetime });
  }

  //Echo string script
  export async function pythonEcho(string_input: string): Promise<string> {
    return await grok.functions.call('CvmTests:PythonEcho', { string_input });
  }

  //df input/output
  export async function pythonEmptyDataframe(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonEmptyDataframe', {});
  }

  //anchors count in html
  export async function pythonAnchorsCount(html: string): Promise<number> {
    return await grok.functions.call('CvmTests:PythonAnchorsCount', { html });
  }

  //df input/output
  export async function pythonException(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonException', { df });
  }

  //image pixel count
  export async function imagePixelCount(fileInput: DG.FileInfo): Promise<number> {
    return await grok.functions.call('CvmTests:ImagePixelCount', { fileInput });
  }

  //graphics output column input
  export async function pythonGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CvmTests:PythonGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function pythonSimpleCached(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CvmTests:PythonSimpleCached', { integer_input, double_input, bool_input, string_input });
  }

  //primitive inputs and outputs test
  export async function pythonSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CvmTests:PythonSimple', { integer_input, double_input, bool_input, string_input });
  }

  //map input/output
  export async function pythonMap(input_map: any, unique_key: string): Promise<any> {
    return await grok.functions.call('CvmTests:PythonMap', { input_map, unique_key });
  }

  //file lines count
  export async function pythonLinesCount(file: DG.FileInfo): Promise<number> {
    return await grok.functions.call('CvmTests:PythonLinesCount', { file });
  }

  //df performance
  export async function pythonSingleDfCached(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonSingleDfCached', { df });
  }

  //df performance
  export async function pythonSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:PythonSingleDf', { df });
  }

  //choices test
  export async function pythonStringChoices(choices: string): Promise<string> {
    return await grok.functions.call('CvmTests:PythonStringChoices', { choices });
  }

  export async function rParamsTest(i: number, d: number, b: boolean, s: string, dt: any, df: DG.DataFrame): Promise<number> {
    return await grok.functions.call('CvmTests:RParamsTest', { i, d, b, s, dt, df });
  }

  //file and blob
  export async function rFileBlobInputOutput(fileInput: DG.FileInfo, blobInput: any): Promise<DG.FileInfo> {
    return await grok.functions.call('CvmTests:RFileBlobInputOutput', { fileInput, blobInput });
  }

  //calc column
  export async function rCalcColumn(x: number): Promise<number> {
    return await grok.functions.call('CvmTests:RCalcColumn', { x });
  }

  //column list input
  export async function rColumnList(df: DG.DataFrame, cols: string[]): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:RColumnList', { df, cols });
  }

  export async function rIntColumn(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:RIntColumn', {});
  }

  //df input/output
  export async function rDataframe(df: DG.DataFrame, dfNumerical: DG.DataFrame, dfCategorical: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:RDataframe', { df, dfNumerical, dfCategorical });
  }

  //datetime input/output
  export async function rDate(input_datetime: any): Promise<any> {
    return await grok.functions.call('CvmTests:RDate', { input_datetime });
  }

  //file lines count
  export async function rLinesCount(file: DG.FileInfo, header: boolean, separator: string, dec: string): Promise<number> {
    return await grok.functions.call('CvmTests:RLinesCount', { file, header, separator, dec });
  }

  //graphics output column input
  export async function rGraphics(df: DG.DataFrame, xName: DG.Column, yName: DG.Column): Promise<any> {
    return await grok.functions.call('CvmTests:RGraphics', { df, xName, yName });
  }

  //primitive inputs and outputs test
  export async function rSimple(integer_input: number, double_input: number, bool_input: boolean, string_input: string): Promise<number> {
    return await grok.functions.call('CvmTests:RSimple', { integer_input, double_input, bool_input, string_input });
  }

  //df performance
  export async function rSingleDf(df: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:RSingleDf', { df });
  }

  export async function rowsTest(): Promise<DG.DataFrame> {
    return await grok.functions.call('CvmTests:RowsTest', {});
  }
}

export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('CvmTests:Info', {});
  }
}
