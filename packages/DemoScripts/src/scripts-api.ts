import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export async function juliaDup(s: string): Promise<string> {
  return await grok.functions.call('DemoScripts:JuliaDup', { s });
}

export async function pcaJulia(table: DG.DataFrame, features: string[], components: number): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:PCAJulia', { table, features, components });
}

export async function scatterPlotJulia(t: DG.DataFrame, xColumnName: DG.Column, yColumnName: DG.Column, colorColumnName: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:ScatterPlotJulia', { t, xColumnName, yColumnName, colorColumnName });
}

export async function tTestJulia(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
  return await grok.functions.call('DemoScripts:TTestJulia', { data, x, y });
}

export async function octaveParamsTest(BOOL_INPUT: boolean, STRING_INPUT: string, DOUBLE_INPUT: number, INT_INPUT: number, table: DG.DataFrame, COLUMN_INPUT: DG.Column, COLUMN_LIST_INPUT: string[]): Promise<boolean> {
  return await grok.functions.call('DemoScripts:OctaveParamsTest', { BOOL_INPUT, STRING_INPUT, DOUBLE_INPUT, INT_INPUT, table, COLUMN_INPUT, COLUMN_LIST_INPUT });
}

export async function bmi(height: number, weight: number): Promise<number> {
  return await grok.functions.call('DemoScripts:BMI', { height, weight });
}

export async function cellImagingSegmentation(file: DG.FileInfo): Promise<number> {
  return await grok.functions.call('DemoScripts:CellImagingSegmentation', { file });
}

export async function datawigTest(): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:DatawigTest', {});
}

export async function exif(file: DG.FileInfo): Promise<any> {
  return await grok.functions.call('DemoScripts:EXIF', { file });
}

export async function funcParamsTest(country: string, vegetable: string, saltiness: number): Promise<void> {
  return await grok.functions.call('DemoScripts:FuncParamsTest', { country, vegetable, saltiness });
}

export async function imageClassification(file: DG.FileInfo): Promise<any> {
  return await grok.functions.call('DemoScripts:ImageClassification', { file });
}

export async function knnPython(data: DG.DataFrame, imputeColumns: string[], dataColumns: string[], neighbours: number): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:KNNPython', { data, imputeColumns, dataColumns, neighbours });
}

export async function pythonDup(s: string): Promise<string> {
  return await grok.functions.call('DemoScripts:PythonDup', { s });
}

export async function scalogramPython(data: DG.DataFrame, signal: DG.Column, sampleRate: number, w0: number, removeDc: boolean): Promise<any> {
  return await grok.functions.call('DemoScripts:ScalogramPython', { data, signal, sampleRate, w0, removeDc });
}

export async function scatterPlotPython(t: DG.DataFrame, xColumnName: DG.Column, yColumnName: DG.Column, colorColumnName: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:ScatterPlotPython', { t, xColumnName, yColumnName, colorColumnName });
}

export async function tTestPython(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
  return await grok.functions.call('DemoScripts:TTestPython', { data, x, y });
}

export async function pythonParamsTest(i: number, d: number, b: boolean, s: string, dt: any, df: DG.DataFrame, col: DG.Column, cols: string[]): Promise<number> {
  return await grok.functions.call('DemoScripts:PythonParamsTest', { i, d, b, s, dt, df, col, cols });
}

export async function ulmo(): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:Ulmo', {});
}

export async function anova(table: DG.DataFrame, categories: DG.Column, variable: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:ANOVA', { table, categories, variable });
}

export async function acf(data: DG.DataFrame, columns: string[]): Promise<any> {
  return await grok.functions.call('DemoScripts:ACF', { data, columns });
}

export async function bsa(height: number, weight: number): Promise<number> {
  return await grok.functions.call('DemoScripts:BSA', { height, weight });
}

export async function contourPlot(data: DG.DataFrame, columns: string[]): Promise<any> {
  return await grok.functions.call('DemoScripts:ContourPlot', { data, columns });
}

export async function financialChart(rates: DG.DataFrame, date: DG.Column, open: DG.Column, high: DG.Column, low: DG.Column, close: DG.Column, volume: DG.Column, adjusted: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:FinancialChart', { rates, date, open, high, low, close, volume, adjusted });
}

export async function fittingDRC(table: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:FittingDRC', { table, x, y });
}

export async function iirFilter(data: DG.DataFrame, signal: DG.Column, order: number, frequency: number, sampleRate: number, type: string): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:IIRFilter', { data, signal, order, frequency, sampleRate, type });
}

export async function knnR(data: DG.DataFrame, columns: string[], neighbours: number): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:KNNR', { data, columns, neighbours });
}

export async function ksTest(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
  return await grok.functions.call('DemoScripts:KSTest', { data, x, y });
}

export async function lda(table: DG.DataFrame, predict: DG.Column, features: string[], perc: number): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:LDA', { table, predict, features, perc });
}

export async function linearRegression(data: DG.DataFrame, x: DG.Column, y: DG.Column, interceptZero: boolean): Promise<number> {
  return await grok.functions.call('DemoScripts:LinearRegression', { data, x, y, interceptZero });
}

export async function listPackages(): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:ListPackages', {});
}

export async function lmer(table: DG.DataFrame, features: string[], random: DG.Column, predict: DG.Column, perc: number): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:LMER', { table, features, random, predict, perc });
}

export async function manova(table: DG.DataFrame, variable1: DG.Column, variable2: DG.Column, variable3: DG.Column): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:MANOVA', { table, variable1, variable2, variable3 });
}

export async function pcaR(T: DG.DataFrame, columns: string[], numComp: number, center: boolean, scale: boolean): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:PCAR', { T, columns, numComp, center, scale });
}

export async function pls(table: DG.DataFrame, predict: DG.Column, features: string[], components: number): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:PLS', { table, predict, features, components });
}

export async function predictiveModelSVM(datasetPredict: DG.DataFrame, predict: DG.Column, dataset: DG.DataFrame, columns: string[], fillMissing: boolean, perc: number): Promise<any> {
  return await grok.functions.call('DemoScripts:PredictiveModelSVM', { datasetPredict, predict, dataset, columns, fillMissing, perc });
}

export async function rdup(s: string): Promise<string> {
  return await grok.functions.call('DemoScripts:RDup', { s });
}

export async function renvSpellingExample(): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:RenvSpellingExample', {});
}

export async function scalogramR(data: DG.DataFrame, signal: DG.Column, sampleRate: number, octaves: number, voices: number, removeDc: boolean): Promise<any> {
  return await grok.functions.call('DemoScripts:ScalogramR', { data, signal, sampleRate, octaves, voices, removeDc });
}

export async function scatterPlotR(t: DG.DataFrame, xColumnName: DG.Column, yColumnName: DG.Column, colorColumnName: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:ScatterPlotR', { t, xColumnName, yColumnName, colorColumnName });
}

export async function sentimentClassification(data: DG.DataFrame, col: DG.Column): Promise<DG.DataFrame> {
  return await grok.functions.call('DemoScripts:SentimentClassification', { data, col });
}

export async function spectrogram(data: DG.DataFrame, signal: DG.Column, sampleRate: number, windowLength: number, timeStep: number, removeDc: boolean): Promise<any> {
  return await grok.functions.call('DemoScripts:Spectrogram', { data, signal, sampleRate, windowLength, timeStep, removeDc });
}

export async function surfacePlot(t: DG.DataFrame, X: DG.Column, Y: DG.Column, Z: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:SurfacePlot', { t, X, Y, Z });
}

export async function tTestR(data: DG.DataFrame, x: DG.Column, y: DG.Column): Promise<number> {
  return await grok.functions.call('DemoScripts:TTestR', { data, x, y });
}

export async function ternaryPlot(t: DG.DataFrame, topColumnName: DG.Column, leftColumnName: DG.Column, rightColumnName: DG.Column, pointSize: number): Promise<any> {
  return await grok.functions.call('DemoScripts:TernaryPlot', { t, topColumnName, leftColumnName, rightColumnName, pointSize });
}

export async function rParamsTest(i: number, d: number, b: boolean, s: string, dt: any, df: DG.DataFrame, col: DG.Column, cols: string[]): Promise<number> {
  return await grok.functions.call('DemoScripts:RParamsTest', { i, d, b, s, dt, df, col, cols });
}

export async function timeSeriesDecomposition(data: DG.DataFrame, dates: DG.Column, observations: DG.Column): Promise<any> {
  return await grok.functions.call('DemoScripts:TimeSeriesDecomposition', { data, dates, observations });
}

export async function arimaForecasting(data: DG.DataFrame, dates: DG.Column, observations: DG.Column, P: number, D: number, Q: number, obsForecast: number): Promise<any> {
  return await grok.functions.call('DemoScripts:ARIMAForecasting', { data, dates, observations, P, D, Q, obsForecast });
}
