var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import * as grok from 'datagrok-api/grok';
export function juliaDup(s) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:JuliaDup', { s });
    });
}
export function pcaJulia(table, features, components) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:PCAJulia', { table, features, components });
    });
}
export function scatterPlotJulia(t, xColumnName, yColumnName, colorColumnName) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ScatterPlotJulia', { t, xColumnName, yColumnName, colorColumnName });
    });
}
export function tTestJulia(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:TTestJulia', { data, x, y });
    });
}
export function octaveParamsTest(BOOL_INPUT, STRING_INPUT, DOUBLE_INPUT, INT_INPUT, table, COLUMN_INPUT, COLUMN_LIST_INPUT) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:OctaveParamsTest', { BOOL_INPUT, STRING_INPUT, DOUBLE_INPUT, INT_INPUT, table, COLUMN_INPUT, COLUMN_LIST_INPUT });
    });
}
export function bmi(height, weight) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:BMI', { height, weight });
    });
}
export function cellImagingSegmentation(file) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:CellImagingSegmentation', { file });
    });
}
export function datawigTest() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:DatawigTest', {});
    });
}
export function exif(file) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:EXIF', { file });
    });
}
export function funcParamsTest(country, vegetable, saltiness) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:FuncParamsTest', { country, vegetable, saltiness });
    });
}
export function imageClassification(file) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ImageClassification', { file });
    });
}
export function knnPython(data, imputeColumns, dataColumns, neighbours) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:KNNPython', { data, imputeColumns, dataColumns, neighbours });
    });
}
export function pythonDup(s) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:PythonDup', { s });
    });
}
export function scalogramPython(data, signal, sampleRate, w0, removeDc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ScalogramPython', { data, signal, sampleRate, w0, removeDc });
    });
}
export function scatterPlotPython(t, xColumnName, yColumnName, colorColumnName) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ScatterPlotPython', { t, xColumnName, yColumnName, colorColumnName });
    });
}
export function tTestPython(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:TTestPython', { data, x, y });
    });
}
export function pythonParamsTest(i, d, b, s, dt, df, col, cols) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:PythonParamsTest', { i, d, b, s, dt, df, col, cols });
    });
}
export function ulmo() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:Ulmo', {});
    });
}
export function anova(table, categories, variable) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ANOVA', { table, categories, variable });
    });
}
export function acf(data, columns) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ACF', { data, columns });
    });
}
export function bsa(height, weight) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:BSA', { height, weight });
    });
}
export function contourPlot(data, columns) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ContourPlot', { data, columns });
    });
}
export function financialChart(rates, date, open, high, low, close, volume, adjusted) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:FinancialChart', { rates, date, open, high, low, close, volume, adjusted });
    });
}
export function fittingDRC(table, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:FittingDRC', { table, x, y });
    });
}
export function iirFilter(data, signal, order, frequency, sampleRate, type) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:IIRFilter', { data, signal, order, frequency, sampleRate, type });
    });
}
export function knnR(data, columns, neighbours) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:KNNR', { data, columns, neighbours });
    });
}
export function ksTest(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:KSTest', { data, x, y });
    });
}
export function lda(table, predict, features, perc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:LDA', { table, predict, features, perc });
    });
}
export function linearRegression(data, x, y, interceptZero) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:LinearRegression', { data, x, y, interceptZero });
    });
}
export function listPackages() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ListPackages', {});
    });
}
export function lmer(table, features, random, predict, perc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:LMER', { table, features, random, predict, perc });
    });
}
export function manova(table, variable1, variable2, variable3) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:MANOVA', { table, variable1, variable2, variable3 });
    });
}
export function pcaR(T, columns, numComp, center, scale) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:PCAR', { T, columns, numComp, center, scale });
    });
}
export function pls(table, predict, features, components) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:PLS', { table, predict, features, components });
    });
}
export function predictiveModelSVM(datasetPredict, predict, dataset, columns, fillMissing, perc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:PredictiveModelSVM', { datasetPredict, predict, dataset, columns, fillMissing, perc });
    });
}
export function rdup(s) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:RDup', { s });
    });
}
export function renvSpellingExample() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:RenvSpellingExample', {});
    });
}
export function scalogramR(data, signal, sampleRate, octaves, voices, removeDc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ScalogramR', { data, signal, sampleRate, octaves, voices, removeDc });
    });
}
export function scatterPlotR(t, xColumnName, yColumnName, colorColumnName) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ScatterPlotR', { t, xColumnName, yColumnName, colorColumnName });
    });
}
export function sentimentClassification(data, col) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:SentimentClassification', { data, col });
    });
}
export function spectrogram(data, signal, sampleRate, windowLength, timeStep, removeDc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:Spectrogram', { data, signal, sampleRate, windowLength, timeStep, removeDc });
    });
}
export function surfacePlot(t, X, Y, Z) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:SurfacePlot', { t, X, Y, Z });
    });
}
export function tTestR(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:TTestR', { data, x, y });
    });
}
export function ternaryPlot(t, topColumnName, leftColumnName, rightColumnName, pointSize) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:TernaryPlot', { t, topColumnName, leftColumnName, rightColumnName, pointSize });
    });
}
export function rParamsTest(i, d, b, s, dt, df, col, cols) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:RParamsTest', { i, d, b, s, dt, df, col, cols });
    });
}
export function timeSeriesDecomposition(data, dates, observations) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:TimeSeriesDecomposition', { data, dates, observations });
    });
}
export function arimaForecasting(data, dates, observations, P, D, Q, obsForecast) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Demo:Scripts:ARIMAForecasting', { data, dates, observations, P, D, Q, obsForecast });
    });
}
