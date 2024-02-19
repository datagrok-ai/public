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
// test: juliaDup("abc") == "abcabc"
export function juliaDup(s) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:JuliaDup', { s });
    });
}
export function pcaJulia(table, features, components) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:PCAJulia', { table, features, components });
    });
}
export function scatterPlotJulia(t, xColumnName, yColumnName, colorColumnName) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ScatterPlotJulia', { t, xColumnName, yColumnName, colorColumnName });
    });
}
export function tTestJulia(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:TTestJulia', { data, x, y });
    });
}
export function octaveParamsTest(BOOL_INPUT, STRING_INPUT, DOUBLE_INPUT, INT_INPUT, table, COLUMN_INPUT, COLUMN_LIST_INPUT) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:OctaveParamsTest', { BOOL_INPUT, STRING_INPUT, DOUBLE_INPUT, INT_INPUT, table, COLUMN_INPUT, COLUMN_LIST_INPUT });
    });
}
export function bmi(height, weight) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:BMI', { height, weight });
    });
}
export function cellImagingSegmentation(file) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:CellImagingSegmentation', { file });
    });
}
export function datawigTest() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:DatawigTest', {});
    });
}
export function exif(file) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:EXIF', { file });
    });
}
export function funcParamsTest(country, vegetable, saltiness) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:FuncParamsTest', { country, vegetable, saltiness });
    });
}
export function imageClassification(file) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ImageClassification', { file });
    });
}
export function knnPython(data, imputeColumns, dataColumns, neighbours) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:KNNPython', { data, imputeColumns, dataColumns, neighbours });
    });
}
export function pythonDup(s) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:PythonDup', { s });
    });
}
export function scalogramPython(data, signal, sampleRate, w0, removeDc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ScalogramPython', { data, signal, sampleRate, w0, removeDc });
    });
}
export function scatterPlotPython(t, xColumnName, yColumnName, colorColumnName) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ScatterPlotPython', { t, xColumnName, yColumnName, colorColumnName });
    });
}
export function tTestPython(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:TTestPython', { data, x, y });
    });
}
export function pythonParamsTest(i, d, b, s, dt, df, col, cols) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:PythonParamsTest', { i, d, b, s, dt, df, col, cols });
    });
}
export function ulmo() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:Ulmo', {});
    });
}
export function anova(table, categories, variable) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ANOVA', { table, categories, variable });
    });
}
export function acf(data, columns) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ACF', { data, columns });
    });
}
export function bsa(height, weight) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:BSA', { height, weight });
    });
}
export function contourPlot(data, columns) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ContourPlot', { data, columns });
    });
}
export function financialChart(rates, date, open, high, low, close, volume, adjusted) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:FinancialChart', { rates, date, open, high, low, close, volume, adjusted });
    });
}
export function fittingDRC(table, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:FittingDRC', { table, x, y });
    });
}
export function iirFilter(data, signal, order, frequency, sampleRate, type) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:IIRFilter', { data, signal, order, frequency, sampleRate, type });
    });
}
export function knnR(data, columns, neighbours) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:KNNR', { data, columns, neighbours });
    });
}
export function ksTest(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:KSTest', { data, x, y });
    });
}
export function lda(table, predict, features, perc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:LDA', { table, predict, features, perc });
    });
}
export function linearRegression(data, x, y, interceptZero) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:LinearRegression', { data, x, y, interceptZero });
    });
}
export function listPackages() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ListPackages', {});
    });
}
export function lmer(table, features, random, predict, perc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:LMER', { table, features, random, predict, perc });
    });
}
export function manova(table, variable1, variable2, variable3) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:MANOVA', { table, variable1, variable2, variable3 });
    });
}
export function pcaR(T, columns, numComp, center, scale) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:PCAR', { T, columns, numComp, center, scale });
    });
}
export function pls(table, predict, features, components) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:PLS', { table, predict, features, components });
    });
}
export function predictiveModelSVM(datasetPredict, predict, dataset, columns, fillMissing, perc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:PredictiveModelSVM', { datasetPredict, predict, dataset, columns, fillMissing, perc });
    });
}
export function rdup(s) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:RDup', { s });
    });
}
export function renvSpellingExample() {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:RenvSpellingExample', {});
    });
}
export function scalogramR(data, signal, sampleRate, octaves, voices, removeDc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ScalogramR', { data, signal, sampleRate, octaves, voices, removeDc });
    });
}
export function scatterPlotR(t, xColumnName, yColumnName, colorColumnName) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ScatterPlotR', { t, xColumnName, yColumnName, colorColumnName });
    });
}
export function sentimentClassification(data, col) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:SentimentClassification', { data, col });
    });
}
export function spectrogram(data, signal, sampleRate, windowLength, timeStep, removeDc) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:Spectrogram', { data, signal, sampleRate, windowLength, timeStep, removeDc });
    });
}
export function surfacePlot(t, X, Y, Z) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:SurfacePlot', { t, X, Y, Z });
    });
}
export function tTestR(data, x, y) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:TTestR', { data, x, y });
    });
}
export function ternaryPlot(t, topColumnName, leftColumnName, rightColumnName, pointSize) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:TernaryPlot', { t, topColumnName, leftColumnName, rightColumnName, pointSize });
    });
}
export function rParamsTest(i, d, b, s, dt, df, col, cols) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:RParamsTest', { i, d, b, s, dt, df, col, cols });
    });
}
export function timeSeriesDecomposition(data, dates, observations) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:TimeSeriesDecomposition', { data, dates, observations });
    });
}
export function arimaForecasting(data, dates, observations, P, D, Q, obsForecast) {
    return __awaiter(this, void 0, void 0, function* () {
        return yield grok.functions.call('Samples:Scripts:ARIMAForecasting', { data, dates, observations, P, D, Q, obsForecast });
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoic2NyaXB0cy1hcGkuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyJzY3JpcHRzLWFwaS50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7Ozs7Ozs7QUFBQSxPQUFPLEtBQUssSUFBSSxNQUFNLG1CQUFtQixDQUFDO0FBSTFDLG9DQUFvQztBQUNwQyxNQUFNLFVBQWdCLFFBQVEsQ0FBQyxDQUFTOztRQUN0QyxPQUFPLE1BQU0sSUFBSSxDQUFDLFNBQVMsQ0FBQyxJQUFJLENBQUMsMEJBQTBCLEVBQUUsRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDO0lBQ3RFLENBQUM7Q0FBQTtBQUVELE1BQU0sVUFBZ0IsUUFBUSxDQUFDLEtBQW1CLEVBQUUsUUFBa0IsRUFBRSxVQUFrQjs7UUFDeEYsT0FBTyxNQUFNLElBQUksQ0FBQyxTQUFTLENBQUMsSUFBSSxDQUFDLDBCQUEwQixFQUFFLEVBQUUsS0FBSyxFQUFFLFFBQVEsRUFBRSxVQUFVLEVBQUUsQ0FBQyxDQUFDO0lBQ2hHLENBQUM7Q0FBQTtBQUVELE1BQU0sVUFBZ0IsZ0JBQWdCLENBQUMsQ0FBZSxFQUFFLFdBQXNCLEVBQUUsV0FBc0IsRUFBRSxlQUEwQjs7UUFDaEksT0FBTyxNQUFNLElBQUksQ0FBQyxTQUFTLENBQUMsSUFBSSxDQUFDLGtDQUFrQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxXQUFXLEVBQUUsZUFBZSxFQUFFLENBQUMsQ0FBQztJQUN6SCxDQUFDO0NBQUE7QUFFRCxNQUFNLFVBQWdCLFVBQVUsQ0FBQyxJQUFrQixFQUFFLENBQVksRUFBRSxDQUFZOztRQUM3RSxPQUFPLE1BQU0sSUFBSSxDQUFDLFNBQVMsQ0FBQyxJQUFJLENBQUMsNEJBQTRCLEVBQUUsRUFBRSxJQUFJLEVBQUUsQ0FBQyxFQUFFLENBQUMsRUFBRSxDQUFDLENBQUM7SUFDakYsQ0FBQztDQUFBO0FBRUQsTUFBTSxVQUFnQixnQkFBZ0IsQ0FBQyxVQUFtQixFQUFFLFlBQW9CLEVBQUUsWUFBb0IsRUFBRSxT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QixXQUFXLENBQUMsQ0FBUyxFQUFFLENBQVMsRUFBRSxDQUFVLEVBQUUsQ0FBUyxFQUFFLEVBQU8sRUFBRSxFQUFnQixFQUFFLEdBQWMsRUFBRSxJQUFjOztRQUN0SSxPQUFPLE1BQU0sSUFBSSxDQUFDLFNBQVMsQ0FBQyxJQUFJLENBQUMsNkJBQTZCLEVBQUUsRUFBRSxDQUFDLEVBQUUsQ0FBQyxFQUFFLENBQUMsRUFBRSxDQUFDLEVBQUUsRUFBRSxFQUFFLEVBQUUsRUFBRSxHQUFHLEVBQUUsSUFBSSxFQUFFLENBQUMsQ0FBQztJQUNyRyxDQUFDO0NBQUE7QUFFRCxNQUFNLFVBQWdCLHVCQUF1QixDQUFDLElBQWtCLEVBQUUsS0FBZ0IsRUFBRSxZQUF1Qjs7UUFDekcsT0FBTyxNQUFNLElBQUksQ0FBQyxTQUFTLENBQUMsSUFBSSxDQUFDLHlDQUF5QyxFQUFFLEVBQUUsSUFBSSxFQUFFLEtBQUssRUFBRSxZQUFZLEVBQUUsQ0FBQyxDQUFDO0lBQzdHLENBQUM7Q0FBQTtBQUVELE1BQU0sVUFBZ0IsZ0JBQWdCLENBQUMsSUFBa0IsRUFBRSxLQUFnQixFQUFFLFlBQXVCLEVBQUUsQ0FBUyxFQUFFLENBQVMsRUFBRSxDQUFTLEVBQUUsV0FBbUI7O1FBQ3hKLE9BQU8sTUFBTSxJQUFJLENBQUMsU0FBUyxDQUFDLElBQUksQ0FBQyxrQ0FBa0MsRUFBRSxFQUFFLElBQUksRUFBRSxLQUFLLEVBQUUsWUFBWSxFQUFFLENBQUMsRUFBRSxDQUFDLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxDQUFDLENBQUM7SUFDNUgsQ0FBQztDQUFBIiwic291cmNlc0NvbnRlbnQiOlsiaW1wb3J0ICogYXMgZ3JvayBmcm9tICdkYXRhZ3Jvay1hcGkvZ3Jvayc7XG5pbXBvcnQgKiBhcyBERyBmcm9tICdkYXRhZ3Jvay1hcGkvZGcnO1xuXG5cbi8vIHRlc3Q6IGp1bGlhRHVwKFwiYWJjXCIpID09IFwiYWJjYWJjXCJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBqdWxpYUR1cChzOiBzdHJpbmcpOiBQcm9taXNlPHN0cmluZz4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkp1bGlhRHVwJywgeyBzIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gcGNhSnVsaWEodGFibGU6IERHLkRhdGFGcmFtZSwgZmVhdHVyZXM6IHN0cmluZ1tdLCBjb21wb25lbnRzOiBudW1iZXIpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlBDQUp1bGlhJywgeyB0YWJsZSwgZmVhdHVyZXMsIGNvbXBvbmVudHMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBzY2F0dGVyUGxvdEp1bGlhKHQ6IERHLkRhdGFGcmFtZSwgeENvbHVtbk5hbWU6IERHLkNvbHVtbiwgeUNvbHVtbk5hbWU6IERHLkNvbHVtbiwgY29sb3JDb2x1bW5OYW1lOiBERy5Db2x1bW4pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlNjYXR0ZXJQbG90SnVsaWEnLCB7IHQsIHhDb2x1bW5OYW1lLCB5Q29sdW1uTmFtZSwgY29sb3JDb2x1bW5OYW1lIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdFRlc3RKdWxpYShkYXRhOiBERy5EYXRhRnJhbWUsIHg6IERHLkNvbHVtbiwgeTogREcuQ29sdW1uKTogUHJvbWlzZTxudW1iZXI+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpUVGVzdEp1bGlhJywgeyBkYXRhLCB4LCB5IH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gb2N0YXZlUGFyYW1zVGVzdChCT09MX0lOUFVUOiBib29sZWFuLCBTVFJJTkdfSU5QVVQ6IHN0cmluZywgRE9VQkxFX0lOUFVUOiBudW1iZXIsIElOVF9JTlBVVDogbnVtYmVyLCB0YWJsZTogREcuRGF0YUZyYW1lLCBDT0xVTU5fSU5QVVQ6IERHLkNvbHVtbiwgQ09MVU1OX0xJU1RfSU5QVVQ6IHN0cmluZ1tdKTogUHJvbWlzZTxib29sZWFuPiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6T2N0YXZlUGFyYW1zVGVzdCcsIHsgQk9PTF9JTlBVVCwgU1RSSU5HX0lOUFVULCBET1VCTEVfSU5QVVQsIElOVF9JTlBVVCwgdGFibGUsIENPTFVNTl9JTlBVVCwgQ09MVU1OX0xJU1RfSU5QVVQgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBibWkoaGVpZ2h0OiBudW1iZXIsIHdlaWdodDogbnVtYmVyKTogUHJvbWlzZTxudW1iZXI+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpCTUknLCB7IGhlaWdodCwgd2VpZ2h0IH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gY2VsbEltYWdpbmdTZWdtZW50YXRpb24oZmlsZTogREcuRmlsZUluZm8pOiBQcm9taXNlPG51bWJlcj4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkNlbGxJbWFnaW5nU2VnbWVudGF0aW9uJywgeyBmaWxlIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gZGF0YXdpZ1Rlc3QoKTogUHJvbWlzZTxERy5EYXRhRnJhbWU+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpEYXRhd2lnVGVzdCcsIHt9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGV4aWYoZmlsZTogREcuRmlsZUluZm8pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkVYSUYnLCB7IGZpbGUgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBmdW5jUGFyYW1zVGVzdChjb3VudHJ5OiBzdHJpbmcsIHZlZ2V0YWJsZTogc3RyaW5nLCBzYWx0aW5lc3M6IG51bWJlcik6IFByb21pc2U8dm9pZD4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkZ1bmNQYXJhbXNUZXN0JywgeyBjb3VudHJ5LCB2ZWdldGFibGUsIHNhbHRpbmVzcyB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGltYWdlQ2xhc3NpZmljYXRpb24oZmlsZTogREcuRmlsZUluZm8pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkltYWdlQ2xhc3NpZmljYXRpb24nLCB7IGZpbGUgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBrbm5QeXRob24oZGF0YTogREcuRGF0YUZyYW1lLCBpbXB1dGVDb2x1bW5zOiBzdHJpbmdbXSwgZGF0YUNvbHVtbnM6IHN0cmluZ1tdLCBuZWlnaGJvdXJzOiBudW1iZXIpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOktOTlB5dGhvbicsIHsgZGF0YSwgaW1wdXRlQ29sdW1ucywgZGF0YUNvbHVtbnMsIG5laWdoYm91cnMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBweXRob25EdXAoczogc3RyaW5nKTogUHJvbWlzZTxzdHJpbmc+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpQeXRob25EdXAnLCB7IHMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBzY2Fsb2dyYW1QeXRob24oZGF0YTogREcuRGF0YUZyYW1lLCBzaWduYWw6IERHLkNvbHVtbiwgc2FtcGxlUmF0ZTogbnVtYmVyLCB3MDogbnVtYmVyLCByZW1vdmVEYzogYm9vbGVhbik6IFByb21pc2U8YW55PiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6U2NhbG9ncmFtUHl0aG9uJywgeyBkYXRhLCBzaWduYWwsIHNhbXBsZVJhdGUsIHcwLCByZW1vdmVEYyB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHNjYXR0ZXJQbG90UHl0aG9uKHQ6IERHLkRhdGFGcmFtZSwgeENvbHVtbk5hbWU6IERHLkNvbHVtbiwgeUNvbHVtbk5hbWU6IERHLkNvbHVtbiwgY29sb3JDb2x1bW5OYW1lOiBERy5Db2x1bW4pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlNjYXR0ZXJQbG90UHl0aG9uJywgeyB0LCB4Q29sdW1uTmFtZSwgeUNvbHVtbk5hbWUsIGNvbG9yQ29sdW1uTmFtZSB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRUZXN0UHl0aG9uKGRhdGE6IERHLkRhdGFGcmFtZSwgeDogREcuQ29sdW1uLCB5OiBERy5Db2x1bW4pOiBQcm9taXNlPG51bWJlcj4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlRUZXN0UHl0aG9uJywgeyBkYXRhLCB4LCB5IH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gcHl0aG9uUGFyYW1zVGVzdChpOiBudW1iZXIsIGQ6IG51bWJlciwgYjogYm9vbGVhbiwgczogc3RyaW5nLCBkdDogYW55LCBkZjogREcuRGF0YUZyYW1lLCBjb2w6IERHLkNvbHVtbiwgY29sczogc3RyaW5nW10pOiBQcm9taXNlPG51bWJlcj4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlB5dGhvblBhcmFtc1Rlc3QnLCB7IGksIGQsIGIsIHMsIGR0LCBkZiwgY29sLCBjb2xzIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdWxtbygpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlVsbW8nLCB7fSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBhbm92YSh0YWJsZTogREcuRGF0YUZyYW1lLCBjYXRlZ29yaWVzOiBERy5Db2x1bW4sIHZhcmlhYmxlOiBERy5Db2x1bW4pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkFOT1ZBJywgeyB0YWJsZSwgY2F0ZWdvcmllcywgdmFyaWFibGUgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBhY2YoZGF0YTogREcuRGF0YUZyYW1lLCBjb2x1bW5zOiBzdHJpbmdbXSk6IFByb21pc2U8YW55PiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6QUNGJywgeyBkYXRhLCBjb2x1bW5zIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gYnNhKGhlaWdodDogbnVtYmVyLCB3ZWlnaHQ6IG51bWJlcik6IFByb21pc2U8bnVtYmVyPiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6QlNBJywgeyBoZWlnaHQsIHdlaWdodCB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGNvbnRvdXJQbG90KGRhdGE6IERHLkRhdGFGcmFtZSwgY29sdW1uczogc3RyaW5nW10pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkNvbnRvdXJQbG90JywgeyBkYXRhLCBjb2x1bW5zIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gZmluYW5jaWFsQ2hhcnQocmF0ZXM6IERHLkRhdGFGcmFtZSwgZGF0ZTogREcuQ29sdW1uLCBvcGVuOiBERy5Db2x1bW4sIGhpZ2g6IERHLkNvbHVtbiwgbG93OiBERy5Db2x1bW4sIGNsb3NlOiBERy5Db2x1bW4sIHZvbHVtZTogREcuQ29sdW1uLCBhZGp1c3RlZDogREcuQ29sdW1uKTogUHJvbWlzZTxhbnk+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpGaW5hbmNpYWxDaGFydCcsIHsgcmF0ZXMsIGRhdGUsIG9wZW4sIGhpZ2gsIGxvdywgY2xvc2UsIHZvbHVtZSwgYWRqdXN0ZWQgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBmaXR0aW5nRFJDKHRhYmxlOiBERy5EYXRhRnJhbWUsIHg6IERHLkNvbHVtbiwgeTogREcuQ29sdW1uKTogUHJvbWlzZTxhbnk+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpGaXR0aW5nRFJDJywgeyB0YWJsZSwgeCwgeSB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGlpckZpbHRlcihkYXRhOiBERy5EYXRhRnJhbWUsIHNpZ25hbDogREcuQ29sdW1uLCBvcmRlcjogbnVtYmVyLCBmcmVxdWVuY3k6IG51bWJlciwgc2FtcGxlUmF0ZTogbnVtYmVyLCB0eXBlOiBzdHJpbmcpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOklJUkZpbHRlcicsIHsgZGF0YSwgc2lnbmFsLCBvcmRlciwgZnJlcXVlbmN5LCBzYW1wbGVSYXRlLCB0eXBlIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24ga25uUihkYXRhOiBERy5EYXRhRnJhbWUsIGNvbHVtbnM6IHN0cmluZ1tdLCBuZWlnaGJvdXJzOiBudW1iZXIpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOktOTlInLCB7IGRhdGEsIGNvbHVtbnMsIG5laWdoYm91cnMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBrc1Rlc3QoZGF0YTogREcuRGF0YUZyYW1lLCB4OiBERy5Db2x1bW4sIHk6IERHLkNvbHVtbik6IFByb21pc2U8bnVtYmVyPiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6S1NUZXN0JywgeyBkYXRhLCB4LCB5IH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gbGRhKHRhYmxlOiBERy5EYXRhRnJhbWUsIHByZWRpY3Q6IERHLkNvbHVtbiwgZmVhdHVyZXM6IHN0cmluZ1tdLCBwZXJjOiBudW1iZXIpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOkxEQScsIHsgdGFibGUsIHByZWRpY3QsIGZlYXR1cmVzLCBwZXJjIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gbGluZWFyUmVncmVzc2lvbihkYXRhOiBERy5EYXRhRnJhbWUsIHg6IERHLkNvbHVtbiwgeTogREcuQ29sdW1uLCBpbnRlcmNlcHRaZXJvOiBib29sZWFuKTogUHJvbWlzZTxudW1iZXI+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpMaW5lYXJSZWdyZXNzaW9uJywgeyBkYXRhLCB4LCB5LCBpbnRlcmNlcHRaZXJvIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gbGlzdFBhY2thZ2VzKCk6IFByb21pc2U8REcuRGF0YUZyYW1lPiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6TGlzdFBhY2thZ2VzJywge30pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gbG1lcih0YWJsZTogREcuRGF0YUZyYW1lLCBmZWF0dXJlczogc3RyaW5nW10sIHJhbmRvbTogREcuQ29sdW1uLCBwcmVkaWN0OiBERy5Db2x1bW4sIHBlcmM6IG51bWJlcik6IFByb21pc2U8REcuRGF0YUZyYW1lPiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6TE1FUicsIHsgdGFibGUsIGZlYXR1cmVzLCByYW5kb20sIHByZWRpY3QsIHBlcmMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBtYW5vdmEodGFibGU6IERHLkRhdGFGcmFtZSwgdmFyaWFibGUxOiBERy5Db2x1bW4sIHZhcmlhYmxlMjogREcuQ29sdW1uLCB2YXJpYWJsZTM6IERHLkNvbHVtbik6IFByb21pc2U8REcuRGF0YUZyYW1lPiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6TUFOT1ZBJywgeyB0YWJsZSwgdmFyaWFibGUxLCB2YXJpYWJsZTIsIHZhcmlhYmxlMyB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHBjYVIoVDogREcuRGF0YUZyYW1lLCBjb2x1bW5zOiBzdHJpbmdbXSwgbnVtQ29tcDogbnVtYmVyLCBjZW50ZXI6IGJvb2xlYW4sIHNjYWxlOiBib29sZWFuKTogUHJvbWlzZTxERy5EYXRhRnJhbWU+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpQQ0FSJywgeyBULCBjb2x1bW5zLCBudW1Db21wLCBjZW50ZXIsIHNjYWxlIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gcGxzKHRhYmxlOiBERy5EYXRhRnJhbWUsIHByZWRpY3Q6IERHLkNvbHVtbiwgZmVhdHVyZXM6IHN0cmluZ1tdLCBjb21wb25lbnRzOiBudW1iZXIpOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlBMUycsIHsgdGFibGUsIHByZWRpY3QsIGZlYXR1cmVzLCBjb21wb25lbnRzIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gcHJlZGljdGl2ZU1vZGVsU1ZNKGRhdGFzZXRQcmVkaWN0OiBERy5EYXRhRnJhbWUsIHByZWRpY3Q6IERHLkNvbHVtbiwgZGF0YXNldDogREcuRGF0YUZyYW1lLCBjb2x1bW5zOiBzdHJpbmdbXSwgZmlsbE1pc3Npbmc6IGJvb2xlYW4sIHBlcmM6IG51bWJlcik6IFByb21pc2U8YW55PiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6UHJlZGljdGl2ZU1vZGVsU1ZNJywgeyBkYXRhc2V0UHJlZGljdCwgcHJlZGljdCwgZGF0YXNldCwgY29sdW1ucywgZmlsbE1pc3NpbmcsIHBlcmMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiByZHVwKHM6IHN0cmluZyk6IFByb21pc2U8c3RyaW5nPiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6UkR1cCcsIHsgcyB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHJlbnZTcGVsbGluZ0V4YW1wbGUoKTogUHJvbWlzZTxERy5EYXRhRnJhbWU+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpSZW52U3BlbGxpbmdFeGFtcGxlJywge30pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gc2NhbG9ncmFtUihkYXRhOiBERy5EYXRhRnJhbWUsIHNpZ25hbDogREcuQ29sdW1uLCBzYW1wbGVSYXRlOiBudW1iZXIsIG9jdGF2ZXM6IG51bWJlciwgdm9pY2VzOiBudW1iZXIsIHJlbW92ZURjOiBib29sZWFuKTogUHJvbWlzZTxhbnk+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpTY2Fsb2dyYW1SJywgeyBkYXRhLCBzaWduYWwsIHNhbXBsZVJhdGUsIG9jdGF2ZXMsIHZvaWNlcywgcmVtb3ZlRGMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBzY2F0dGVyUGxvdFIodDogREcuRGF0YUZyYW1lLCB4Q29sdW1uTmFtZTogREcuQ29sdW1uLCB5Q29sdW1uTmFtZTogREcuQ29sdW1uLCBjb2xvckNvbHVtbk5hbWU6IERHLkNvbHVtbik6IFByb21pc2U8YW55PiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6U2NhdHRlclBsb3RSJywgeyB0LCB4Q29sdW1uTmFtZSwgeUNvbHVtbk5hbWUsIGNvbG9yQ29sdW1uTmFtZSB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHNlbnRpbWVudENsYXNzaWZpY2F0aW9uKGRhdGE6IERHLkRhdGFGcmFtZSwgY29sOiBERy5Db2x1bW4pOiBQcm9taXNlPERHLkRhdGFGcmFtZT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlNlbnRpbWVudENsYXNzaWZpY2F0aW9uJywgeyBkYXRhLCBjb2wgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBzcGVjdHJvZ3JhbShkYXRhOiBERy5EYXRhRnJhbWUsIHNpZ25hbDogREcuQ29sdW1uLCBzYW1wbGVSYXRlOiBudW1iZXIsIHdpbmRvd0xlbmd0aDogbnVtYmVyLCB0aW1lU3RlcDogbnVtYmVyLCByZW1vdmVEYzogYm9vbGVhbik6IFByb21pc2U8YW55PiB7XG4gIHJldHVybiBhd2FpdCBncm9rLmZ1bmN0aW9ucy5jYWxsKCdTYW1wbGVzOlNjcmlwdHM6U3BlY3Ryb2dyYW0nLCB7IGRhdGEsIHNpZ25hbCwgc2FtcGxlUmF0ZSwgd2luZG93TGVuZ3RoLCB0aW1lU3RlcCwgcmVtb3ZlRGMgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBzdXJmYWNlUGxvdCh0OiBERy5EYXRhRnJhbWUsIFg6IERHLkNvbHVtbiwgWTogREcuQ29sdW1uLCBaOiBERy5Db2x1bW4pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlN1cmZhY2VQbG90JywgeyB0LCBYLCBZLCBaIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdFRlc3RSKGRhdGE6IERHLkRhdGFGcmFtZSwgeDogREcuQ29sdW1uLCB5OiBERy5Db2x1bW4pOiBQcm9taXNlPG51bWJlcj4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlRUZXN0UicsIHsgZGF0YSwgeCwgeSB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRlcm5hcnlQbG90KHQ6IERHLkRhdGFGcmFtZSwgdG9wQ29sdW1uTmFtZTogREcuQ29sdW1uLCBsZWZ0Q29sdW1uTmFtZTogREcuQ29sdW1uLCByaWdodENvbHVtbk5hbWU6IERHLkNvbHVtbiwgcG9pbnRTaXplOiBudW1iZXIpOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlRlcm5hcnlQbG90JywgeyB0LCB0b3BDb2x1bW5OYW1lLCBsZWZ0Q29sdW1uTmFtZSwgcmlnaHRDb2x1bW5OYW1lLCBwb2ludFNpemUgfSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiByUGFyYW1zVGVzdChpOiBudW1iZXIsIGQ6IG51bWJlciwgYjogYm9vbGVhbiwgczogc3RyaW5nLCBkdDogYW55LCBkZjogREcuRGF0YUZyYW1lLCBjb2w6IERHLkNvbHVtbiwgY29sczogc3RyaW5nW10pOiBQcm9taXNlPG51bWJlcj4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlJQYXJhbXNUZXN0JywgeyBpLCBkLCBiLCBzLCBkdCwgZGYsIGNvbCwgY29scyB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRpbWVTZXJpZXNEZWNvbXBvc2l0aW9uKGRhdGE6IERHLkRhdGFGcmFtZSwgZGF0ZXM6IERHLkNvbHVtbiwgb2JzZXJ2YXRpb25zOiBERy5Db2x1bW4pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuY2FsbCgnU2FtcGxlczpTY3JpcHRzOlRpbWVTZXJpZXNEZWNvbXBvc2l0aW9uJywgeyBkYXRhLCBkYXRlcywgb2JzZXJ2YXRpb25zIH0pO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gYXJpbWFGb3JlY2FzdGluZyhkYXRhOiBERy5EYXRhRnJhbWUsIGRhdGVzOiBERy5Db2x1bW4sIG9ic2VydmF0aW9uczogREcuQ29sdW1uLCBQOiBudW1iZXIsIEQ6IG51bWJlciwgUTogbnVtYmVyLCBvYnNGb3JlY2FzdDogbnVtYmVyKTogUHJvbWlzZTxhbnk+IHtcbiAgcmV0dXJuIGF3YWl0IGdyb2suZnVuY3Rpb25zLmNhbGwoJ1NhbXBsZXM6U2NyaXB0czpBUklNQUZvcmVjYXN0aW5nJywgeyBkYXRhLCBkYXRlcywgb2JzZXJ2YXRpb25zLCBQLCBELCBRLCBvYnNGb3JlY2FzdCB9KTtcbn1cbiJdfQ==