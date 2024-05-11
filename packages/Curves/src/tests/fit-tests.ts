import {
	getSeriesFitFunction,
	getCurve,
	fitSeries,
	getSeriesConfidenceInterval,
	getSeriesStatistics,
	getPointsArrays,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {sigmoid, FIT_FUNCTION_SIGMOID, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {calculateBoxPlotStatistics} from '@datagrok-libraries/statistics/src/box-plot-statistics';
import {category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';


const sigmoidSeries: IFitSeries = {"fitLineColor": "#1f77b4", "pointColor": "#1f77b4", "showCurveConfidenceInterval": false, "points": [{"x": 0.10000000149011612, "y": 1.7412786483764648}, {"x": 0.6000000238418579, "y": 1.8561450242996216}, {"x": 1.100000023841858, "y": 1.6065685749053955}, {"x": 1.600000023841858,"y": 1.70476496219635}, {"x": 2.0999999046325684, "y": 1.5737264156341553}, {"x": 2.5999999046325684, "y": 1.6007002592086792}, {"x": 3.0999999046325684, "y": 1.6796687841415405}, {"x": 3.5999999046325684, "y": 1.656104326248169}, {"x": 4.099999904632568, "y": 1.782997488975525}, {"x": 4.599999904632568, "y": 1.530208945274353}, {"x": 5.099999904632568, "y": 1.1572397947311401}, {"x": 5.599999904632568, "y": 0.8691964745521545}, {"x": 6.099999904632568, "y": 0.3228665590286255}, {"x": 6.599999904632568,"y": 0.2990703880786896}, {"x": 7.099999904632568, "y": 0.23361243307590485}], "fitFunction": "sigmoid", "clickToToggle": false, "showFitLine": true, "showPoints": "points", "parameters": [1.6914372095641517, 1.1536998642628853, 5.410173358224149, 0.2089689354045083]};
const polynomialSeries: IFitSeries = {"fitLineColor": "#2ca02c", "pointColor": "#2ca02c", "showCurveConfidenceInterval": true, "fitFunction": {"name": "Polynomial", "function": "([p1, p2, p3, p4], x) => p1 * x * x * x + p2 * x * x + p3 * x + p4", "getInitialParameters": "(xs, ys) => [0.1, -1, 4, 4]", "parameterNames": ["Slope", "Intercept", "Parameter3", "Parameter4"]}, "points": [{"x": 0.10000000149011612, "y": 7.365033149719238}, {"x": 0.6000000238418579, "y": 6.595034599304199}, {"x": 1.100000023841858, "y": 7.05179500579834}, {"x": 1.600000023841858, "y": 7.251461982727051}, {"x": 2.0999999046325684, "y": 7.775498867034912}, {"x": 2.5999999046325684, "y": 7.748039722442627}, {"x": 3.0999999046325684, "y": 6.8391804695129395}, {"x": 3.5999999046325684, "y": 7.570991516113281}, {"x": 4.099999904632568, "y": 6.387666702270508}, {"x": 4.599999904632568, "y": 6.464827537536621}, {"x": 5.099999904632568, "y": 3.932436227798462}, {"x": 5.599999904632568, "y": 1.4741199016571045}, {"x": 6.099999904632568, "y": 0.47604307532310486}, {"x": 6.599999904632568, "y": 0.5836760401725769}, {"x": 7.099999904632568, "y": 1.0317600965499878}], "clickToToggle": false, "showFitLine": true, "showPoints": "points", "parameters": [0.07070206940832963, -1.138666933634074, 3.947509109407932, 3.9947960440506685]};
const sigmoidFitFunc = getSeriesFitFunction(sigmoidSeries);
const polynomialFitFunc = getSeriesFitFunction(polynomialSeries);


category('fit', () => {
	test('getSeriesFitFunction', async () => {
		const {xs: sigmoidXs, ys: sigmoidYs} = getPointsArrays(sigmoidSeries.points);
		const {xs: polynomialXs, ys: polynomialYs} = getPointsArrays(polynomialSeries.points);

		expect(sigmoidFitFunc.name, FIT_FUNCTION_SIGMOID);
		expect(polynomialFitFunc.name, 'Polynomial');
		expectArray(sigmoidFitFunc.parameterNames, ['Top', 'Bottom', 'Slope', 'IC50']);
		expectArray(polynomialFitFunc.parameterNames, ['Slope', 'Intercept', 'Parameter3', 'Parameter4']);
		expectArray(sigmoidFitFunc.getInitialParameters(sigmoidXs, sigmoidYs), [1.8561450242996216, 1, 5.099999904632568, 0.23361243307590485]);
		expectArray(polynomialFitFunc.getInitialParameters(polynomialXs, polynomialYs), [0.1, -1, 4, 4]);
		const params = new Float32Array(sigmoidSeries.parameters?.length!);
		params.set(sigmoidSeries.parameters!);
		expect(sigmoidFitFunc.y(params, 1.1), 1.6914214213007113);
		expect(polynomialFitFunc.y(polynomialFitFunc.getInitialParameters(polynomialXs, polynomialYs), 1.1), 7.3231);
	});

	test('getCurve', async () => {
		const sigmoidCurve = getCurve(sigmoidSeries, sigmoidFitFunc);
		const polynomialCurve = getCurve(polynomialSeries, polynomialFitFunc);

		expect(sigmoidCurve(2.5), 1.6907865884456907);
		expect(sigmoidCurve(5.8858), 0.5356637992990292);
		expect(polynomialCurve(1.56), 7.6502648725700535);
		expect(polynomialCurve(3.99876), 6.09328601677405);
	});

	test('fitSeries', async () => {
		const start = Date.now();
		const sigmoidFitSeries = fitSeries(sigmoidSeries, sigmoidFitFunc);
		const polynomialFitSeries = fitSeries(polynomialSeries, polynomialFitFunc);
		const stop = Date.now();

		expect(sigmoidFitSeries.fittedCurve(2.5), 1.690742266220136);
		expect(polynomialFitSeries.fittedCurve(3.99876), 6.09328601677405);
		expectArray(sigmoidFitSeries.parameters, [1.691390784119389, 1.1540819399621425, 5.4104428146970704, 0.20886765646057662]);
		expectArray(polynomialFitSeries.parameters, [0.07070206940832963, -1.138666933634074, 3.947509109407932, 3.9947960440506685]);
		return `${stop.valueOf() - start.valueOf()} ms`;
	});

	test('getSeriesConfidenceInterval', async () => {
		const sigmoidSeriesConfidenceIntervals = getSeriesConfidenceInterval(sigmoidSeries, sigmoidFitFunc, true);
		const polynomialSeriesConfidenceIntervals = getSeriesConfidenceInterval(polynomialSeries, polynomialFitFunc, true);

		expect(sigmoidSeriesConfidenceIntervals.confidenceTop(2.578), 1.8540854906189528);
		expect(polynomialSeriesConfidenceIntervals.confidenceTop(2.3342), 10.035230015978911);
		expect(sigmoidSeriesConfidenceIntervals.confidenceBottom(3.987), 1.494930602732917);
		expect(polynomialSeriesConfidenceIntervals.confidenceBottom(3.8796), 4.168651824724921);
	});

	test('getSeriesStatistics', async () => {
		const sigmoidSeriesStatistics = getSeriesStatistics(sigmoidSeries, sigmoidFitFunc);
		const polynomialSeriesStatistics = getSeriesStatistics(polynomialSeries, polynomialFitFunc);

		expect(sigmoidSeriesStatistics.auc, 9.33542235488813);
		expect(polynomialSeriesStatistics.auc, 36.52554692824738);
		expect(sigmoidSeriesStatistics.rSquared, 0.9781915962446398);
		expect(polynomialSeriesStatistics.rSquared, 0.8460557736673078);
	});

	test('calculateBoxPlotStatistics', async () => {
		const values = [0.7654603719711304, 0.8199243545532227, 0.8257747292518616, 0.9558155536651611, 0.9596694707870483];
		const boxPlotStats = calculateBoxPlotStatistics(values);

		expect(boxPlotStats.q1, 0.8199243545532227);
		expect(boxPlotStats.q2, 0.8257747292518616);
		expect(boxPlotStats.q3, 0.9558155536651611);
		expect(boxPlotStats.lowerAdjacentValue, 0.7654603719711304);
		expect(boxPlotStats.upperAdjacentValue, 0.9558155536651611);
	});

	test('sigmoid', async () => {
		const parameters = [2.612687371569039, -1.421800778574126, 5.16688669906018, 0.341788492203575];
		const params = new Float32Array(4);
		params.set(parameters);
		expect(sigmoid(params, 0.5), 0.34178901799627626);
		expect(sigmoid(params, 2.5), 0.3421552218929382);
		expect(sigmoid(params, 5.678976), 2.254881739737238);
	});
});
