import {
	getSeriesFitFunction,
	getCurve,
	fitSeries,
	getSeriesConfidenceInterval,
	getSeriesStatistics,
	getPointsArrays,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {
	sigmoid,
	FIT_FUNCTION_SIGMOID,
	IFitSeries,
	FIT_FUNCTION_LINEAR, FIT_FUNCTION_LOG_LINEAR, FIT_FUNCTION_EXPONENTIAL, linear, logLinear, exponential
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {calculateBoxPlotStatistics} from '@datagrok-libraries/statistics/src/box-plot-statistics';
import {category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';


const sigmoidSeries: IFitSeries = {"fitLineColor": "#1f77b4", "pointColor": "#1f77b4", "showCurveConfidenceInterval": false, "points": [{"x": 0.10000000149011612, "y": 1.7412786483764648}, {"x": 0.6000000238418579, "y": 1.8561450242996216}, {"x": 1.100000023841858, "y": 1.6065685749053955}, {"x": 1.600000023841858,"y": 1.70476496219635}, {"x": 2.0999999046325684, "y": 1.5737264156341553}, {"x": 2.5999999046325684, "y": 1.6007002592086792}, {"x": 3.0999999046325684, "y": 1.6796687841415405}, {"x": 3.5999999046325684, "y": 1.656104326248169}, {"x": 4.099999904632568, "y": 1.782997488975525}, {"x": 4.599999904632568, "y": 1.530208945274353}, {"x": 5.099999904632568, "y": 1.1572397947311401}, {"x": 5.599999904632568, "y": 0.8691964745521545}, {"x": 6.099999904632568, "y": 0.3228665590286255}, {"x": 6.599999904632568,"y": 0.2990703880786896}, {"x": 7.099999904632568, "y": 0.23361243307590485}], "fitFunction": "sigmoid", "clickToToggle": false, "showFitLine": true, "showPoints": "points", "parameters": [1.6914372095641517, 1.1536998642628853, 5.410173358224149, 0.2089689354045083]};
const linearSeries: IFitSeries = {"fitLineColor": "#1f77b4", "pointColor": "#1f77b4", "showCurveConfidenceInterval": false, "points": [{"x": 0.10000000149011612, "y": 1.5412786483764648}, {"x": 0.6000000238418579, "y": 1.9261450242996216}, {"x": 1.100000023841858, "y": 1.6065685749053955}, {"x": 1.600000023841858,"y": 1.70476496219635}, {"x": 2.0999999046325684, "y": 1.5737264156341553}, {"x": 2.5999999046325684, "y": 1.6007002592086792}, {"x": 3.0999999046325684, "y": 1.6796687841415405}, {"x": 3.5999999046325684, "y": 1.656104326248169}, {"x": 4.099999904632568, "y": 1.782997488975525}, {"x": 4.599999904632568, "y": 1.530208945274353}, {"x": 5.099999904632568, "y": 1.1572397947311401}, {"x": 5.599999904632568, "y": 0.8691964745521545}, {"x": 6.099999904632568, "y": 0.3228665590286255}, {"x": 6.599999904632568,"y": 0.2990703880786896}, {"x": 7.099999904632568, "y": 0.23361243307590485}], "fitFunction": "linear", "clickToToggle": false, "showFitLine": true, "showPoints": "points", "parameters": [-5.45641324, 10.456486765]};
const logLinearSeries: IFitSeries = {"fitLineColor": "#1f77b4", "pointColor": "#1f77b4", "showCurveConfidenceInterval": false, "points": [{"x": 0.10000000149011612, "y": 1.789712786483764648}, {"x": 0.6000000238418579, "y": 1.210161450242996216}, {"x": 1.100000023841858, "y": 1.7895685749053955}, {"x": 1.600000023841858,"y": 1.82476496219635}, {"x": 2.0999999046325684, "y": 1.5737264156341553}, {"x": 2.5999999046325684, "y": 1.6007002592086792}, {"x": 3.0999999046325684, "y": 1.6796687841415405}, {"x": 3.5999999046325684, "y": 1.656104326248169}, {"x": 4.099999904632568, "y": 1.782997488975525}, {"x": 4.599999904632568, "y": 1.530208945274353}, {"x": 5.099999904632568, "y": 1.1572397947311401}, {"x": 5.599999904632568, "y": 0.8691964745521545}, {"x": 6.099999904632568, "y": 0.3228665590286255}, {"x": 6.599999904632568,"y": 0.2990703880786896}, {"x": 7.099999904632568, "y": 0.23361243307590485}], "fitFunction": "log-linear", "clickToToggle": false, "showFitLine": true, "showPoints": "points", "parameters": [-8.45641324, 90.4456456486765]};
const exponentialSeries: IFitSeries = {"fitLineColor": "#1f77b4", "pointColor": "#1f77b4", "showCurveConfidenceInterval": false, "points": [{"x": 0.10000000149011612, "y": 1.789712786483764648}, {"x": 0.6000000238418579, "y": 1.210161450242996216}, {"x": 1.100000023841858, "y": 1.7895685749053955}, {"x": 1.600000023841858,"y": 1.82476496219635}, {"x": 2.0999999046325684, "y": 1.5737264156341553}, {"x": 2.5999999046325684, "y": 1.6007002592086792}, {"x": 3.0999999046325684, "y": 1.6796687841415405}, {"x": 3.5999999046325684, "y": 1.656104326248169}, {"x": 4.099999904632568, "y": 1.782997488975525}, {"x": 4.599999904632568, "y": 1.530208945274353}, {"x": 5.099999904632568, "y": 1.1572397947311401}, {"x": 5.599999904632568, "y": 0.8691964745521545}, {"x": 6.099999904632568, "y": 0.3228665590286255}, {"x": 6.599999904632568,"y": 0.2990703880786896}, {"x": 7.099999904632568, "y": 0.23361243307590485}], "fitFunction": "exponential", "clickToToggle": false, "showFitLine": true, "showPoints": "points", "parameters": [88.45641324, -10.4456456486765]};
const polynomialSeries: IFitSeries = {"fitLineColor": "#2ca02c", "pointColor": "#2ca02c", "showCurveConfidenceInterval": true, "fitFunction": {"name": "Polynomial", "function": "([p1, p2, p3, p4], x) => p1 * x * x * x + p2 * x * x + p3 * x + p4", "getInitialParameters": "(xs, ys) => [0.1, -1, 4, 4]", "parameterNames": ["Slope", "Intercept", "Parameter3", "Parameter4"]}, "points": [{"x": 0.10000000149011612, "y": 7.365033149719238}, {"x": 0.6000000238418579, "y": 6.595034599304199}, {"x": 1.100000023841858, "y": 7.05179500579834}, {"x": 1.600000023841858, "y": 7.251461982727051}, {"x": 2.0999999046325684, "y": 7.775498867034912}, {"x": 2.5999999046325684, "y": 7.748039722442627}, {"x": 3.0999999046325684, "y": 6.8391804695129395}, {"x": 3.5999999046325684, "y": 7.570991516113281}, {"x": 4.099999904632568, "y": 6.387666702270508}, {"x": 4.599999904632568, "y": 6.464827537536621}, {"x": 5.099999904632568, "y": 3.932436227798462}, {"x": 5.599999904632568, "y": 1.4741199016571045}, {"x": 6.099999904632568, "y": 0.47604307532310486}, {"x": 6.599999904632568, "y": 0.5836760401725769}, {"x": 7.099999904632568, "y": 1.0317600965499878}], "clickToToggle": false, "showFitLine": true, "showPoints": "points", "parameters": [0.07070206940832963, -1.138666933634074, 3.947509109407932, 3.9947960440506685]};
const sigmoidFitFunc = getSeriesFitFunction(sigmoidSeries);
const linearFitFunc = getSeriesFitFunction(linearSeries);
const logLinearFitFunc = getSeriesFitFunction(logLinearSeries);
const exponentialFitFunc = getSeriesFitFunction(exponentialSeries);
const polynomialFitFunc = getSeriesFitFunction(polynomialSeries);


category('fit', () => {
	test('getSeriesFitFunction', async () => {
		const {xs: sigmoidXs, ys: sigmoidYs} = getPointsArrays(sigmoidSeries.points);
		const {xs: linearXs, ys: linearYs} = getPointsArrays(linearSeries.points);
		const {xs: logLinearXs, ys: logLinearYs} = getPointsArrays(logLinearSeries.points);
		const {xs: exponentialXs, ys: exponentialYs} = getPointsArrays(exponentialSeries.points);
		const {xs: polynomialXs, ys: polynomialYs} = getPointsArrays(polynomialSeries.points);

		expect(sigmoidFitFunc.name, FIT_FUNCTION_SIGMOID);
		expect(linearFitFunc.name, FIT_FUNCTION_LINEAR);
		expect(logLinearFitFunc.name, FIT_FUNCTION_LOG_LINEAR);
		expect(exponentialFitFunc.name, FIT_FUNCTION_EXPONENTIAL);
		expect(polynomialFitFunc.name, 'Polynomial');
		expectArray(sigmoidFitFunc.parameterNames, ['Top', 'Bottom', 'Slope', 'IC50']);
		expectArray(linearFitFunc.parameterNames, ['Slope', 'Intercept']);
		expectArray(logLinearFitFunc.parameterNames, ['Slope', 'Intercept']);
		expectArray(exponentialFitFunc.parameterNames, ['Mantissa', 'Power']);
		expectArray(polynomialFitFunc.parameterNames, ['Slope', 'Intercept', 'Parameter3', 'Parameter4']);
		expectArray(sigmoidFitFunc.getInitialParameters(sigmoidXs, sigmoidYs), [1.8561450242996216, 1, 5.099999904632568, 0.23361243307590485]);
		expectArray(linearFitFunc.getInitialParameters(linearXs, linearYs), [-0.1868094652891159, 1.5599596500396729]);
		expectArray(logLinearFitFunc.getInitialParameters(logLinearXs, logLinearYs), [-5, 100]);
		expectArray(exponentialFitFunc.getInitialParameters(exponentialXs, exponentialYs), [100, -2]);
		expectArray(polynomialFitFunc.getInitialParameters(polynomialXs, polynomialYs), [0.1, -1, 4, 4]);
		const params = new Float32Array(sigmoidSeries.parameters?.length!);
		params.set(sigmoidSeries.parameters!);
		expect(sigmoidFitFunc.y(params, 1.1), 1.6914214561555851);
		expect(linearFitFunc.y(linearFitFunc.getInitialParameters(linearXs, linearYs), 1.1), 1.3544692382216454);
		expect(logLinearFitFunc.y(logLinearFitFunc.getInitialParameters(logLinearXs, logLinearYs), 1.1), 96.29031327635312);
		expect(exponentialFitFunc.y(exponentialFitFunc.getInitialParameters(exponentialXs, exponentialYs), 1.1), 11.080315836233387);
		expect(polynomialFitFunc.y(polynomialFitFunc.getInitialParameters(polynomialXs, polynomialYs), 1.1), 7.3231);
	});

	test('getCurve', async () => {
		const sigmoidCurve = getCurve(sigmoidSeries, sigmoidFitFunc);
		const linearCurve = getCurve(linearSeries, linearFitFunc);
		const logLinearCurve = getCurve(logLinearSeries, logLinearFitFunc);
		const exponentialCurve = getCurve(exponentialSeries, exponentialFitFunc);
		const polynomialCurve = getCurve(polynomialSeries, polynomialFitFunc);

		expect(sigmoidCurve(2.5), 1.690786623428712);
		expect(sigmoidCurve(5.8858), 0.5356638447565519);
		expect(linearCurve(2.5), -3.18454647064209);
		expect(linearCurve(5.8858), -21.658870516967774);
		expect(logLinearCurve(2.5), 79.85176680360948);
		expect(logLinearCurve(5.8858), 74.12932588432805);
		expect(exponentialCurve(0.01), 79.68277758636006);
		expect(exponentialCurve(0.00001345), 88.4439865528934);
		expect(polynomialCurve(1.56), 7.650264639782428);
		expect(polynomialCurve(3.99876), 6.093284856699922);
	});

	test('fitSeries', async () => {
		const start = Date.now();
		const sigmoidFitSeries = fitSeries(sigmoidSeries, sigmoidFitFunc);
		const linearFitSeries = fitSeries(linearSeries, linearFitFunc);
		const logLinearFitSeries = fitSeries(logLinearSeries, logLinearFitFunc);
		const exponentialFitSeries = fitSeries(exponentialSeries, exponentialFitFunc);
		const polynomialFitSeries = fitSeries(polynomialSeries, polynomialFitFunc);
		const stop = Date.now();

		expect(sigmoidFitSeries.fittedCurve(2.5), 1.7048618582775712);
		expect(linearFitSeries.fittedCurve(2.5), 1.5337481647729874);
		expect(logLinearFitSeries.fittedCurve(2.5), 1.4879121793209378);
		expect(exponentialFitSeries.fittedCurve(0.00001345), 2.739771889293024);
		expect(polynomialFitSeries.fittedCurve(3.99876), 6.449473383870631);
		expectArray(sigmoidFitSeries.parameters, [1.7049071788787842, 1.569833517074585, 5.36868143081665, 0.2605934143066406]);
		expectArray(linearFitSeries.parameters, [-0.21358928084373474, 2.067721366882324]);
		expectArray(logLinearFitSeries.parameters, [-1.1107122898101807, 2.879371404647827]);
		expectArray(exponentialFitSeries.parameters, [2.739849090576172, -2.0949888229370117]);
		expectArray(polynomialFitSeries.parameters, [0.05000000074505806, -0.887523889541626, 2.861894130706787, 6]);
		return `${stop.valueOf() - start.valueOf()} ms`;
	});

	test('fitSeries sigmoid benchmark', async () => {
		const start = Date.now();
		fitSeries(sigmoidSeries, sigmoidFitFunc);
		const stop = Date.now();
		return `${stop.valueOf() - start.valueOf()} ms`;
	}, {benchmark: true});

	test('fitSeries linear benchmark', async () => {
		const start = Date.now();
		fitSeries(linearSeries, linearFitFunc);
		const stop = Date.now();
		return `${stop.valueOf() - start.valueOf()} ms`;
	}, {benchmark: true});

	test('fitSeries log-linear benchmark', async () => {
		const start = Date.now();
		fitSeries(logLinearSeries, logLinearFitFunc);
		const stop = Date.now();
		return `${stop.valueOf() - start.valueOf()} ms`;
	}, {benchmark: true});

	test('fitSeries exponential benchmark', async () => {
		const start = Date.now();
		fitSeries(exponentialSeries, exponentialFitFunc);
		const stop = Date.now();
		return `${stop.valueOf() - start.valueOf()} ms`;
	}, {benchmark: true});

	test('fitSeries polynomial benchmark', async () => {
		const start = Date.now();
		fitSeries(polynomialSeries, polynomialFitFunc);
		const stop = Date.now();
		return `${stop.valueOf() - start.valueOf()} ms`;
	}, {benchmark: true});

	test('getSeriesConfidenceInterval', async () => {
		const sigmoidSeriesConfidenceIntervals = getSeriesConfidenceInterval(sigmoidSeries, sigmoidFitFunc, true);
		const linearSeriesConfidenceIntervals = getSeriesConfidenceInterval(linearSeries, linearFitFunc, true);
		const logLinearSeriesConfidenceIntervals = getSeriesConfidenceInterval(logLinearSeries, logLinearFitFunc, true);
		const exponentialSeriesConfidenceIntervals = getSeriesConfidenceInterval(exponentialSeries, exponentialFitFunc, true);
		const polynomialSeriesConfidenceIntervals = getSeriesConfidenceInterval(polynomialSeries, polynomialFitFunc, true);

		expect(sigmoidSeriesConfidenceIntervals.confidenceTop(2.578), 1.8540855269320198);
		expect(linearSeriesConfidenceIntervals.confidenceTop(2.578), 26.645777275332037);
		expect(logLinearSeriesConfidenceIntervals.confidenceTop(2.578), 231.75768109324696);
		expect(exponentialSeriesConfidenceIntervals.confidenceTop(0.00001345), 103.51168317113999);
		expect(polynomialSeriesConfidenceIntervals.confidenceTop(2.3342), 10.03522960294218);
		expect(sigmoidSeriesConfidenceIntervals.confidenceBottom(3.987), 1.4949306416468298);
		expect(linearSeriesConfidenceIntervals.confidenceBottom(3.987), -41.55415698266847);
		expect(logLinearSeriesConfidenceIntervals.confidenceBottom(3.987), -75.23471085303483);
		expect(exponentialSeriesConfidenceIntervals.confidenceBottom(0.01), 64.61508096811347);
		expect(polynomialSeriesConfidenceIntervals.confidenceBottom(3.8796), 4.1686506915250385);
	});

	test('getSeriesStatistics', async () => {
		const sigmoidSeriesStatistics = getSeriesStatistics(sigmoidSeries, sigmoidFitFunc);
		// TODO: statistics tests on linear, log linear and exponential series
		const polynomialSeriesStatistics = getSeriesStatistics(polynomialSeries, polynomialFitFunc);

		expect(sigmoidSeriesStatistics.auc, 9.335422628533761);
		expect(polynomialSeriesStatistics.auc, 36.525538394763096);
		expect(sigmoidSeriesStatistics.rSquared, 0.9781915962461156);
		expect(polynomialSeriesStatistics.rSquared, 0.846055768956569);
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
		expect(sigmoid(params, 0.5), 0.341789026340178);
		expect(sigmoid(params, 2.5), 0.34215523020636207);
		expect(sigmoid(params, 5.678976), 2.2548815999597016);
	});

	test('linear', async () => {
		const parameters = [-5.12346546, 6.24654654];
		const params = new Float32Array(2);
		params.set(parameters);
		expect(linear(params, 0.5), 3.684813976287842);
		expect(linear(params, 2.5), -6.562117099761963);
		expect(linear(params, 5.678976), -22.849491081970214);
	});

	test('log-linear', async () => {
		const parameters = [-5.12346546, 91.5465454654];
		const params = new Float32Array(2);
		params.set(parameters);
		expect(logLinear(params, 0.5), 89.46916042777144);
		expect(logLinear(params, 2.5), 85.12805903963536);
		expect(logLinear(params, 5.678976), 81.81726682791052);
	});

	test('exponential', async () => {
		const parameters = [95.12346546, -11.5465454654];
		const params = new Float32Array(2);
		params.set(parameters);
		expect(exponential(params, 0.0005), 94.57587501550823);
		expect(exponential(params, 0.5), 0.2957925783034655);
		expect(exponential(params, 0.98746), 0.0010630898492656174);
	});
});
