import * as DG from 'datagrok-api/dg';
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


const sigmoidSeries: IFitSeries = {"fitLineColor": "#1f77b4", "pointColor": "#1f77b4", "showCurveConfidenceInterval": false, "points": [{"x": 0.10000000149011612, "y": 1.7412786483764648}, {"x": 0.6000000238418579, "y": 1.8561450242996216}, {"x": 1.100000023841858, "y": 1.6065685749053955}, {"x": 1.600000023841858,"y": 1.70476496219635}, {"x": 2.0999999046325684, "y": 1.5737264156341553}, {"x": 2.5999999046325684, "y": 1.6007002592086792}, {"x": 3.0999999046325684, "y": 1.6796687841415405}, {"x": 3.5999999046325684, "y": 1.656104326248169}, {"x": 4.099999904632568, "y": 1.782997488975525}, {"x": 4.599999904632568, "y": 1.530208945274353}, {"x": 5.099999904632568, "y": 1.1572397947311401}, {"x": 5.599999904632568, "y": 0.8691964745521545}, {"x": 6.099999904632568, "y": 0.3228665590286255}, {"x": 6.599999904632568,"y": 0.2990703880786896}, {"x": 7.099999904632568, "y": 0.23361243307590485}], "fitFunction": "sigmoid", "clickToToggle": false, "showFitLine": true, "showPoints": "points", "showBoxPlot": true, "parameters": [1.6914372095641517, 1.1536998642628853, 5.410173358224149, 0.2089689354045083]};
const polynomialSeries: IFitSeries = {"fitLineColor": "#2ca02c", "pointColor": "#2ca02c", "showCurveConfidenceInterval": true, "fitFunction": {"name": "Polynomial", "function": "([p1, p2, p3, p4], x) => p1 * x * x * x + p2 * x * x + p3 * x + p4", "getInitialParameters": "(xs, ys) => [0.1, -1, 4, 4]", "parameterNames": ["Slope", "Intercept", "Parameter3", "Parameter4"]}, "points": [{"x": 0.10000000149011612, "y": 7.365033149719238}, {"x": 0.6000000238418579, "y": 6.595034599304199}, {"x": 1.100000023841858, "y": 7.05179500579834}, {"x": 1.600000023841858, "y": 7.251461982727051}, {"x": 2.0999999046325684, "y": 7.775498867034912}, {"x": 2.5999999046325684, "y": 7.748039722442627}, {"x": 3.0999999046325684, "y": 6.8391804695129395}, {"x": 3.5999999046325684, "y": 7.570991516113281}, {"x": 4.099999904632568, "y": 6.387666702270508}, {"x": 4.599999904632568, "y": 6.464827537536621}, {"x": 5.099999904632568, "y": 3.932436227798462}, {"x": 5.599999904632568, "y": 1.4741199016571045}, {"x": 6.099999904632568, "y": 0.47604307532310486}, {"x": 6.599999904632568, "y": 0.5836760401725769}, {"x": 7.099999904632568, "y": 1.0317600965499878}], "clickToToggle": false, "showFitLine": true, "showPoints": "points", "showBoxPlot": true, "parameters": [0.07070206940832963, -1.138666933634074, 3.947509109407932, 3.9947960440506685]};
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
		expectArray(sigmoidFitFunc.getInitialParameters(sigmoidXs, sigmoidYs), [1.8561450242996216, 1.2, 5.099999904632568, 0.23361243307590485]);
		expectArray(polynomialFitFunc.getInitialParameters(polynomialXs, polynomialYs), [0.1, -1, 4, 4]);
		expect(sigmoidFitFunc.y(sigmoidSeries.parameters!, 1.1), 1.6914214213007113);
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
		let sigmoidSeriesTest = sigmoidSeries, sigmoidFitFuncTest = sigmoidFitFunc;
		let polynomialSeriesTest = polynomialSeries, polynomialFitFuncTest = polynomialFitFunc;

		if (DG.Test.isInBenchmark) {
			sigmoidSeriesTest = {"fitLineColor":"#2ca02c","pointColor":"#2ca02c","showCurveConfidenceInterval":true,"fitFunction": "sigmoid", "clickToToggle": false, "showFitLine": true, "showPoints": "points", "showBoxPlot": true,"points":[{"x":0.10000000149011612,"y":1.8080450296401978},{"x":0.30000001192092896,"y":1.7073274850845337},{"x":0.5,"y":1.8684450387954712},{"x":0.699999988079071,"y":1.8083453178405762},{"x":0.8999999761581421,"y":1.871852993965149},{"x":1.100000023841858,"y":2.0314087867736816},{"x":1.2999999523162842,"y":1.8348426818847656},{"x":1.5,"y":1.716798186302185},{"x":1.7000000476837158,"y":1.9205678701400757},{"x":1.899999976158142,"y":2.040295124053955},{"x":2.0999999046325684,"y":1.915327787399292},{"x":2.299999952316284,"y":1.7997850179672241},{"x":2.5,"y":1.8200938701629639},{"x":2.700000047683716,"y":1.8472325801849365},{"x":2.9000000953674316,"y":1.892165184020996},{"x":3.0999999046325684,"y":1.6691789627075195},{"x":3.299999952316284,"y":1.7463487386703491},{"x":3.5,"y":1.856266975402832},{"x":3.700000047683716,"y":1.8447580337524414},{"x":3.9000000953674316,"y":2.044322967529297},{"x":4.099999904632568,"y":1.946313500404358},{"x":4.300000190734863,"y":1.7538775205612183},{"x":4.5,"y":1.7557778358459473},{"x":4.699999809265137,"y":1.743997573852539},{"x":4.900000095367432,"y":1.8651076555252075},{"x":5.099999904632568,"y":1.9675723314285278},{"x":5.300000190734863,"y":1.9454524517059326},{"x":5.5,"y":1.8253188133239746},{"x":5.699999809265137,"y":2.0103540420532227},{"x":5.900000095367432,"y":2.0355498790740967},{"x":6.099999904632568,"y":1.9354662895202637},{"x":6.300000190734863,"y":1.7890597581863403},{"x":6.5,"y":2.1066558361053467},{"x":6.699999809265137,"y":2.1674845218658447},{"x":6.900000095367432,"y":1.8644132614135742},{"x":7.099999904632568,"y":2.1932010650634766},{"x":7.300000190734863,"y":1.9899929761886597},{"x":7.5,"y":2.0531647205352783},{"x":7.699999809265137,"y":2.2773845195770264},{"x":7.900000095367432,"y":1.98076331615448},{"x":8.100000381469727,"y":2.2184560298919678},{"x":8.300000190734863,"y":2.323612689971924},{"x":8.5,"y":2.2961959838867188},{"x":8.699999809265137,"y":2.2926039695739746},{"x":8.899999618530273,"y":2.043086051940918},{"x":9.100000381469727,"y":2.1240429878234863},{"x":9.300000190734863,"y":1.9789459705352783},{"x":9.5,"y":2.192298173904419},{"x":9.699999809265137,"y":2.0946543216705322},{"x":0,"y":0.11970949918031693}]};
			sigmoidFitFuncTest = getSeriesFitFunction(sigmoidSeriesTest);
			polynomialSeriesTest = {"fitLineColor":"#2ca02c","pointColor":"#2ca02c","showCurveConfidenceInterval":true,"fitFunction": {"name": "Polynomial", "function": "([p1, p2, p3, p4], x) => p1 * x * x * x + p2 * x * x + p3 * x + p4", "getInitialParameters": "(xs, ys) => [0.1, -1, 4, 4]", "parameterNames": ["Slope", "Intercept", "Parameter3", "Parameter4"]}, "clickToToggle": false, "showFitLine": true, "showPoints": "points", "showBoxPlot": true,"points":[{"x":0.10000000149011612,"y":1.8080450296401978},{"x":0.30000001192092896,"y":1.7073274850845337},{"x":0.5,"y":1.8684450387954712},{"x":0.699999988079071,"y":1.8083453178405762},{"x":0.8999999761581421,"y":1.871852993965149},{"x":1.100000023841858,"y":2.0314087867736816},{"x":1.2999999523162842,"y":1.8348426818847656},{"x":1.5,"y":1.716798186302185},{"x":1.7000000476837158,"y":1.9205678701400757},{"x":1.899999976158142,"y":2.040295124053955},{"x":2.0999999046325684,"y":1.915327787399292},{"x":2.299999952316284,"y":1.7997850179672241},{"x":2.5,"y":1.8200938701629639},{"x":2.700000047683716,"y":1.8472325801849365},{"x":2.9000000953674316,"y":1.892165184020996},{"x":3.0999999046325684,"y":1.6691789627075195},{"x":3.299999952316284,"y":1.7463487386703491},{"x":3.5,"y":1.856266975402832},{"x":3.700000047683716,"y":1.8447580337524414},{"x":3.9000000953674316,"y":2.044322967529297},{"x":4.099999904632568,"y":1.946313500404358},{"x":4.300000190734863,"y":1.7538775205612183},{"x":4.5,"y":1.7557778358459473},{"x":4.699999809265137,"y":1.743997573852539},{"x":4.900000095367432,"y":1.8651076555252075},{"x":5.099999904632568,"y":1.9675723314285278},{"x":5.300000190734863,"y":1.9454524517059326},{"x":5.5,"y":1.8253188133239746},{"x":5.699999809265137,"y":2.0103540420532227},{"x":5.900000095367432,"y":2.0355498790740967},{"x":6.099999904632568,"y":1.9354662895202637},{"x":6.300000190734863,"y":1.7890597581863403},{"x":6.5,"y":2.1066558361053467},{"x":6.699999809265137,"y":2.1674845218658447},{"x":6.900000095367432,"y":1.8644132614135742},{"x":7.099999904632568,"y":2.1932010650634766},{"x":7.300000190734863,"y":1.9899929761886597},{"x":7.5,"y":2.0531647205352783},{"x":7.699999809265137,"y":2.2773845195770264},{"x":7.900000095367432,"y":1.98076331615448},{"x":8.100000381469727,"y":2.2184560298919678},{"x":8.300000190734863,"y":2.323612689971924},{"x":8.5,"y":2.2961959838867188},{"x":8.699999809265137,"y":2.2926039695739746},{"x":8.899999618530273,"y":2.043086051940918},{"x":9.100000381469727,"y":2.1240429878234863},{"x":9.300000190734863,"y":1.9789459705352783},{"x":9.5,"y":2.192298173904419},{"x":9.699999809265137,"y":2.0946543216705322},{"x":0,"y":0.11970949918031693}]};
			polynomialFitFuncTest = getSeriesFitFunction(polynomialSeriesTest);
		}

		const sigmoidFitSeries = fitSeries(sigmoidSeriesTest, sigmoidFitFuncTest);
		const polynomialFitSeries = fitSeries(polynomialSeriesTest, polynomialFitFuncTest);

		expect(sigmoidFitSeries.fittedCurve(2.5), DG.Test.isInBenchmark ? 1.779457823734546 : 1.6907865884456907);
		expect(polynomialFitSeries.fittedCurve(3.99876), DG.Test.isInBenchmark ? 4.166024343916158 : 6.09328601677405);
		expectArray(sigmoidFitSeries.parameters, DG.Test.isInBenchmark ? [0.10961644628818075, 0.028443471965179174, 5.226265783814711, 3.775729673923325] :
			[1.6914372095641517, 1.1536998642628853, 5.410173358224149, 0.2089689354045083]);
		expectArray(polynomialFitSeries.parameters, DG.Test.isInBenchmark ? [0.10211958221951015, -1.3588606666390697, 3.8567403757829686, 3.942561068097394] :
			[0.07070206940832963, -1.138666933634074, 3.947509109407932, 3.9947960440506685]);
	});

	test('getSeriesConfidenceInterval', async () => {
		const sigmoidSeriesConfidenceIntervals = getSeriesConfidenceInterval(sigmoidSeries, sigmoidFitFunc, true);
		const polynomialSeriesConfidenceIntervals = getSeriesConfidenceInterval(polynomialSeries, polynomialFitFunc, true);

		expect(sigmoidSeriesConfidenceIntervals.confidenceTop(2.578), 1.7380288296703128);
		expect(polynomialSeriesConfidenceIntervals.confidenceTop(2.3342), 8.522118220066847);
		expect(sigmoidSeriesConfidenceIntervals.confidenceBottom(3.987), 1.610987263681557);
		expect(polynomialSeriesConfidenceIntervals.confidenceBottom(3.8796), 5.681763620636985);
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
		expect(sigmoid(parameters, 0.5), 0.34178901799627626);
		expect(sigmoid(parameters, 2.5), 0.3421552218929382);
		expect(sigmoid(parameters, 5.678976), 2.254881739737238);
	});
});
