/* eslint-disable max-len */
import {
  getSeriesFitFunction,
  getCurve,
  getSeriesConfidenceInterval,
  getSeriesStatistics,
  getPointsArrays,
  getSeriesFit,
  toDataSpace,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {getDataPoints} from '@datagrok-libraries/statistics/src/fit/fit-points';
import {
  fitSeries,
  fitFunctions,
  getFitFunction,
  getStatistic,
  getStatisticProperty,
  isFit,
  FitTypeMap,
  FitSeries,
  SigmoidFit,
  LinearFit,
  LogLinearFit,
  ExponentialFit,
  FourPLDoseResponseFit,
  JSFunctionFit,
} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {
  sigmoid,
  FIT_FUNCTION_SIGMOID,
  IFitSeries,
  FIT_FUNCTION_LINEAR, FIT_FUNCTION_LOG_LINEAR, FIT_FUNCTION_EXPONENTIAL, FIT_FUNCTION_4PL_DOSE_RESPONSE,
  FIT_FUNCTION_4PL_REGRESSION, linear, logLinear, exponential
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {calculateBoxPlotStatistics} from '@datagrok-libraries/statistics/src/box-plot-statistics';
import {category, test, expect, expectArray, expectFloat} from '@datagrok-libraries/test/src/test';
import {droplineFraction} from '../fit/render-utils';


const sigmoidSeries: IFitSeries = {'fitLineColor': '#1f77b4', 'pointColor': '#1f77b4', 'showCurveConfidenceInterval': false, 'points': [{'x': 0.10000000149011612, 'y': 1.7412786483764648}, {'x': 0.6000000238418579, 'y': 1.8561450242996216}, {'x': 1.100000023841858, 'y': 1.6065685749053955}, {'x': 1.600000023841858, 'y': 1.70476496219635}, {'x': 2.0999999046325684, 'y': 1.5737264156341553}, {'x': 2.5999999046325684, 'y': 1.6007002592086792}, {'x': 3.0999999046325684, 'y': 1.6796687841415405}, {'x': 3.5999999046325684, 'y': 1.656104326248169}, {'x': 4.099999904632568, 'y': 1.782997488975525}, {'x': 4.599999904632568, 'y': 1.530208945274353}, {'x': 5.099999904632568, 'y': 1.1572397947311401}, {'x': 5.599999904632568, 'y': 0.8691964745521545}, {'x': 6.099999904632568, 'y': 0.3228665590286255}, {'x': 6.599999904632568, 'y': 0.2990703880786896}, {'x': 7.099999904632568, 'y': 0.23361243307590485}], 'fitFunction': 'sigmoid', 'clickToToggle': false, 'showFitLine': true, 'showPoints': 'points', 'parameters': [1.6914372095641517, 1.1536998642628853, 5.410173358224149, 0.2089689354045083]};
const linearSeries: IFitSeries = {'fitLineColor': '#1f77b4', 'pointColor': '#1f77b4', 'showCurveConfidenceInterval': false, 'points': [{'x': 0.10000000149011612, 'y': 1.5412786483764648}, {'x': 0.6000000238418579, 'y': 1.9261450242996216}, {'x': 1.100000023841858, 'y': 1.6065685749053955}, {'x': 1.600000023841858, 'y': 1.70476496219635}, {'x': 2.0999999046325684, 'y': 1.5737264156341553}, {'x': 2.5999999046325684, 'y': 1.6007002592086792}, {'x': 3.0999999046325684, 'y': 1.6796687841415405}, {'x': 3.5999999046325684, 'y': 1.656104326248169}, {'x': 4.099999904632568, 'y': 1.782997488975525}, {'x': 4.599999904632568, 'y': 1.530208945274353}, {'x': 5.099999904632568, 'y': 1.1572397947311401}, {'x': 5.599999904632568, 'y': 0.8691964745521545}, {'x': 6.099999904632568, 'y': 0.3228665590286255}, {'x': 6.599999904632568, 'y': 0.2990703880786896}, {'x': 7.099999904632568, 'y': 0.23361243307590485}], 'fitFunction': 'linear', 'clickToToggle': false, 'showFitLine': true, 'showPoints': 'points', 'parameters': [-5.45641324, 10.456486765]};
const logLinearSeries: IFitSeries = {'fitLineColor': '#1f77b4', 'pointColor': '#1f77b4', 'showCurveConfidenceInterval': false, 'points': [{'x': 0.10000000149011612, 'y': 1.789712786483764648}, {'x': 0.6000000238418579, 'y': 1.210161450242996216}, {'x': 1.100000023841858, 'y': 1.7895685749053955}, {'x': 1.600000023841858, 'y': 1.82476496219635}, {'x': 2.0999999046325684, 'y': 1.5737264156341553}, {'x': 2.5999999046325684, 'y': 1.6007002592086792}, {'x': 3.0999999046325684, 'y': 1.6796687841415405}, {'x': 3.5999999046325684, 'y': 1.656104326248169}, {'x': 4.099999904632568, 'y': 1.782997488975525}, {'x': 4.599999904632568, 'y': 1.530208945274353}, {'x': 5.099999904632568, 'y': 1.1572397947311401}, {'x': 5.599999904632568, 'y': 0.8691964745521545}, {'x': 6.099999904632568, 'y': 0.3228665590286255}, {'x': 6.599999904632568, 'y': 0.2990703880786896}, {'x': 7.099999904632568, 'y': 0.23361243307590485}], 'fitFunction': 'log-linear', 'clickToToggle': false, 'showFitLine': true, 'showPoints': 'points', 'parameters': [-8.45641324, 90.4456456486765]};
const exponentialSeries: IFitSeries = {'fitLineColor': '#1f77b4', 'pointColor': '#1f77b4', 'showCurveConfidenceInterval': false, 'points': [{'x': 0.10000000149011612, 'y': 1.789712786483764648}, {'x': 0.6000000238418579, 'y': 1.210161450242996216}, {'x': 1.100000023841858, 'y': 1.7895685749053955}, {'x': 1.600000023841858, 'y': 1.82476496219635}, {'x': 2.0999999046325684, 'y': 1.5737264156341553}, {'x': 2.5999999046325684, 'y': 1.6007002592086792}, {'x': 3.0999999046325684, 'y': 1.6796687841415405}, {'x': 3.5999999046325684, 'y': 1.656104326248169}, {'x': 4.099999904632568, 'y': 1.782997488975525}, {'x': 4.599999904632568, 'y': 1.530208945274353}, {'x': 5.099999904632568, 'y': 1.1572397947311401}, {'x': 5.599999904632568, 'y': 0.8691964745521545}, {'x': 6.099999904632568, 'y': 0.3228665590286255}, {'x': 6.599999904632568, 'y': 0.2990703880786896}, {'x': 7.099999904632568, 'y': 0.23361243307590485}], 'fitFunction': 'exponential', 'clickToToggle': false, 'showFitLine': true, 'showPoints': 'points', 'parameters': [88.45641324, -10.4456456486765]};
const polynomialSeries: IFitSeries = {'fitLineColor': '#2ca02c', 'pointColor': '#2ca02c', 'showCurveConfidenceInterval': true, 'fitFunction': {'name': 'Polynomial', 'function': '([p1, p2, p3, p4], x) => p1 * x * x * x + p2 * x * x + p3 * x + p4', 'getInitialParameters': '(xs, ys) => [0.1, -1, 4, 4]', 'parameterNames': ['Slope', 'Intercept', 'Parameter3', 'Parameter4']}, 'points': [{'x': 0.10000000149011612, 'y': 7.365033149719238}, {'x': 0.6000000238418579, 'y': 6.595034599304199}, {'x': 1.100000023841858, 'y': 7.05179500579834}, {'x': 1.600000023841858, 'y': 7.251461982727051}, {'x': 2.0999999046325684, 'y': 7.775498867034912}, {'x': 2.5999999046325684, 'y': 7.748039722442627}, {'x': 3.0999999046325684, 'y': 6.8391804695129395}, {'x': 3.5999999046325684, 'y': 7.570991516113281}, {'x': 4.099999904632568, 'y': 6.387666702270508}, {'x': 4.599999904632568, 'y': 6.464827537536621}, {'x': 5.099999904632568, 'y': 3.932436227798462}, {'x': 5.599999904632568, 'y': 1.4741199016571045}, {'x': 6.099999904632568, 'y': 0.47604307532310486}, {'x': 6.599999904632568, 'y': 0.5836760401725769}, {'x': 7.099999904632568, 'y': 1.0317600965499878}], 'clickToToggle': false, 'showFitLine': true, 'showPoints': 'points', 'parameters': [0.07070206940832963, -1.138666933634074, 3.947509109407932, 3.9947960440506685]};
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
    expectArray(sigmoidFitFunc.parameterNames, ['Top', 'Slope', 'IC50', 'Bottom']);
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

    expect(sigmoidFitSeries.fittedCurve(2.5), 1.6903875076701773);
    expect(linearFitSeries.fittedCurve(2.5), 1.5348530188202858);
    expect(logLinearFitSeries.fittedCurve(2.5), 1.3696734183366701);
    expect(exponentialFitSeries.fittedCurve(0.00001345), 2.005337264075174);
    expect(polynomialFitSeries.fittedCurve(3.99876), 5.911350894323636);
    expectArray(sigmoidFitSeries.parameters, [1.690911054611206, 1.1880898475646973, 5.402130603790283, 0.22163867950439453]);
    expectArray(linearFitSeries.parameters, [-0.21431826055049896, 2.070648670196533]);
    expectArray(logLinearFitSeries.parameters, [-0.6303421854972839, 2.1593427658081055]);
    expectArray(exponentialFitSeries.parameters, [2.005340814590454, -0.13163800537586212]);
    expectArray(polynomialFitSeries.parameters, [0.06199381500482559, -0.9554166793823242, 2.872941017150879, 5.736424446105957]);
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

  test('fillParams', async () => {
    const sigmoidCurve = fitSeries(sigmoidSeries, sigmoidFitFunc);
    const sigmoidFit = sigmoidFitFunc.fillParams(sigmoidCurve, sigmoidSeries) as SigmoidFit;
    expect(sigmoidFit.name, FIT_FUNCTION_SIGMOID);
    expect(sigmoidFit.top, sigmoidCurve.parameters[0]);
    expect(sigmoidFit.slope, sigmoidCurve.parameters[1]);
    expect(sigmoidFit.ic50, sigmoidCurve.parameters[2]);
    expect(sigmoidFit.bottom, sigmoidCurve.parameters[3]);
    expect(sigmoidFit.interceptY, sigmoidCurve.fittedCurve(sigmoidCurve.parameters[2]));

    const linearCurve = fitSeries(linearSeries, linearFitFunc);
    const linearFit = linearFitFunc.fillParams(linearCurve, linearSeries) as LinearFit;
    expect(linearFit.slope, linearCurve.parameters[0]);
    expect(linearFit.intercept, linearCurve.parameters[1]);

    const logLinearCurve = fitSeries(logLinearSeries, logLinearFitFunc);
    const logLinearFit = logLinearFitFunc.fillParams(logLinearCurve, logLinearSeries) as LogLinearFit;
    expect(logLinearFit.name, FIT_FUNCTION_LOG_LINEAR);
    expect(logLinearFit.slope, logLinearCurve.parameters[0]);
    expect(logLinearFit.intercept, logLinearCurve.parameters[1]);

    const exponentialCurve = fitSeries(exponentialSeries, exponentialFitFunc);
    const exponentialFit = exponentialFitFunc.fillParams(exponentialCurve, exponentialSeries) as ExponentialFit;
    expect(exponentialFit.mantissa, exponentialCurve.parameters[0]);
    expect(exponentialFit.power, exponentialCurve.parameters[1]);
  });

  test('fillParams dose-response', async () => {
    const fitFunc = fitFunctions[FIT_FUNCTION_4PL_DOSE_RESPONSE];
    const curve = fitSeries(sigmoidSeries, fitFunc);
    const fit = fitFunc.fillParams(curve, sigmoidSeries) as FourPLDoseResponseFit;

    expect(fit.name, FIT_FUNCTION_4PL_DOSE_RESPONSE);
    expect(fit.ic50, curve.parameters[2]);
    expect(fit.ec50, curve.parameters[2]);
    expect(fit.top, curve.parameters[0]);
    expect(fit.bottom, curve.parameters[3]);
  });

  test('fillParams custom function', async () => {
    const curve = fitSeries(polynomialSeries, polynomialFitFunc);
    const fit = polynomialFitFunc.fillParams(curve, polynomialSeries) as JSFunctionFit;

    expectArray(fit.parameters, curve.parameters);
    for (let i = 0; i < polynomialFitFunc.parameterNames.length; i++)
      expect(fit[polynomialFitFunc.parameterNames[i]], curve.parameters[i]);
  });

  test('fillParams data points', async () => {
    const curve = fitSeries(sigmoidSeries, sigmoidFitFunc);
    const defaultFit = sigmoidFitFunc.fillParams(curve, sigmoidSeries);
    const explicitFit = sigmoidFitFunc.fillParams(curve, sigmoidSeries,
      getDataPoints(sigmoidSeries, undefined, false));
    const logFit = sigmoidFitFunc.fillParams(curve, sigmoidSeries, undefined, {logX: true, logY: false});

    expect(explicitFit.rSquared, defaultFit.rSquared);
    expect(explicitFit.auc, defaultFit.auc);
    expect(logFit.rSquared === defaultFit.rSquared, false);
  });

  test('statisticsProperties', async () => {
    expectArray(sigmoidFitFunc.statisticsProperties.map((p) => p.name),
      ['rSquared', 'auc', 'top', 'slope', 'ic50', 'bottom', 'interceptY', 'maxY', 'minY', 'pIC50']);
    expectArray(linearFitFunc.statisticsProperties.map((p) => p.name),
      ['rSquared', 'auc', 'slope', 'intercept']);
    expectArray(logLinearFitFunc.statisticsProperties.map((p) => p.name),
      ['rSquared', 'auc', 'slope', 'intercept']);
    expectArray(exponentialFitFunc.statisticsProperties.map((p) => p.name),
      ['rSquared', 'auc', 'mantissa', 'power']);
    expectArray(fitFunctions[FIT_FUNCTION_4PL_DOSE_RESPONSE].statisticsProperties.map((p) => p.name),
      ['rSquared', 'auc', 'top', 'slope', 'ic50', 'bottom', 'interceptY', 'maxY', 'minY', 'pIC50']);
    expectArray(polynomialFitFunc.statisticsProperties.map((p) => p.name),
      ['rSquared', 'auc', 'Slope', 'Intercept', 'Parameter3', 'Parameter4']);
  });

  test('statisticsProperties friendly names', async () => {
    // parameter labels come from the fit function's own parameterNames, so the two cannot drift
    expectArray(sigmoidFitFunc.statisticsProperties.map((p) => p.friendlyName),
      ['R²', 'AUC', 'Top', 'Slope', 'IC50', 'Bottom', 'Y at IC50', 'Max Y', 'Min Y', 'pIC50']);
    expectArray(fitFunctions[FIT_FUNCTION_4PL_REGRESSION].statisticsProperties.map((p) => p.friendlyName),
      ['R²', 'AUC', 'Top', 'Slope', 'EC50', 'Bottom', 'Y at EC50', 'Max Y', 'Min Y']);
    expectArray(fitFunctions[FIT_FUNCTION_4PL_DOSE_RESPONSE].statisticsProperties.map((p) => p.friendlyName),
      ['R²', 'AUC', 'Max', 'Hill', 'IC50', 'Min', 'Y at IC50', 'Max Y', 'Min Y', 'pIC50']);
    expectArray(linearFitFunc.statisticsProperties.map((p) => p.friendlyName),
      ['R²', 'AUC', 'Slope', 'Intercept']);
  });

  test('getStatistic legacy names', async () => {
    // getFitFunction is typed, so these need no casts - the fit type is known at compile time
    const sigmoidFit = getSeriesFit(sigmoidSeries, getFitFunction('sigmoid'));
    expect(getStatistic(sigmoidFit, 'interceptX'), sigmoidFit.ic50);
    // the property panel's "+" button passes the current names, extracted columns from saved projects
    // pass the legacy ones - both have to resolve off the same fit
    expect(getStatistic(sigmoidFit, 'ic50'), sigmoidFit.ic50);
    expect(getStatistic(sigmoidFit, 'maxY'), sigmoidFit.maxY);
    expect(getStatistic(sigmoidFit, 'minY'), sigmoidFit.minY);
    expect(getStatistic(sigmoidFit, 'top'), sigmoidFit.top);
    expect(getStatistic(sigmoidFit, 'bottom'), sigmoidFit.bottom);
    expect(getStatistic(sigmoidFit, 'slope'), sigmoidFit.slope);

    const regressionFit = getSeriesFit({...sigmoidSeries, parameters: undefined}, getFitFunction('4pl-regression'));
    expect(getStatistic(regressionFit, 'interceptX'), regressionFit.ec50);

    const doseResponseFit = getSeriesFit({...sigmoidSeries, parameters: undefined},
      getFitFunction('4pl-dose-response'));
    expect(getStatistic(doseResponseFit, 'interceptX'), doseResponseFit.ic50);

    // legacy names a 2-parameter fit lacks fall back to the slot the pre-typed API read positionally,
    // and the descriptor has to agree - a value with no descriptor renders nothing on the plot
    const linearFit = getSeriesFit(linearSeries, getFitFunction('linear'));
    expect(getStatistic(linearFit, 'top'), linearFit.parameters[0]);
    expect(getStatisticProperty(getFitFunction('linear'), 'top') !== undefined, true);
    expect(getStatisticProperty(getFitFunction('linear'), 'interceptX') === undefined, true);
    expect(getStatistic(linearFit, 'interceptX') === undefined);
    expect(getStatistic(linearFit, 'interceptY') === undefined);
    expect(getStatistic(linearFit, 'bottom') === undefined);
    expect(getStatistic(linearFit, 'slope'), linearFit.slope);

    // the pre-typed API reported interceptY for every fit function, custom ones included
    const jsFit = getSeriesFit(polynomialSeries, getSeriesFitFunction(polynomialSeries));
    expect(typeof getStatistic(jsFit, 'interceptY'), 'number');
    // and a value with no descriptor is skipped by the plot, so both have to resolve
    expect(getStatisticProperty(getSeriesFitFunction(polynomialSeries), 'interceptY') !== undefined, true);

    // never resolves to a non-numeric field of the fit
    expect(getStatistic(sigmoidFit, 'series') === undefined);
    expect(getStatistic(sigmoidFit, 'parameters') === undefined);
    expect(getStatistic(sigmoidFit, 'name') === undefined);
  });

  test('FitSeries options', async () => {
    // this only compiles if FitSeries still exposes every IFitSeriesOptions member
    const series = new FitSeries([{x: 1, y: 2}]);
    series.name = 'a';
    series.fitFunction = FIT_FUNCTION_SIGMOID;
    series.parameters = [1, 2, 3, 4];
    series.parameterBounds = [{min: 0, max: 1}];
    series.markerType = 'circle';
    series.outlierMarkerType = 'outlier';
    series.lineStyle = 'dashed';
    series.pointColor = '#111111';
    series.fitLineColor = '#222222';
    series.confidenceIntervalColor = '#333333';
    series.outlierColor = '#444444';
    series.connectDots = false;
    series.showFitLine = true;
    series.showPoints = 'points';
    series.showOutliers = true;
    series.showCurveConfidenceInterval = false;
    series.errorModel = 'constant';
    series.clickToToggle = true;
    series.labels = {batch: 'A1'};
    series.droplines = ['IC50'];
    series.columnName = 'curve';
    series.auxLegendName = 'aux';

    expect(series.name, 'a');
    expect(series.showPoints, 'points');
    expect(series.droplines![0], 'IC50');
    expect(series.auxLegendName, 'aux');

    // a series stays JSON-serializable - an own `fit` property used to make this circular
    const json = JSON.parse(JSON.stringify(series));
    expect(json.name, 'a');
    expect(json.columnName, 'curve');
  });

  test('maxY and minY', async () => {
    const sigmoidFit = getSeriesFit(sigmoidSeries, getFitFunction('sigmoid'));
    expect(sigmoidFit.maxY, Math.max(sigmoidFit.top, sigmoidFit.bottom));
    expect(sigmoidFit.minY, Math.min(sigmoidFit.top, sigmoidFit.bottom));
    expect(sigmoidFit.maxY >= sigmoidFit.minY, true);
  });

  test('toDataSpace', async () => {
    const logX = {logX: true, logY: false};
    const fitSpace = getSeriesFit(sigmoidSeries, getFitFunction('sigmoid'), undefined, logX);
    const rawIc50 = fitSpace.ic50;
    // toDataSpace mutates in place, so capture before converting
    const untouchedTop = fitSpace.top;
    const dataSpace = toDataSpace(fitSpace, logX);

    // the inflection point is fitted in log space and must be reported as a concentration
    expect(dataSpace.ic50, Math.pow(10, rawIc50));
    // pIC50 is derived from the concentration, so only after the conversion
    expect(dataSpace.pIC50, -Math.log10(Math.pow(10, rawIc50)));
    expect(dataSpace.top, untouchedTop);

    // the y asymptotes are reported in the space they were fitted in. A stored bottom of 0 has no
    // finite logarithm, so the forward map has to skip it while the inverse would still return 1 -
    // the two could never be exact inverses, so neither converts
    const logY = {logX: false, logY: true};
    const yFit = getSeriesFit(sigmoidSeries, getFitFunction('sigmoid'), undefined, logY);
    const rawTop = yFit.top;
    const rawIc50Y = yFit.ic50;
    toDataSpace(yFit, logY);
    expect(yFit.top, rawTop);
    expect(yFit.ic50, rawIc50Y); // x untouched when logX is off

    // no log axes means no conversion at all
    const plain = getSeriesFit(sigmoidSeries, getFitFunction('sigmoid'));
    const plainIc50 = plain.ic50;
    toDataSpace(plain, {logX: false, logY: false});
    expect(plain.ic50, plainIc50);
  });

  test('isFit narrowing', async () => {
    // a runtime-dispatched fit function yields a plain Fit; isFit narrows it to the concrete type
    const fit = getSeriesFit(sigmoidSeries, getSeriesFitFunction(sigmoidSeries));
    expect(isFit(fit, 'sigmoid'));
    expect(isFit(fit, 'linear') === false);
    if (isFit(fit, 'sigmoid'))
      expect(fit.ic50, getStatistic(fit, 'interceptX'));
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

  test('ICxx is the x at that fraction of the curve', async () => {
    // travelled fraction, measured from the low-x asymptote - so the same call means IC90 on a
    // descending curve and EC90 on an ascending one
    const at = (fitFunc: keyof FitTypeMap, params: number[], fraction: number,
      lowX: number, highX: number): number => {
      const f = getFitFunction(fitFunc)!;
      const array = Float32Array.from(params);
      const start = f.y(array, lowX);
      const end = f.y(array, highX);
      return (f.y(array, f.inverse(array, fraction)!) - start) / (end - start);
    };

    for (const fitFunc of ['sigmoid', '4pl-dose-response'] as (keyof FitTypeMap)[]) {
      for (const slope of [1.2, -1.2]) {
        for (const fraction of [0.1, 0.5, 0.9]) {
          expectFloat(at(fitFunc, [100, slope, -6.5, 5], fraction, -1e6, 1e6), fraction, 0.001,
            `${fitFunc}, slope ${slope}, IC${fraction * 100}`);
        }
      }
    }
    // the power parameterization, where x is a concentration rather than its logarithm
    for (const slope of [1.2, -1.2]) {
      for (const fraction of [0.1, 0.5, 0.9])
        expectFloat(at('4pl-regression', [100, slope, 3e-7, 5], fraction, 1e-30, 1e30), fraction, 0.001);
    }

    // IC50 is the inflection point, which is what the droplines drew before there was a family
    const sigmoidFunc = getFitFunction('sigmoid')!;
    expect(sigmoidFunc.inverse(Float32Array.from([100, 1.2, -6.5, 5]), 0.5), -6.5);

    // nothing to travel between, so nothing to name
    for (const fitFunc of ['linear', 'log-linear', 'exponential'] as (keyof FitTypeMap)[])
      expect(getFitFunction(fitFunc)!.inverse(Float32Array.from([1, 2]), 0.9) === undefined, true, fitFunc);
  });

  test('a dropline name asks for a fraction', async () => {
    expect(droplineFraction('IC50'), 0.5);
    expect(droplineFraction('IC90'), 0.9);
    expect(droplineFraction('EC50'), 0.5);
    expect(droplineFraction('ic55'), 0.55);
    expect(droplineFraction('IC12.5'), 0.125);
    // the asymptotes themselves are never reached, and these are not concentrations at all
    for (const name of ['IC0', 'IC100', 'AUC', ''])
      expect(droplineFraction(name) === undefined, true, name);
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

  test('the combined error model fits to fixed values', async () => {
    // noisy on purpose: a perfect fit leaves no residuals for the error model to weigh
    const series = {fitFunction: 'sigmoid', errorModel: 'combined',
      points: [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]
        .map((x, i) => ({x: x, y: (5 + 95 / (1 + Math.pow(10, Math.log10(x) + 6.5))) * (1 + 0.08 * Math.sin(i * 2.4))}))};
    const logOptions = {logX: true, logY: false};
    const fit = fitSeries(series as any, getSeriesFitFunction(series as any), undefined, logOptions);
    expectArray(Array.from(fit.parameters).map((p) => Number(p.toFixed(6))),
      [100.771515, 0.984725, -6.517356, 5.065719]);
  });
});
