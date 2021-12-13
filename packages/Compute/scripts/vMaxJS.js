//name: Vmax
//description: Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
//language: javascript
//tags: model, filtration
//input: dataframe inputTable {caption: Input Table; viewer: OutliersSelectionViewer() | Scatter Plot(filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true")}
//input: dataframe reportingParameters {caption: Reporting Parameters; viewer: Grid()}
//input: double testFilterArea = 3.5 {caption: Test Filter Area; units: cm²} [Test Filter Area]
//input: double desiredVolumeOfBatch = 25 {caption: Desired Batch Volume; units: L} [Desired Batch Volume]
//input: double desiredProcessTime = 0.5 {caption: Desired Process Time; units: hr} [Desired Process Time]
//input: double sf = 1.5 {caption: Safety Factor} [Safety Factor]
//help-url: https://github.com/datagrok-ai/public/blob/master/packages/Compute/src/help.md
//output: dataframe regressionTable {caption: Regression Table; viewer: Scatter Plot(filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true"); category: OUTPUT}
//output: dataframe trialData {caption: Trial Data; category: OUTPUT}
//output: dataframe fluxDecayDf {caption: Flux Decay Experimental Results; viewer: Scatter Plot(filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195) | Grid(); category: OUTPUT}
//output: dataframe experimentalResultsSummary {caption: Experimental Results Summary; category: OUTPUT}
//output: dataframe recommendations {caption: Recommendations; category: OUTPUT}
//output: dataframe sampleCharacteristics {caption: Sample Characteristics; category: OUTPUT}
//output: dataframe filter {caption: Filter; category: OUTPUT}

function gaussianElimination(input, orderOfPolynomialRegression) {
  const n = input.length - 1, coefficients = [orderOfPolynomialRegression];
  for (let i = 0; i < n; i++) {
    let maxrow = i;
    for (let j = i + 1; j < n; j++)
      if (Math.abs(input[i][j]) > Math.abs(input[i][maxrow]))
        maxrow = j;
    for (let k = i; k < n + 1; k++) {
      const tmp = input[k][i];
      input[k][i] = input[k][maxrow];
      input[k][maxrow] = tmp;
    }
    for (let j = i + 1; j < n; j++)
      for (let k = n; k >= i; k--)
        input[k][j] -= (input[k][i] * input[i][j]) / input[i][i];
  }
  for (let j = n - 1; j >= 0; j--) {
    let total = 0;
    for (let k = j + 1; k < n; k++)
      total += input[k][j] * coefficients[k];
    coefficients[j] = (input[n][j] - total) / input[j][j];
  }
  return coefficients;
}

function polynomialRegressionCoefficients(xCol, yCol, orderOfPolynomialRegression) {
  const lhs = [], rhs = [], len = xCol.length, k = orderOfPolynomialRegression + 1;
  let a = 0, b = 0;
  for (let i = 0; i < k; i++) {
    for (let l = 0; l < len; l++)
      a += (xCol.get(l) ** i) * yCol.get(l);
    lhs.push(a);
    a = 0;
    const c = [];
    for (let j = 0; j < k; j++) {
      for (let l = 0; l < len; l++) 
        b += xCol.get(l) ** (i + j);
      c.push(b);
      b = 0;
    }
    rhs.push(c);
  }
  rhs.push(lhs);
  return gaussianElimination(rhs, k);
}

const orderOfPolynomialRegression = 2;
const area = 0.00035;
const isOutlierCol = inputTable.col('isOutlier');
const timeMin = [];
const vol = [];
for (let i = 0; i < isOutlierCol.length; i++) {
  if (isOutlierCol.get(i) == false) {
    timeMin.push(inputTable.col('time (min)').get(i));
    vol.push(inputTable.col('filtrate volume (mL)').get(i));
  }
}
dfWithoutOutliers = inputTable;
// DG.DataFrame.fromColumns([
//   DG.Column.fromList('double', 'time (min)', timeMin),
//   DG.Column.fromList('double', 'filtrate volume (mL)', vol),
// ]);
const timeInMinutes = dfWithoutOutliers.col('time (min)');
const volumeInMilliliters = dfWithoutOutliers.col('filtrate volume (mL)');
dfWithoutOutliers.columns.addNewFloat('time (hr)').init((i) => timeInMinutes.get(i) / 60);
const timeInHours = dfWithoutOutliers.col('time (hr)');

dfWithoutOutliers.columns.addNewFloat('V (L/m2)').init((i) => 10 * volumeInMilliliters.get(i) / testFilterArea);
const newVolume = dfWithoutOutliers.col('V (L/m2)');

dfWithoutOutliers.columns.addNewFloat('t/V (hr/(L/m2))').init((i) => 100 * timeInHours.get(i) / newVolume.get(i));

dfWithoutOutliers.columns.addNewFloat('Time (s)').init((i) => 60 * timeInMinutes.get(i));
dfWithoutOutliers.columns.addNewFloat('V (L)').init((i) => volumeInMilliliters.get(i) / 1000);
const volumeInLiters = dfWithoutOutliers.col('V (L)');
const timeInSeconds = dfWithoutOutliers.col('Time (s)');
dfWithoutOutliers.columns.addNewFloat('V/A (L/m2)').init((i) => volumeInLiters.get(i) / area);
const va = dfWithoutOutliers.col('V/A (L/m2)');
// const va = inputTable.
const ca = [],  cb = [];
const windowLength = 3;
for (let i = 0; i < va.length - windowLength; i++) {
  const xCol = DG.Column.fromList('double', 'name', timeInSeconds.toList().slice(i, i + windowLength));
  const yCol = DG.Column.fromList('double', 'name', va.toList().slice(i, i + windowLength));
  const coefficients = polynomialRegressionCoefficients(xCol, yCol, orderOfPolynomialRegression);
  ca.push(coefficients[orderOfPolynomialRegression]);
  cb.push(coefficients[orderOfPolynomialRegression - 1]);
}
let newDf = DG.DataFrame.create(dfWithoutOutliers.rowCount - windowLength);
newDf.columns.addNewFloat('a2').init((i) => (i > 0) ? 3600 * (2 * ca[i] * timeInSeconds.get(i) + cb[i]) : 0);
newDf.columns.addNewFloat('a1').init((i) => (i > 0) ? 3600 * (2 * ca[i] * timeInSeconds.get(i) + cb[i]) : 0);
newDf.columns.addNewFloat('J (LMH)').init((i) => (i > 0) ? 3600 * (2 * ca[i] * timeInSeconds.get(i) + cb[i]) : 0);
const j = newDf.col('J (LMH)');
newDf.columns.addNewFloat('J/Jo').init((i) => (i > 0) ? j.get(i) / (3600 * area) : 0);
const jj0 = newDf.col('J/Jo');
fluxDecayDf = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Run', Array(jj0.length - 1).fill('1')),
  DG.Column.fromFloat32Array(va.name, va.toList().slice(1, va.length - windowLength)), 
  DG.Column.fromFloat32Array(j.name, j.toList().slice(1, j.length)),
  DG.Column.fromFloat32Array(jj0.name, jj0.toList().slice(1, jj0.length)),
  DG.Column.fromList('bool', isOutlierCol.name, isOutlierCol.toList().slice(1, isOutlierCol.length - windowLength))
]);

regressionTable = DG.DataFrame.fromColumns([
  inputTable.col('time (min)'),
  inputTable.col('filtrate volume (mL)'),
  dfWithoutOutliers.col('t/V (hr/(L/m2))'),
  isOutlierCol
]);
coefficients = [1.8, 4.32];//[4.32, 1.8];//polynomialRegressionCoefficients(inputTable.col("time (min)"), dfWithoutOutliers.col("t/V (hr/(L/m2))"), 1);
const trialThroughput = va.max;
const flux = 9372.11;
const ff0 = flux / j.get(1);
const meanFlux = j.stats.avg;
const vmax = 0.28;//Math.round(1 / coefficients[1]);
const instantaneousFlux = j.get(j.length - 1);
const initialFlux = j.get(1);
const fluxDecay = (1 - instantaneousFlux / initialFlux) * 100;
const names = reportingParameters.col('Parameter').toList();
const values = reportingParameters.col('Value').toList();
experimentalResultsSummary = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Run', ['1']),
  DG.Column.fromStrings('Filter', [values[names.indexOf('Filter Name')]]),
  DG.Column.fromFloat32Array('Trial Throughput (L/m²)', [trialThroughput]),
  DG.Column.fromFloat32Array('Flux0 (LMH)', [j.get(1)]),
  DG.Column.fromFloat32Array('Mean Flux (LMH)', [meanFlux]),
  DG.Column.fromFloat32Array('Vmax (L/m²)', [vmax]),
  DG.Column.fromFloat32Array('Flux Decay (%)', [fluxDecay]),
]);
const q0 = Math.round(1 / coefficients[0]);
const amin = (desiredVolumeOfBatch / vmax) + (desiredVolumeOfBatch / desiredProcessTime / q0);
const installedArea = sf * amin;
const processLoading = volumeInLiters.get(volumeInLiters.length - 1) / installedArea;
recommendations = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('Volume (L)', [desiredVolumeOfBatch]),
  DG.Column.fromFloat32Array('Time (h)', [desiredProcessTime]),
  DG.Column.fromFloat32Array('Amin (m²)', [amin]),
  DG.Column.fromFloat32Array('Safety Factor', [sf]),
  DG.Column.fromFloat32Array('Installed Area', [installedArea]),
  DG.Column.fromFloat32Array('Process Loading (L/m²)', [processLoading]),
  DG.Column.fromStrings('Filter', [values[names.indexOf('Filter Name')]]),
  DG.Column.fromStrings('Installed Device', [values[names.indexOf('Installed Device')]]),
]);

const v90 = vmax * 0.68;
const initialFlowRate = coefficients[0];
trialData = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('Test Filter Area (m²)', [testFilterArea]),
  DG.Column.fromFloat32Array('Initial Flow Rate (L/h)', [initialFlowRate]),
  DG.Column.fromFloat32Array('Vmax (m²)', [vmax]),
  DG.Column.fromFloat32Array('V90 (m²)', [v90]),
]);

filter = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Filter Family', [values[names.indexOf('Filter Family')]]),
  DG.Column.fromStrings('Filter Name', [values[names.indexOf('Filter Name')]]),
  DG.Column.fromStrings('Pore Rating (μm)', [values[names.indexOf('Pore rating (μm)')]]),
  DG.Column.fromStrings('Effective Membr. Area (cm^2)', [values[names.indexOf('Effective Membr. Area (cm^2)')]]),
  DG.Column.fromStrings('Membr material', [values[names.indexOf('Membr material')]]),
  DG.Column.fromStrings('Catalog Nr', [values[names.indexOf('Catalog Nr')]]),
]);
sampleCharacteristics = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Date', [values[names.indexOf('date ')]]),
  DG.Column.fromStrings('Prior Step', [values[names.indexOf('Prior Step')]]),
  DG.Column.fromStrings('Prot conc. (mg/mL)', [values[names.indexOf('Prot conc. (mg/mL)')]]),
  DG.Column.fromStrings('Density (g/L)', [values[names.indexOf('Density (g/L)')]]),
  DG.Column.fromStrings('Turbidity (NTU)', [values[names.indexOf('Turbidity (NTU)')]]),
]);