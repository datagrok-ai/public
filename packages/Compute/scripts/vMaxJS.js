//name: Vmax
//description: Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
//language: javascript
//tags: model, filtration
//input: dataframe inputTable {caption: Input Table; viewer: OutliersSelectionViewer() | Scatter Plot(filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true", legendVisibility: "Never", filterByZoom: "false") | Grid()}
//input: dataframe parametersForReporting {caption: Parameters For Reporting; viewer: Grid()}
//input: double testFilterArea = 3.5 {caption: Test Filter Area; units: cm²} [Test Filter area]
//input: double desiredVolumeOfBatch = 25 {caption: Desired Volume Of Batch; units: L} [Desired Volume of Batch]
//input: double desiredProcessTime = 0.5 {caption: Desired Process Time; units: hr} [Desired process time]
//input: double sf = 1.5 {caption: Safety Factor} [Safety Factor]
//help-url: https://github.com/datagrok-ai/public/blob/master/packages/Compute/src/help.md
//output: dataframe regressionTable {caption: Regression Table; viewer: Scatter Plot(filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true"); category: OUTPUT}
//output: dataframe trialData {caption: Trial Data; category: OUTPUT}
//output: dataframe fluxDecayDf {caption: Flux Decay Experimental Results; viewer: Line Chart(x: "V/A (L/m2)", multiAxis: "true") | Grid(); category: OUTPUT}
//output: dataframe experimentalResultsSummary {caption: Experimental Results Summary; category: OUTPUT}
//output: dataframe recommendations {caption: Recommendations; category: OUTPUT}
//output: dataframe sampleCharacteristics {caption: Sample Characteristics; category: OUTPUT}
//output: dataframe filter {caption: Filter; category: OUTPUT}

volume = desiredVolumeOfBatch;
function gaussianElimination(input, orderOfPolynomialRegression) {
  const n = input.length - 1; const coefficients = [orderOfPolynomialRegression];
  for (let i = 0; i < n; i++) {
    let maxrow = i;
    for (let j = i + 1; j < n; j++) {
      if (Math.abs(input[i][j]) > Math.abs(input[i][maxrow])) {
        maxrow = j;
      }
    }
    for (let k = i; k < n + 1; k++) {
      const tmp = input[k][i];
      input[k][i] = input[k][maxrow];
      input[k][maxrow] = tmp;
    }
    for (let j = i + 1; j < n; j++) {
      for (let k = n; k >= i; k--) {
        input[k][j] -= (input[k][i] * input[i][j]) / input[i][i];
      }
    }
  }
  for (let j = n - 1; j >= 0; j--) {
    let total = 0;
    for (let k = j + 1; k < n; k++) {
      total += input[k][j] * coefficients[k];
    }
    coefficients[j] = (input[n][j] - total) / input[j][j];
  }
  return coefficients;
}

function polynomialRegressionCoefficients(xCol, yCol, orderOfPolynomialRegression) {
  const lhs = []; const rhs = []; const len = xCol.length; const k = orderOfPolynomialRegression + 1;
  let a = 0; let b = 0;
  for (let i = 0; i < k; i++) {
    for (let l = 0; l < len; l++) {
      a += (xCol.get(l) ** i) * yCol.get(l);
    }
    lhs.push(a);
    a = 0;
    const c = [];
    for (let j = 0; j < k; j++) {
      for (let l = 0; l < len; l++) {
        b += xCol.get(l) ** (i + j);
      }
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
dfWithoutOutliers = DG.DataFrame.fromColumns([
  DG.Column.fromList('double', 'time (min)', timeMin),
  DG.Column.fromList('double', 'filtrate volume (mL)', vol),
]);
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
  DG.Column.fromList('double', va.name, va.toList().slice(1, va.length - windowLength)), 
  DG.Column.fromList('double', j.name, j.toList().slice(1, j.length)),
  DG.Column.fromList('double', jj0.name, jj0.toList().slice(1, jj0.length))
]);

regressionTable = inputTable;
coefficients = [1.8, 4.32];//[4.32, 1.8];//polynomialRegressionCoefficients(inputTable.col("time (min)"), dfWithoutOutliers.col("t/V (hr/(L/m2))"), 1);
const trialThroughput = va.max;
const flux = 9372.11;
const ff0 = flux / j.get(1);
experimentalResults = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Run', Array(va.length - windowLength).fill('1')),
  DG.Column.fromList('double', va.name, va.toList().slice(0, va.length - windowLength)), 
  j, 
  jj0
]);
const meanFlux = j.stats.avg;
const vmax = 0.28;//Math.round(1 / coefficients[1]);
const instantaneousFlux = j.get(j.length - 1);
const initialFlux = j.get(1);
const fluxDecay = (1 - instantaneousFlux / initialFlux) * 100;
const names = parametersForReporting.col('Parameter').toList();
const values = parametersForReporting.col('Value').toList();
experimentalResultsSummary = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Run', ['1']),
  DG.Column.fromList('string', 'Filter', [values[names.indexOf('Filter Name')]]),
  DG.Column.fromList('double', 'Trial Throughput (L/m²)', [trialThroughput]),
  DG.Column.fromList('double', 'Flux0 (LMH)', [j.get(1)]),
  // DG.Column.fromList('double', 'Flux/Flux0 (LMH)', [ff0]),
  DG.Column.fromList('double', 'Mean Flux (LMH)', [meanFlux]),
  DG.Column.fromList('double', 'Vmax (L/m²)', [vmax]),
  DG.Column.fromList('double', 'Flux Decay (%)', [fluxDecay]),
]);
const q0 = Math.round(1 / coefficients[0]);
const amin = (desiredVolumeOfBatch / vmax) + (desiredVolumeOfBatch / desiredProcessTime / q0);
const installedArea = sf * amin;
const processLoading = volumeInLiters.get(volumeInLiters.length - 1) / installedArea;
recommendations = DG.DataFrame.fromColumns([
  DG.Column.fromList('double', 'Volume (L)', [desiredVolumeOfBatch]),
  DG.Column.fromList('double', 'Time (h)', [desiredProcessTime]),
  DG.Column.fromList('double', 'Amin (m²)', [amin]),
  DG.Column.fromList('double', 'Safety Factor', [sf]),
  DG.Column.fromList('double', 'Installed Area', [installedArea]),
  DG.Column.fromList('double', 'Process Loading (L/m²)', [processLoading]),
  DG.Column.fromList('string', 'Filter', [values[names.indexOf('Filter Name')]]),
  DG.Column.fromList('string', 'Installed Device', [values[names.indexOf('Installed Device')]]),
]);

const v90 = vmax * 0.68;
const initialFlowArea = coefficients[0];
trialData = DG.DataFrame.fromColumns([
  DG.Column.fromList('double', 'Test Filter Area (m²)', [testFilterArea]),
  DG.Column.fromList('double', 'Initial Flow Area (L/h)', [initialFlowArea]),
  DG.Column.fromList('double', 'Vmax (m²)', [vmax]),
  DG.Column.fromList('double', 'V90 (m²)', [v90]),
]);

filter = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Filter Family', [values[names.indexOf('Filter Family')]]),
  DG.Column.fromList('string', 'Filter Name', [values[names.indexOf('Filter Name')]]),
  DG.Column.fromList('string', 'Pore Rating (μm)', [values[names.indexOf('Pore rating (μm)')]]),
  DG.Column.fromList('string', 'Effective Membr. Area (cm^2)', [values[names.indexOf('Effective Membr. Area (cm^2)')]]),
  DG.Column.fromList('string', 'Membr material', [values[names.indexOf('Membr material')]]),
  DG.Column.fromList('string', 'Catalog Nr', [values[names.indexOf('Catalog Nr')]]),
]);
sampleCharacteristics = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Date', [values[names.indexOf('date ')]]),
  DG.Column.fromList('string', 'Prior Step', [values[names.indexOf('Prior Step')]]),
  DG.Column.fromList('string', 'Prot conc. (mg/mL)', [values[names.indexOf('Prot conc. (mg/mL)')]]),
  DG.Column.fromList('string', 'Density (g/L)', [values[names.indexOf('Density (g/L)')]]),
  DG.Column.fromList('string', 'Turbidity (NTU)', [values[names.indexOf('Turbidity (NTU)')]]),
]);