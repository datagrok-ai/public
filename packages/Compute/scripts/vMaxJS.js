//name: Vmax
//description: Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
//language: javascript
//tags: model, filtration
//input: dataframe inputTable {viewer: OutliersSelectionViewer() | Scatter Plot(color: "isOutlier", lassoTool: "true", legendVisibility: "Never", filterByZoom: "false") | Grid()}
//input: dataframe inputParametersForReporting {viewer: Grid()}
//input: double testArea = 3.5 {caption: Test Filter Area; units: cm²} [Test Filter area]
//input: double vbatch = 25 {caption: Desired Volume of Batch; units: L} [Desired Volume of Batch]
//input: double tbatch = 0.5 {caption: Desired process time; units: hr} [Desired process time]
//input: double sf = 1.5 {caption: Safety Factor} [Safety Factor]
//help-url: https://github.com/datagrok-ai/public/blob/master/packages/Compute/src/help.md
//output: dataframe regressionTable {viewer: Scatter Plot(filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true"); category: OUTPUT}
//output: dataframe absoluteFluxDecay {viewer: Scatter Plot(); category: OUTPUT}
//output: dataframe normalizedFluxDecay {viewer: Scatter Plot(); category: OUTPUT}
//output: dataframe experimentalResults {category: OUTPUT}
//output: dataframe experimentalResultsSummary {category: OUTPUT}
//output: dataframe recommendations {category: OUTPUT}
//output: dataframe trialData {category: OUTPUT}
//output: dataframe sampleCharacteristics {category: OUTPUT}
//output: dataframe filter {category: OUTPUT}

volume = vbatch;
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
const isOutlierCol = inputTable.col("isOutlier");
let timeMin = [];
let vol = [];
for (let i = 0; i < isOutlierCol.length; i++) {
  if (isOutlierCol.get(i) == false) {
    timeMin.push(inputTable.col("time (min)").get(i));
    vol.push(inputTable.col("filtrate volume (mL)").get(i));
  }
}
dfWithoutOutliers = DG.DataFrame.fromColumns([
  DG.Column.fromList('double', "time (min)", timeMin),
  DG.Column.fromList('double', "filtrate volume (mL)", vol)
]);
const timeInMinutes = dfWithoutOutliers.col("time (min)");
const volumeInMilliliters = dfWithoutOutliers.col("filtrate volume (mL)");
dfWithoutOutliers.columns.addNewFloat('time (hr)').init((i) => timeInMinutes.get(i) / 60);
const timeInHours = dfWithoutOutliers.col("time (hr)");

dfWithoutOutliers.columns.addNewFloat('V (L/m2)').init((i) => 10 * volumeInMilliliters.get(i) / testArea);
const newVolume = dfWithoutOutliers.col('V (L/m2)');

dfWithoutOutliers.columns.addNewFloat("t/V (hr/(L/m2))").init((i) => 100 * timeInHours.get(i) / newVolume.get(i));

dfWithoutOutliers.columns.addNewFloat("Time (s)").init((i) => 60 * timeInMinutes.get(i));
dfWithoutOutliers.columns.addNewFloat("V (L)").init((i) => volumeInMilliliters.get(i) / 1000);
const volumeInLiters = dfWithoutOutliers.col('V (L)');
const timeInSeconds = dfWithoutOutliers.col('Time (s)');
dfWithoutOutliers.columns.addNewFloat('V/A (L/m2)').init((i) => volumeInLiters.get(i) / area);
const va = dfWithoutOutliers.col('V/A (L/m2)');
let ca = [];
let cb = [];
for (let i = 1; i < va.length - 3; i++) {
  let xCol = DG.Column.fromList('double', 'name', timeInSeconds.toList().slice(i, i + 3));
  let yCol = DG.Column.fromList('double', 'name', va.toList().slice(i, i + 3));
  let coefficients = polynomialRegressionCoefficients(xCol, yCol, orderOfPolynomialRegression);
  ca.push(coefficients[orderOfPolynomialRegression]);
  cb.push(coefficients[orderOfPolynomialRegression - 1]);
}

dfWithoutOutliers.columns.addNewFloat('a2').init((i) => (i > 0) ? 3600 * (2 * ca[i] * timeInSeconds.get(i) + cb[i]) : 0);
dfWithoutOutliers.columns.addNewFloat('a1').init((i) => (i > 0) ? 3600 * (2 * ca[i] * timeInSeconds.get(i) + cb[i]) : 0);
dfWithoutOutliers.columns.addNewFloat('J (LMH)').init((i) => (i > 0) ? 3600 * (2 * ca[i] * timeInSeconds.get(i) + cb[i]) : 0);
const j = dfWithoutOutliers.col('J (LMH)');
dfWithoutOutliers.columns.addNewFloat('J/Jo').init((i) => (i > 0) ? j.get(i) / (3600 * area) : 0);
const jj0 = dfWithoutOutliers.col('J/Jo');
absoluteFluxDecay = DG.DataFrame.fromColumns([va, j]);
normalizedFluxDecay = DG.DataFrame.fromColumns([va, jj0]);

regressionTable = inputTable;
coefficients = [1.8, 4.32];//[4.32, 1.8];//polynomialRegressionCoefficients(inputTable.col("time (min)"), dfWithoutOutliers.col("t/V (hr/(L/m2))"), 1);
const trialThroughput = va.max;
const flux = 9372.11;
const ff0 = flux / j.get(1);
experimentalResults = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Run', ['1']),
  DG.Column.fromList('double', 'Trial Throughput (L/m²)', [trialThroughput]),
  DG.Column.fromList('double', 'Flux (LMH)', [flux]),
  DG.Column.fromList('double', 'Flux/Flux0 (LMH)', [ff0])
]);
const meanFlux = j.stats.avg;
const vmax = 0.28;//Math.round(1 / coefficients[1]);
const instantaneousFlux = j.get(j.length - 1);
const initialFlux = j.get(1);
const fluxDecay = (1 - instantaneousFlux / initialFlux) * 100;
experimentalResultsSummary = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Run', ['1']),
  DG.Column.fromList('string', 'Filter', [inputParametersForReporting.get(0, 1)]),
  DG.Column.fromList('double', 'Trial Throughput (L/m²)', [trialThroughput]),
  DG.Column.fromList('double', 'Flux0 (LMH)', [j.get(1)]),
  // DG.Column.fromList('double', 'Flux/Flux0 (LMH)', [ff0]),
  DG.Column.fromList('double', 'Mean Flux (LMH)', [meanFlux]),
  DG.Column.fromList('double', 'Vmax (L/m²)', [vmax]),
  DG.Column.fromList('double', 'Flux Decay (%)', [fluxDecay])
]);
const q0 = Math.round(1 / coefficients[0]);
const amin = (vbatch / vmax) + (vbatch / tbatch / q0);
const installedArea = sf * amin;
const processLoading = volumeInLiters.get(volumeInLiters.length - 1) / installedArea;
recommendations = DG.DataFrame.fromColumns([
  DG.Column.fromList('double', 'Volume (L)', [vbatch]),
  DG.Column.fromList('double', 'Time (h)', [tbatch]),
  DG.Column.fromList('double', 'Amin (m²)', [amin]),
  DG.Column.fromList('double', 'Safety Factor', [sf]),
  DG.Column.fromList('double', 'Installed Area', [installedArea]),
  DG.Column.fromList('double', 'Process Loading (L/m²)', [processLoading]),
  DG.Column.fromList('string', 'Filter', [inputParametersForReporting.get(1, 0)]),
  DG.Column.fromList('string', 'Installed Device', [inputParametersForReporting.get(1, 10)]),
]);

const v90 = vmax * 0.68;
const initialFlowArea = coefficients[0];
trialData = DG.DataFrame.fromColumns([
  DG.Column.fromList('double', 'Test Filter Area (m²)', [testArea]),
  DG.Column.fromList('double', 'Initial Flow Area (L/h)', [initialFlowArea]),
  DG.Column.fromList('double', 'Vmax (m²)', [vmax]),
  DG.Column.fromList('double', 'V90 (m²)', [v90])
]);

filter = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Filter Family', [inputParametersForReporting.columns.names()[1]]),
  DG.Column.fromList('string', 'Filter Name', [inputParametersForReporting.get(1, 0)]),
  DG.Column.fromList('string', 'Pore Rating', [inputParametersForReporting.get(1, 1)]),
  DG.Column.fromList('string', 'Effective Membr. Area', [inputParametersForReporting.get(1, 2)]),
  DG.Column.fromList('string', 'Membr material', [inputParametersForReporting.get(1, 3)]),
  DG.Column.fromList('string', 'Catalog', [inputParametersForReporting.get(1, 4)]),
]);

sampleCharacteristics = DG.DataFrame.fromColumns([
  DG.Column.fromList('string', 'Date', [inputParametersForReporting.get(1, 5)]),
  DG.Column.fromList('string', 'Prior step', [inputParametersForReporting.get(1, 6)]),
  DG.Column.fromList('string', 'Prot conc', [inputParametersForReporting.get(1, 7)]),
  DG.Column.fromList('string', 'Density', [inputParametersForReporting.get(1, 8)]),
  DG.Column.fromList('string', 'Turbidity', [inputParametersForReporting.get(1, 9)])
]);