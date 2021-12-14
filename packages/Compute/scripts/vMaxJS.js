//name: Vmax
//description: Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
//language: javascript
//tags: model, filtration
//input: dataframe inputTable {caption: Input table; viewer: OutliersSelectionViewer(block: 25) | Scatter Plot(y: "t/V (hr/(L/m²))", block: 75, filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true")}
//input: dataframe reportingParameters {caption: Reporting parameters; viewer: Grid()}
//input: double testFilterArea = 3.5 {caption: Test filter area; units: cm²} [Test filter area]
//input: double desiredVolumeOfBatch = 25 {caption: Desired batch volume; units: L} [Desired batch volume]
//input: double desiredProcessTime = 0.5 {caption: Desired process time; units: hr} [Desired process time]
//input: double sf = 1.5 {caption: Safety factor} [Safety factor]
//help-url: https://github.com/datagrok-ai/public/blob/master/packages/Compute/src/help.md
//output: dataframe kineticFiltration {caption: Kinetic filtration; viewer: Scatter Plot(block: 75, filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true"); category: OUTPUT}
//output: dataframe trialData {caption: Trial data; viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); category: OUTPUT}
//output: dataframe fluxDecayDf {caption: Flux decay experimental results; viewer: Scatter Plot(description: "Absolute flux decay", y: "J (LMH)", block: 50, filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195) | Scatter Plot(description: "Relative flux decay", y: "J/Jo", block: 50, filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195); category: OUTPUT}
//output: dataframe experimentalResultsSummary {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Experimental results summary; category: OUTPUT}
//output: dataframe recommendations {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Recommendations; category: OUTPUT}
//output: dataframe sampleCharacteristics {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Sample characteristics; category: OUTPUT}
//output: dataframe filter {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Filter; category: OUTPUT}

// inputTable = grok.shell.tables[0];
// reportingParameters = grok.shell.tables[1];
// testFilterArea = 3.5;
// desiredVolumeOfBatch = 25;
// desiredProcessTime = 0.5;
// sf = 1.5;

const expectedInputColNames = {
  "time": "time (min)",
  "volume": "filtrate volume (mL)",
  "isOutlier": "isOutlier"
};

const expectedReportParamsNames = {
  1: "Filter family",
  2: "Filter name",
  3: "Pore rating (μm)",
  4: "Effective membr. area (cm^2)",
  5: "Membr material",
  6: "Catalog nr",
  7: "Date",
  8: "Prior step",
  9: "Prot conc. (mg/mL)",
  10: "Density (g/L)",
  11: "Turbidity (NTU)",
  12: "Installed device",
  13: "Project name",
  14: "Filtration stage",
  15: "ELN reference",
  16: "Bioreactor ID",
  17: "Bioreactor PCV (%)",
  18: "Bioreactor turbidity (NTU)",
  19: "Filter feed turbidity (NTU)",
  20: "Filter feed PCV (%)",
  21: "Filter feed pH",
  22: "Filter feed LDH (U/L)",
  23: "Filter feed conductivity (mS/cm)"
};

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
const timeInMinutesCol = inputTable.col(expectedInputColNames['time']);
const volumeInMillilitersCol = inputTable.col(expectedInputColNames['volume']);
const isOutlierCol = inputTable.col(expectedInputColNames["isOutlier"]);

let df = DG.DataFrame.fromColumns([timeInMinutesCol, volumeInMillilitersCol, isOutlierCol]);
df.columns.addNewFloat('time (hr)').init((i) => timeInMinutesCol.get(i) / 60);
const timeInHours = df.col('time (hr)');
df.columns.addNewFloat('V (L/m2)').init((i) => 10 * volumeInMillilitersCol.get(i) / testFilterArea);
const newVolume = df.col('V (L/m2)');
df.columns.addNewFloat('t/V (hr/(L/m2))').init((i) => 100 * timeInHours.get(i) / newVolume.get(i));
df.columns.addNewFloat('Time (s)').init((i) => 60 * timeInMinutesCol.get(i));
df.columns.addNewFloat('V (L)').init((i) => volumeInMillilitersCol.get(i) / 1000);
const volumeInLiters = df.col('V (L)');
const timeInSeconds = df.col('Time (s)');
df.columns.addNewFloat('V/A (L/m2)').init((i) => volumeInLiters.get(i) / area);
const va = df.col('V/A (L/m2)');

const slopes = [], intercepts = [];
const windowLength = 3;
for (let i = 0; i < inputTable.rowCount - windowLength; i++) {
  const xCol = DG.Column.fromFloat32Array('name', timeInSeconds.toList().slice(i, i + windowLength));
  const yCol = DG.Column.fromFloat32Array('name', va.toList().slice(i, i + windowLength));
  const coefficients = polynomialRegressionCoefficients(xCol, yCol, orderOfPolynomialRegression);
  slopes.push(coefficients[orderOfPolynomialRegression]);
  intercepts.push(coefficients[orderOfPolynomialRegression - 1]);
}

let newDf = DG.DataFrame.create(inputTable.rowCount - windowLength);
newDf.columns.addNewFloat('a2').init((i) => (i > 0) ? 3600 * (2 * slopes[i] * timeInSeconds.get(i) + intercepts[i]) : 0);
newDf.columns.addNewFloat('a1').init((i) => (i > 0) ? 3600 * (2 * slopes[i] * timeInSeconds.get(i) + intercepts[i]) : 0);
newDf.columns.addNewFloat('J (LMH)').init((i) => (i > 0) ? 3600 * (2 * slopes[i] * timeInSeconds.get(i) + intercepts[i]) : 0);
const j = newDf.col('J (LMH)');
newDf.columns.addNewFloat('J/Jo').init((i) => (i > 0) ? j.get(i) / (3600 * area) : 0);
let jj0 = newDf.col('J/Jo');
jj0 = DG.Column.fromFloat32Array('J/Jo', jj0.toList().map((e) => 100 * e / jj0.max));
fluxDecayDf = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Run', Array(jj0.length - 1).fill('1')),
  DG.Column.fromFloat32Array(va.name, va.toList().slice(1, inputTable.rowCount - windowLength)),
  DG.Column.fromFloat32Array(j.name, j.toList().slice(1, j.length)),
  DG.Column.fromFloat32Array(jj0.name, jj0.toList().slice(1, jj0.length)),
  DG.Column.fromList('bool', isOutlierCol.name, isOutlierCol.toList().slice(1, isOutlierCol.length - windowLength))
]);

kineticFiltration = DG.DataFrame.fromColumns([
  timeInMinutesCol,
  inputTable.col('t/V (hr/(L/m²))'),
  volumeInMillilitersCol,
  isOutlierCol
]);

const timeMin = [], vol = [], isOutlier = [];
const isOutlierColList = isOutlierCol.toList();
for (let i = 0; i < inputTable.rowCount - windowLength; i++) {
  if (isOutlierCol.get(i) == false) {
    timeMin.push(timeInMinutesCol.get(i));
    vol.push(volumeInMillilitersCol.get(i));
    isOutlier.push(true in isOutlierColList.slice(i, i + windowLength));
  }
}

coefficients = polynomialRegressionCoefficients(DG.Column.fromFloat32Array("time (min)", timeMin), DG.Column.fromFloat32Array("t/V (hr/(L/m2))", vol), 1);
const initialFlowRate = coefficients[0];
const trialThroughput = Math.round(1000 * va.max) / 1000;
const flux = initialFlowRate / testFilterArea;
const initialFlux = Math.round(1000 * j.get(1)) / 1000;
const ff0 = flux / initialFlux;
const meanFlux = Math.round(1000 * j.stats.avg) / 1000;
const vmax = Math.round(1000 * coefficients[1]) / 1000;
const v90 = Math.round(1000 * vmax * 0.68) / 1000;
const instantaneousFlux = j.get(j.length - 1);
const fluxDecay = Math.round(1000 * (1 - instantaneousFlux / initialFlux) * 100) / 1000;
const names = reportingParameters.col('Parameter name').toList();
const values = reportingParameters.col('Parameter value').toList();
const q0 = 1 / coefficients[0];
const amin = Math.round(1000 * ((desiredVolumeOfBatch / vmax) + (desiredVolumeOfBatch / desiredProcessTime / q0))) / 1000;
const installedArea = Math.round(1000 * sf * amin) / 1000;
const processLoading = Math.round(1000 * volumeInLiters.get(volumeInLiters.length - 1) / installedArea) / 1000;

experimentalResultsSummary = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Name', ['Run', 'Filter', 'Trial Throughput (L/m²)', 'Flux0 (LMH)', 'Mean Flux (LMH)', 'Vmax (L/m²)', 'Flux Decay (%)']),
  DG.Column.fromStrings('Value', ['1', values[names.indexOf(expectedReportParamsNames[1])], String(trialThroughput), String(initialFlux), String(meanFlux), String(vmax), String(fluxDecay)])
]);

recommendations = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Name', ['Volume (L)', 'Time (h)', 'Amin (m²)', 'Safety Factor', 'Installed Area', 'Process Loading (L/m²)', expectedReportParamsNames[1], expectedReportParamsNames[12]]),
  DG.Column.fromStrings('Value', [String(desiredVolumeOfBatch), String(desiredProcessTime), String(amin), String(sf), String(installedArea), String(processLoading), values[names.indexOf(expectedReportParamsNames[1])], values[names.indexOf(expectedReportParamsNames[12])]])
]);

trialData = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Name', ['Test Filter Area (m²)', 'Initial Flow Rate (L/h)', 'Vmax (m²)', 'V90 (m²)']),
  DG.Column.fromStrings('Value', [String(Math.round(1000 * testFilterArea) / 1000), String(Math.round(1000 * initialFlowRate) / 1000), String(vmax), String(v90)])
]);

filter = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Name', [expectedReportParamsNames[1], expectedReportParamsNames[2], expectedReportParamsNames[3], expectedReportParamsNames[4], expectedReportParamsNames[5], expectedReportParamsNames[6]]),
  DG.Column.fromStrings('Value', [values[names.indexOf(expectedReportParamsNames[1])], values[names.indexOf(expectedReportParamsNames[2])], values[names.indexOf(expectedReportParamsNames[3])], values[names.indexOf(expectedReportParamsNames[4])], values[names.indexOf(expectedReportParamsNames[5])], values[names.indexOf(expectedReportParamsNames[6])]])
]);

sampleCharacteristics = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Name', [expectedReportParamsNames[7], expectedReportParamsNames[8], expectedReportParamsNames[9], expectedReportParamsNames[10], expectedReportParamsNames[11]]),
  DG.Column.fromStrings('Value', [values[names.indexOf(expectedReportParamsNames[7])].slice(0, 10), values[names.indexOf(expectedReportParamsNames[8])], values[names.indexOf(expectedReportParamsNames[9])], values[names.indexOf(expectedReportParamsNames[10])], values[names.indexOf(expectedReportParamsNames[11])]])
]);