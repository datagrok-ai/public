//name: Vmax
//description: Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
//language: javascript
//tags: filtration, scale-up
//meta.department: BTDS
//meta.status: Upstream
//input: dataframe inputTable {caption: Input table; viewer: OutliersSelectionViewer(block: 25) | Scatter Plot(y: "t/V (hr/(L/m²))", block: 75, filter: "!${isOutlier}", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true", markers: "isOutlier", lassoTool: "true", legendVisibility: "Never", filterByZoom: "false")}
//input: dataframe reportingParameters {caption: Reporting parameters; viewer: Grid()}
//input: double testFilterArea = 3.5 {caption: Test filter area; units: cm²} [Test filter area]
//input: double desiredVolumeOfBatch = 25 {caption: Desired batch volume; units: L} [Desired batch volume]
//input: double desiredProcessTime = 0.5 {caption: Desired process time; units: hr} [Desired process time]
//input: double sf = 1.5 {caption: Safety factor} [Safety factor]
//help-url: https://github.com/datagrok-ai/public/blob/master/packages/Compute/src/help.md
//output: dataframe filtrationKinetics {caption: Filtration kinetics; viewer: Scatter Plot(block: 75, filter: "!${isOutlier}", markers: "isOutlier", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, showRegressionLine: "true", legendVisibility: "Never"); category: OUTPUT}
//output: dataframe trialData {caption: Trial data; viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); category: OUTPUT}
//output: dataframe fluxDecayDf {caption: Flux decay experimental results; viewer: Scatter Plot(description: "Absolute flux decay", y: "Flux", block: 50, filter: "!${isOutlier}", markers: "isOutlier", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, legendVisibility: "Never") | Scatter Plot(description: "Relative flux decay", y: "Rel Flux", block: 50, filter: "!${isOutlier}", markers: "isOutlier", showFilteredOutPoints: "true",  filteredOutRowsColor: 4293991195, legendVisibility: "Never"); category: OUTPUT}
//output: dataframe experimentalResultsSummary {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Experimental results summary; category: OUTPUT}
//output: dataframe recommendations {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Recommendations; category: OUTPUT}
//output: dataframe sampleCharacteristics {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Sample characteristics; category: OUTPUT}
//output: dataframe filter {viewer: Grid(block: 25, showColumnLabels: false, showRowHeader: false); caption: Filter; category: OUTPUT}
//editor: Compute:ComputationView
// inputTable = grok.shell.tables[0]; reportingParameters = grok.shell.tables[1];
// testFilterArea = 3.5; desiredVolumeOfBatch = 25; desiredProcessTime = 0.5; sf = 1.5;

const COL_NAMES = {
  INPUT_TIME_IN_MINUTES: "time (min)",
  INPUT_VOLUME_IN_MILLILITERS: "filtrate volume (mL)",
  INPUT_IS_OUTLIER: "isOutlier",
  TIME_IN_HOURS: "Time [hr]",
  TRIAL_THROUGHPUT: "Throughput [L/m²]",
  T_V: "Time/volume/area [hr/L/m²]",
  FLOW_RATE: "Flow rate",
  FLUX: "Flux",
  RELATIVE_FLUX: "Rel Flux"
};

const PARAM_NAMES = {
  FILTER_FAMILY: "Filter family",
  FILTER_NAME: "Filter name",
  PORE_RATING: "Pore rating (μm)",
  EFFECTIVE_MEMBR_AREA: "Effective membr. area (cm^2)",
  MEMBR_MATERIAL: "Membr material",
  CATALOG_NR: "Catalog nr",
  DATE: "Date",
  PRIOR_STEP: "Prior step",
  PROT_CONC: "Prot conc. (mg/mL)",
  DENSITY: "Density (g/L)",
  TURBIDITY: "Turbidity (NTU)",
  INSTALLED_DEVICE: "Installed device",
  PROJECT_NAME: "Project name",  // PARAMETERS BELOW ARE NOT USED
  FILTRATION_STAGE: "Filtration stage",
  ELN_REFERENCE: "ELN reference",
  BIOREACTOR_ID: "Bioreactor ID",
  BIOREACTOR_PCV: "Bioreactor PCV (%)",
  BIOREACTOR_TURBIDITY: "Bioreactor turbidity (NTU)",
  FILTER_FEED_TURBIDITY: "Filter feed turbidity (NTU)",
  FILTER_FEED_PCV: "Filter feed PCV (%)",
  FILTER_FEED_PH: "Filter feed pH",
  FILTER_FEED_LDH: "Filter feed LDH (U/L)",
  FILTER_FEED_CONDUCTIVITY: "Filter feed conductivity (mS/cm)"
};

function linearRegression(x, y) {
  let sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
  for (let i = 0; i < y.length; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_xy += x[i] * y[i];
    sum_xx += x[i] * x[i];
  }
  const slope = (y.length * sum_xy - sum_x * sum_y) / (y.length * sum_xx - sum_x * sum_x);
  const intercept = (sum_y - slope * sum_x) / y.length;
  return { slope: slope, intercept: intercept };
}

const timeInMinutes = inputTable.col(COL_NAMES.INPUT_TIME_IN_MINUTES);
const volumeInMilliliters = inputTable.col(COL_NAMES.INPUT_VOLUME_IN_MILLILITERS);
const isOutlier = inputTable.col(COL_NAMES.INPUT_IS_OUTLIER);

let df = DG.DataFrame.fromColumns([timeInMinutes, volumeInMilliliters, isOutlier]);

df.columns.addNewFloat(COL_NAMES.TIME_IN_HOURS).init((i) => timeInMinutes.get(i) / 60);
const timeInHours = df.col(COL_NAMES.TIME_IN_HOURS);

df.columns.addNewFloat(COL_NAMES.TRIAL_THROUGHPUT).init((i) => 10 * volumeInMilliliters.get(i) / testFilterArea);
const trialThroughputCol = df.col(COL_NAMES.TRIAL_THROUGHPUT);

df.columns.addNewFloat(COL_NAMES.T_V).init((i) => timeInHours.get(i) / trialThroughputCol.get(i));
const tv = df.col(COL_NAMES.T_V);

df.columns.addNewFloat(COL_NAMES.FLOW_RATE).init((i) => (i > 0) ?
  0.06 * (volumeInMilliliters.get(i) - volumeInMilliliters.get(i - 1)) / (timeInMinutes.get(i) - timeInMinutes.get(i - 1)) : 0
);
const flowRate = df.col(COL_NAMES.FLOW_RATE);

df.columns.addNewFloat(COL_NAMES.FLUX).init((i) => 10000 * flowRate.get(i) / testFilterArea);
const flux = df.col(COL_NAMES.FLUX);

df.columns.addNewFloat(COL_NAMES.RELATIVE_FLUX).init((i) => flux.get(i) / flux.get(1));
const relativeFlux = df.col(COL_NAMES.RELATIVE_FLUX);

fluxDecayDf = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array(COL_NAMES.TRIAL_THROUGHPUT, df.col(COL_NAMES.TRIAL_THROUGHPUT).toList().slice(1)),
  DG.Column.fromFloat32Array(COL_NAMES.FLUX, df.col(COL_NAMES.FLUX).toList().slice(1)),
  DG.Column.fromFloat32Array(COL_NAMES.RELATIVE_FLUX, df.col(COL_NAMES.RELATIVE_FLUX).toList().slice(1)),
  DG.Column.fromList(DG.COLUMN_TYPE.BOOL, isOutlier.name, isOutlier.toList().slice(1))
]);

filtrationKinetics = DG.DataFrame.fromColumns([
  timeInHours,
  df.col(COL_NAMES.T_V),
  volumeInMilliliters,
  isOutlier
]);

df.col(COL_NAMES.INPUT_IS_OUTLIER).markers.assign('true', DG.MARKER_TYPE.OUTLIER);
filtrationKinetics.col(COL_NAMES.INPUT_IS_OUTLIER).markers.assign('true', DG.MARKER_TYPE.OUTLIER);

for (let i = df.rowCount - 1; i > -1; i--)
  if (isOutlier.get(i))
    df.rows.removeAt(i, 1, false);

const coefficients = linearRegression(df.col(COL_NAMES.TIME_IN_HOURS).toList(), df.col(COL_NAMES.T_V).toList());
const vmax = 1 / coefficients.slope;
const j0 = 1 / coefficients.intercept;
const flux0 = flux.get(1);
const amin = desiredVolumeOfBatch / vmax + desiredVolumeOfBatch / desiredProcessTime / j0;
const installedArea = sf * amin;
const processLoading = desiredVolumeOfBatch / installedArea;
const trialThroughput = df.col(COL_NAMES.TRIAL_THROUGHPUT).get(df.rowCount - 1);
const fluxDecay = (1 - df.col(COL_NAMES.RELATIVE_FLUX).get(df.rowCount - 1)) * 100;
const initialFlowRate = flowRate.get(1);
const meanFlux = flux.stats.avg * df.rowCount / (df.rowCount - 1);  // exclude zero from calculation

const names = reportingParameters.col('Parameter name').toList();
const values = reportingParameters.col('Parameter value').toList();

experimentalResultsSummary = DG.DataFrame.fromObjects([
  { Name: 'Run', Value: '1' },
  { Name: 'Filter', Value: values[names.indexOf(PARAM_NAMES.FILTER_FAMILY)] },
  { Name: 'Trial Throughput [L/m²]', Value: Math.round(1000 * trialThroughput) / 1000 },
  { Name: 'Flux0 [LMH]', Value: Math.round(1000 * flux0) / 1000 },
  { Name: 'Mean Flux [LMH]', Value: Math.round(1000 * meanFlux) / 1000 },
  { Name: 'Vmax [L/m²]', Value: Math.round(1000 * vmax) / 1000 },
  { Name: 'Flux Decay [%]', Value: Math.round(1000 * fluxDecay) / 1000 }
]);

recommendations = DG.DataFrame.fromObjects([
  { Name: 'Volume [L]', Value: desiredVolumeOfBatch },
  { Name: 'Time [h]', Value: desiredProcessTime },
  { Name: 'Amin [m²]', Value: Math.round(1000 * amin) / 1000 },
  { Name: 'Safety Factor', Value: sf },
  { Name: 'Installed Area', Value: Math.round(1000 * installedArea) / 1000 },
  { Name: 'Process Loading [L/m²]', Value: Math.round(1000 * processLoading) / 1000 },
  { Name: PARAM_NAMES.FILTER_FAMILY, Value: values[names.indexOf(PARAM_NAMES.FILTER_FAMILY)] },
  { Name: PARAM_NAMES.INSTALLED_DEVICE, Value: values[names.indexOf(PARAM_NAMES.INSTALLED_DEVICE)] }
]);

trialData = DG.DataFrame.fromObjects([
  { Name: 'Test Filter Area [cm²]', Value: testFilterArea },
  { Name: 'Initial Flow Rate [L/h]', Value: Math.round(1000 * initialFlowRate) / 1000 },
  { Name: 'Vmax [L/m²]', Value: Math.round(1000 * vmax) / 1000 }
]);

filter = DG.DataFrame.fromObjects([
  { Name: PARAM_NAMES.FILTER_FAMILY, Value: values[names.indexOf(PARAM_NAMES.FILTER_FAMILY)] },
  { Name: PARAM_NAMES.FILTER_NAME, Value: values[names.indexOf(PARAM_NAMES.FILTER_NAME)] },
  { Name: PARAM_NAMES.PORE_RATING, Value: values[names.indexOf(PARAM_NAMES.PORE_RATING)] },
  { Name: PARAM_NAMES.EFFECTIVE_MEMBR_AREA, Value: values[names.indexOf(PARAM_NAMES.EFFECTIVE_MEMBR_AREA)] },
  { Name: PARAM_NAMES.MEMBR_MATERIAL, Value: values[names.indexOf(PARAM_NAMES.MEMBR_MATERIAL)] },
  { Name: PARAM_NAMES.CATALOG_NR, Value: values[names.indexOf(PARAM_NAMES.CATALOG_NR)] }
]);

sampleCharacteristics = DG.DataFrame.fromObjects([
  { Name: PARAM_NAMES.DATE, Value: values[names.indexOf(PARAM_NAMES.DATE)].slice(0, 10) },
  { Name: PARAM_NAMES.PRIOR_STEP, Value: values[names.indexOf(PARAM_NAMES.PRIOR_STEP)] },
  { Name: PARAM_NAMES.PROT_CONC, Value: values[names.indexOf(PARAM_NAMES.PROT_CONC)] },
  { Name: PARAM_NAMES.DENSITY, Value: values[names.indexOf(PARAM_NAMES.DENSITY)].slice(0, 5) },
  { Name: PARAM_NAMES.TURBIDITY, Value: values[names.indexOf(PARAM_NAMES.TURBIDITY)] }
]);