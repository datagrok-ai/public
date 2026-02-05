// Constants and type definitions for probabilistic scoring (pMPO)
// Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/

/** Minimum number of samples required to compute pMPO */
export const MIN_SAMPLES_COUNT = 10;

export const PMPO_NON_APPLICABLE = 'pMPO is not applicable';
export const PMPO_COMPUTE_FAILED = 'Failed to compute pMPO parameters';

/** Basic statistics for desired and non-desired compounds */
export type BasicStats = {
  desAvg: number,
  desStd: number,
  nonDesAvg: number,
  nonDesStd: number,
  min: number,
  max: number,
};

/** Descriptor statistics including basic stats, t-statistics and p-value */
export type DescriptorStatistics = BasicStats & {
  desLen: number,
  nonSesLen: number,
  tstat: number,
  pValue: number,
};

/** Cutoff parameters for the basic functions of the pMPO model */
export type Cutoff = {
  cutoff: number,
  cutoffDesired: number,
  cutoffNotDesired: number,
};

/** Generalized Sigmoid parameters for the desirability functions of the pMPO model */
export type SigmoidParams = {
  pX0: number,
  b: number,
  c: number,
};

/** pMPO parameters including basic stats, cutoffs, sigmoid params, z-score, weight, intersections */
export type PmpoParams = BasicStats & Cutoff & SigmoidParams & {
  zScore: number,
  weight: number,
  intersections: number[],
  x0: number,
  xBound: number,
  inflection: number,
};

export type CorrelationTriple = [string, string, number];

const DESIRED = 'desired';
const NON_DESIRED = 'non-desired';
const MEAN = 'Mean';
const STD = 'Std';
const T_STAT = 't-statistics';
export const P_VAL = 'p-value';
const MEAN_DES = `${MEAN}(${DESIRED})`;
const MEAN_NON_DES = `${MEAN}(${NON_DESIRED})`;
const STD_DES = `${STD}(${DESIRED})`;
const STD_NON_DES = `${STD}(${NON_DESIRED})`;

/** Map of statistic field names to their display titles  */
export const STAT_TO_TITLE_MAP = new Map([
  ['desAvg', MEAN_DES],
  ['desStd', STD_DES],
  ['nonDesAvg', MEAN_NON_DES],
  ['nonDesStd', STD_NON_DES],
  ['tstat', T_STAT],
  ['pValue', P_VAL],
]);

export const DESCR_TITLE = 'Descriptor';
export const DESCR_TABLE_TITLE = DESCR_TITLE + ' Statistics';
export const SELECTED_TITLE = 'Selected';
export const WEIGHT_TITLE = 'Weight';
export const SCORES_TITLE = 'pMPO score';
export const DESIRABILITY_COL_NAME = 'Desirability Curve';

/** Default p-value threshold for filtering descriptors */
export const P_VAL_TRES_DEFAULT = 0.01;

/** Minimum p-value threshold for filtering descriptors */
export const P_VAL_TRES_MIN = 0.01;

/** Maximum p-value threshold for filtering descriptors */
export const P_VAL_TRES_MAX = 1;

/** Default R-squared threshold for filtering correlated descriptors */
export const R2_DEFAULT = 0.53;

/** Minimum R-squared threshold for filtering correlated descriptors */
export const R2_MIN = 0.01;

/** Maximum R-squared threshold for filtering correlated descriptors */
export const R2_MAX = 1.0;

/** Default q-cutoff for descriptors in the pMPO model */
export const Q_CUTOFF_DEFAULT = 0.05;

/** Minimum q-cutoff for descriptors in the pMPO model */
export const Q_CUTOFF_MIN = 0.01;

/** Maximum q-cutoff for descriptors in the pMPO model */
export const Q_CUTOFF_MAX = 1;

/** Default setting for using sigmoid correction in pMPO */
export const USE_SIGMOID_DEFAULT = true;

export const FORMAT = '0.00';

/** Colors used for selected and skipped descriptors */
export enum COLORS {
  SELECTED = 'rgb(26, 146, 26)',
  SKIPPED = 'rgb(208, 57, 67)',
};

export const TINY = 1e-8;

/** Folder path for storing pMPO models */
export const FOLDER = 'System:AppData/Chem/mpo';

/** Desirability profile properties type */
export type DesirabilityProfileProperties = Record<string, {
  line: [number, number][],
  weight: number,
  min?: number,
  max?: number,
}>;

export const STAT_GRID_HEIGHT = 75;
export const DESIRABILITY_COLUMN_WIDTH = 305;

const POSITIVE_BASIC_RANGE_SIGMA_COEFFS = [0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5];

/** Basic range sigma coefficients for desirability profile points */
export const BASIC_RANGE_SIGMA_COEFFS = POSITIVE_BASIC_RANGE_SIGMA_COEFFS
  .slice(1)
  .map((v) => -v)
  .reverse()
  .concat(POSITIVE_BASIC_RANGE_SIGMA_COEFFS);

const EXTRA_RANGE_SIGMA_COEFFS = [0.12, 0.37, 0.63, 0.75, 0.88, 1.25, 1.75, 2.25, 2.75];
const EXTENDED_POSITIVE_RANGE_SIGMA_COEFFS = POSITIVE_BASIC_RANGE_SIGMA_COEFFS.concat(EXTRA_RANGE_SIGMA_COEFFS).sort();

/** Extended range sigma coefficients for desirability profile points */
export const EXTENDED_RANGE_SIGMA_COEFFS = EXTENDED_POSITIVE_RANGE_SIGMA_COEFFS
  .slice(1)
  .map((v) => -v)
  .reverse()
  .concat(EXTENDED_POSITIVE_RANGE_SIGMA_COEFFS);

/** Confusion matrix type */
export type ConfusionMatrix = {
  TP: number,
  TN: number,
  FP: number,
  FN: number,
};

// Titles for ROC curve columns
export const TPR_TITLE = 'TPR (Sensitivity)';
export const FPR_TITLE = 'FPR (1 - Specificity)';
export const THRESHOLD = 'Threshold';

// Number of points in ROC curve
const ROC_POINTS = 100;
export const ROC_TRESHOLDS_COUNT = ROC_POINTS + 1;

/** ROC curve thresholds from 0.0 to 1.0 */
export const ROC_TRESHOLDS = new Float32Array(Array.from({length: ROC_TRESHOLDS_COUNT}, (_, i) => i / ROC_POINTS));

/** Sample dataframe for pMPO training: https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/ */
export const SOURCE_PATH = 'System:AppData/Eda/drugs-props-train.csv';

/** Scores of the sample dataframe computed using https://github.com/Merck/pmpo */
export const SCORES_PATH = 'System:AppData/Eda/drugs-props-train-scores.csv';

/** Name of the synthetic drug used in the sample dataframe */
export const SYNTHETIC_DRUG_NAME = 'Synthetic drug';

/** pMPO model evaluation result type */
export type ModelEvaluationResult = {
  auc: number,
  threshold: number,
  tpr: Float32Array,
  fpr: Float32Array,
};

/** Maximum number of rows for which auto-tuning is applicable */
export const AUTO_TUNE_MAX_APPLICABLE_ROWS = 1000;

/** Minimum number of rows to show a warning about auto-tuning */
export const AUTO_TUNE_WARNING_MIN_ROWS = 3000;

/** Default settings for optimization in pMPO parameter tuning */
export const DEFAULT_OPTIMIZATION_SETTINGS = new Map<string, number>([
  ['tolerance', 0.0001],
  ['maxIter', 50],
  ['nonZeroParam', 0.0001],
  ['initialScale', 0.02],
  ['scaleReflaction', 1],
  ['scaleExpansion', 2],
  ['scaleContraction', -0.5],
]);

/** Optimal point type for pMPO parameter tuning */
export type OptimalPoint = {
  pValTresh: number,
  r2Tresh: number,
  qCutoff: number,
  success: boolean,
};

/** Minimum bounds for pMPO parameters during optimization */
export const LOW_PARAMS_BOUNDS = new Float32Array([P_VAL_TRES_MIN, R2_MIN, Q_CUTOFF_MIN]);

/** Maximum bounds for pMPO parameters during optimization */
export const HIGH_PARAMS_BOUNDS = new Float32Array([P_VAL_TRES_MAX, R2_MAX, Q_CUTOFF_MAX]);
