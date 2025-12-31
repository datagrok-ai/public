// Probabilistic MPO definitions & constants
export const MIN_SAMPLES_COUNT = 10;
export const PMPO_NON_APPLICABLE = 'pMPO is not applicable';
export const PMPO_COMPUTE_FAILED = 'Failed to compute pMPO parameters';

export type BasicStats = {
  desAvg: number,
  desStd: number,
  nonDesAvg: number,
  nonDesStd: number,
};

export type DescriptorStatistics = BasicStats & {
  desLen: number,
  nonSesLen: number,
  tstat: number,
  pValue: number,
};

export type Cutoff = {
  cutoff: number,
  cutoffDesired: number,
  cutoffNotDesired: number,
};

export type SigmodParams = {
  pX0: number,
  b: number,
  c: number,
};

export type PmpoParams = BasicStats & Cutoff & SigmodParams & {
  zScore: number,
  weight: number,
  intersections: number[],
  x0: number,
  xBound: number,
};

const DESIRED = 'desired';
const NON_DESIRED = 'non-desired';
const MEAN = 'mean';
const STD = 'std';
const T_STAT = 't-statistics';
export const P_VAL = 'p-value';
const Z_SCORE = 'z-score';
const MEAN_DES = `${MEAN}(${DESIRED})`;
const MEAN_NON_DES = `${MEAN}(${NON_DESIRED})`;
const STD_DES = `${STD}(${DESIRED})`;
const STD_NON_DES = `${STD}(${NON_DESIRED})`;

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
export const WEIGHT_TITLE = 'Weights';
export const WEIGHT_TABLE_TITLE = DESCR_TITLE + ' ' + WEIGHT_TITLE;
export const SCORES_TITLE = 'pMPO score';

export const P_VAL_TRES_MIN = 0.01;
export const R2_MIN = 0.01;
export const Q_CUTOFF_MIN = 0.01;

export enum COLORS {
  SELECTED = '#84d984ff',
  SKIPPED = '#d9979bff',
};

export const TINY = 1e-8;

export const FOLDER = 'System:AppData/Chem/mpo';

export type DesirabilityProfileProperties = Record<string, {
  line: [number, number][],
  weight: number,
  min?: number,
  max?: number,
}>;

export const SCORE_SCALE = 100;
