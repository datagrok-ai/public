// Probabilistic MPO definitions & constants
export const MIN_SAMPLES_COUNT = 10;
export const PMPO_NON_APPLICABLE = 'pMPO is not applicable';

export type DescriptorStatistics = {
  desAvg: number,
  desStd: number,
  desLen: number,
  nonDesAvg: number,
  nonDesStd: number,
  nonSesLen: number,
  tstat: number,
  pValue: number,
  zScore: number,
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
  ['zScore', Z_SCORE],
]);

export const DESCR_TITLE = 'Descriptor';
export const DESCR_TABLE_TITLE = DESCR_TITLE + ' Statistics';
export const SELECTED = 'Selected';

export const P_VAL_TRES_MIN = 1e-8;
export const R2_MIN = 0.01;

export enum COLORS {
  SELECTED = '#84d984ff',
  SKIPPED = '#d9979bff',
};
