/**
 * Public entry point for the stats module.
 *
 * Re-exports a clean TypeScript API. jstat is used internally as the
 * mathematical backend (special functions, distributions, linear algebra)
 * and is not exposed here.
 */

export * from './types';

// Distribution wrappers (typed)
export * as distributions from './distributions';

// Tests
export {welchTTest} from './tests/welch-t';
export {mannWhitneyU} from './tests/mann-whitney';
export {spearman, severityTrend} from './tests/spearman';
export {hedgesG} from './tests/hedges-g';
export {fisherExact2x2} from './tests/fisher-exact';
export {jonckheere} from './tests/jonckheere';
export {welchPairwise} from './tests/welch-pairwise';
export type {TreatedGroup as WelchTreatedGroup, WelchPairwiseResult} from './tests/welch-pairwise';
export {dunnettPairwise} from './tests/dunnett';
export type {TreatedGroup as DunnettTreatedGroup, DunnettPairwise} from './tests/dunnett';
export {
  cochranArmitage,
  cochranArmitageBasic,
  thresholdTest,
} from './tests/cochran-armitage';
export type {
  CASettings, CAResult, Variance,
  ThresholdSettings, ThresholdStep,
} from './tests/cochran-armitage';
export {
  pavaIncreasing, pavaDecreasing, williamsTest,
} from './tests/williams';
export type {
  WilliamsSettings, WilliamsStepResult, WilliamsResult, WilliamsDirection,
} from './tests/williams';
export {lookup1971, lookup1972} from './tests/williams-tables';
export {runAncova} from './tests/ancova';
export type {
  AncovaSettings, AncovaResult,
  AdjustedMean, AncovaPairwise, AncovaSlope,
  AncovaSlopeHomogeneity, AncovaEffectDecomposition,
} from './tests/ancova';

// Multiple comparison
export {bonferroniCorrect} from './multiple-comparison/bonferroni';
