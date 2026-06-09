/**
 * Public entry point for the NCA computation core.
 *
 * Re-exports the type contract and the implemented core functions.
 * Concrete numerical functions are added incrementally by Phase 1 tasks.
 */

export * from './types';
export {mulberry32, deriveWorkerSeeds} from './prng';
export {applyBlqStrategy} from './blq';
export {findCmax} from './cmax';
export {lambdaZBestFit, lambdaZManual} from './lambda-z';
export {estimateC0, insertC0, C0_DEFAULT_METHODS} from './c0';
export type {C0Method, C0Options} from './c0';
export {
  halfLifeFromLambdaZ,
  clearance,
  volumeTerminal,
  pctExtrapolated,
  meanResidenceTime,
  volumeSteadyState,
  pctExtrapolatedAumc,
  tlag,
} from './derived';
export {computeNca} from './compute-nca';
export {
  aucLinearNaive,
  aucLogLinearNaive,
  aucLinearUpLogDownNaive,
  aucLinearCompensated,
  aucLogLinearCompensated,
  aucLinearUpLogDownCompensated,
  aucExtrapolateToInfinity,
  neumaierSum,
} from './auc';
export {
  aumcLinearNaive,
  aumcLogLinearNaive,
  aumcLinearUpLogDownNaive,
  aumcLinearCompensated,
  aumcLogLinearCompensated,
  aumcLinearUpLogDownCompensated,
  aumcExtrapolateToInfinity,
} from './aumc';
