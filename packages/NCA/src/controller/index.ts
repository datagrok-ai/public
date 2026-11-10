/**
 * Controller barrel — single entry point for downstream consumers (Phase 5
 * column-mapping dialog, Phase 6 Datagrok function, Phase 7 tests).
 *
 * @see ./nca-controller.ts
 */

export {runNcaPipeline, NcaController} from './nca-controller';
export type {
  NcaRunOutcome, PerProfileResult, RunNcaPipelineOptions, OverrideAction,
} from './nca-controller';
