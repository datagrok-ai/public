// Thin re-export shim for backwards-compat with files that still import
// from `./utils/prolif-panel`. The actual implementation lives in
// `@datagrok-libraries/bio/src/prolif/prolif-panel` — moved there per
// PR review (https://github.com/datagrok-ai/public/pull/3783, comment #11).
export {
  hasNonWaterHetatm,
  detectNonWaterHetatmInstances,
  makeProlifWidget,
  renderInteractionBreakdown,
  interactionsColForDiagram,
  getPlHtmlForRow,
  runPlBatch,
  PL_DIAGRAM_SEM_TYPE,
  PROLIF_SKIP_RESNAMES,
  type ProlifBatchCtx,
  type PlScriptArgs,
  type RunPlBatchOptions,
} from '@datagrok-libraries/bio/src/prolif/prolif-panel';
