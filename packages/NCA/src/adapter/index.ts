/**
 * Adapter layer — Phase 3 surface.
 *
 * Phase 4 adds: NcaController, audit-log wiring, route-from-string mapping
 * carry-over for the controller→writer pipeline.
 *
 * Phase 5 (review hygiene): re-export the sci-comp public unit/rules types
 * downstream consumers (views) need, so the views layer never imports
 * `@datagrok-libraries/sci-comp` directly — single boundary point.
 */

import type {nca as ncaNs} from '@datagrok-libraries/sci-comp';

export type TimeUnit = ncaNs.TimeUnit;
export type ConcentrationUnit = ncaNs.ConcentrationUnit;
export type DoseUnit = ncaNs.DoseUnit;
export type NcaRules = ncaNs.NcaRules;
export type RouteCode = ncaNs.RouteCode;
export type ParameterWarning = ncaNs.ParameterWarning;

export type {
  ColumnRawType,
  NumericRawArray,
  IntegerRawArray,
  BooleanRawArray,
  StringRawArray,
  ColumnSpec,
  ColumnContract,
  ValidationCode,
  ValidationError,
} from './column-contract';

export {
  DEFAULT_RULES, DEFAULT_BLQ_STRATEGY, withBlq, withLambdaZ,
  DEFAULT_ACCEPTANCE_RSQUARED, ACCEPTANCE_RSQUARED_MIN, ACCEPTANCE_RSQUARED_MAX,
} from './defaults';

// The sanctioned DG-column boundary helpers (nca-dg-native-missing-values F1/F2).
// `views/` must read a result column through one of the READ accessors, never
// through a bare `Column.get()` + `Number.isFinite` — `get()` returns the raw
// FLOAT_NULL (2.68e-34) even on a null cell, and that value is finite and positive.
//
// `f64Column` (the WRITE side) is exported too. An earlier revision withheld it on
// the premise that "only `adapter/` builds result frames" — the Arm A inter-rater
// pass (2026-08-03) falsified that: `nca-mapping-pane.ts`'s `appendTlastColumn`
// adds `TLAST` to the RESULTS frame, and `computeAccumulationRatios` reads it back
// as a result column. Withholding the helper is what let that site keep writing a
// raw NaN. Any module that builds a result column must be able to reach it.
export {numAt, finiteMinMax, toPlatformNull, log10OrNull, f64Column, ensureBoolBlq, DERIVED_BLQ_COL}
  from './adapter-utils';
export {scanNumericSentinels, describeSentinelEntries} from './adapter-utils';
export type {SentinelScanEntry} from './adapter-utils';

export {AuditLog, getAuditLog, setAuditLog} from './audit-log';
export type {AuditEntry, AuditStatus} from './audit-log';

export {
  validateObservations,
  validateUnits,
  validateScalarsForDose,
  validateScalarsForMolar,
  validateMultiAnalyte,
  validateLloqSanity,
  validateSparseDesign,
} from './validators';

export {buildProfileIndices} from './profile-index';
export type {ProfileIndex} from './profile-index';

export {extractProfileInputs, routeStringToCode, ROUTE_STRING_TO_CODE} from './profile-inputs';
export type {ProfileScalars} from './profile-inputs';

export {
  normalizeRouteString,
  normalizeRouteColumn,
  normalizeRouteToken,
  routeScalarColumnConflict,
  BARE_IV_DEFAULTED_TO_BOLUS,
  ROUTE_COLUMN_MAP_THRESHOLD,
  DERIVED_ROUTE_CODE_COL,
} from './route-normalize';
export type {RouteNormalization, RouteColumnResult} from './route-normalize';

export {computeBioavailability, isIntravascular} from './bioavailability';
export type {BioavailabilityOptions, BioavailabilityResult, BioavailabilityWarning} from './bioavailability';

export {
  computeAccumulationRatios,
  accumulationIndexTheoretical,
  computeExposureMargin,
  computeSexFold,
} from './toxicokinetics';
export type {
  TkWarning,
  TkWarningCode,
  AccumulationOptions,
  AccumulationResult,
  ExposureMarginOptions,
  ExposureMarginResult,
  SexFoldOptions,
  SexFoldResult,
} from './toxicokinetics';

export {
  computeDoseProportionality,
  shortDoseProportionalityNote,
  computeCriticalRegion,
  classifyAgainstCR,
  minNAdvisory,
  DEFAULT_THETA_LOW,
  DEFAULT_THETA_HIGH,
  DEFAULT_CI_LEVEL,
  DEFAULT_MIN_N_SEARCH_CAP,
  HETEROSCEDASTICITY_SD_RATIO,
  EXTREME_DOSE_DROP_SXX_RATIO,
  MIN_DOSE_LEVELS,
  MIN_ANIMALS_PER_LEVEL,
  SEX_INTERACTION_ALPHA,
} from './dose-proportionality';
export type {
  AucBasis,
  VerdictMode,
  DoseProportionalityVerdict,
  DoseProportionalityWarningCode,
  DoseProportionalityWarning,
  DoseProportionalityOptions,
  DoseProportionalityResult,
  DoseProportionalityFit,
  SexHomogeneity,
  DnAucAnova,
  MinNAdvisory,
} from './dose-proportionality';

export {
  impliedClPerKgWarning,
  CL_PLAUSIBILITY_FLOOR_LHKG,
  CL_PLAUSIBILITY_CEILING_LHKG,
} from './cl-plausibility';

export {
  computeSparseNca, shortSparseNote, isSparseDesign, computeSparseForRecipe, sparseCompositeUnits,
  computeSparseStrata, sparseCompositeKey, sparseDoseGroups, SPARSE_INFUSION_BLOCKED_REASON,
} from './sparse-nca';
export type {
  SparseComputeMapping, SparseUnits, SparseStrataSpec, SparseStrata, SparseStratumCell, SparseBlockedStratum,
  SparseDoseGroup,
} from './sparse-nca';
export type {
  AucCiMethod,
  SparseNcaOptions,
  SparseNcaWarningCode,
  SparseNcaWarning,
  SparseBootstrapCi,
  SparseGroupResult,
  SparseNcaResult,
} from './sparse-nca';

// λz acceptance predicate — pure core of the ACCEPTANCE verdict (nca-workflow-stages AD-4).
export {
  evaluateAcceptance, criterionReasonText, fmtThreshold, SPAN_RATIO_POOR,
  ACCEPT_ACCEPTED, ACCEPT_FLAGGED, ACCEPT_NA,
  criterionShortLabel, ACCEPT_REASONS_COL, ACCEPTANCE_COL, FLAG_COLS,
  REASONS_SEPARATOR, reasonsCellValue, tagReasonPillsColumn,
  ACCEPTANCE_FLAG_CODES,
} from './acceptance';
export type {
  AcceptanceReading, AcceptanceInputs, AcceptanceThresholds, AcceptanceVerdict,
  AcceptanceFlagCode, FlagState,
} from './acceptance';

// nca-workflow-stages: exclusion feed, review units, reasons, fingerprints, stage (S1).
export {
  LZ_DEPENDENT_PARAMS, isLzDependentParam, SPARSE_LZ_DEPENDENT_COLS, isSparseLzDependentCol,
  excludedCompositeRows, reportableResults, REVIEW_COL, REVIEW_EXCLUDED,
  LZ_REVIEW_INFO_COL, REVIEW_TO_REVIEW, reviewStateOf, upsertStringColumn, writeReviewColumns,
} from './reportable';
export {buildSummaryCache, buildCurvesByCompound, countDoseClusters, CROSS_DOSE_POOL_CODE} from './summary-cache';
export type {ObservationsCtx, CurvesByCompound} from './summary-cache';
export type {ReportableResults, ReviewInfo} from './reportable';
export {
  classifyReviewUnits, pendingUnits, profileMeasurability, sparseUnitKey, SPARSE_UNIT_PREFIX, sharedFlagCodes,
  denseProfilesViable,
} from './review-units';
export type {
  ReviewUnit, ReviewNeed, UnitCategory, ProfileUnitInput, CompositeUnitInput, ClassifyInput,
  Measurability,
} from './review-units';
export {
  judgementReasons, criterionReasons, formatReasons, reasonsValid, ACTIONS_BY_CATEGORY,
  OTHER_REASON_CODE, CRITERION_REASON_PREFIX, REVIEW_ACTION_LABEL, REVIEW_STATE_LABEL,
} from './review-reasons';
export type {ReviewCategory, ReasonDef, ReasonExportCategory} from './review-reasons';
export {
  computeRecipeFingerprint, computeDecisionFingerprint, computeRecipeHash, RECIPE_HASH_VERSION,
} from './fingerprint';
export type {DecisionFingerprintInput} from './fingerprint';
export {
  deriveStage, formatStageBadge, formatNotReviewed, formatStudyAnalysesSummary, formatNextHint,
  formatPendingCount, STAGE_LABEL, pipelineFlagCount, formatPipelineFlags, PER_ANIMAL_SPARSE_NOTE,
  formatSupersedes, formatSupersededBy,
} from './workflow-stage';
export type {StageInput, StageResult, StageUnit} from './workflow-stage';

export {writeResultsToDataFrame, patchResultRow, ROUTE_CODE_TO_STRING, FIXED_RESULT_COLUMNS} from './results-writer';
export type {ResultRowValues} from './results-writer';

// λz manual point-selection (VR-11): the service-column contract shared by the
// four exclusion surfaces (AC-3.1) + the pure translation/floor core (P-B1c).
export {
  isLzServiceColumn, formatOverrideTimestamp,
  LZ_SERVICE_COLUMNS, LZ_ACTION_COL, LZ_USER_COL, LZ_TIME_COL,
  buildAugmentedProfile, buildManualPoints, pointsUsedToRows, subjectHasManualOverride,
  applyMinPointsFloor, LAMBDAZ_BELOW_MIN_POINTS,
  isLambdaZNonContiguous, LAMBDAZ_NONCONTIGUOUS,
} from './lambdaz-manual-points';
export type {LzAction, AugmentedProfile} from './lambdaz-manual-points';

export {
  aggregateNcaResults, resolveAggregateGroupKeys, DEFAULT_AGG_PARAMS,
  blankMixedBasisFamily, IV_EV_INCOMPARABLE_FAMILY,
} from './results-aggregator';
export type {AggregateOptions} from './results-aggregator';

export {compositeGroupKey, KEY_SEP} from './composite-group-key';

export {appendGroupingColumns, groupByCandidates} from './grouping-columns';
export {passthroughCandidate, recipeRoleNames, isConstantWithinSubject} from './passthrough-gate';

export {buildMeanSdProfile} from './mean-sd-profile';
export type {MeanSdProfileSpec} from './mean-sd-profile';

export {
  detectMapping, isBoolLikeNumeric, detectAnalyteColumns, detectEventDrivenColumns,
  isCompoundIdColumnName, detectMixedCompoundGroups, mixedCompoundMessage, findColByPattern,
} from './column-autodetect';
export type {MixedCompoundGroup} from './column-autodetect';
export type {DetectedMapping, AnalyteDetection, EventDrivenDetection} from './column-autodetect';

export {
  normalizeEventDrivenObservations,
  expandDoseLadder,
  assertDoseLadder,
  SUBJECT_COL,
  TAD_COL,
  CONC_COL,
  DOSE_COL,
  OCCASION_COL,
  DERIVED_INFUSION_DURATION_COL,
  NONMEM_SUBJECT_FIRST_ROW_NOT_DOSE,
  NONMEM_OCCASION_KEY_COLLISION,
  NONMEM_INFUSION_RATE_INFERRED,
  NONMEM_ROUTE_INFERRED,
  EVID_RESET_NEW_EPISODE,
  EVID_OTHER_IGNORED,
  NONMEM_STEADY_STATE_ASSUMED,
  sparseDoseSource,
} from './event-driven';
export type {
  EventContract,
  NormalizeEventOptions,
  NormalizeEventResult,
  InferredRoute,
  ClassificationCounts,
} from './event-driven';

export {poolCompositeByCell, STUDY_DAY_COL, CELL_BW_COL} from './composite-pool';
export type {CompositePoolColumns, CompositePoolResult} from './composite-pool';

export {
  defaultSexAssignment, uniqueSexValues, applySexMap, sexMapChangesValues, SEX_SYNONYMS,
} from './sex-mapping';
export type {SexCode} from './sex-mapping';

export {fnv1a, recipeOf, makeAnalysisEntry} from './study-schema';
export {INSTALLED_ENGINE_VERSION} from './engine-version';
export type {
  SourceMode, RecipeColumns, SamplingType, RecipeModes, RecipeUnits, AnalysisRecipe,
  SourceKind, SourcePointer, SignatureKind, SourceSignature,
  MeanCv, CurveSeries, CompoundCurve, CompoundWarning, GroupKeyValue,
  KeyValueCount, KeyInventoryEntry, PassthroughColumnSpec, PassthroughEntry,
  AnalysisTierSummary, CompoundTierRow, SummaryCache,
  AnalysisEntry, StudyAccessFlags, NcaStudy,
  WorkflowStage, StageCache, ConfirmedRecipe, LzOverride, ReviewAction, ReviewReason,
  SupersededBy, ReviewDecisionRecord, ReviewLedger, AnalysisOrigin,
  StandaloneOrigin, PipelineOrigin, PipelineFlag, MakeAnalysisEntryArgs,
} from './study-schema';

export {
  buildStudyGrid, studyGridRowTiers, normalizeGroupKey, dedupGroupValuesByCanonical,
  buildCurveFrame, buildOverlayCurveFrame,
  ANALYSIS_COL, NO_GROUP, ROUTE_COL, CURVE_COL, FLAGS_COL, ANALYSIS_ID_COL, COMPOUND_ID_COL,
  PARAM_TAG, STAT_TAG, PASSTHROUGH_TAG, PASSTHROUGH_STRUCT, PASSTHROUGH_TEXT, RECOMPUTE_FLAG, STALE_ENGINE_FLAG,
  STUDY_EXCLUDED_N_TEMP, studyCellExclusion,
  BASIS_DEPENDENT_PARAMS, STUDY_BASIS_STALE_TEMP, studyCellBasisStale, STUDY_CURVE_PENDING_TEMP, studyCurvePending,
} from './study-grid';

export {
  STANDALONE_ROUTE_SEGMENT, NCA_APP_ROOT, NCA_STUDIES_CHANGED_EVENT, newId, uniqueLabel, isStale, isEntryStale,
  isBasisStale,
  computeSnapshotHash, snapshotSignature, snapshotDiff, SourceReplacedError, ContentChangedError,
  KeepsChangingError, SnapshotFidelityError, getNcaRepository, setNcaRepository,
  registerNcaRepositoryFactory, isCurvePending, pendingCurveIds, withCurves,
} from './nca-repository';
export type {
  NcaRepository, StudyAccess, AnalysisHeader, StudyHeader, LoadedStudy,
  LoadedAnalysis, WorkflowState, NewAnalysisInput, SaveAnalysisInput, SummaryWrite, MoveTarget, MoveResult,
  PipelineBase, PipelineWrite, WrittenSnapshot,
} from './nca-repository';
export {EmsNcaRepository} from './ems-nca-repository';

export {
  scaleProfileInputs,
  CANONICAL_UNIT_TAGS,
  TIME_TO_HOUR,
  CONC_TO_MG_L,
  DOSE_TO_MG,
  BW_TO_KG,
  convertConcentrationFactor,
  BW_HEURISTIC_THRESHOLD_KG,
  INFUSION_DURATION_HEURISTIC_MARGIN,
  CONC_PLAUSIBILITY_CEILING_MG_L,
  CONC_GUARD_MIN_POSITIVE_COUNT,
  normaliseClVolumes,
  compositeDoseWeight,
  isPerKgDose,
  NO_BODY_WEIGHT,
  L_PER_H_TO_ML_PER_MIN,
  RESULTS_BASIS,
  CL_PER_KG_UNIT,
  VOLUME_PER_KG_UNIT,
  CL_ABS_UNIT,
  VOLUME_ABS_UNIT,
  BW_UNIT,
} from './unit-conversion';
export type {
  BodyWeightUnit, ScaledProfileInputs, UnitScalars, BodyWeightContext, EngineClVolumes, NormalisedClVolumes,
  CompositeMember, CompositeDoseWeight,
} from './unit-conversion';
