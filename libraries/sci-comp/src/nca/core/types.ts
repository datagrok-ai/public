/**
 * Public type contracts for the NCA computation core.
 *
 * Only types that are part of the pure-math contract live here. Types tied
 * to the Datagrok DataFrame layer (`ColumnContract`, `ValidationError`,
 * `ProfileIndex`, ...) live in the `packages/NCA/` package, not in this
 * library. See `docs/nca_core_interface_v2_1.md` for the architectural
 * rationale.
 */

// ──────────────────────────────────────────────────────────────────────────
// Units
// ──────────────────────────────────────────────────────────────────────────

/** Concentration unit string literal. */
export type ConcentrationUnit =
  | 'ng/mL' | 'µg/mL' | 'µg/L' | 'mg/L' | 'µmol/L' | 'nmol/L';

/** Time unit string literal. */
export type TimeUnit = 'h' | 'min' | 's' | 'd';

/** Dose unit string literal. */
export type DoseUnit = 'mg' | 'µg' | 'ng' | 'mg/kg' | 'µg/kg';

/** Volume unit string literal (terminal volume of distribution). */
export type VolumeUnit = 'L' | 'mL' | 'L/kg' | 'mL/kg';

/** Clearance unit string literal. */
export type ClearanceUnit = 'L/h' | 'mL/min' | 'mL/min/kg' | 'L/h/kg';

// ──────────────────────────────────────────────────────────────────────────
// Routes
// ──────────────────────────────────────────────────────────────────────────

/** Human-readable route label. Used in ProfileKey / metadata. */
export type Route = 'IV-bolus' | 'IV-infusion' | 'PO' | 'SC' | 'IM' | 'other';

/**
 * Integer encoding of {@link Route} as carried inside the core. The mapping
 * between string and code is owned by the adapter; the core only sees codes.
 *
 * Numerical assignments (used by `computeNca`):
 *   0 = IV-bolus, 1 = IV-infusion, 2 = PO, 3 = SC, 4 = IM, 5 = other.
 */
export type RouteCode = 0 | 1 | 2 | 3 | 4 | 5;

/** Numeric route constants for use inside the core. */
export const ROUTE_IV_BOLUS: RouteCode = 0;
export const ROUTE_IV_INFUSION: RouteCode = 1;
export const ROUTE_PO: RouteCode = 2;
export const ROUTE_SC: RouteCode = 3;
export const ROUTE_IM: RouteCode = 4;
export const ROUTE_OTHER: RouteCode = 5;

// ──────────────────────────────────────────────────────────────────────────
// Profile inputs — what the core receives for a single subject's profile
// ──────────────────────────────────────────────────────────────────────────

/**
 * One profile's data plus its dosing context. The arrays are sorted by time.
 * The adapter is responsible for copying the right slice out of a DataFrame
 * into freshly allocated `Float64Array` / `Uint8Array` instances.
 */
export interface ProfileInputs {
  /** Time points in the units of {@link timeUnits}, sorted ascending. */
  readonly time: Float64Array;
  /** Observed concentrations in the units of {@link concentrationUnits}. */
  readonly conc: Float64Array;
  /** 1 = below LOQ at this time point, 0 = measurable. Same length as `time`. */
  readonly blqMask: Uint8Array;
  /** Lower limit of quantification — per-row array, or a scalar applied to all rows. */
  readonly lloq: Float64Array | number;

  /** Administered dose in the units of {@link doseUnits}. */
  readonly dose: number;
  readonly doseUnits: DoseUnit;
  readonly concentrationUnits: ConcentrationUnit;
  readonly timeUnits: TimeUnit;
  readonly route: RouteCode;
  /** Infusion duration in {@link timeUnits} for IV-infusion route, otherwise null. */
  readonly infusionDuration: number | null;

  /** Body weight in kg, when available. Used for weight-normalised parameters. */
  readonly bodyWeight: number | null;
}

// ──────────────────────────────────────────────────────────────────────────
// Rules — solver configuration for one profile
// ──────────────────────────────────────────────────────────────────────────

/** Trapezoidal AUC integration scheme. */
export type AucMethod = 'linear' | 'log-linear' | 'linear-up-log-down';

/** Action taken on a single BLQ data point during pre-processing. */
export type BlqRule = 'set-zero' | 'set-half-lloq' | 'exclude' | 'missing';

/**
 * BLQ-handling configuration broken out by phase of the concentration profile.
 * Conventions follow PKNCA (see `scripts/generate-nca-fixtures.R`).
 */
export interface BlqStrategy {
  /** BLQ points before the first measurable observation. */
  readonly preFirstMeasurable: BlqRule;
  /** BLQ points sitting between two measurable observations. */
  readonly embedded: BlqRule;
  /** First BLQ point after the last measurable observation. */
  readonly afterLast: BlqRule;
  /** Second and subsequent consecutive BLQ points after the last measurable. */
  readonly consecutiveAfterLast: BlqRule;
}

/**
 * Strategy for selecting the terminal-phase data points for the lambda_z
 * (terminal slope) regression.
 */
export interface LambdaZStrategy {
  readonly mode: 'auto-best-fit' | 'manual-points' | 'manual-time-range';
  /** Minimum number of points to accept a candidate fit. */
  readonly minPoints: number;
  /** Minimum adjusted R² required to accept the auto best-fit result. */
  readonly minRSquared: number;
  /** When true, the Cmax point is excluded from the lambda_z candidate window. */
  readonly excludeCmax: boolean;
  /**
   * PKNCA-style tie-breaking factor for adj-R². When two candidate subsets
   * have adj-R² within this distance, the one with more points wins.
   * Implemented as an additive score `adjRSquared + adjRSquaredFactor * n`.
   * PKNCA default is `1e-4`. Set to `0` to disable tie-breaking and use
   * strict max adj-R² (then the smaller subset wins on ties).
   */
  readonly adjRSquaredFactor?: number;
  /** Indices of selected points (manual-points mode). */
  readonly manualPoints?: Int32Array;
  /** Inclusive index range of selected points (manual-time-range mode). */
  readonly manualRange?: { readonly startIdx: number; readonly endIdx: number };
}

/**
 * Effective rules for one profile. The hierarchical inheritance
 * (study → group → profile) is resolved by the adapter; the core only sees
 * the resolved object.
 */
export interface NcaRules {
  readonly aucMethod: AucMethod;
  readonly blq: BlqStrategy;
  readonly lambdaZ: LambdaZStrategy;
  /** Soft warning threshold for % AUC extrapolated. */
  readonly extrapWarnPct: number;
  /** Hard error threshold for % AUC extrapolated. */
  readonly extrapErrorPct: number;
  /** When true, AUC sums use Neumaier-compensated summation. */
  readonly compensatedSummation: boolean;
}

// ──────────────────────────────────────────────────────────────────────────
// Sub-results returned by individual core functions
// ──────────────────────────────────────────────────────────────────────────

/** Outcome of {@link findCmax}. */
export interface CmaxResult {
  readonly cmax: number;
  readonly tmax: number;
  /** Index of the Cmax point in the post-BLQ time array. */
  readonly cmaxIdx: number;
}

/** Outcome of {@link lambdaZBestFit} or {@link lambdaZManual}. */
export interface LambdaZResult {
  /** Terminal-phase rate constant (1/time-unit). */
  readonly lambdaZ: number;
  /** Intercept of the log-linear fit, in concentration-unit space. */
  readonly intercept: number;
  readonly rSquared: number;
  readonly adjRSquared: number;
  /** Indices of points used for the fit (in the original time array). */
  readonly pointsUsed: Int32Array;
  readonly tStart: number;
  readonly tEnd: number;
}

/** Outcome of {@link applyBlqStrategy}. */
export interface BlqProcessingResult {
  /** Concentration array after BLQ rules have been applied. */
  readonly conc: Float64Array;
  /** Indices of points fully excluded from downstream computation. */
  readonly excluded: Int32Array;
}

// ──────────────────────────────────────────────────────────────────────────
// Final compute result
// ──────────────────────────────────────────────────────────────────────────

/**
 * Numeric NCA parameters for one profile.
 *
 * `cl` and `vz` follow the PKNCA universal-key convention: their physical
 * interpretation (`CL` vs `CL/F`, `Vz` vs `Vz/F`) follows from the
 * {@link ProfileInputs.route} of the profile, not from the field name.
 */
export interface ParameterValues {
  readonly cmax: number;
  readonly tmax: number;
  readonly aucLast: number;
  readonly aucInf: number;
  readonly pctExtrap: number;
  readonly lambdaZ: number;
  readonly halfLife: number;
  readonly cl: number;
  readonly vz: number;
}

/**
 * Diagnostic warning emitted by `computeNca` when a parameter sits in a
 * dubious regime (high extrapolation, low R², high BLQ fraction, etc.).
 */
export interface ParameterWarning {
  readonly code:
    | 'AUC_EXTRAP_HIGH'
    | 'LAMBDAZ_LOW_R2'
    | 'LAMBDAZ_FEW_POINTS'
    | 'BLQ_HIGH_FRACTION'
    | string;
  readonly severity: 'info' | 'warning' | 'error';
  readonly message: string;
}

/**
 * Provenance — everything an auditor needs to understand how the parameter
 * values were obtained from the inputs.
 */
export interface ProfileProvenance {
  readonly lambdaZ: LambdaZResult | null;
  readonly blqApplied: BlqProcessingResult;
  readonly aucMethod: AucMethod;
  readonly compensated: boolean;
  readonly warnings: ReadonlyArray<ParameterWarning>;
}

/** Top-level result of {@link computeNca}. */
export interface ComputeResult {
  readonly values: ParameterValues;
  readonly provenance: ProfileProvenance;
  /**
   * - `'ok'`      — every requested parameter computed successfully.
   * - `'partial'` — Cmax/AUClast available but lambda_z (and dependents:
   *                 AUCinf, t½, CL, Vz) could not be obtained.
   * - `'failed'`  — Cmax not available (e.g. all observations BLQ).
   */
  readonly status: 'ok' | 'partial' | 'failed';
}
