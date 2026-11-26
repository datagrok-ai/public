/*
 * The control-growth priors of the engine: files/engine_priors_v3.json, exported offline from the curated atlas
 * (aggregates only: files/SOURCES.md; v3 = v2's MC38, CT26, B16F10 verbatim plus 4T1, LLC, Panc02, Hepa1-6,
 * EMT6, MDA-MB-231, A549, PDX-NSCLC, PDX-PDAC, the same fields; R7 freeze: EMT6 in, SKOV3 out). Reduced Gompertz per
 * cell line:
 *   log V(t) = log V0 + (alpha / beta) (1 - exp(-beta t))
 * A new study's (log alpha, log beta) come from the labs' bivariate normal (study_prior); the MAP prior of log alpha
 * is the one computed offline (per-paper values are not distributed).
 */
import type {ReferenceRow} from './exposure';
import {_package} from './package';

export const PRIORS_FILE = 'engine_priors_v3.json';

export interface MixtureComponent {kind: string, w: number, a: number, b: number}

/** One fit of a line's growth spreads, as the fields it puts in the line's prior. */
export interface PriorFit {sd_within: {log_alpha: number, log_beta: number},
  study_prior: {sd_log_alpha: number, sd_log_beta: number}, distance_mean6: number, best?: boolean}

export interface LinePrior {
  papers: number;
  animals: number;
  log_alpha: {mixture: MixtureComponent[], map_mean: number, map_sd: number, map_q025: number, map_q975: number,
    tau: number, ess: number};
  log_beta: number;
  /** log of the volume (mm3) alpha refers to: d log V / dt = alpha - beta (log V - log_v_ref) */
  log_v_ref: number;
  sd_within: {log_alpha: number, log_beta: number, log_v0: number};
  /** v2: a new study's (log alpha, log beta), the labs' joint distribution (deconvolution of the atlas's controls,
   *  its size calibrated) */
  study_prior?: {mean_log_alpha: number, sd_log_alpha: number, mean_log_beta: number, sd_log_beta: number,
    rho: number, size_factor?: number};
  /** v2: the observation error, (c0 + c1 / sqrt(V)) times a Student t with df degrees of freedom (k already in c0 and
   *  c1), and a shift per study and day shared by its animals (sd_session; cancels between arms) */
  noise?: {c0: number, c1: number, k: number, df: number, sd_session: number};
  /** v2: the digital twin's own within-study model (Crux stats.tumor_twin), with
   *  study_prior, noise, log_v_ref and start_sd: between-animal sds, process noise, the first measurement's error
   *  factor, the start's t df, the slow animals' share and log alpha shift, the error's skew, and its validation */
  twin?: {sd_log_alpha: number, sd_log_beta: number, sd_process: number, first_factor: number, start_df: number,
    slow_share: number, slow_shift: number, skew: number,
    validation: {loo_coverage_90: number, measurements: number, below: number, above: number},
    /** The probability levels at which the predictive of an animal's last value, built after its session (Crux
     *  stats.tumor_monitor's last_flag), covers 90% and 99% of published controls, and how often a correct value
     *  fell below and above each on studies left out (an average per study; false: not confirmed). Not applied by
     *  the package: Crux takes one pair of levels for that flag and for the ranges shown before measuring
     *  (monitor.ts checkChance) */
    conformal?: {levels: {[level: string]: {lo: number, hi: number}}, papers: number, points: number,
      pooled: boolean, validated: {[level: string]: {miss_lo: number, miss_hi: number, cover: number}} | false,
      source: string},
    /** How often a correct value fell below and above the ranges Monitor shows before measuring ("90": next
     *  expected, "99": measure again), at Crux's own levels on the atlas's published controls: an average per study */
    next_ranges?: {levels: {[level: string]: {lo: number, hi: number}},
      observed: {[level: string]: {miss_lo: number, miss_hi: number}}, papers: number, points: number,
      pooled: boolean, source: string}};
  /** v3: `identification` on a line whose growth spreads the atlas's controls do not pin down: the fits that
   *  describe them equally well (calibrate_v2.py's starts within the rule), the one this prior carries marked */
  calibration?: {identification?: {status: string, rule: string, why: string, solutions?: PriorFit[]}};
  /** sd of log V between a study's animals at randomisation, measured in the atlas (measurement_error.py) */
  start_sd: number;
  sd_residual_log_v: number;
  take: {nonprogressing_map: number, nonprogressing: number, animals: number, definition: string};
  facet_shifts_log_rate: {facet: string, shift_log_lambda: number, se: number, ci95: number[], papers: number}[];
  map_settings: {[k: string]: string | number};
}

/** The caliper's measurement error, measured in the atlas: sd(log V) = c0 + c1 / sqrt(V mm3). */
export interface MeasurementError {model: string, c0: number, c1: number, triples: number, overall_sd: number}

export interface EnginePriors {
  version: string;
  measurement_error: MeasurementError;
  model: string;
  equation: string;
  provenance: {[k: string]: string};
  lines: {[line: string]: LinePrior};
  /** v3: every line's tumor model, for picking one (class, search by tumor type or organ) */
  models?: {[line: string]: TumorModel};
}

/** A line's tumor model in the atlas: its class (syngeneic: mouse line in immunocompetent mice; cdx: human cell line,
 *  pdx: patient tumor, both in immunodeficient mice), the tumor type and organ it models, the host, the atlas tier its
 *  priors come from (Full / Own spreads / Thin) and the papers and animals behind them. */
export interface TumorModel {
  class: 'syngeneic' | 'cdx' | 'pdx';
  tumour_type: string;
  organ: string;
  host_class: string;
  tier: string;
  papers: number;
  animals: number;
}

/** The atlas's tumor models by line, in the file's order. */
export async function tumorModels(): Promise<{[line: string]: TumorModel}> {
  return (await enginePriors()).models ?? {};
}

/** A package file's text, read once; a failed read is not kept, so the next call tries again. */
function cachedText(file: string): () => Promise<string> {
  let read: Promise<string> | null = null;
  return () => read ??= _package.files.readAsText(file).catch((e) => {
    read = null;
    console.error(`PMX Studio: ${file}`, e);
    throw new Error('PMX Studio\'s reference data could not be loaded. Check your connection and try again.');
  });
}

/** A package file's JSON, read once; a failed read is not kept, so the next call tries again. */
function cachedJson<T>(file: string, text = cachedText(file)): () => Promise<T> {
  let read: Promise<T> | null = null;
  return () => read ??= text().then((t) => JSON.parse(t) as T).catch((e) => {
    read = null;
    throw e;
  });
}

/** The atlas file's text: read once for its priors and its hash (plan-file.ts's provenance). */
export const priorsText = cachedText(PRIORS_FILE);
export const enginePriors = cachedJson<EnginePriors>(PRIORS_FILE, priorsText);

/** The key of a tumor model planned from its own history (history-models.ts), not from the atlas. */
export const HISTORY_PREFIX = 'new:pdx:';
/** The models planned from their history in this session, by key: their readable name and their growth prior. */
type HistoryModel = {name: string, prior: () => Promise<LinePrior>};
const historyModels = new Map<string, HistoryModel>();
/** Registers a model planned from its history (history-models.ts registerHistory): linePrior then gives its prior. */
export function setHistoryModel(key: string, name: string, prior: () => Promise<LinePrior>): void {
  historyModels.set(key, {name, prior});
}
/** The readable name of a model planned from its history; null for an atlas model or one not loaded. */
export const historyModelName = (line: string | null | undefined): string | null =>
  historyModels.get(line ?? '')?.name ?? null;
/** How the user sees a tumor model: its readable name for one planned from its history, else the line. */
export const modelName = (line: string): string => historyModels.get(line)?.name ?? line;
/** A tumor model's host class: the atlas's; a model planned from its history is a PDX, in immunodeficient mice;
 *  null when not known. */
export const hostClassOf = (line: string | null | undefined,
  models: {[line: string]: TumorModel} | null | undefined): string | null =>
  line?.startsWith(HISTORY_PREFIX) ? 'immunodeficient' : models?.[line ?? '']?.host_class ?? null;
/** The keys of the models planned from their history in this session. */
export const historyKeys = (): string[] => [...historyModels.keys()];

export async function linePrior(line: string): Promise<LinePrior> {
  if (line.startsWith(HISTORY_PREFIX)) {
    const h = historyModels.get(line);
    if (!h) {
      throw new Error('This tumor model comes from its own history, which is not loaded: open the plan file ' +
        '(.pmxplan) it was planned in.');
    }
    return h.prior();
  }
  const all = await enginePriors();
  const p = all.lines[line];
  if (!p) {
    throw new Error(`PMX Studio has no data for the tumor model "${line}": pick one of ` +
      `${Object.keys(all.lines).join(', ')}.`);
  }
  return p;
}

/** The PDX reference a model's own history is weighed against (model-history.ts): files/pdx_history_v1.json,
 *  exported offline from the atlas corpus (aggregates; its provenance inside). In ln of a control group's growth
 *  rate (1/day): the class mean and the sds between models, between studies of one model and between animals; the
 *  window volume (mm³) that rate refers to; the sd a different lab adds; the animals assumed for group means. */
export interface HistoryAnchor {mu: number, tau_model: number, tau_study: number, sigma_animal: number,
  window_volume_mm3: number, tau_lab: number, n_for_group_means: number, provenance: {[k: string]: string}}
export const pdxHistoryAnchor = cachedJson<HistoryAnchor>('pdx_history_v1.json');

/** A tumor's density, mg per mm³ of caliper volume: files/tumour_density_v0.json, exported offline from the atlas
 *  papers that weighed the excised tumors. Per line the
 *  papers' median, with the papers, animals and range behind it; spread_log the sds of log density between labs,
 *  between groups of one study and between animals (Crux stats.tumor_body_weight's lab_sd, group_sd, within_sd). */
export interface LineDensity {density: number, papers: number, animals: number, range: [number, number]}
export interface TumourDensity {lines: {[line: string]: LineDensity}, all: LineDensity,
  spread_log: {between_labs: number, between_groups: number, between_animals: number, caliper: number}}

export const tumourDensity = cachedJson<TumourDensity>('tumour_density_v0.json');

/** How much wider the engine's bands are drawn in the lines a validation measured them in:
 *  files/band_calibration_v1.json (its provenance inside). Per line, the factor of each half-width around the median
 *  in log scale; a line without an entry, or a factor null, keeps the engine's band. */
export interface BandFactors {end80: number | null, volume90: number | null}
export interface BandCalibration {version: number, provenance: string, lines: {[line: string]: BandFactors}}
export const bandCalibration = cachedJson<BandCalibration>('band_calibration_v1.json');

/** PMX Studio's reference table of patients' exposure: files/human_exposure_v1.json (exposure.ts's
 *  ReferenceRow; its provenance inside). */
export const humanReference = cachedJson<{version: string, provenance: string, drugs: ReferenceRow[]}>(
  'human_exposure_v1.json');

/** The line's factors, null where there are none or the file cannot be read: shared by Plan, Analyze and Monitor.
 *  A model planned from its history keeps the engine's bands (factor 1): what its validation measured. */
export const bandFactors = (line: string): Promise<BandFactors | null> => line.startsWith(HISTORY_PREFIX) ?
  Promise.resolve(null) : bandCalibration().then((c) => c.lines[line] ?? null).catch(() => null);

/** A quantile `x` moved `s` times as far from the median, in log scale (`floor`: the least any of them may be). */
export function widened(median: number, x: number, s: number, floor = 0): number {
  const [m, q] = [Math.max(median, floor), Math.max(x, floor)];
  return Math.max(floor, Math.exp(Math.log(m) + s * (Math.log(q) - Math.log(m))));
}
