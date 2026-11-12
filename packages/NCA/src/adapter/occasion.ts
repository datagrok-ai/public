/**
 * Occasion (study day / period), occasion design and animal resolution — the ONE place
 * every derived analysis (dose proportionality, F%, sex comparison, accumulation, the study
 * grid's %F, the Pipeline) learns which column is the occasion and which profiles belong to
 * one animal. Topic tk-study-day-stratification, F1 (F1a occasion + design, F1b animal).
 *
 * Before this module there were two name detectors that had drifted apart (the view's
 * `/^study_?day$/i` and the Pipeline's `period|visitdy|…`) and neither matched NONMEM `DAY`
 * (guardrail "one function, every caller"). Both are replaced by {@link isOccasionName} /
 * {@link detectOccasionColumn}.
 *
 * Nothing here changes a per-profile NCA value: the occasion and the animal are read from a
 * frame that already holds one row per profile (the results frame, or the Pipeline's result
 * keys), so the profile key is untouched.
 *
 * Scientific contract (research tk-study-day-stratification, ICH S3A §2, §4.3, Note 3):
 * - The occasion is a reporting stratum: DP, the sex comparison and accumulation are reported
 *   per occasion whenever an occasion with ≥ 2 values resolves, whatever the design, except
 *   `escalation` (one dose per occasion — no within-occasion dose range).
 * - Only F% pairing depends on the design: `repeat-dose` → same-occasion pairing; `crossover`
 *   → across periods (today's behaviour); `unknown` → today's behaviour plus a warning.
 * - A proposed design never changes F% pairing until the analyst confirms it (R2 F5: a
 *   within-subject DDI study on "Day 1 / Day 8" satisfies the repeat-dose evidence too).
 */

import type * as DG from 'datagrok-api/dg';

import {f64Column, numAt} from './adapter-utils';

/** Results column holding AUC0–τ (time after dose). */
export const AUC_TAU_COL = 'AUC_TAU';
/** TLAST tolerance (h) for the AUClast ≡ AUC0–τ identity — the TK aggregator's own default. */
export const AUC_TAU_TLAST_TOL_H = 0.01;

/**
 * Write `AUC_TAU` (AUC over [0, τ], time after dose) on a results frame. Where the profile's
 * last measurable sample is at τ (TLAST within {@link AUC_TAU_TLAST_TOL_H}) AUC0–τ IS the
 * engine's AUClast — no new arithmetic. Otherwise the value comes from `intervalAuc(row)`
 * (sci-comp `nca.aucInterval`, F4b; PKNCA `aucint.last`) when supplied and the profile was
 * SAMPLED through τ: TLAST > τ, or `lastSampleH(row)` (the last sampling time, BLQ samples
 * included) ≥ τ — a BLQ trough at τ leaves TLAST < τ although τ was spanned, and the area to τ
 * then ends at TLAST. Sampling that stops before τ → missing (τ not spanned → the labelled
 * AUClast fallback). τ unknown (NaN) → no column is written and the function returns false.
 */
export function appendAucTauColumn(
  results: DG.DataFrame, tauH: number, intervalAuc?: (row: number) => number,
  lastSampleH?: (row: number) => number,
): boolean {
  if (!(tauH > 0)) {
    if (results.col(AUC_TAU_COL) !== null) results.columns.remove(AUC_TAU_COL);
    return false;
  }
  const last = results.col('AUC_LAST');
  const tlast = results.col('TLAST');
  const out = new Float64Array(results.rowCount);
  for (let i = 0; i < results.rowCount; i++) {
    const tl = tlast !== null ? numAt(tlast, i) : NaN;
    const sampledThroughTau = (Number.isFinite(tl) && tl > tauH) ||
      (lastSampleH !== undefined && lastSampleH(i) >= tauH - AUC_TAU_TLAST_TOL_H);
    if (Number.isFinite(tl) && Math.abs(tl - tauH) <= AUC_TAU_TLAST_TOL_H)
      out[i] = last !== null ? numAt(last, i) : NaN;
    else if (sampledThroughTau && intervalAuc !== undefined) out[i] = intervalAuc(i);
    else out[i] = NaN;
  }
  const col = f64Column(AUC_TAU_COL, out);
  const lastUnits = last?.getTag('units');
  if (lastUnits !== null && lastUnits !== undefined && lastUnits !== '') col.setTag('units', lastUnits);
  col.setTag('label', 'AUC0–τ');
  if (results.col(AUC_TAU_COL) !== null) results.columns.remove(AUC_TAU_COL);
  results.columns.add(col);
  return true;
}

/** The occasion names (every name in exactly one class), each class in detection-priority
 *  order (the first match wins). */
const DAY_LIKE_ORDER = ['studyday', 'day', 'nomdy', 'visitdy', 'visitday', 'week'];
const PERIOD_LIKE_ORDER = ['visit', 'period', 'occasion', 'occ'];
const DAY_LIKE_NAMES: ReadonlySet<string> = new Set(DAY_LIKE_ORDER);
const PERIOD_LIKE_NAMES: ReadonlySet<string> = new Set(PERIOD_LIKE_ORDER);

/** Column-name normalisation shared with the Pipeline mapping (case, spaces, `_`, `-`). */
export const normOccasionName = (name: string): string => name.toLowerCase().replace(/[\s_-]+/g, '');

/** True for a column named like an occasion (study day, week, visit, period, occasion). */
export function isOccasionName(name: string): boolean {
  const n = normOccasionName(name);
  return DAY_LIKE_NAMES.has(n) || PERIOD_LIKE_NAMES.has(n);
}

/** True for a day-like occasion name (`StudyDay`, `DAY`, `Week`, …): repeat-dose by default. */
export function isDayLikeOccasionName(name: string): boolean {
  return DAY_LIKE_NAMES.has(normOccasionName(name));
}

/** Detection priority: day-like names before period-like ones. */
const NAME_PRIORITY = [...DAY_LIKE_ORDER, ...PERIOD_LIKE_ORDER];

/**
 * The occasion column among `profileConstant` (columns constant within every profile — a
 * `Day` that varies inside a profile is a time axis and is never an occasion, architect S3).
 * Returns null when no name matches. `mapped` (the analyst's pick) wins when it is one of the
 * candidates; `''` is the analyst's explicit "(none)" (the days are pooled).
 */
export function detectOccasionColumn(
  profileConstant: ReadonlyArray<string>, mapped?: string | null,
): string | null {
  if (mapped === '') return null;
  if (mapped !== undefined && mapped !== null && profileConstant.includes(mapped))
    return mapped;
  let best: string | null = null;
  let bestRank = Infinity;
  for (const name of profileConstant) {
    const rank = NAME_PRIORITY.indexOf(normOccasionName(name));
    if (rank >= 0 && rank < bestRank) {best = name; bestRank = rank;}
  }
  return best;
}

/** Prefixes dropped before an occasion token is read as a number (`D28`, `Day 28`, `W4`). */
const TOKEN_PREFIX = /^(day|d|wk|w|visit|v)[\s_-]*/i;

/** The occasion value as a number after dropping a `D` / `Day` / `W` / `Wk` / `Visit`
 *  prefix, or NaN when it does not parse. */
export function occasionNumber(value: string): number {
  const s = value.trim().replace(TOKEN_PREFIX, '');
  if (s === '' || !/^[-+]?\d+(\.\d+)?$/.test(s)) return NaN;
  return Number(s);
}

/** Occasion levels in display order: numeric when every value parses, else first appearance.
 *  The shown label is always the raw value. */
export function orderOccasionLevels(values: ReadonlyArray<string>): {levels: string[]; numeric: boolean} {
  const seen: string[] = [];
  for (const v of values) if (v !== '' && !seen.includes(v)) seen.push(v);
  const numeric = seen.length > 0 && seen.every((v) => Number.isFinite(occasionNumber(v)));
  if (numeric) seen.sort((a, b) => occasionNumber(a) - occasionNumber(b));
  return {levels: seen, numeric};
}

// ── F1b — animal ─────────────────────────────────────────────────────────────────

/** One row per profile, as the resolution reads it. `occasion` is '' when absent. */
export interface ProfileRow {
  readonly subject: string;
  readonly occasion: string;
  readonly sex: string;
  readonly dose: number;
  readonly route: string;
  /** The mapped animal column's value, when an animal column is mapped. */
  readonly animal?: string;
}

/** Column names of a one-row-per-profile frame (the results frame, a Pipeline result). */
export interface ProfileFrameCols {
  readonly subject: string;
  readonly occasion: string | null;
  readonly sex: string | null;
  readonly dose: string | null;
  readonly route: string | null;
  readonly animal?: string | null;
}

/** Read the {@link ProfileRow}s of a one-row-per-profile frame. A missing column reads as ''
 *  (strings) or NaN (dose). */
export function profileRowsOf(frame: DG.DataFrame, c: ProfileFrameCols): ProfileRow[] {
  const col = (n: string | null | undefined): DG.Column | null => n == null ? null : frame.col(n);
  const subj = col(c.subject);
  const occ = col(c.occasion);
  const sex = col(c.sex);
  const dose = col(c.dose);
  const route = col(c.route);
  const animal = col(c.animal);
  const str = (k: DG.Column | null, i: number): string => k === null || k.isNone(i) ? '' : k.getString(i);
  const out: ProfileRow[] = [];
  for (let i = 0; i < frame.rowCount; i++) {
    out.push({
      subject: str(subj, i), occasion: str(occ, i), sex: str(sex, i),
      dose: dose === null ? NaN : numAt(dose, i), route: str(route, i),
      ...(animal !== null ? {animal: str(animal, i)} : {}),
    });
  }
  return out;
}

export type AnimalRule = 'mapped' | 'repeated-subject' | 'day-coded-id';

export interface AnimalResolution {
  /** Animal per profile row, or null when no rule applies (N then counts profiles). */
  readonly animals: ReadonlyArray<string> | null;
  readonly rule: AnimalRule | null;
  /** Animals failing verification (two profiles on one occasion, sex or dose changing). */
  readonly failing: ReadonlyArray<{animal: string; reason: string}>;
}

/**
 * Rule (iii): strip a trailing separator + token from a day-coded subject id when the token,
 * after dropping a `D` / `Day` / `W` / `Wk` prefix, equals the profile's occasion as a number.
 * `T001_D28` on occasion 28 → `T001`; `T001_D1` on occasion 2 → null (AC-1.4).
 */
export function stripDayCode(subject: string, occasion: string): string | null {
  const m = /^(.+?)[\s_.-]+([A-Za-z]*\s*\d+(?:\.\d+)?)$/.exec(subject.trim());
  if (m === null) return null;
  const tok = occasionNumber(m[2]);
  const occ = occasionNumber(occasion);
  if (!Number.isFinite(tok) || !Number.isFinite(occ) || tok !== occ) return null;
  return m[1];
}

/** Resolve the animal per profile: (i) mapped column, (ii) the subject id repeats across
 *  occasions, (iii) a day-coded subject id. */
export function resolveAnimals(rows: ReadonlyArray<ProfileRow>, hasOccasion: boolean): AnimalResolution {
  const none: AnimalResolution = {animals: null, rule: null, failing: []};
  if (rows.length === 0) return none;
  if (rows.every((r) => r.animal !== undefined && r.animal !== ''))
    return verifyAnimals(rows, rows.map((r) => r.animal!), 'mapped');
  if (!hasOccasion) return none;
  // (ii) the same subject id on more than one occasion.
  const occBySubject = new Map<string, Set<string>>();
  for (const r of rows) {
    const s = occBySubject.get(r.subject) ?? new Set<string>();
    s.add(r.occasion);
    occBySubject.set(r.subject, s);
  }
  if ([...occBySubject.values()].some((s) => s.size > 1))
    return verifyAnimals(rows, rows.map((r) => r.subject), 'repeated-subject');
  // (iii) day-coded ids: every profile must strip, and at least one animal must span 2 occasions.
  const stripped = rows.map((r) => stripDayCode(r.subject, r.occasion));
  if (stripped.every((a): a is string => a !== null)) {
    const occByAnimal = new Map<string, Set<string>>();
    stripped.forEach((a, i) => occByAnimal.set(a, (occByAnimal.get(a) ?? new Set()).add(rows[i].occasion)));
    if ([...occByAnimal.values()].some((s) => s.size > 1))
      return verifyAnimals(rows, stripped, 'day-coded-id');
  }
  return none;
}

/** Verification: each occasion at most once per animal; sex constant; dose constant unless
 *  the dose changes on every occasion of the animal (an escalation is legal, judged by the
 *  design). Failing animals are listed; the key stays `animal` for the rest. */
function verifyAnimals(rows: ReadonlyArray<ProfileRow>, animals: string[], rule: AnimalRule): AnimalResolution {
  const byAnimal = new Map<string, number[]>();
  animals.forEach((a, i) => byAnimal.set(a, [...(byAnimal.get(a) ?? []), i]));
  const failing: {animal: string; reason: string}[] = [];
  for (const [a, idx] of byAnimal) {
    const occ = idx.map((i) => rows[i].occasion);
    if (new Set(occ).size !== occ.length) {failing.push({animal: a, reason: 'two profiles on one occasion'}); continue;}
    if (new Set(idx.map((i) => rows[i].sex)).size > 1) failing.push({animal: a, reason: 'sex changes'});
  }
  return {animals, rule, failing};
}

// ── F1a — design ─────────────────────────────────────────────────────────────────

export type OccasionDesign = 'repeat-dose' | 'crossover' | 'escalation' | 'unknown';
export const OCCASION_DESIGNS: ReadonlyArray<OccasionDesign> = ['repeat-dose', 'crossover', 'escalation', 'unknown'];

export interface DesignProposal {
  readonly design: OccasionDesign;
  /** The evidence shown beside the proposal, e.g. "same animal, route and dose constant on 1 and 28". */
  readonly evidence: string;
}

/**
 * Propose the occasion design from within-animal variation (never a silent decision — the
 * pane shows the proposal with its evidence and the analyst confirms it). With a verified
 * animal: route varies across an animal's occasions → `crossover`; dose varies, route
 * constant → `escalation`; both constant → `repeat-dose` for a day-like name, `unknown` for a
 * period-like one. No animal → `unknown`.
 */
export function classifyOccasionDesign(
  rows: ReadonlyArray<ProfileRow>, animal: AnimalResolution, occasionCol: string, levels: ReadonlyArray<string>,
): DesignProposal {
  const days = levels.join(', ');
  if (animal.animals === null)
    return {design: 'unknown', evidence: `no animal key: the profiles on ${days} cannot be linked to one animal`};
  const bad = new Set(animal.failing.map((f) => f.animal));
  const byAnimal = new Map<string, number[]>();
  animal.animals.forEach((a, i) => {if (!bad.has(a)) byAnimal.set(a, [...(byAnimal.get(a) ?? []), i]);});
  let routeVaries = false;
  let doseVaries = false;
  let linked = 0;
  for (const idx of byAnimal.values()) {
    if (idx.length < 2) continue;
    linked++;
    if (new Set(idx.map((i) => rows[i].route)).size > 1) routeVaries = true;
    if (new Set(idx.map((i) => rows[i].dose)).size > 1) doseVaries = true;
  }
  if (linked === 0)
    return {design: 'unknown', evidence: `no animal is observed on more than one of ${days}`};
  if (routeVaries) return {design: 'crossover', evidence: `same animal, route changes between ${days}`};
  if (doseVaries) return {design: 'escalation', evidence: `same animal, dose changes between ${days}, route constant`};
  const dayLike = isDayLikeOccasionName(occasionCol);
  return {
    design: dayLike ? 'repeat-dose' : 'unknown',
    evidence: `same animal, route and dose constant on ${days}` +
      (dayLike ? '' : ` ('${occasionCol}' may be a crossover period — confirm the design)`),
  };
}

// ── The resolution every consumer reads ──────────────────────────────────────────

/** The analyst's choices (recipe `derived`), all optional: absent = automatic. */
export interface DerivedChoices {
  readonly occasionCol?: string;
  readonly animalCol?: string;
  readonly occasionDesign?: OccasionDesign;
  readonly tauH?: number;
}

export interface OccasionResolution {
  /** The occasion column, or null. */
  readonly col: string | null;
  readonly levels: ReadonlyArray<string>;
  readonly numericOrder: boolean;
  readonly animal: AnimalResolution;
  readonly proposal: DesignProposal;
  /** The design the analyst confirmed, or null while it is only proposed. */
  readonly confirmed: OccasionDesign | null;
  /** True when DP, the sex comparison and accumulation report per occasion. */
  readonly stratify: boolean;
  /** F% pairing: `same-occasion` only for a CONFIRMED repeat-dose design (or a Pipeline run,
   *  which has nobody to confirm — see `resolveOccasion`'s `autoConfirm`). */
  readonly fPctPairing: 'same-occasion' | 'across-occasions';
}

/**
 * Resolve occasion, animal and design for one frame of profile rows.
 * `autoConfirm` (Pipeline only): adopt the proposal as the confirmed design — the Pipeline has
 * no analyst; a period-like name still proposes `unknown`, so crossover pairing is kept there.
 */
export function resolveOccasion(
  occasionCol: string | null, rows: ReadonlyArray<ProfileRow>, choices: DerivedChoices = {},
  autoConfirm = false,
): OccasionResolution {
  const {levels, numeric} = orderOccasionLevels(rows.map((r) => r.occasion));
  const hasOcc = occasionCol !== null && levels.length >= 2;
  const animal = resolveAnimals(rows, hasOcc);
  const proposal: DesignProposal = hasOcc ?
    classifyOccasionDesign(rows, animal, occasionCol!, levels) :
    {design: 'unknown', evidence: ''};
  const confirmed = choices.occasionDesign ?? (autoConfirm && hasOcc ? proposal.design : null);
  const effective = confirmed ?? proposal.design;
  // Dose constant within animal unless the design is an escalation (F1b step 2).
  let verified = animal;
  if (animal.animals !== null && effective !== 'escalation') {
    const doses = new Map<string, Set<number>>();
    animal.animals.forEach((a, i) => doses.set(a, (doses.get(a) ?? new Set()).add(rows[i].dose)));
    const already = new Set(animal.failing.map((f) => f.animal));
    const extra = [...doses].filter(([a, d]) => d.size > 1 && !already.has(a))
      .map(([a]) => ({animal: a, reason: 'dose changes'}));
    if (extra.length > 0) verified = {...animal, failing: [...animal.failing, ...extra]};
  }
  return {
    col: occasionCol, levels, numericOrder: numeric, animal: verified, proposal, confirmed,
    stratify: hasOcc && effective !== 'escalation',
    fPctPairing: hasOcc && confirmed === 'repeat-dose' ? 'same-occasion' : 'across-occasions',
  };
}

/** Notices the resolution raises (shown over the results / F% panels). Text lives here (S8). */
export function occasionNotices(r: OccasionResolution): {key: string; severity: 'info' | 'warning'; text: string}[] {
  const out: {key: string; severity: 'info' | 'warning'; text: string}[] = [];
  if (r.col === null || r.levels.length < 2) return out;
  if (r.animal.animals === null) {
    out.push({key: 'OCCASION_NO_ANIMAL', severity: 'info',
      text: `Per-${r.col} results count profiles, not animals: no animal key links the profiles on ` +
        `${r.levels.join(', ')}. Pick the Animal column to count animals.`});
  }
  if (r.animal.failing.length > 0) {
    out.push({key: 'OCCASION_FAILING_ANIMALS', severity: 'warning',
      text: `${r.animal.failing.length} animal(s) excluded from per-animal outputs: ` +
        r.animal.failing.map((f) => `${f.animal} (${f.reason})`).join(', ') + '.'});
  }
  if ((r.confirmed ?? r.proposal.design) === 'unknown') {
    out.push({key: 'OCCASION_DESIGN_UNKNOWN', severity: 'warning',
      text: `The design of '${r.col}' (${r.levels.join(', ')}) could not be classified ` +
        `(${r.proposal.evidence}). F% pairs IV and extravascular profiles across occasions, as for a ` +
        'crossover. Confirm the design to pair within each occasion.'});
  }
  if (r.confirmed === null && r.proposal.design !== 'unknown') {
    out.push({key: 'OCCASION_DESIGN_PROPOSED', severity: 'info',
      text: `'${r.col}' looks like ${r.proposal.design === 'escalation' ? 'an' : 'a'} ${r.proposal.design} ` +
        `design (${r.proposal.evidence}). ` +
        'F% keeps cross-occasion pairing until you confirm the design.'});
  }
  if ((r.confirmed ?? r.proposal.design) === 'escalation') {
    out.push({key: 'OCCASION_ESCALATION', severity: 'warning',
      text: 'Dose changes within animal across occasions: dose proportionality pools the occasions; ' +
        'within-animal dependence not modelled; CI anti-conservative.'});
  }
  return out;
}

// ── AC-1.9 — may an occasion-named column become a Pipeline profile key? ─────────

/**
 * True when an occasion-named column splits each subject into genuine profiles: it is not the
 * time column, each (subject, value) block holds ≥ 3 distinct times (the λz minimum), and time
 * restarts per block (every block's minimum time is within 10 % of τ — of the median block
 * span when τ is unknown — of the smallest block minimum). A per-row nominal `Day` on one
 * daily-sampled profile fails the restart test, so it never splits a profile.
 */
export function splitsProfiles(
  subjects: ReadonlyArray<string>, values: ReadonlyArray<string>, times: ArrayLike<number>,
  tauH?: number,
): boolean {
  const blocks = new Map<string, number[]>();
  for (let i = 0; i < subjects.length; i++) {
    if (!Number.isFinite(times[i])) continue;
    const k = `${subjects[i]}\u0001${values[i]}`;
    const b = blocks.get(k);
    if (b === undefined) blocks.set(k, [times[i]]); else b.push(times[i]);
  }
  if (blocks.size === 0) return false;
  const mins: number[] = [];
  const spans: number[] = [];
  for (const t of blocks.values()) {
    if (new Set(t).size < 3) return false;
    const mn = Math.min(...t);
    mins.push(mn);
    spans.push(Math.max(...t) - mn);
  }
  spans.sort((a, b) => a - b);
  const scale = tauH !== undefined && tauH > 0 ? tauH : spans[Math.floor(spans.length / 2)];
  const lo = Math.min(...mins);
  return mins.every((m) => m - lo <= 0.1 * scale);
}
