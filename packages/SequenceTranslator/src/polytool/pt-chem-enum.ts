/* eslint-disable max-len */
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/**
 * PolyTool Chemical Enumeration — pure logic module.
 *
 * Responsibilities:
 *   • Recognize and normalize the many spellings of R-group labels used in SMILES
 *     (`[1*]`, `[*:1]`, `[R1]`, `[R:1]`, `[*1]`, including multi-digit R numbers).
 *   • Validate cores (≥1 R-label) and R-groups (exactly one R-label).
 *   • Join a core SMILES with one R-group per R-number via SMILES-concatenation
 *     using shared ring-closure digits, then canonicalize through RDKit.
 *   • Enumerate across multiple cores and R-group lists in Zip or Cartesian mode.
 *   • Enforce a hard result cap.
 *
 * No Datagrok or UI dependencies — safe to unit-test with just an RDKit module.
 */

// ─── Constants ──────────────────────────────────────────────────────────────

export const CHEM_ENUM_MAX_RESULTS = 1_000_000;

export const ChemEnumModes = {
  Zip: 'Zip',
  Cartesian: 'Cartesian',
} as const;
export type ChemEnumMode = typeof ChemEnumModes[keyof typeof ChemEnumModes];

// ─── Types ──────────────────────────────────────────────────────────────────

export interface ChemEnumCore {
  /** Normalized SMILES — all R-labels rewritten to `[*:N]` form. */
  smiles: string;
  /** SMILES as supplied by the user (pre-normalization). */
  originalSmiles: string;
  /** Unique, sorted ascending, R numbers present in the core. */
  rNumbers: number[];
  /** Human-readable source, e.g. `"Drawn #1"` or `"Cores[3]"`. */
  id: string;
  /** If set, the core is invalid and should be excluded from enumeration. */
  error?: string;
}

export interface ChemEnumRGroup {
  /** Normalized SMILES with its single R-label remapped to the target R number. */
  smiles: string;
  /** SMILES as supplied (pre-normalization and pre-remap). */
  originalSmiles: string;
  /** Target R number this group fills in. */
  rNumber: number;
  /** R number as originally written in `originalSmiles` (pre-remap). */
  sourceRNumber?: number;
  id: string;
  error?: string;
}

export interface ChemEnumParams {
  cores: ChemEnumCore[];
  /** Key = R number (1-based). Value = list of R-groups for that slot. */
  rGroups: Map<number, ChemEnumRGroup[]>;
  mode: ChemEnumMode;
  /** Overrides `CHEM_ENUM_MAX_RESULTS` for tests. */
  maxResults?: number;
}

export interface ChemEnumResult {
  /** Canonical SMILES of the assembled molecule. */
  smiles: string;
  /** `originalSmiles` of the core used. */
  coreSmiles: string;
  /** R-number → `originalSmiles` of the R-group used at that position. */
  rGroupSmilesByNum: Map<number, string>;
}

export interface ChemEnumValidation {
  /** Overall validity flag. */
  ok: boolean;
  /** Free-form top-level messages (e.g. zip-length mismatch, cap exceeded). */
  errors: string[];
  /** Predicted total result count (capped at MAX+1 to signal "too many"). */
  predictedCount: number;
  /** True when `predictedCount > maxResults`. */
  overCap: boolean;
}

// ─── R-label recognition ────────────────────────────────────────────────────

/**
 * Matches all supported R-label spellings as a single bracketed atom:
 *   [N*]   [*:N]   [*N]   [RN]   [R:N]
 * Capture groups 1..5 hold the numeric portion (exactly one is non-empty).
 */
const R_LABEL_SOURCE = String.raw`\[(?:(\d+)\*|\*:(\d+)|\*(\d+)|R(\d+)|R:(\d+))\]`;

export function rLabelRegex(): RegExp {
  return new RegExp(R_LABEL_SOURCE, 'g');
}

function pickNum(groups: string[]): number | null {
  for (const g of groups) if (g !== undefined) return parseInt(g, 10);
  return null;
}

/** Replaces every supported R-label spelling with the canonical `[*:N]` form. */
export function normalizeRLabels(smi: string): string {
  return smi.replace(rLabelRegex(), (_m, g1, g2, g3, g4, g5) => {
    const n = pickNum([g1, g2, g3, g4, g5]);
    return n === null ? _m : `[*:${n}]`;
  });
}

/** Returns R numbers found in the SMILES, sorted ascending, deduplicated. */
export function extractRNumbers(smi: string): number[] {
  const seen = new Set<number>();
  for (const m of smi.matchAll(rLabelRegex())) {
    const n = pickNum([m[1], m[2], m[3], m[4], m[5]]);
    if (n !== null) seen.add(n);
  }
  return [...seen].sort((a, b) => a - b);
}

/** Returns every R-label occurrence with its source number (order-preserving). */
export function findRLabels(smi: string): { source: number, match: string, index: number }[] {
  const out: { source: number, match: string, index: number }[] = [];
  for (const m of smi.matchAll(rLabelRegex())) {
    const n = pickNum([m[1], m[2], m[3], m[4], m[5]]);
    if (n !== null) out.push({source: n, match: m[0], index: m.index!});
  }
  return out;
}

/**
 * Rewrites every R-label in `smi` to `[*:newN]`, regardless of source number
 * or spelling. Intended for single-R groups being remapped into their assigned slot.
 */
export function remapSingleRLabel(smi: string, newN: number): string {
  return smi.replace(rLabelRegex(), () => `[*:${newN}]`);
}

// ─── Core / R-group construction ────────────────────────────────────────────

export function makeCore(originalSmiles: string, id: string, rdkit?: RDModule): ChemEnumCore {
  const trimmed = (originalSmiles ?? '').trim();
  if (trimmed === '')
    return {smiles: '', originalSmiles, rNumbers: [], id, error: 'Empty SMILES'};

  const normalized = normalizeRLabels(trimmed);
  const rNumbers = extractRNumbers(normalized);

  if (rNumbers.length === 0)
    return {smiles: normalized, originalSmiles, rNumbers, id, error: 'Core must contain at least one R group'};

  if (rdkit) {
    const err = tryParse(normalized, rdkit);
    if (err) return {smiles: normalized, originalSmiles, rNumbers, id, error: `Invalid SMILES: ${err}`};
  }

  return {smiles: normalized, originalSmiles, rNumbers, id};
}

export function makeRGroup(
  originalSmiles: string, targetRNumber: number, id: string, rdkit?: RDModule,
): ChemEnumRGroup {
  const trimmed = (originalSmiles ?? '').trim();
  if (trimmed === '')
    return {smiles: '', originalSmiles, rNumber: targetRNumber, id, error: 'Empty SMILES'};

  const normalized = normalizeRLabels(trimmed);
  const rNumbers = extractRNumbers(normalized);

  if (rNumbers.length === 0) {
    return {
      smiles: normalized, originalSmiles, rNumber: targetRNumber, id,
      error: 'R-group must contain exactly one R label (found none)'};
  }
  if (rNumbers.length > 1) {
    return {
      smiles: normalized, originalSmiles, rNumber: targetRNumber, id,
      sourceRNumber: rNumbers[0],
      error: `R-group must contain exactly one R label (found ${rNumbers.length}: ${rNumbers.map((n) => 'R' + n).join(', ')})`};
  }

  const sourceRNumber = rNumbers[0];
  const remapped = remapSingleRLabel(normalized, targetRNumber);

  if (rdkit) {
    const err = tryParse(remapped, rdkit);
    if (err) {
      return {
        smiles: remapped, originalSmiles, rNumber: targetRNumber, id, sourceRNumber,
        error: `Invalid SMILES: ${err}`};
    }
  }

  return {smiles: remapped, originalSmiles, rNumber: targetRNumber, id, sourceRNumber};
}

function tryParse(smi: string, rdkit: RDModule): string | null {
  let mol: RDMol | null = null;
  try {
    mol = rdkit.get_mol(smi);
    if (!mol || !mol.is_valid()) return 'failed to parse';
    return null;
  } catch (err: any) {
    return (err?.message ?? String(err)).toString().slice(0, 120);
  } finally {
    mol?.delete();
  }
}

// ─── Count + validation ─────────────────────────────────────────────────────

/** Results per core: depends on mode and the R-numbers the core uses. */
export function countForCore(
  core: ChemEnumCore, rGroups: Map<number, ChemEnumRGroup[]>, mode: ChemEnumMode,
): { count: number, uncovered: number[] } {
  const uncovered: number[] = [];
  const counts: number[] = [];
  for (const n of core.rNumbers) {
    const list = rGroups.get(n);
    if (!list || list.length === 0) { uncovered.push(n); continue; }
    counts.push(list.filter((g) => !g.error).length);
  }
  if (uncovered.length > 0) return {count: 0, uncovered};
  if (counts.some((c) => c === 0)) return {count: 0, uncovered};

  if (mode === ChemEnumModes.Zip) {
    if (counts.length === 0) return {count: 0, uncovered};
    const first = counts[0];
    return {count: counts.every((c) => c === first) ? first : -1, uncovered};
  }
  // Cartesian
  return {count: counts.reduce((p, c) => p * c, 1), uncovered};
}

/**
 * Quickly validates the overall enumeration parameters, predicting the total
 * result count and flagging structural issues (invalid inputs, uncovered R-numbers,
 * zip-length mismatches, cap exceedance).
 */
export function validateParams(params: ChemEnumParams): ChemEnumValidation {
  const max = params.maxResults ?? CHEM_ENUM_MAX_RESULTS;
  const errors: string[] = [];
  let total = 0;

  const validCores = params.cores.filter((c) => !c.error);
  if (validCores.length === 0) errors.push('No valid cores provided.');

  // Collect all R-numbers used by any valid core
  const usedRs = new Set<number>();
  for (const c of validCores) for (const n of c.rNumbers) usedRs.add(n);

  // Per-R-number R-group validity counts
  const rgCounts = new Map<number, number>();
  for (const n of usedRs) {
    const list = params.rGroups.get(n) ?? [];
    const valid = list.filter((g) => !g.error).length;
    rgCounts.set(n, valid);
    if (valid === 0) errors.push(`No valid R-group provided for R${n}.`);
  }

  // Zip mode: all used R-group lists must share a length > 0
  if (params.mode === ChemEnumModes.Zip) {
    const lens = [...rgCounts.values()].filter((v) => v > 0);
    if (lens.length > 1 && !lens.every((v) => v === lens[0]))
      errors.push(`Zip mode requires every R-group list to have the same number of entries. Got ${[...rgCounts.entries()].map(([n, v]) => `R${n}=${v}`).join(', ')}.`);
  }

  for (const c of validCores) {
    const {count, uncovered} = countForCore(c, params.rGroups, params.mode);
    if (uncovered.length > 0) {
      errors.push(`Core "${c.id}" references uncovered R-number${uncovered.length > 1 ? 's' : ''}: ${uncovered.map((n) => 'R' + n).join(', ')}.`);
      continue;
    }
    if (count < 0) continue; // already covered by the global zip-mismatch message
    total += count;
    if (total > max) break;
  }

  const overCap = total > max;
  if (overCap)
    errors.push(`Too many combinations (> ${max.toLocaleString()}). Reduce the number of cores or R-groups, or switch to Zip mode.`);

  return {ok: errors.length === 0, errors, predictedCount: total, overCap};
}

// ─── SMILES assembly ────────────────────────────────────────────────────────

/**
 * Pure-string join: core + one R-group per R-number share ring-closure digits
 * across a disconnected SMILES. Returns a SMILES that RDKit can parse, but
 * **not canonicalized** — call `Chem:convertNotation` or {@link assembleMolecule}
 * for canonical output. No RDKit calls — safe to run on millions of rows without
 * blocking the main thread.
 *
 *   core: `C[*:1]N[*:2]`, R1=`O[*:1]`, R2=`S[*:2]`
 *   → `C%50N%51.O%50.S%51` (uncanonical but valid)
 */
export function buildJoinedSmiles(
  coreSmiles: string,
  rgSmilesByNum: Map<number, string>,
): string | null {
  if (rgSmilesByNum.size === 0) return null;

  const coreFixed = moveStartRLabelToBranch(coreSmiles);
  const rgsFixed = new Map<number, string>();
  for (const [k, s] of rgSmilesByNum) rgsFixed.set(k, moveStartRLabelToBranch(s));

  const allPieces = [coreFixed, ...rgsFixed.values()];
  const digits = pickFreeRingDigits(allPieces, rgSmilesByNum.size);
  if (digits.length < rgSmilesByNum.size) return null;

  const digitByNum = new Map<number, string>();
  let i = 0;
  for (const k of rgSmilesByNum.keys()) digitByNum.set(k, formatRingDigit(digits[i++]));

  let assembledCore = coreFixed;
  for (const [k, d] of digitByNum)
    assembledCore = substituteRLabelWithRingDigit(assembledCore, k, d);

  const assembledRgs: string[] = [];
  for (const [k, s] of rgsFixed)
    assembledRgs.push(substituteRLabelWithRingDigit(s, k, digitByNum.get(k)!));

  return [assembledCore, ...assembledRgs].join('.');
}

/**
 * Joins a core with one R-group per R-number and canonicalizes via RDKit.
 * Per-molecule sync RDKit call — **do not use in bulk**; prefer {@link buildJoinedSmiles}
 * + a batched `Chem:convertNotation` over the whole column.
 */
export function assembleMolecule(
  coreSmiles: string,
  rgSmilesByNum: Map<number, string>,
  rdkit: RDModule,
): string | null {
  const joined = buildJoinedSmiles(coreSmiles, rgSmilesByNum);
  if (!joined) return null;
  let mol: RDMol | null = null;
  try {
    mol = rdkit.get_mol(joined);
    if (!mol || !mol.is_valid()) return null;
    return mol.get_smiles();
  } catch {
    return null;
  } finally {
    mol?.delete();
  }
}

/**
 * `[*:N]X…` or `[*:N]=X…` at SMILES start becomes `X([*:N])…` / `X(=[*:N])…`
 * so every R-label is preceded by an atom — required for ring-digit substitution.
 */
export function moveStartRLabelToBranch(smi: string): string {
  const m = smi.match(/^(\[\*:\d+\])([-=#:/\\])?(\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnops])(.*)$/s);
  if (!m) return smi;
  const [, rlab, bond, atom, rest] = m;
  return `${atom}(${bond ?? ''}${rlab})${rest}`;
}

/** Replaces `[*:n]` and `([*:n])` in `smi` with a ring-closure token. */
export function substituteRLabelWithRingDigit(smi: string, n: number, digitToken: string): string {
  const target = `[*:${n}]`;
  // Collapse a lone-branch form first so `(` / `)` don't linger: X([*:n]) → X<digit>
  const branchForm = new RegExp(`\\(\\s*\\[\\*:${n}\\]\\s*\\)`, 'g');
  smi = smi.replace(branchForm, digitToken);
  return smi.split(target).join(digitToken);
}

/**
 * Picks `count` ring-closure digits not already in use in any of the pieces,
 * formatted as bare digits when possible and `%NN` otherwise.
 */
export function pickFreeRingDigits(pieces: string[], count: number): number[] {
  const used = new Set<number>();
  for (const p of pieces) {
    const stripped = p.replace(/\[[^\]]*\]/g, ''); // atoms are bracketed — ignore their digits
    for (const m of stripped.matchAll(/%(\d{2})/g))
      used.add(parseInt(m[1], 10));
    for (const ch of stripped) {
      const v = ch.charCodeAt(0) - 48;
      if (v >= 0 && v <= 9) used.add(v);
    }
  }
  const free: number[] = [];
  for (let d = 1; d < 100 && free.length < count; d++)
    if (!used.has(d)) free.push(d);
  return free;
}

export function formatRingDigit(n: number): string {
  if (n < 0 || n > 99) throw new Error(`Ring digit out of range: ${n}`);
  return n <= 9 ? `${n}` : `%${n.toString().padStart(2, '0')}`;
}

// ─── Enumeration ────────────────────────────────────────────────────────────

/** Enumerates R-group assignments per core, yielding up to `params.maxResults`. */
export function* iterateAssignments(params: ChemEnumParams): Generator<{core: ChemEnumCore, assignment: Map<number, ChemEnumRGroup>}> {
  const max = params.maxResults ?? CHEM_ENUM_MAX_RESULTS;
  let produced = 0;

  for (const core of params.cores) {
    if (core.error) continue;

    const rNums = core.rNumbers;
    const lists: ChemEnumRGroup[][] = [];
    let uncovered = false;
    for (const n of rNums) {
      const list = (params.rGroups.get(n) ?? []).filter((g) => !g.error);
      if (list.length === 0) { uncovered = true; break; }
      lists.push(list);
    }
    if (uncovered) continue;

    if (params.mode === ChemEnumModes.Zip) {
      if (lists.length === 0) continue;
      const N = lists[0].length;
      if (!lists.every((l) => l.length === N)) continue;
      for (let i = 0; i < N; i++) {
        const assignment = new Map<number, ChemEnumRGroup>();
        for (let j = 0; j < rNums.length; j++) assignment.set(rNums[j], lists[j][i]);
        yield {core, assignment};
        if (++produced >= max) return;
      }
    } else {
      // Cartesian — odometer iteration
      const idx = new Array<number>(lists.length).fill(0);
      while (true) {
        const assignment = new Map<number, ChemEnumRGroup>();
        for (let j = 0; j < rNums.length; j++) assignment.set(rNums[j], lists[j][idx[j]]);
        yield {core, assignment};
        if (++produced >= max) return;

        let k = idx.length - 1;
        while (k >= 0) {
          idx[k]++;
          if (idx[k] < lists[k].length) break;
          idx[k] = 0;
          k--;
        }
        if (k < 0) break;
      }
    }
  }
}

/**
 * Runs the full enumeration. Returns `null` when validation fails (errors available
 * via {@link validateParams}). Silently skips assignments that fail to assemble.
 */
export function enumerate(params: ChemEnumParams, rdkit: RDModule): ChemEnumResult[] | null {
  const v = validateParams(params);
  if (!v.ok) return null;

  const out: ChemEnumResult[] = [];
  for (const {core, assignment} of iterateAssignments(params)) {
    const rgSmiByNum = new Map<number, string>();
    for (const [n, rg] of assignment)
      rgSmiByNum.set(n, rg.smiles);
    const smi = assembleMolecule(core.smiles, rgSmiByNum, rdkit);
    if (!smi) continue;

    const originalRgs = new Map<number, string>();
    for (const [n, rg] of assignment)
      originalRgs.set(n, rg.originalSmiles);
    out.push({smiles: smi, coreSmiles: core.originalSmiles, rGroupSmilesByNum: originalRgs});
  }
  return out;
}

/**
 * Reservoir-samples up to `sampleSize` results for a live preview.
 * Total iteration count is capped at `params.maxResults` for safety.
 */
export function enumerateSample(
  params: ChemEnumParams, rdkit: RDModule, sampleSize: number, rand: () => number = Math.random,
): ChemEnumResult[] {
  const reservoir: ChemEnumResult[] = [];
  let seen = 0;
  for (const {core, assignment} of iterateAssignments(params)) {
    const rgSmiByNum = new Map<number, string>();
    for (const [n, rg] of assignment) rgSmiByNum.set(n, rg.smiles);
    const smi = assembleMolecule(core.smiles, rgSmiByNum, rdkit);
    if (!smi) continue;

    const originalRgs = new Map<number, string>();
    for (const [n, rg] of assignment) originalRgs.set(n, rg.originalSmiles);
    const item: ChemEnumResult = {smiles: smi, coreSmiles: core.originalSmiles, rGroupSmilesByNum: originalRgs};

    if (reservoir.length < sampleSize) {
      reservoir.push(item);
    } else {
      const j = Math.floor(rand() * (seen + 1));
      if (j < sampleSize) reservoir[j] = item;
    }
    seen++;
  }
  return reservoir;
}

/**
 * No-RDKit enumeration — returns *uncanonical* joined SMILES per assignment.
 * Intended as the first stage of a bulk pipeline: collect these into a column
 * and canonicalize with one parallel `Chem:convertNotation` call instead of
 * per-row sync RDKit work.
 */
export function enumerateRaw(params: ChemEnumParams): ChemEnumResult[] | null {
  const v = validateParams(params);
  if (!v.ok) return null;

  const out: ChemEnumResult[] = [];
  for (const {core, assignment} of iterateAssignments(params)) {
    const rgSmiByNum = new Map<number, string>();
    for (const [n, rg] of assignment) rgSmiByNum.set(n, rg.smiles);
    const smi = buildJoinedSmiles(core.smiles, rgSmiByNum);
    if (!smi) continue;

    const originalRgs = new Map<number, string>();
    for (const [n, rg] of assignment) originalRgs.set(n, rg.originalSmiles);
    out.push({smiles: smi, coreSmiles: core.originalSmiles, rGroupSmilesByNum: originalRgs});
  }
  return out;
}

/** Reservoir-sample with the no-RDKit join. Output SMILES are uncanonical but parseable. */
export function enumerateSampleRaw(
  params: ChemEnumParams, sampleSize: number, rand: () => number = Math.random,
): ChemEnumResult[] {
  const reservoir: ChemEnumResult[] = [];
  let seen = 0;
  for (const {core, assignment} of iterateAssignments(params)) {
    const rgSmiByNum = new Map<number, string>();
    for (const [n, rg] of assignment) rgSmiByNum.set(n, rg.smiles);
    const smi = buildJoinedSmiles(core.smiles, rgSmiByNum);
    if (!smi) continue;

    const originalRgs = new Map<number, string>();
    for (const [n, rg] of assignment) originalRgs.set(n, rg.originalSmiles);
    const item: ChemEnumResult = {smiles: smi, coreSmiles: core.originalSmiles, rGroupSmilesByNum: originalRgs};

    if (reservoir.length < sampleSize) {
      reservoir.push(item);
    } else {
      const j = Math.floor(rand() * (seen + 1));
      if (j < sampleSize) reservoir[j] = item;
    }
    seen++;
  }
  return reservoir;
}
