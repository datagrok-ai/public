// Synced from Rusty-HMMER web/anarci/schemes.ts by web/sync-datagrok.mjs — do not edit here.
/* eslint-disable max-len */
/** Numbering schemes of ANARCI 2024.05.21 (BSD-3-Clause, Copyright (C) 2016
 * Oxford Protein Informatics Group; notice in licenses/ANARCI.txt).
 *
 * A line-by-line port of `anarci/schemes.py` plus `number_sequence_from_alignment`
 * and `validate_numbering` from `anarci/anarci.py`. Each Python function has a
 * camelCase counterpart. The 128 HMM match states are IMGT positions; each scheme
 * maps them through a state string (`X` own position, `I` insertion), a region
 * string and per-region offsets (`_number_regions`), then renumbers the regions
 * whose insertions or gaps the alignment cannot place.
 *
 * Python semantics that differ from JavaScript are made explicit: `at`/`setAt`
 * index like Python (negative from the end, `IndexError` past it), `sorted`
 * orders tuples element-wise, and missing dictionary keys raise `KeyError`.
 * Python's `assert` raises `AnarciAssertion`; other Python exceptions raise a
 * `PythonError` named after the Python exception type. Sequences are indexed
 * by code point, as Python indexes `str`. */

import {AnarciAssertion, check} from './types.ts';
import type {ChainType, NumberedDomain, Numbering, Position, Scheme, State, StateEntry, StateVector} from './types.ts';

/** A Python exception other than `AssertionError` (`IndexError`, `KeyError`,
 * `ValueError`, `TypeError`); `name: message` is Python's `type: str(e)`. */
export class PythonError extends Error {
  constructor(name: string, message: string) {
    super(message);
    this.name = name;
  }
}

/** A scheme function's `(numbering, start index, end index)`. Python gives None
 * (here null) for the indices when the state vector holds no residue. */
type SchemeDomain = [Numbering, number | null, number | null];

// ---------------------------------------------------------------------------
// Python helpers

/** Python `items[index]`: negative indices count from the end. */
function at<T>(items: ArrayLike<T>, index: number, kind = 'list'): T {
  const i = index < 0 ? index + items.length : index;
  if (!(i >= 0 && i < items.length)) throw new PythonError('IndexError', `${kind} index out of range`);
  return items[i] as T;
}

/** Python `items[index] = value`. */
function setAt<T>(items: T[], index: number, value: T): void {
  const i = index < 0 ? index + items.length : index;
  if (!(i >= 0 && i < items.length)) throw new PythonError('IndexError', 'list assignment index out of range');
  items[i] = value;
}

/** Python `range(stop)` or `range(start, stop, step)`. */
function range(start: number, stop?: number, step = 1): number[] {
  if (stop === undefined) [start, stop] = [0, start];
  const out: number[] = [];
  if (step > 0) for (let i = start; i < stop; i += step) out.push(i);
  else for (let i = start; i > stop; i += step) out.push(i);
  return out;
}

/** Python `zip(a, b)`: stops at the shorter list. */
function zip<A, B>(a: readonly A[], b: readonly B[]): [A, B][] {
  return range(Math.min(a.length, b.length)).map((i): [A, B] => [a[i] as A, b[i] as B]);
}

/** Python `repr` of a str, as it appears in a `KeyError` message. */
function repr(value: string): string {
  const quote = value.includes('\'') && !value.includes('"') ? '"' : '\'';
  return quote + value.replaceAll('\\', '\\\\') + quote;
}

/** Python `dict[key]` for a str key. */
function lookup<T>(dict: Readonly<Record<string, T>>, key: string): T {
  if (!Object.hasOwn(dict, key)) throw new PythonError('KeyError', repr(key));
  return dict[key] as T;
}

/** Python subscripting of a value that may be None. */
function notNone<T>(value: T | null | undefined): T {
  if (value === null || value === undefined)
    throw new PythonError('TypeError', '\'NoneType\' object is not subscriptable');

  return value;
}

/** Python tuple ordering of `(number, insertion code)`; str compares by code point. */
function comparePositions(a: Position, b: Position): number {
  if (a[0] !== b[0]) return a[0] < b[0] ? -1 : 1;
  return a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0;
}

/** Python `sorted` of positions (Array#sort is stable, like `sorted`). */
function sortedPositions(positions: readonly Position[]): Position[] {
  return [...positions].sort(comparePositions);
}

/** Python `sorted` of ints (numeric, not Array#sort's default string order). */
function sortedNumbers(numbers: readonly number[]): number[] {
  return [...numbers].sort((a, b) => a - b);
}

/** `(n, ' ')` for each n: positions without insertion code. */
function plain(numbers: readonly number[]): Position[] {
  return numbers.map((n): Position => [n, ' ']);
}

/** `[(number, alphabet[i]) for i in range(count)]`: `count` insertions on `number`. */
function insertionsOn(number: number, count: number): Position[] {
  return range(count).map((i): Position => [number, at(ALPHABET, i)]);
}

/** Python `list.index` of a position. */
function indexOfPosition(positions: readonly Position[], position: Position): number {
  const i = positions.findIndex((p) => p[0] === position[0] && p[1] === position[1]);
  if (i < 0) throw new PythonError('ValueError', `(${position[0]}, ${repr(position[1])}) is not in list`);
  return i;
}

/** `[(annotations[i], region[i][1]) for i in range(len(region))]`. */
function annotate(annotations: readonly Position[], region: Numbering): Numbering {
  return region.map(([, residue], i): [Position, string] => [at(annotations, i), residue]);
}

/** `[((annotations[i], " "), region[i][1]) for i in range(len(region))]` (Wolfguy). */
function annotateNumbers(annotations: readonly number[], region: Numbering): Numbering {
  return region.map(([, residue], i): [Position, string] => [[at(annotations, i), ' '], residue]);
}

/** Python `sum(lists, [])`. */
function concatAll(numbering: readonly Numbering[]): Numbering {
  return ([] as Numbering).concat(...numbering);
}

/** Python `str[index]` by code point; None raises TypeError. */
function residueAt(residues: readonly string[], index: number | null): string {
  if (index === null) throw new PythonError('TypeError', 'string indices must be integers, not \'NoneType\'');
  return at(residues, index, 'string');
}

// ---------------------------------------------------------------------------
// Tables

/** Alphabet used for insertion codes; the last (-1th) entry is a blank for no insertion. */
const ALPHABET: readonly string[] = [
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
  'X', 'Y', 'Z', 'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN', 'OO', 'PP', 'QQ',
  'RR', 'SS', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY', 'ZZ', ' ',
];

/** ANARCI's `blosum62` dict as `<first><second>:<score>` (one-letter keys, both
 * orders only where ANARCI lists them; `_get_wolfguy_L1` tries both). */
const BLOSUM62_TABLE = `
BN:3 WL:-2 GG:6 XS:0 XD:-1 KG:-2 SE:0 XM:-1 YE:-2 WR:-3 IR:-3 XZ:-1 HE:0 VM:1 NR:0 ID:-3 FD:-3 WC:-2
NA:-2 WQ:-2 LQ:-2 SN:1 ZK:1 VN:-3 QN:0 MK:-1 VH:-3 GE:-2 SL:-2 PR:-2 DA:-2 SC:-1 ED:2 YG:-3 WP:-4
XX:-1 ZL:-3 QA:-1 VY:-1 WA:-3 GD:-1 XP:-2 KD:-1 TN:0 YF:3 WW:11 ZM:-1 LD:-4 MR:-1 YK:-2 FE:-3 ME:-2
SS:4 XC:-2 YL:-1 HR:0 PP:7 KC:-3 SA:1 PI:-3 QQ:5 LI:2 PF:-4 BA:-2 ZN:0 MQ:0 VI:3 QC:-3 IH:-3 ZD:1
ZP:-1 YW:2 TG:-2 BP:-2 PA:-1 CD:-3 YH:2 XV:-1 BB:4 ZF:-3 ML:2 FG:-3 SM:-1 MG:-3 ZQ:3 SQ:0 XA:0 VT:0
WF:1 SH:-1 XN:-1 BQ:0 KA:-1 IQ:-3 XW:-2 NN:6 WT:-2 PD:-1 BC:-3 IC:-1 VK:-2 XY:-1 KR:2 ZR:0 WE:-3
TE:-1 BR:-1 LR:-2 QR:1 XF:-1 TS:1 BD:4 ZA:-1 MN:-2 VD:-3 FA:-2 XE:-1 FH:-1 MA:-1 KQ:1 ZS:0 XG:-1
VV:4 WD:-4 XH:-1 SF:-2 XL:-1 BS:0 SG:0 PM:-2 YM:-1 HD:-1 BE:1 ZB:1 IE:-3 VE:-2 XT:0 XR:-1 RR:5 ZT:-1
YD:-3 VW:-3 FL:0 TC:-1 XQ:-1 BT:-1 KN:0 TH:-2 YI:-1 FQ:-3 TI:-1 TQ:-1 PL:-3 RA:-1 BF:-3 ZC:-3 MH:-2
VF:-1 FC:-2 LL:4 MC:-1 CR:-3 DD:6 ER:0 VP:-2 SD:0 EE:5 WG:-2 PC:-3 FR:-3 BG:-1 CC:9 IG:-4 VG:-3
WK:-3 GN:0 IN:-3 ZV:-2 AA:4 VQ:-2 FK:-3 TA:0 BV:-3 KL:-2 LN:-3 YN:-2 FF:6 LG:-4 BH:0 ZE:4 QD:0 XB:-1
ZW:-3 SK:0 XK:-1 VR:-3 KE:1 IA:-1 PH:-2 BW:-4 KK:5 HC:-3 EN:0 YQ:-1 HH:8 BI:-3 CA:0 II:4 VA:0 WI:-3
TF:-2 VS:-2 TT:5 FM:0 LE:-3 MM:5 ZG:-2 DR:-2 MD:-3 WH:-2 GC:-3 SR:-1 SI:-2 PQ:-1 YA:-2 XI:-1 EA:-1
BY:-3 KI:-3 HA:-2 PG:-2 FN:-3 HN:1 BK:0 VC:-1 TL:-1 PK:-1 WS:-3 TD:-1 TM:-1 PN:-2 KH:-1 TR:-1 YR:-2
LC:-1 BL:-4 ZY:-2 WN:-4 GA:0 SP:-1 EQ:2 CN:-3 HQ:0 DN:1 YC:-2 LH:-3 EC:-4 ZH:0 HG:-2 PE:-1 YS:-2
GR:-2 BM:-3 ZZ:4 WM:-1 YT:-2 YP:-3 YY:7 TK:-1 ZI:-3 TP:-1 VL:1 FI:0 GQ:-2 LA:-1 MI:1`;

/** BLOSUM62 keyed by `first + '\t' + second`. */
const BLOSUM62: ReadonlyMap<string, number> = new Map(BLOSUM62_TABLE.trim().split(/\s+/).map((entry): [string, number] => {
  const [pair = '', score = ''] = entry.split(':');
  return [`${pair[0]}\t${pair[1]}`, Number(score)];
}));

/** `blosum62[(a, b)]`, or undefined when the key is missing. */
function blosum62(a: string, b: string): number | undefined {
  return BLOSUM62.get(`${a}\t${b}`);
}

// ---------------------------------------------------------------------------
// Alignment smoothing and region numbering

/** Insertion patterns enforced at the ends of the framework regions (`smooth_insertions`). */
const ENFORCED_PATTERNS: readonly (readonly State[])[] = [
  [[25, 'm'], [26, 'm'], [27, 'm'], [28, 'i']],
  [[38, 'i'], [38, 'm'], [39, 'm'], [40, 'm']],
  [[54, 'm'], [55, 'm'], [56, 'm'], [57, 'i']],
  [[65, 'i'], [65, 'm'], [66, 'm'], [67, 'm']],
  [[103, 'm'], [104, 'm'], [105, 'm'], [106, 'i']],
  [[117, 'i'], [117, 'm'], [118, 'm'], [119, 'm']],
];

/** `smooth_insertions`: moves HMMER insertions at the ends of the framework
 * regions into the CDRs, and N-terminal insertions into missing N-terminal
 * states. A buffer still open when the vector ends is dropped, as in Python. */
export function smoothInsertions(stateVector: StateVector): StateVector {
  let stateBuffer: StateEntry[] = [];
  const sv: StateVector = [];
  // Python leaves `reg` unbound until a state is buffered; it is read only when the buffer is not empty.
  let reg = -1;
  for (const entry of stateVector) {
    const stateId = entry[0][0];
    if (stateId < 23) { // Everything before the cysteine at 23.
      stateBuffer.push(entry);
      reg = -1;
    } else if (25 <= stateId && stateId < 28) {
      stateBuffer.push(entry);
      reg = 0;
    } else if (37 < stateId && stateId <= 40) {
      stateBuffer.push(entry);
      reg = 1;
    } else if (54 <= stateId && stateId < 57) {
      stateBuffer.push(entry);
      reg = 2;
    } else if (64 < stateId && stateId <= 67) {
      stateBuffer.push(entry);
      reg = 3;
    } else if (103 <= stateId && stateId < 106) {
      stateBuffer.push(entry);
      reg = 4;
    } else if (116 < stateId && stateId <= 119) {
      stateBuffer.push(entry);
      reg = 5;
    } else if (stateBuffer.length !== 0) { // Add the buffer and reset.
      const nins = stateBuffer.filter((s) => s[0][1] === 'i').length;
      if (nins > 0) {
        if (reg === -1) { // FW1: only adjust with at least as many N-terminal deletions as insertions.
          let ntDels = at(stateBuffer, 0)[0][0] - 1; // Missing states.
          for (const [[, type], bufferSi] of stateBuffer) { // Explicit deletion states.
            if (type === 'd' || bufferSi === null) ntDels += 1;
            else break; // First residue found.
          }
          if (ntDels >= nins) { // Likely misalignment.
            let newStates: State[] = stateBuffer.filter(([s]) => s[1] === 'm').map(([s]) => s);
            const first = at(newStates, 0)[0];
            stateBuffer = stateBuffer.filter((s) => s[0][1] !== 'd');
            const add = stateBuffer.length - newStates.length;
            check(add >= 0, 'Implementation logic error');
            newStates = [...range(first - add, first).map((id): State => [id, 'm']), ...newStates];
            check(newStates.length === stateBuffer.length, 'Implementation logic error');
            for (let i = 0; i < stateBuffer.length; i++) sv.push([at(newStates, i), at(stateBuffer, i)[1]]);
          } else
            sv.push(...stateBuffer); // Let the alignment place the insertions.
        } else {
          stateBuffer = stateBuffer.filter((s) => s[0][1] !== 'd');
          const pattern = at(ENFORCED_PATTERNS, reg);
          const length = stateBuffer.length;
          let newStates: State[];
          if (reg % 2) { // N-terminal end of a framework region.
            newStates = [...Array<State>(Math.max(0, length - 3)).fill(at(pattern, 0)),
              ...pattern.slice(Math.max(4 - length, 1))];
          } else { // C-terminal end of a framework region.
            newStates = [...pattern.slice(0, 3), ...Array<State>(Math.max(0, length - 3)).fill(at(pattern, 2))];
          }
          for (let i = 0; i < length; i++) sv.push([at(newStates, i), at(stateBuffer, i)[1]]);
        }
      } else { // All match or deletion states.
        sv.push(...stateBuffer);
      }
      sv.push(entry);
      stateBuffer = [];
    } else
      sv.push(entry);
  }
  return sv;
}

/** `_number_regions`: numbers the aligned part of the sequence region by region.
 * `rels` (offset of the scheme number from the IMGT state at the start of each
 * region) is updated in place, as Python updates the caller's dict. */
function numberRegions(sequence: string, stateVector: StateVector, stateString: string, regionString: string,
  regionIndexDict: Readonly<Record<string, number>>, rels: number[], nRegions: number,
  excludeDeletions: readonly number[]): [Numbering[], number | null, number | null] {
  const residues = Array.from(sequence);
  const smoothed = smoothInsertions(stateVector);
  const regions: Numbering[] = range(nRegions).map(() => []);

  let insertion = -1; // -1 is the blank insertion code.
  let previousStateId = 1;
  let previousStateType = 'd';
  let startIndex: number | null = null;
  let endIndex: number | null = null;
  let region: number | null = null;

  for (const [[stateId, stateType], si] of smoothed) {
    // A new region never starts with an insertion.
    if (stateType !== 'i' || region === null) region = lookup(regionIndexDict, at(regionString, stateId - 1, 'string'));

    if (stateType === 'm') {
      if (at(stateString, stateId - 1, 'string') === 'I') { // Treated as an insertion in this scheme.
        if (previousStateType !== 'd') insertion += 1; // Unless a deletion precedes it.
        rels[region] = at(rels, region) - 1;
      } else
        insertion = -1;

      at(regions, region).push([[stateId + at(rels, region), at(ALPHABET, insertion)], residueAt(residues, si)]);
      previousStateId = stateId;
      if (startIndex === null) startIndex = si;
      endIndex = si;
      previousStateType = stateType;
    } else if (stateType === 'i') {
      insertion += 1;
      at(regions, region).push([[previousStateId + at(rels, region), at(ALPHABET, insertion)], residueAt(residues, si)]);
      if (startIndex === null) startIndex = si;
      endIndex = si;
      previousStateType = stateType;
    } else { // A deletion.
      previousStateType = stateType;
      if (at(stateString, stateId - 1, 'string') === 'I') { // Irrelevant to the scheme.
        rels[region] = at(rels, region) - 1;
        continue;
      }
      insertion = -1;
      previousStateId = stateId;
    }

    // Reset the insertion index where the region will be renumbered anyway.
    if (insertion >= 25 && excludeDeletions.includes(region)) insertion = 0;
    check(insertion < 25, 'Too many insertions for numbering scheme to handle');
  }
  return [regions, startIndex, endIndex];
}

// ---------------------------------------------------------------------------
// IMGT

/** `number_imgt`: IMGT for all chain types, with CDR1-3 renumbered symmetrically. */
export function numberImgt(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX';
  const regionString = '11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777';
  const regionIndexDict = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6};
  const rels = [0, 0, 0, 0, 0, 0, 0, 0];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 7, [1, 3, 5]);

  const numbering: Numbering[] = [at(regions, 0), [], at(regions, 2), [], at(regions, 4), [], at(regions, 6)];

  /** Renumbers one CDR region symmetrically between `start` and `end` (exclusive). */
  const renumber = (index: number, maxlength: number, start: number, end: number): void => {
    const cdrSeq = at(regions, index).map((x) => x[1]).filter((a) => a !== '-');
    const target = at(numbering, index);
    let si = 0;
    let previousState = start - 1;
    for (const ann of getImgtCdr(cdrSeq.length, maxlength, start, end)) {
      if (!ann) {
        target.push([[previousState + 1, ' '], '-']);
        previousState += 1;
      } else {
        target.push([ann, at(cdrSeq, si, 'string')]);
        previousState = ann[0];
        si += 1;
      }
    }
  };

  renumber(1, 12, 27, 39); // CDR1: 27 (inc.) to 39 (exc.), maximum length 12.
  renumber(3, 10, 56, 66); // CDR2: 56 (inc.) to 66 (exc.), maximum length 10.
  // FW3 insertions are placed by the alignment.
  // CDR3: 105 (inc.) to 118 (exc.), insertions on 111 and 112 symmetrically.
  const cdr3Length = at(regions, 5).filter((x) => x[1] !== '-').length;
  if (cdr3Length > 117) return [[], startIndex, endIndex]; // Too many insertions.
  renumber(5, 13, 105, 118);

  return [gapMissing(numbering), startIndex, endIndex];
}

/** `get_imgt_cdr`: symmetric annotations of a CDR of `length` residues; None
 * (null) marks a gap. Insertions go on both sides of the centre. */
export function getImgtCdr(length: number, maxlength: number, start: number, end: number): (Position | null)[] {
  const annotations: (Position | null)[] = Array<Position | null>(Math.max(length, maxlength)).fill(null);
  if (length === 0) return annotations;
  if (length === 1) {
    setAt(annotations, 0, [start, ' ']);
    return annotations;
  }

  let front = 0;
  let back = -1; // Python negative index.
  const az = ALPHABET.slice(0, -1);
  const za = [...az].reverse();

  for (let i = 0; i < Math.min(length, maxlength); i++) {
    if (i % 2) {
      setAt(annotations, back, [end + back, ' ']);
      back -= 1;
    } else {
      setAt(annotations, front, [start + front, ' ']);
      front += 1;
    }
  }

  // Add insertions around the centre point.
  const centrepoint = range(annotations.length).filter((i) => annotations[i] === null);
  if (centrepoint.length === 0) return annotations;

  const centreLeft = notNone(at(annotations, Math.min(...centrepoint) - 1))[0];
  const centreRight = notNone(at(annotations, Math.max(...centrepoint) + 1))[0];

  const half = Math.floor(maxlength / 2);
  const [frontfactor, backfactor] = maxlength % 2 ? [half + 1, half] : [half, half];

  for (let i = 0; i < Math.max(0, length - maxlength); i++) {
    if (!(i % 2)) {
      setAt(annotations, back, [centreRight, at(za, back + backfactor)]); // Negative index into `za`.
      back -= 1;
    } else {
      setAt(annotations, front, [centreLeft, at(az, front - frontfactor)]);
      front += 1;
    }
  }
  return annotations;
}

// ---------------------------------------------------------------------------
// AHo

/** AHo CDR1 (25-42) gap order per chain type. */
const AHO_CDR1_DELETIONS: Readonly<Record<string, readonly number[]>> = {
  L: [28, 36, 35, 37, 34, 38, 27, 29, 33, 39, 32, 40, 26, 30, 25, 31, 41, 42],
  K: [28, 27, 36, 35, 37, 34, 38, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  H: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  A: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  B: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  D: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  G: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
};

/** `number_aho`: AHo numbering; gap order in CDR1/CDR2 depends on the chain type. */
export function numberAho(stateVector: StateVector, sequence: string, chainType: ChainType | string): SchemeDomain {
  const stateString = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX';
  const regionString = 'BBBBBBBBBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFHHHHHHHHHHHHHHHHIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK';
  const regionIndexDict = Object.fromEntries(zip([...'ABCDEFGHIJK'], range(11)));
  const rels = [0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 21];
  const [regions, startIndex, regionsEnd] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 11, [1, 3, 4, 5, 7, 9]);
  let endIndex = regionsEnd;

  const numbering: Numbering[] = [at(regions, 0), at(regions, 1), at(regions, 2), [], at(regions, 4), [],
    at(regions, 6), [], at(regions, 8), at(regions, 9), at(regions, 10)];

  // FW1: move the indel onto 8. The first recognised residue sets the expected
  // stretch length, so N-terminal deletions are not placed at 8.
  const fw1 = at(regions, 1);
  if (fw1.length > 0) {
    const start = at(fw1, 0)[0][0];
    const stretchLen = 10 - (start - 1);
    let annotations: Position[];
    if (fw1.length > stretchLen) { // Insertions: place on 8.
      annotations = [...plain(range(start, 9)), ...insertionsOn(8, fw1.length - stretchLen), [9, ' '], [10, ' ']];
    } else {
      const orderedDeletions = plain([8, ...range(start, 11).filter((p) => p !== 8)]);
      annotations = sortedPositions(orderedDeletions.slice(Math.max(stretchLen - fw1.length, 0)));
    }
    numbering[1] = annotate(annotations, fw1);
  }

  /** Gaps region `index` in `orderedDeletions` order; insertions go linearly on
   * `insertAt`, whose index must be `expected`. False when there are too many insertions. */
  const regap = (index: number, orderedDeletions: readonly number[], insertAt: number, expected: number): boolean => {
    const region = at(regions, index);
    const size = orderedDeletions.length;
    let annotations = plain(sortedNumbers(orderedDeletions.slice(Math.max(size - region.length, 0))));
    const insertions = Math.max(region.length - size, 0);
    if (insertions > 26) return false;
    if (insertions > 0) {
      const insertat = indexOfPosition(annotations, [insertAt, ' ']) + 1;
      check(insertat === expected, 'AHo numbering failed');
      annotations = [...annotations.slice(0, insertat), ...insertionsOn(insertAt, insertions),
        ...annotations.slice(insertat)];
    }
    numbering[index] = annotate(annotations, region);
    return true;
  };

  // CDR1 (25-42): gap order by chain type; insertions on 36.
  if (!regap(3, lookup(AHO_CDR1_DELETIONS, chainType), 36, 12)) return [[], startIndex, endIndex];

  // CDR2 (58-77): gaps symmetric about 63 (VA also at 74 and 73); insertions on 63.
  const cdr2Deletions = chainType === 'A' ?
    [74, 73, 63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 75, 76, 77] :
    [63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77];
  if (!regap(5, cdr2Deletions, 63, 6)) return [[], startIndex, endIndex];

  // FW3 (78-93): deletions on 86 then 85; insertions on 85.
  if (!regap(7, [86, 85, 87, 84, 88, 83, 89, 82, 90, 81, 91, 80, 92, 79, 93, 78], 85, 8))
    return [[], startIndex, endIndex];


  // CDR3 (107-138): deletions symmetric about 123; insertions on 123.
  const cdr3Deletions = [123, 124, 122, 125, 121, 126, 120, 127, 119, 128, 118, 129, 117, 130, 116, 131, 115, 132,
    114, 133, 113, 134, 112, 135, 111, 136, 110, 137, 109, 138, 108, 107];
  if (!regap(9, cdr3Deletions, 123, 17)) return [[], startIndex, endIndex];

  // AHo has one more position than IMGT for light chains: when the last state is
  // 148 and a residue follows, number it 149.
  const result = gapMissing(numbering);
  if (result.length > 0) {
    const [lastPosition, lastResidue] = at(result, -1);
    if (lastPosition[0] === 148 && lastPosition[1] === ' ' && lastResidue !== '-') {
      const residues = Array.from(sequence);
      const next = notNone(endIndex) + 1;
      if (next < residues.length) {
        result.push([[149, ' '], at(residues, next, 'string')]);
        endIndex = next;
      }
    }
  }
  return [result, startIndex, endIndex];
}

// ---------------------------------------------------------------------------
// Chothia

/** Heavy-chain FW1 renumbering shared by Chothia, Kabat and Martin: all
 * insertions found by the HMM go on `insertAt`, followed by `after`. */
function renumberHeavyFw1(region: Numbering, insertAt: number, after: readonly number[]): Numbering {
  const insertions = region.filter((x) => x[0][1] !== ' ').length;
  if (!insertions) return region;
  const start = at(region, 0)[0][0]; // The starting number found by the HMM.
  return annotate([...plain(range(start, insertAt + 1)), ...insertionsOn(insertAt, insertions), ...plain(after)],
    region);
}

/** Heavy-chain CDR2 (50-57) shared by Chothia, Kabat and Martin: deletions in the
 * order 52, 51, 50, 53, ...; insertions on 52. */
function renumberHeavyCdr2(region: Numbering): Numbering {
  const length = region.length;
  const insertions = Math.max(length - 8, 0);
  const annotations = [
    ...plain([50, 51, 52]).slice(0, Math.max(0, length - 5)),
    ...insertionsOn(52, insertions),
    ...plain([53, 54, 55, 56, 57]).slice(Math.abs(Math.min(0, length - 5))),
  ];
  return annotate(annotations, region);
}

/** Chothia/Martin heavy CDR1 (23-33): insertions on 31. */
function renumberChothiaHeavyCdr1(region: Numbering): Numbering {
  const length = region.length;
  const insertions = Math.max(length - 11, 0); // Pulled back to the cysteine.
  const annotations = insertions ?
    [...plain(range(23, 32)), ...insertionsOn(31, insertions), ...plain([32, 33])] :
    // Python slices: `[:length - 2]` counts from the end when negative (Array#slice does the same).
    [...plain(range(23, 32)).slice(0, length - 2), ...plain([32, 33]).slice(0, length)];
  return annotate(annotations, region);
}

/** Light-chain CDR2 (51-54) shared by Chothia, Kabat and Martin: insertions on
 * 52; the alignment places deletions. */
function renumberLightCdr2(region: Numbering): Numbering {
  const insertions = Math.max(region.length - 4, 0);
  if (insertions > 0) return annotate([...plain([51, 52]), ...insertionsOn(52, insertions), ...plain([53, 54])], region);
  return region;
}

/** `number_chothia_heavy`. */
export function numberChothiaHeavy(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX';
  const regionString = '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888';
  const regionIndexDict = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7};
  const rels = [0, -1, -1, -5, -5, -8, -12, -15];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 8, [0, 2, 4, 6]);

  const numbering: Numbering[] = [[], at(regions, 1), [], at(regions, 3), [], at(regions, 5), [], at(regions, 7)];
  numbering[0] = renumberHeavyFw1(at(regions, 0), 6, [7, 8, 9]); // Insertions on 6.
  numbering[2] = renumberChothiaHeavyCdr1(at(regions, 2)); // Insertions on 31.
  numbering[4] = renumberHeavyCdr2(at(regions, 4)); // Insertions on 52.
  // FW3: the alignment places insertions (82A-C are ordinary positions).
  // CDR3 (93-102): insertions on 100.
  const cdr3 = at(regions, 6);
  if (cdr3.length > 36) return [[], startIndex, endIndex]; // Too many insertions.
  numbering[6] = annotate(getCdr3Annotations(cdr3.length, 'chothia', 'heavy') as Position[], cdr3);

  return [gapMissing(numbering), startIndex, endIndex];
}

/** `number_chothia_light` (also Martin light). */
export function numberChothiaLight(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX';
  const regionString = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777';
  const regionIndexDict = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6};
  const rels = [0, 0, -6, -6, -13, -16, -20];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 7, [1, 3, 4, 5]);

  const numbering: Numbering[] = [at(regions, 0), [], at(regions, 2), [], at(regions, 4), [], at(regions, 6)];

  // CDR1 (24-34): insertions on 30, deletions forward from 31.
  const cdr1 = at(regions, 1);
  numbering[1] = annotate([
    ...plain([24, 25, 26, 27, 28, 29, 30]).slice(0, Math.max(0, cdr1.length)),
    ...insertionsOn(30, Math.max(cdr1.length - 11, 0)),
    ...plain([31, 32, 33, 34]).slice(Math.abs(Math.min(0, cdr1.length - 11))),
  ], cdr1);

  numbering[3] = renumberLightCdr2(at(regions, 3)); // CDR2: insertions on 52.

  // FW3: insertions on 68, first deletion on 68; otherwise the alignment places them.
  const fw3 = at(regions, 4);
  const fw3Insertions = Math.max(fw3.length - 34, 0);
  if (fw3Insertions > 0)
    numbering[4] = annotate([...plain(range(55, 69)), ...insertionsOn(68, fw3Insertions), ...plain(range(69, 89))], fw3);
  else if (fw3.length === 33)
    numbering[4] = annotate([...plain(range(55, 68)), ...plain(range(69, 89))], fw3);
  else
    numbering[4] = fw3;


  // CDR3 (89-97): insertions on 95.
  const cdr3 = at(regions, 5);
  if (cdr3.length > 35) return [[], startIndex, endIndex]; // Too many insertions.
  numbering[5] = annotate(getCdr3Annotations(cdr3.length, 'chothia', 'light') as Position[], cdr3);

  return [gapMissing(numbering), startIndex, endIndex];
}

// ---------------------------------------------------------------------------
// Kabat

/** `number_kabat_heavy`. */
export function numberKabatHeavy(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX';
  const regionString = '11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888';
  const regionIndexDict = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7};
  const rels = [0, -1, -1, -5, -5, -8, -12, -15];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 8, [2, 4, 6]);

  const numbering: Numbering[] = [[], at(regions, 1), [], at(regions, 3), [], at(regions, 5), [], at(regions, 7)];
  numbering[0] = renumberHeavyFw1(at(regions, 0), 6, [7, 8, 9]); // Insertions on 6.

  // CDR1 (23-35): insertions on 35, deletions from 35 backwards.
  const cdr1 = at(regions, 2);
  numbering[2] = annotate([...plain(range(23, 36)).slice(0, cdr1.length),
    ...insertionsOn(35, Math.max(0, cdr1.length - 13))], cdr1);

  numbering[4] = renumberHeavyCdr2(at(regions, 4)); // Insertions on 52.
  // FW3: the alignment places insertions.
  // CDR3 (93-102): insertions on 100 (as Chothia).
  const cdr3 = at(regions, 6);
  if (cdr3.length > 36) return [[], startIndex, endIndex]; // Too many insertions.
  numbering[6] = annotate(getCdr3Annotations(cdr3.length, 'kabat', 'heavy') as Position[], cdr3);

  return [gapMissing(numbering), startIndex, endIndex];
}

/** `number_kabat_light`. */
export function numberKabatLight(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX';
  const regionString = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777';
  const regionIndexDict = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6};
  const rels = [0, 0, -6, -6, -13, -16, -20];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 7, [1, 3, 5]);

  const numbering: Numbering[] = [at(regions, 0), [], at(regions, 2), [], at(regions, 4), [], at(regions, 6)];

  // CDR1 (24-34): insertions on 27, deletions forward from 28.
  const cdr1 = at(regions, 1);
  numbering[1] = annotate([
    ...plain([24, 25, 26, 27]).slice(0, Math.max(0, cdr1.length)),
    ...insertionsOn(27, Math.max(cdr1.length - 11, 0)),
    ...plain([28, 29, 30, 31, 32, 33, 34]).slice(Math.abs(Math.min(0, cdr1.length - 11))),
  ], cdr1);

  numbering[3] = renumberLightCdr2(at(regions, 3)); // CDR2: insertions on 52.
  // FW3: all insertions are placed by the alignment.
  // CDR3 (89-97): insertions on 95.
  const cdr3 = at(regions, 5);
  if (cdr3.length > 35) return [[], startIndex, endIndex]; // Too many insertions.
  numbering[5] = annotate(getCdr3Annotations(cdr3.length, 'kabat', 'light') as Position[], cdr3);

  return [gapMissing(numbering), startIndex, endIndex];
}

// ---------------------------------------------------------------------------
// Martin (extended Chothia)

/** `number_martin_heavy`. */
export function numberMartinHeavy(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX';
  const regionString = '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888';
  const regionIndexDict = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7};
  const rels = [0, -1, -1, -5, -5, -8, -12, -15];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 8, [2, 4, 5, 6]);

  const numbering: Numbering[] = [[], at(regions, 1), [], at(regions, 3), [], at(regions, 5), [], at(regions, 7)];
  numbering[0] = renumberHeavyFw1(at(regions, 0), 8, [9]); // Insertions on 8.
  numbering[2] = renumberChothiaHeavyCdr1(at(regions, 2)); // Insertions on 31.
  numbering[4] = renumberHeavyCdr2(at(regions, 4)); // Insertions on 52.

  // FW3 (58-92): all insertions on 72; the alignment places gaps.
  const fw3 = at(regions, 5);
  const fw3Insertions = Math.max(fw3.length - 35, 0);
  if (fw3Insertions > 0)
    numbering[5] = annotate([...plain(range(58, 73)), ...insertionsOn(72, fw3Insertions), ...plain(range(73, 93))], fw3);
  else {
    // Upstream writes `_numbering[4] = _regions[4]` here, replacing the CDR2
    // renumbering above with the alignment's numbering; kept for parity.
    numbering[4] = at(regions, 4);
  }

  // CDR3 (93-102): insertions on 100.
  const cdr3 = at(regions, 6);
  if (cdr3.length > 36) return [[], startIndex, endIndex]; // Too many insertions.
  numbering[6] = annotate(getCdr3Annotations(cdr3.length, 'chothia', 'heavy') as Position[], cdr3);

  return [gapMissing(numbering), startIndex, endIndex];
}

/** `number_martin_light`: identical to Chothia light. */
export function numberMartinLight(stateVector: StateVector, sequence: string): SchemeDomain {
  return numberChothiaLight(stateVector, sequence);
}

// ---------------------------------------------------------------------------
// Wolfguy

/** `number_wolfguy_heavy`: heavy chains numbered 101-499; no gap filling. */
export function numberWolfguyHeavy(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX';
  const regionString = '11111111111111111111111111222222222222223333333333333344444444444444444444555555555555555555555555555555666666666666677777777777';
  const regionIndexDict = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6};
  const rels = [100, 124, 160, 196, 226, 244, 283];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 7, [1, 3, 5]);

  const numbering: Numbering[] = [at(regions, 0), [], at(regions, 2), [], at(regions, 4), [], at(regions, 6)];

  // CDRH1 (151-199): delete symmetrically about 175/176, right first.
  let orderedDeletions = [151];
  for (const [p1, p2] of zip(range(152, 176), range(199, 175, -1))) orderedDeletions.push(p1, p2);
  const cdr1 = at(regions, 1);
  numbering[1] = annotateNumbers(sortedNumbers(orderedDeletions.slice(0, cdr1.length)), cdr1);

  // CDRH2 (251-299): delete symmetrically about 271, right first; then right from 290.
  orderedDeletions = [251];
  for (const [p1, p2] of zip(range(252, 271), range(290, 271, -1))) orderedDeletions.push(p1, p2);
  orderedDeletions.push(271);
  orderedDeletions = [...range(299, 290, -1), ...orderedDeletions];
  const cdr2 = at(regions, 3);
  numbering[3] = annotateNumbers(sortedNumbers(orderedDeletions.slice(0, cdr2.length)), cdr2);

  // CDRH3 (331, 332, 351-399): delete symmetrically about 374, right first.
  orderedDeletions = [];
  for (const [p1, p2] of zip(range(356, 374), range(391, 373, -1))) orderedDeletions.push(p1, p2);
  orderedDeletions = [354, 394, 355, 393, 392, ...orderedDeletions];
  orderedDeletions = [331, 332, 399, 398, 351, 352, 397, 353, 396, 395, ...orderedDeletions];
  const cdr3 = at(regions, 5);
  if (cdr3.length > orderedDeletions.length) return [[], startIndex, endIndex]; // Too many insertions.
  numbering[5] = annotateNumbers(sortedNumbers(orderedDeletions.slice(0, cdr3.length)), cdr3);

  return [concatAll(numbering), startIndex, endIndex];
}

/** `number_wolfguy_light`: light chains numbered 501-799; no gap filling. */
export function numberWolfguyLight(stateVector: StateVector, sequence: string): SchemeDomain {
  const stateString = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX';
  const regionString = '1111111AAABBBBBBBBBBBBB222222222222222223333333333333334444444444444455555555555666677777777777777777777888888888888899999999999';
  const regionIndexDict = {'1': 0, 'A': 1, 'B': 2, '2': 3, '3': 4, '4': 5, '5': 6, '6': 7, '7': 8, '8': 9, '9': 10};
  const rels = [500, 500, 500, 527, 560, 595, 631, 630, 630, 646, 683];
  const [regions, startIndex, endIndex] = numberRegions(sequence, stateVector, stateString, regionString,
    regionIndexDict, rels, 11, [1, 3, 5, 7, 9]);

  const numbering: Numbering[] = [at(regions, 0), [], at(regions, 2), [], at(regions, 4), [], at(regions, 6), [],
    at(regions, 8), [], at(regions, 10)];

  // FW1 gaps go on 508 (instead of the IMGT 510 equivalent).
  const fw1 = at(regions, 1);
  numbering[1] = annotate(sortedPositions([...plain([510, 509, 508]).slice(0, fw1.length),
    ...ALPHABET.slice(0, Math.max(0, fw1.length - 3)).map((a): Position => [508, a])]), fw1);

  // CDRL1 (551-599): by canonical class, predicted from length and sequence.
  const cdr1 = at(regions, 3);
  numbering[3] = annotateNumbers(getWolfguyL1(cdr1, cdr1.length), cdr1);

  // CDRL2 (651-699): delete about 673, then right from 694; 651 last.
  let orderedDeletions: number[] = [];
  for (const [p1, p2] of zip(range(652, 673), range(694, 672, -1))) orderedDeletions.push(p2, p1);
  orderedDeletions = [651, ...range(699, 694, -1), ...orderedDeletions, 673];
  const cdr2 = at(regions, 5);
  numbering[5] = annotateNumbers(sortedNumbers(orderedDeletions.slice(0, cdr2.length)), cdr2);

  // FW3 indel on 713/714 (IMGT places it differently).
  const fw3 = at(regions, 7);
  const insertions = Math.max(0, fw3.length - 4);
  numbering[7] = annotate([...plain([711, 712, 713, 714]).slice(0, fw3.length),
    ...ALPHABET.slice(0, insertions).map((a): Position => [714, a])], fw3);

  // CDRL3 (751-799): delete symmetrically about 775, right first.
  orderedDeletions = [];
  for (const [p1, p2] of zip(range(751, 775), range(799, 775, -1))) orderedDeletions.push(p1, p2);
  orderedDeletions.push(775);
  const cdr3 = at(regions, 9);
  if (cdr3.length > orderedDeletions.length) return [[], startIndex, endIndex]; // Too many insertions.
  numbering[9] = annotateNumbers(sortedNumbers(orderedDeletions.slice(0, cdr3.length)), cdr3);

  return [concatAll(numbering), startIndex, endIndex];
}

/** Wolfguy CDRL1 canonical forms by length: name, consensus, positions. */
const WOLFGUY_L1: Readonly<Record<number, readonly [string, string, readonly number[]][]>> = {
  9: [['9', 'XXXXXXXXX', [551, 552, 554, 556, 563, 572, 597, 598, 599]]],
  10: [['10', 'XXXXXXXXXX', [551, 552, 553, 556, 561, 562, 571, 597, 598, 599]]],
  11: [['11a', 'RASQDISSYLA', [551, 552, 553, 556, 561, 562, 571, 596, 597, 598, 599]],
    ['11b', 'GGNNIGSKSVH', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]],
    ['11b.2', 'SGDQLPKKYAY', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]]],
  12: [['12a', 'TLSSQHSTYTIE', [551, 552, 553, 554, 555, 556, 561, 563, 572, 597, 598, 599]],
    ['12b', 'TASSSVSSSYLH', [551, 552, 553, 556, 561, 562, 571, 595, 596, 597, 598, 599]],
    ['12c', 'RASQSVxNNYLA', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]],
    ['12d', 'rSShSIrSrrVh', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]]],
  13: [['13a', 'SGSSSNIGNNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 597, 598, 599]],
    ['13b', 'TRSSGSLANYYVQ', [551, 552, 553, 554, 556, 561, 562, 563, 571, 572, 597, 598, 599]]],
  14: [['14a', 'RSSTGAVTTSNYAN', [551, 552, 553, 554, 555, 561, 562, 563, 564, 571, 572, 597, 598, 599]],
    ['14b', 'TGTSSDVGGYNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 596, 597, 598, 599]]],
  15: [['15', 'XXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 594, 595, 596, 597, 598, 599]]],
  16: [['16', 'XXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 594, 595, 596, 597, 598, 599]]],
  17: [['17', 'XXXXXXXXXXXXXXXXX',
    [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 584, 594, 595, 596, 597, 598, 599]]],
};

/** `_get_wolfguy_L1`: positions of the best-scoring (BLOSUM62, first maximum)
 * canonical form of this length, else symmetric about the loop's middle. */
function getWolfguyL1(seq: Numbering, length: number): number[] {
  if (Object.hasOwn(WOLFGUY_L1, length)) {
    let currMax: [readonly [string, string, readonly number[]] | null, number] = [null, -10000];
    for (const canonical of WOLFGUY_L1[length] ?? []) {
      let subScore = 0;
      for (let i = 0; i < length; i++) {
        const a = at(seq, i)[1].toUpperCase();
        const b = at(canonical[1], i, 'string').toUpperCase();
        // Python: try blosum62[(a, b)], on KeyError blosum62[(b, a)] (which may raise).
        const score = blosum62(a, b) ?? blosum62(b, a);
        if (score === undefined) throw new PythonError('KeyError', `(${repr(b)}, ${repr(a)})`);
        subScore += score;
      }
      if (subScore > currMax[1]) currMax = [canonical, subScore];
    }
    return [...notNone(currMax[0])[2]];
  }
  const orderedDeletions: number[] = [];
  for (const [p1, p2] of zip(range(551, 575), range(599, 575, -1))) orderedDeletions.push(p2, p1);
  orderedDeletions.push(575);
  return sortedNumbers(orderedDeletions.slice(0, length));
}

// ---------------------------------------------------------------------------
// Gaps and CDR3 annotations

/** `gap_missing`: concatenates the regions and fills every skipped number with
 * a gap (all schemes except Wolfguy are numbered continuously). */
export function gapMissing(numbering: readonly Numbering[]): Numbering {
  const num: Numbering = [[[0, ' '], '-']];
  for (const [p, a] of concatAll(numbering)) {
    const last = (num[num.length - 1] as [Position, string])[0][0];
    if (p[0] > last + 1) for (let i = last + 1; i < p[0]; i++) num.push([[i, ' '], '-']);
    num.push([p, a]);
  }
  return num.slice(1);
}

/** `get_cdr3_annotations` (deprecated upstream; used by Chothia, Kabat and Martin). */
export function getCdr3Annotations(length: number, scheme = 'imgt', chainType = ''): (Position | null)[] {
  const az = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  const za = 'ZYXWVUTSRQPONMLKJIHGFEDCBA';

  if (scheme === 'imgt') {
    const start = 105;
    const end = 118; // Start inclusive, end exclusive.
    const annotations: (Position | null)[] = Array<Position | null>(Math.max(length, 13)).fill(null);
    let front = 0;
    let back = -1;
    check(length - 13 < 50, 'Too many insertions for numbering scheme to handle');
    for (let i = 0; i < Math.min(length, 13); i++) {
      if (i % 2) {
        setAt(annotations, back, [end + back, ' ']);
        back -= 1;
      } else {
        setAt(annotations, front, [start + front, ' ']);
        front += 1;
      }
    }
    for (let i = 0; i < Math.max(0, length - 13); i++) { // Insertions onto 111 and 112 in turn.
      if (i % 2) {
        setAt(annotations, back, [112, at(za, back + 6, 'string')]);
        back -= 1;
      } else {
        setAt(annotations, front, [111, at(az, front - 7, 'string')]);
        front += 1;
      }
    }
    return annotations;
  }
  if ((scheme === 'chothia' || scheme === 'kabat') && chainType === 'heavy') { // Number forwards from 93.
    const insertions = Math.max(length - 10, 0);
    check(insertions < 27, 'Too many insertions for numbering scheme to handle');
    const orderedDeletions = plain([100, 99, 98, 97, 96, 95, 101, 102, 94, 93]);
    return sortedPositions([...orderedDeletions.slice(Math.max(0, 10 - length)),
      ...Array.from(az.slice(0, insertions), (a): Position => [100, a])]);
  }
  if ((scheme === 'chothia' || scheme === 'kabat') && chainType === 'light') { // Number forwards from 89.
    const insertions = Math.max(length - 9, 0);
    check(insertions < 27, 'Too many insertions for numbering scheme to handle');
    const orderedDeletions = plain([95, 94, 93, 92, 91, 96, 97, 90, 89]);
    return sortedPositions([...orderedDeletions.slice(Math.max(0, 9 - length)),
      ...Array.from(az.slice(0, insertions), (a): Position => [95, a])]);
  }
  throw new AnarciAssertion('Unimplemented scheme');
}

// ---------------------------------------------------------------------------
// anarci.py

/** `number_sequence_from_alignment`: numbers one domain's state vector with a
 * scheme. `chainType` selects heavy/light functions (Python `chain_type in "KL"`
 * is a substring test) and the AHo CDR1 gap order. */
export function numberSequenceFromAlignment(stateVector: StateVector, sequence: string, scheme: Scheme | string,
  chainType: ChainType | string): NumberedDomain {
  const lowered = scheme.toLowerCase();
  const unimplemented = (): AnarciAssertion =>
    new AnarciAssertion(`Unimplemented numbering scheme ${lowered} for chain ${chainType}`);
  const isLight = (): boolean => 'KL'.includes(chainType);
  let result: SchemeDomain;
  if (lowered === 'imgt')
    result = numberImgt(stateVector, sequence);
  else if (lowered === 'chothia') {
    if (chainType === 'H') result = numberChothiaHeavy(stateVector, sequence);
    else if (isLight()) result = numberChothiaLight(stateVector, sequence);
    else throw unimplemented();
  } else if (lowered === 'kabat') {
    if (chainType === 'H') result = numberKabatHeavy(stateVector, sequence);
    else if (isLight()) result = numberKabatLight(stateVector, sequence);
    else throw unimplemented();
  } else if (lowered === 'martin') {
    if (chainType === 'H') result = numberMartinHeavy(stateVector, sequence);
    else if (isLight()) result = numberMartinLight(stateVector, sequence);
    else throw unimplemented();
  } else if (lowered === 'aho')
    result = numberAho(stateVector, sequence, chainType);
  else if (lowered === 'wolfguy') {
    if (chainType === 'H') result = numberWolfguyHeavy(stateVector, sequence);
    else if (isLight()) result = numberWolfguyLight(stateVector, sequence);
    else throw unimplemented();
  } else
    throw unimplemented();

  // Start/end are null (Python None) only for a state vector without residues.
  return result as NumberedDomain;
}

/** `validate_numbering`: the numbers must not decrease and the numbered residues
 * must be a contiguous segment of the sequence. */
export function validateNumbering(domain: NumberedDomain, name: string, sequence: string): NumberedDomain {
  const [numbering, start, end] = domain;
  let last = -1;
  let nseq = '';
  for (const [[index], a] of numbering) {
    check(index >= last, `Numbering was found to decrease along the sequence ${name}. Please report.`);
    last = index;
    nseq += a.replaceAll('-', '');
  }
  check(sequence.replaceAll('-', '').includes(nseq),
    `The algorithm did not number a contiguous segment for sequence ${name}. Please report`);
  return [numbering, start, end];
}
