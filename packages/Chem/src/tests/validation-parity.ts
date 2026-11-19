/**
 * `validateMolecule`'s messages in Crux mode (RDKit 2026.03's lines, Crux's error log) against RDKit mode's (Chem's
 * RDKit.js 2024.09): equal after undoing exactly the differences between the two RDKit versions that crux-core's
 * validation corpus shows (11,998 strings; lane ENGINE's differences.md of run19, 13 kinds, numbered as there) and
 * nothing else, plus Crux's own departures from RDKit 2026.03 that are still open ({@link CRUX_DEPARTURES}). When RDKit
 * goes, the RDKit leg goes and so does this file (engine-utils.ts' header).
 */

/** A message's lines, its "[HH:MM:SS] " prefixes stripped (both sides). */
function lines(message: string | null): string[] {
  return (message ?? '').replace(/^\[\d\d:\d\d:\d\d\] /gm, '').split('\n');
}

const POSITION = /^(SMILES|SMARTS) Parse Error: check for mistakes around position (\d+):$/;
const EXTRA_OPEN = /^(SMILES|SMARTS) Parse Error: extra open parentheses while parsing: /;
const WHILE_PARSING = /^(SMILES|SMARTS) Parse Error: .* while parsing: /;
const FAILED = /^(SMILES|SMARTS) Parse Error: Failed parsing (?:SMILES|SMARTS) '(.*)' for input: '(.*)'$/;
const TWICE = /^(Can't kekulize mol\. {2}Unkekulized atoms: [\d ]+|non-ring atom \d+ marked aromatic)$/;
const CHIRAL = /^SMILES Parse Error: Invalid chiral specification on atom \d+ for input: /;
const VALENCE = /^Explicit valence for atom # (\d+) (\w+), \d+, is greater than permitted$/;
const INVARIANT_FILE = /^Violation occurred on line (\d+) in file \/scratch\/[^/]+\/src\/rdkit\/Code\/GraphMol\/(\S+)$/;

/** 5: the elements 2026.03 refuses a bond to (or an H on) where 2024.09 read the molecule: He Ne Ar Kr; F- Cl- Br- O2-
 * H-; Be2+ B3+ Al3+. */
const VALENCE_2026 = new Set(['He', 'Ne', 'Ar', 'Kr', 'F', 'Cl', 'Br', 'O', 'H', 'Be', 'B', 'Al']);

/** 8, 9, 10: an invariant's "Violation occurred" line, RDKit's file at its line in 2024.09's build -> in 2026.03.6's
 * (the files' lines the corpus shows; PeriodicTable.h's did not move). */
const INVARIANT_LINES: {[file: string]: {[line2024: number]: number}} = {
  'ROMol.cpp': {209: 202, 231: 224, 416: 409, 418: 411, 422: 415},
  'PeriodicTable.h': {93: 93, 177: 177},
  'Atom.cpp': {311: 323},
};

/** 10: the Atom.cpp precondition of a valence of 128 to 255, 2024.09's text -> 2026.03's. */
const ATOM_PRECONDITION: {[line2024: string]: string} = {
  'getExplicitValence() called without call to calcExplicitValence()':
    'getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()',
  'Failed Expression: d_explicitValence > -1': 'Failed Expression: (which == ValenceType::IMPLICIT || d_explicitValence > -1)',
};

/** How a listed string departs: `message`, validateMolecule's answer differs from RDKit 2026.03's (excused in
 * {@link compareValidation}); `read`, the answer is RDKit's, only which of getMolSafe's reads gives a molecule differs. */
export interface CruxDeparture {
  departs: 'message' | 'read';
  what: string;
}

/**
 * Crux's departures from RDKit 2026.03 still open in the vendored chem-wasm: the input -> how it departs. None: at
 * crux-core e881b8e98 all of the corpus' 11,998 strings log and read as RDKit-js 2026.03.6 (the owner's answer 4 of
 * 2026-11-19 fixed 25 of 7eef2854b's 26; the dative SMARTS bond of `[NH3]-->[Pt+2]`, which get_qmol refused, reads as
 * RDKit's since crux-core's CXSMILES fields). The list is exactly what remains: the tests fail on a listed string that
 * no longer departs, and on a string that departs unlisted.
 */
export const CRUX_DEPARTURES = new Map<string, CruxDeparture>([]);

/** A 2024.09 message (RDKit mode) in 2026.03's form: RDKit's source paths, line numbers and Atom.cpp's precondition
 * text in an invariant's record (8, 9, 10: the builds' checkouts and RDKit's sources differ). */
export function normalized2024(message: string | null): string {
  return lines(message).map((l) => {
    const file = INVARIANT_FILE.exec(l);
    const line = file ? INVARIANT_LINES[file[2].replace(/^.*\//, '')]?.[+file[1]] : undefined;
    if (file && line !== undefined)
      return `Violation occurred on line ${line} in file /src/rdkit/Code/GraphMol/${file[2]}`;
    return ATOM_PRECONDITION[l] ?? l;
  }).join('\n');
}

/** RDKit 2026.03's printSyntaxErrorMessage: the input around a grammar error's position `n` (41 characters at most,
 * 20 before it) and the `~~~^` pointer under it, the two lines after "check for mistakes around position N:". */
function positionLines(input: string, n: number): [string, string] {
  const at = n - 1;
  const window = at >= 20 && at + 20 < input.length ? input.substr(at - 20, 41) :
    at >= 20 ? input.substring(at - 20) : input.substring(0, Math.min(input.length, 41));
  return [window, `${'~'.repeat(n >= 20 ? 20 : n - 1)}^`];
}

/** The position `ls[i]` names when it is a "check for mistakes around position N:" line under a "... while parsing: X"
 * line whose window and pointer follow (1), else null (the lines are kept: they differ). */
function position(ls: string[], i: number): number | null {
  const n = POSITION.exec(ls[i]);
  const input = i > 0 ? / while parsing: (.*)$/.exec(ls[i - 1])?.[1] : undefined;
  if (!n || input === undefined || +n[2] < 1)
    return null;
  const [window, pointer] = positionLines(input, +n[2]);
  return ls[i + 1] === window && ls[i + 2] === pointer ? +n[2] : null;
}

/**
 * A 2026.03 message (Crux mode) in 2024.09's form: (1) the "check for mistakes around position N:" line under a grammar
 * error, with its two lines (the input around the position, the `~~~^` pointer: checked), goes; (2) branches left open
 * are reported "extra open parentheses for input: 'X'" alone (2026.03: "... while parsing: X" for each, at its `(`,
 * left to right, then "Failed parsing ... for input: 'X'"), and not at all after another grammar error of the same
 * parse; (3, 4) "Can't kekulize mol." and "non-ring atom N marked aromatic", logged twice by each read, once. Kinds
 * 1, 3 and 4 hold for every such line of 2026.03 (the corpus: 14,237 grammar errors, 969 doubled lines): a line without them is marked, so
 * it differs from 2024.09's.
 */
export function normalized2026(message: string | null): string {
  const ls = lines(message);
  const out: string[] = [];
  for (let i = 0; i < ls.length;) {
    if (position(ls, i) !== null) {
      i += 3;
      continue;
    }
    const extra = EXTRA_OPEN.exec(ls[i]);
    if (extra) {
      // each open branch once, with its position at its `(`, left to right
      const input = ls[i].slice(extra[0].length);
      let j = i;
      let previous = 0;
      for (let n: number | null; j < ls.length && ls[j] === ls[i] && (n = position(ls, j + 1)) !== null &&
        n > previous && input[n - 1] === '('; j += 4)
        previous = n!;
      const last = out[out.length - 1];
      if (j > i && last !== undefined && last.startsWith(extra[1]) && WHILE_PARSING.test(last) &&
        last.endsWith(` while parsing: ${input}`)) {
        i = j;
        continue;
      }
      const failed = j > i && j < ls.length ? FAILED.exec(ls[j]) : null;
      if (failed && failed[1] === extra[1] && failed[2] === input) {
        out.push(`${extra[1]} Parse Error: extra open parentheses for input: '${failed[3]}'`);
        i = j + 1;
        continue;
      }
    }
    if (WHILE_PARSING.test(ls[i]) && position(ls, i + 1) === null) {
      out.push(`${ls[i]} [no position under it]`);
      i++;
      continue;
    }
    if (TWICE.test(ls[i])) {
      out.push(ls[i + 1] === ls[i] ? ls[i] : `${ls[i]} [logged once]`);
      i += ls[i + 1] === ls[i] ? 2 : 1;
      continue;
    }
    out.push(ls[i]);
    i++;
  }
  return out.join('\n');
}

/** A message's non-empty lines. */
const filled = (message: string) => message.split('\n').filter((l) => l !== '');

/** 12: CXSmilesOps.cpp's precondition record, as each of the three reads logs it. */
const SGROUP_RECORD = '\n\n****\nPre-condition Violation\nparse_data_sgroup_attr: first >= last\nViolation occurred on line 742 in ' +
  'file /src/rdkit/Code/GraphMol/SmilesParse/CXSmilesOps.cpp\nFailed Expression: first < last\n****\n\n';

/** The three reads' grammar error (normalised: no position) for an input `x` as RDKit names it in "while parsing". */
function grammarError(error: string, x: string, failed: (parser: string) => string): string[] {
  return ['SMILES', 'SMILES', 'SMARTS'].flatMap((p) => [`${p} Parse Error: ${error} while parsing: ${x}`,
    `${p} Parse Error: Failed parsing ${p} ${failed(p)}`]);
}

/**
 * Where the two versions differ in what they accept, so that no text normalisation applies, the kind of difference
 * (null for none), on the normalised messages, each in the exact shape the corpus shows: (5) 2026.03 refuses a valence
 * 2024.09 let pass (both get_mol reads log it, get_qmol reads the query); (13) so such an atom fails before the one
 * 2024.09 named; (6) the SMARTS read refuses a chirality permutation out of range too (2024.09 read the query); (7) a
 * chiral permutation 0 (`@TH0`) is a grammar error of the three reads; (11) so is a NUL inside the
 * input (the reads see the text before it); (12) a data S-group with too few fields logs a precondition at each read,
 * whose reads failed with no message in 2024.09 ('unknown error').
 */
export function acceptanceDifference(rdkit2024: string, crux2026: string): string | null {
  const [a, b] = [filled(rdkit2024), filled(crux2026)];
  const valences = (ls: string[]) => ls.length === 2 && ls[0] === ls[1] ? VALENCE.exec(ls[0]) : null;
  const [va, vb] = [valences(a), valences(b)];
  if (vb && VALENCE_2026.has(vb[2])) {
    if (a.length === 0)
      return `5: a valence 2026.03 refuses (${vb[2]})`;
    if (va && +va[1] > +vb[1])
      return `13: an atom whose valence 2026.03 refuses fails first (${vb[2]})`;
  }
  if (a.length === 2 && a[0] === a[1] && CHIRAL.test(a[0]) && b.length === 3 && b[0] === a[0] && b[1] === a[0] &&
    b[2] === a[0].replace(/^SMILES/, 'SMARTS'))
    return '6: a chirality permutation out of range, the SMARTS read too';
  const x = / while parsing: (.*)$/.exec(b[0] ?? '')?.[1] ?? '';
  if (a.length === 0 && b.join('\n') ===
    grammarError('chiral permutation cannot be zero', x, () => `'${x}' for input: '${x}'`).join('\n'))
    return '7: a chiral permutation 0';
  const nul = /^SMILES Parse Error: Failed parsing SMILES '.* for input: '(.*\u0000.*)'$/.exec(b[1] ?? '')?.[1];
  if (a.length === 0 && nul !== undefined && b.join('\n') ===
    grammarError('syntax error', x, () => `'${x} for input: '${nul}'`).join('\n'))
    return '11: a NUL inside the input';
  if (rdkit2024 === 'unknown error' && crux2026 === SGROUP_RECORD.repeat(3))
    return '12: a data S-group with too few fields';
  return null;
}

/** How a string's two messages compare: 'same', 'normalised', an acceptance difference, a listed departure of Crux
 * from RDKit 2026.03 in its message (`molecule` in {@link CRUX_DEPARTURES}, `departs: 'message'`), or null (they
 * differ). Equal messages are 'same' only
 * where 2026.03 and 2024.09 write the same (a grammar error or a doubled line in 2024.09's form differs). */
export function compareValidation(rdkit2024: string | null, crux2026: string | null, molecule?: string): string | null {
  const [a, b] = [normalized2024(rdkit2024), normalized2026(crux2026)];
  if (a === b)
    return (rdkit2024 ?? '') === (crux2026 ?? '') ? 'same' : 'normalised';
  const kind = acceptanceDifference(a, b);
  if (kind !== null)
    return kind;
  return molecule !== undefined && CRUX_DEPARTURES.get(molecule)?.departs === 'message' ?
    'a departure of Crux, listed' : null;
}
