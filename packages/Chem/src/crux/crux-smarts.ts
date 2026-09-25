import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {elementsTable} from '../constants';

const ELEMENTS = new Set(elementsTable);
const AROMATIC_SYMBOLS = new Set(['b', 'c', 'n', 'o', 'p', 's', 'se', 'as', 'te', 'si', 'ge', 'sb', 'bi']);
const ORGANIC_SUBSET = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']);
// how RDKit writes a non-query (molecule) atom: [isotope? (#Z | symbol) chirality? H-count? charge? map?]
// eslint-disable-next-line max-len
const NON_QUERY_ATOM = /^(\d*)(#\d+|[A-Z][a-z]?|[a-z]{1,2})(?:@@?(?:TH[12]|AL[12]|SP[1-3]|TB\d{1,2}|OH\d{1,2})?)?(?:H\d*)?([+-]\d*|\+{2,}|-{2,})?(?::\d+)?$/;

/**
 * SMARTS for crux that matches the way Chem's RDKit search matches `queryMol` (as built by getQueryMolSafe),
 * or null when crux cannot express it.
 *
 * RDKit matches molecule (non-query) atoms by element, charge, isotope and radicals only, and the search runs
 * without chirality, while `get_smarts` also writes H counts and stereo. Those are dropped here; radicals,
 * isotopes and primitives crux does not evaluate like RDKit (v, x, h, Rn, z, ^) make the query unsupported.
 */
export function getCruxSmarts(queryMol: RDMol): string | null {
  return hasRadicals(queryMol) ? null : normalizeSmarts(queryMol.get_smarts());
}

function hasRadicals(mol: RDMol): boolean {
  try {
    return /"nRad":[1-9]/.test(mol.get_json());
  } catch {
    return false;
  }
}

function normalizeSmarts(smarts: string): string | null {
  let out = '';
  for (let i = 0; i < smarts.length;) {
    const ch = smarts[i];
    if (ch === '[') {
      const end = matchingClose(smarts, i, '[', ']');
      const atom = end < 0 ? null : normalizeAtom(smarts.slice(i + 1, end));
      if (atom === null)
        return null;
      out += `[${atom}]`;
      i = end + 1;
    } else if (ch === '/' || ch === '\\') {
      out += '-';
      i += smarts[i + 1] === '?' ? 2 : 1;
    } else {
      // crux parses the any-aromatic/any-aliphatic atoms only in brackets
      out += ch === 'a' || ch === 'A' ? `[${ch}]` : ch;
      i++;
    }
  }
  return out;
}

function normalizeAtom(content: string): string | null {
  const m = content.match(NON_QUERY_ATOM);
  const symbol = m?.[2] ?? '';
  if (!m || (symbol[0] !== '#' && !ELEMENTS.has(symbol) && !AROMATIC_SYMBOLS.has(symbol)))
    return normalizeQueryAtom(content);
  if (m[1])
    return null;
  // RDKit matches a molecule atom regardless of aromaticity; lowercase symbols keep their aromatic meaning
  const element = symbol[0] === '#' || ORGANIC_SUBSET.has(symbol) || AROMATIC_SYMBOLS.has(symbol) ?
    symbol : `#${elementsTable.indexOf(symbol) + 1}`;
  return element + (m[3] ?? '');
}

function normalizeQueryAtom(c: string): string | null {
  let out = '';
  for (let i = 0; i < c.length;) {
    const ch = c[i];
    const pair = c.slice(i, i + 2);
    let end = i + 1;
    if ('&;,!*'.includes(ch))
      out += ch;
    else if (ch === '$') {
      end = c[i + 1] === '(' ? matchingClose(c, i + 1, '(', ')') : -1;
      const inner = end < 0 ? null : normalizeSmarts(c.slice(i + 2, end));
      if (inner === null)
        return null;
      out += `$(${inner})`;
      end++;
    } else if (pair.length === 2 && (ELEMENTS.has(pair) || AROMATIC_SYMBOLS.has(pair))) {
      out += pair;
      end = i + 2;
    } else if (ch === '#' || ch === 'H' || ch === 'D' || ch === 'X' || ch === 'r') {
      end = readDigits(c, i + 1);
      if (ch === '#' && end === i + 1)
        return null;
      out += c.slice(i, end);
    } else if (ch === '@') {
      if (c[end] === '@')
        end++;
      end += c.slice(end).match(/^(TH|AL|SP|TB|OH)\d*/)?.[0].length ?? 0;
      if (c[end] === '?')
        end++;
    } else if (ch === '+' || ch === '-') {
      if (c[end] === ch)
        while (c[end] === ch) end++;
      else
        end = readDigits(c, end);
      out += c.slice(i, end);
    } else if (ch === ':')
      end = readDigits(c, end);
    else if (ch === 'R') {
      end = readDigits(c, i + 1);
      const count = c.slice(i + 1, end);
      if (count !== '' && count !== '0')
        return null;
      out += c.slice(i, end);
    } else if (ch === 'A' || ch === 'a' || ELEMENTS.has(ch) || AROMATIC_SYMBOLS.has(ch))
      out += ch;
    else
      return null;
    i = end;
  }
  return out;
}

function matchingClose(s: string, from: number, open: string, close: string): number {
  let depth = 0;
  for (let i = from; i < s.length; i++) {
    if (s[i] === open)
      depth++;
    else if (s[i] === close && --depth === 0)
      return i;
  }
  return -1;
}

function readDigits(s: string, from: number): number {
  let i = from;
  while (i < s.length && s[i] >= '0' && s[i] <= '9') i++;
  return i;
}
