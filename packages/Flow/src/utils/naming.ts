/** TS mirror of core's `propertyNameToFriendly` (prop_gen), so caption-less params read the same as
 *  in `ui.input.forProperty`. One DELIBERATE deviation: an ALL-CAPS word is kept as-is
 *  ('MW' → 'MW', 'maxMW' → 'Max MW') where core folds to 'Mw'. */

const CONJUNCTIONS = new Set(['and', 'or', 'than', 'if', 'but', 'so', 'as', 'that']);
const DELIMITERS = new Set([' ', ',', '.', '-', '!', '"', '\'', ';']);

function isUpperAt(s: string, i: number): boolean {
  const c = s.charCodeAt(i);
  return c >= 65 && c <= 90;
}

/** prop_gen `splitCamelCase`: 'RDKitMol' → RD, Kit, Mol. */
export function splitCamelCase(s: string): string[] {
  const parts: string[] = [];
  for (let start = 0, i = 1; i <= s.length; i++) {
    const lowerToUpper = i < s.length && isUpperAt(s, i) && !isUpperAt(s, i - 1);
    if (i === s.length || lowerToUpper || (isUpperAt(s, i) && i < s.length - 1 && !isUpperAt(s, i + 1))) {
      parts.push(s.substring(start, i));
      start = i;
    }
  }
  return parts;
}

/** prop_gen `camelCaseToWords` with its default options. */
export function camelCaseToWords(s: string): string {
  if (s.toUpperCase() === s) return s;
  if (s.includes(' ')) return s; // already a friendly name
  let out = '';
  let prev = '';
  for (let word of splitCamelCase(s)) {
    if (CONJUNCTIONS.has(word.toLowerCase().trim())) word = word.toLowerCase();
    if (out.length === 0)
      out = word === '' ? word : word[0].toUpperCase() + word.substring(1);
    else {
      if (!prev.endsWith(' ')) out += ' ';
      out += word;
    }
    prev = word;
  }
  return out;
}

/** prop_gen `capitalizeWords`, except an all-caps acronym word is kept exactly as-is. */
export function capitalizeWords(s: string): string {
  let out = '';
  let word = '';
  const flush = (): void => {
    if (word === '') return;
    const isAcronym = word === word.toUpperCase() && word !== word.toLowerCase();
    out += isAcronym ? word : word[0].toUpperCase() + word.substring(1).toLowerCase();
    word = '';
  };
  for (const ch of s) {
    if (DELIMITERS.has(ch)) {
      flush();
      out += ch;
    } else
      word += ch;
  }
  flush();
  return out;
}

/** core `propertyNameToFriendly`: 'maxNumOfSomething' → 'Max Num Of Something'. */
export function propertyNameToFriendly(n: string): string {
  return capitalizeWords(camelCaseToWords(String(n ?? '').replace(/\./g, ' ')));
}
