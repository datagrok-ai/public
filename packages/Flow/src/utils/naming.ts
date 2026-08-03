/** Human-readable fallback names for raw property identifiers — a TS mirror
 *  of core's `propertyNameToFriendly`
 *  (`capitalizeWords(camelCaseToWords(n.replaceAll('.', ' ')))`, helpers in
 *  `core/shared/prop_gen/lib/prop_gen_annotation.dart`), so a caption-less
 *  param reads the same on Flow nodes / panel rows as it does in the
 *  platform's own forms (`ui.input.forProperty`).
 *
 *  One DELIBERATE deviation: core's `capitalizeWords` lowercases a word's
 *  tail, mangling acronyms ('MW' → 'Mw', 'HBA' → 'Hba'); ours keeps an
 *  ALL-CAPS word as-is ('MW' → 'MW', 'maxMW' → 'Max MW'). */

const CONJUNCTIONS = new Set(['and', 'or', 'than', 'if', 'but', 'so', 'as', 'that']);
const DELIMITERS = new Set([' ', ',', '.', '-', '!', '"', '\'', ';']);

function isUpperAt(s: string, i: number): boolean {
  const c = s.charCodeAt(i);
  return c >= 65 && c <= 90;
}

/** prop_gen `splitCamelCase`: boundaries at lower→upper and before the last
 *  upper of an acronym run followed by lowercase ('RDKitMol' → RD, Kit, Mol). */
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

/** prop_gen `camelCaseToWords` with its default options (capitalizeFirst,
 *  no capitalizeNext / capitalizeConjunctions, ' ' separator). */
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

/** prop_gen `capitalizeWords` (first char of every delimiter-separated word
 *  uppercased, the rest lowercased) — EXCEPT an all-caps word (an acronym:
 *  MW, HBA) is kept exactly as-is instead of being folded to 'Mw'/'Hba'. */
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

/** core `propertyNameToFriendly`: 'maxNumOfSomething' → 'Max Num Of Something',
 *  'ratio.split' → 'Ratio Split'. */
export function propertyNameToFriendly(n: string): string {
  return capitalizeWords(camelCaseToWords(String(n ?? '').replace(/\./g, ' ')));
}
