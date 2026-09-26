/** A value as the text a field shows: null and undefined are empty, everything else is String(). */
export function text(value: unknown): string {
  return value === null || value === undefined ? '' : String(value);
}

/** What a required field with nothing in it says — matched by text where a form words the same
 * refusal as "<caption> is required". */
export const EMPTY = 'Value can\'t be empty';

export function isEmpty(value: unknown): boolean {
  return value === null || value === undefined || value === '';
}

/** A count and the word it agrees with: `plural(2, 'row', 'rows')` is "2 rows". The two forms are
 * whole phrases, so a verb agrees with them too — `plural(1, 'row has errors', 'rows have errors')`. */
export function plural(n: number, one: string, many: string): string {
  return `${n.toLocaleString()} ${n === 1 ? one : many}`;
}

const CONJUNCTIONS = new Set(['and', 'or', 'than', 'if', 'but', 'so', 'as', 'that']);

/** prop_gen's `camelCaseToWords` (`prop_gen_annotation.dart:89-112`): all-caps and already-spaced
 * names pass through, humps split, first word capitalized, conjunctions lowercased. */
export function camelCaseToWords(name: string): string {
  if (name === name.toUpperCase() || name.includes(' '))
    return name;
  const words = name.match(/[A-Z]+(?![a-z])|[A-Z]?[^A-Z]+/g) ?? [name];
  return words
    .map((w, i) => {
      const word = CONJUNCTIONS.has(w.toLowerCase()) ? w.toLowerCase() : w;
      return i === 0 ? word.charAt(0).toUpperCase() + word.slice(1) : word;
    })
    .join(' ');
}
