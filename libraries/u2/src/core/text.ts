/** A value as the text a field shows: null and undefined are empty, everything else is String(). */
export function text(value: unknown): string {
  return value === null || value === undefined ? '' : String(value);
}

/** A count and the word it agrees with: `plural(2, 'row', 'rows')` is "2 rows". The two forms are
 * whole phrases, so a verb agrees with them too — `plural(1, 'row has errors', 'rows have errors')`. */
export function plural(n: number, one: string, many: string): string {
  return `${n.toLocaleString()} ${n === 1 ? one : many}`;
}
