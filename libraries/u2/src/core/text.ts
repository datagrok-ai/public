/** A value as the text a field shows: null and undefined are empty, everything else is String(). */
export function text(value: unknown): string {
  return value === null || value === undefined ? '' : String(value);
}
