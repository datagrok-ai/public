// (a leaf: the dosing schedules format their totals with it, and import nothing)

/** `x` with every object's keys in sorted order, so that its JSON text depends on its content alone. */
export const canonical = (x: unknown): unknown => Array.isArray(x) ? x.map(canonical) : x && typeof x === 'object' ?
  Object.fromEntries(Object.keys(x).sort().map((k) => [k, canonical((x as {[k: string]: unknown})[k])])) : x;
/** Whether two values hold the same content, whatever the order of their keys (a document read back from its row's
 *  json has them reordered). */
export const sameJson = (a: unknown, b: unknown): boolean =>
  JSON.stringify(canonical(a)) === JSON.stringify(canonical(b));

/** A measured number as it is read: three significant figures (or `digits`), thousands separated ("12,300"), with
 *  no float noise (a float column's 0.3 reads back as 0.30000001192092896: "0.3"); not a number: "—". Every exposure,
 *  AUC, Cmax, concentration, dose and fraction the app prints goes through it. */
export function sig(x: number, digits = 3): string {
  if (typeof x !== 'number' || !Number.isFinite(x))
    return '—';
  const r = Number(x.toPrecision(digits));
  return Math.abs(r) >= 1000 ? r.toLocaleString('en-US', {maximumFractionDigits: 20}) : String(r);
}

/** A value as entered (a dose): every digit kept (3.125, 112.5), only a 32-bit column's float noise off, thousands
 *  separated. sig rounds what was measured or derived; never what the user typed. */
export const entered = (x: number) => sig(x, 7);

/** The mean; NaN for no values (every caller checks the count first, or shows NaN as no value). */
export const mean = (xs: number[]): number => xs.reduce((t, x) => t + x, 0) / xs.length;
/** The median (the mean of the middle two for an even count); NaN for no values. */
export const median = (xs: number[]): number => {
  if (xs.length === 0)
    return NaN;
  const s = [...xs].sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
};

/** Math.max / Math.min over a list of any length: spread into a call, a table's 120,000 rows overflowed the stack
 *  (QA 0.1.3 round 4 technician R4-1). Empty: -Infinity / Infinity; a NaN in it: NaN, as Math's. */
export function maxOf(xs: Iterable<number>): number {
  let m = -Infinity;
  for (const x of xs) {
    if (Number.isNaN(x))
      return NaN;
    if (x > m)
      m = x;
  }
  return m;
}
export function minOf(xs: Iterable<number>): number {
  let m = Infinity;
  for (const x of xs) {
    if (Number.isNaN(x))
      return NaN;
    if (x < m)
      m = x;
  }
  return m;
}
/** A number's text with a true minus wherever it is signed ("-64%" → "−64%"; "−0" never shown: callers round
 *  first). Every signed figure the app prints goes through it (QA 0.1.3 round 4, analyst U3). */
export const minus = (s: string) => s.replace(/-(?=\d)/g, '−');
