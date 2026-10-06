/* Relative time spans — `-1w`, `2d`, `now` — shared by the filter grammar, the model validation and
   `DateInput.relative`. A resolved span is a `Date` tagged with its text (non-enumerable, so JSON and
   deep equality never see it). */
import {FilterError} from './filter/model.js';

const SPAN = /^(-?\d+)([hdwmy])$/;
const SPAN_MS: Record<string, number> = {h: 3600e3, d: 86400e3, w: 7 * 86400e3, m: 31 * 86400e3, y: 365 * 86400e3};

export function isSpanText(text: string): boolean {
  return text === 'now' || SPAN.test(text);
}

/** `now` plus a signed span — `'-1w'`, `'2d'`, `'now'` — with `m` = 31 days, `y` = 365. */
export function resolveSpan(span: string, now: Date): Date {
  if (span === 'now')
    return new Date(now.getTime());
  const m = SPAN.exec(span);
  if (!m)
    throw new FilterError(`Bad time span "${span}"`);
  return new Date(now.getTime() + Number(m[1]) * SPAN_MS[m[2]]);
}

export function markSpan(date: Date, span: string): Date {
  Object.defineProperty(date, 'span', {value: span, enumerable: false, configurable: true});
  return date;
}

/** The span a resolved date came from (`-1w`, `now`), or undefined for a plain date. */
export function spanOf(date: Date): string | undefined {
  return (date as Date & {span?: string}).span;
}
